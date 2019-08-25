#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <time.h>
#include <zlib.h>
#include <pthread.h>
#include <fcntl.h>
#include <sys/stat.h>

#undef   DEBUG_CHECK
#undef   DEBUG_OUT
#undef   DEBUG_LINES
#undef   DEBUG_AUTO

#include "gene_core.h"

static char *Usage = "[-vsg] [-T<int(4)>] <forward:fast[aq][.gz]> [<reverse:fast[aq][.gz]>";

#define IO_BLOCK 10000000ll

#define INT_MAXLEN 10

int    VERBOSE;  //  Verbose mode?
int    QVS_OUT;  //  Include QV strings
int    GROUP;    //  Group lanes of reads
int    NTHREADS; //  # of threads to use

static int MAXOUT = 0;
static int EOR = '\n';   //  End of record marker

static int NotDNA[256] =
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1,
    1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1,
    1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  };

static int NotQV[256] =
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1
  };


/*******************************************************************************************
 *
 *  Routines to open and close a possibly VGPzip'd fasta/q file
 *
 ********************************************************************************************/

  //  File object + all associated structures for parallel processing

typedef struct
                  // Setup when opened
  { char   *path;   //  Full path name
    char   *root;   //  Ascii root name
    int64   fsize;  //  Size of file in bytes
    int     fastq;  //  Is file a fastq? (else fasta)
    int     zipd;   //  Is file VGPzip'd (has paired .vzi file)
    int     recon;  //  Is file an uncompressed regular zip-file?
    int64   zsize;  //  Size of zip index (in blocks)
    int64  *zoffs;  //  zoffs[i] = offset to compressed block i

                  // Used for block IO
    uint8  *ibufs;  //  NTHREAD buffers for compressed input
    uint8  *zbufs;  //  NTHREAD buffers for decompressed input

                  // Filled in by Check_FastQ in 1st pass
    int64  *zrecd;  //  zrecd[i] = record number at start of block i
    int    *zstat;  //  zstat[i] = automata state at start of block i
    int    *zaux;   //  zaux[i]  = length of string item spanning block i boundary (if a string)
    int    *zgrp;   //  zgrp[i]  = group number spanning block i boundary (if -g)
    int    *zgct;   //  zgct[i]  = record witin group spanning block i boundary (if -g)
    int    *ngrp;   //  ngrp[n]  = # of groups seen by thread n
    int    *rpt;    //  rpt[n]   = 1st block containing data for thread n
    int64 **gcnt;   //  gcnt[n][i] = # of entries in i'th group processed by thread n
  } Fast_Input;

  // Open a possibly VGPzip'd fasta/q file for parallel reading

static void Open_Fastq(char *arg, Fast_Input *input)
{ static char *suffix[5] = { ".fastq.gz", ".fasta.gz", ".fastq", ".fasta", ".gz" };
  static char *sufidx[5] = { ".fastq.vzi", ".fasta.vzi", "", "", ".vzi" };

  struct stat stats;
  char  *pwd, *root, *path;
  int    fid, i;
  int64  fsize;
  int    fastq, zipd, recon;
  int64  zsize = 0;
  int64 *zoffs = NULL;

  pwd = PathTo(arg);
  OPEN2(arg,pwd,root,fid,suffix,5)
  if (fid < 0)
    { fprintf(stderr,"%s: Cannot open %s as an .fast[aq] file\n",Prog_Name,arg);
      exit (1);
    }
  path = Strdup(Catenate(pwd,"/",root,suffix[i]),"Allocating full path name");
  if (i == 4)
    { fastq = (strcmp(".fastq",root+(strlen(root)-6)) == 0);
      zipd  = 1;
    }
  else
    { fastq = (i%2 == 0);
      zipd  = (i < 2);
    }

  if (fstat(fid, &stats) == -1)
    { fprintf(stderr,"%s: Cannot get stats for %s\n",Prog_Name,arg);
      exit (1);
    }
  fsize = stats.st_size;

#ifdef DEBUG_CHECK
  printf("\n%s is a %s file, fastq = %d, zipd = %d, fsize = %lld\n",arg,suffix[i],fastq,zipd,fsize);
#endif

  recon = 0;
  if (zipd)
    { int idx;

      idx = open(Catenate(pwd,"/",root,sufidx[i]),O_RDONLY);
      free(pwd);
      if (idx < 0)
        { if (VERBOSE)
            fprintf(stderr,"  File %s not VGPzip'd, decompressing\n",arg);
          system(Catenate("gunzip -k ",path,"",""));
          path[strlen(path)-3] = '\0';
          zipd  = 0;
          recon = 1;
          if (lstat(path,&stats) == -1)
            { fprintf(stderr,"%s: Cannot get stats for %s\n",Prog_Name,path);
              exit (1);
            }
          fsize = stats.st_size;
          zsize = (fsize-1)/IO_BLOCK+1;
        }
      else
        { read(idx,&zsize,sizeof(int64));
          zoffs = (int64 *) Malloc(sizeof(int64)*(zsize+2),"Allocating VGPzip index");
          if (zoffs == NULL)
            exit (1);
          read(idx,zoffs+1,sizeof(int64)*zsize);
          zoffs[0] = 0;
          zoffs[zsize+1] = zoffs[zsize];
          close(idx);
        }
    }
  else
    { zsize = (fsize-1)/IO_BLOCK+1;
      free(pwd);
    }

  close(fid);

  input->path  = path;
  input->root  = root;
  input->fsize = fsize;
  input->fastq = fastq;
  input->zipd  = zipd;
  input->recon = recon;
  input->zsize = zsize;
  input->zoffs = zoffs;
  input->zrecd = (int64 *) Malloc(sizeof(int64)*(zsize+1),"Allocating block info");
  input->zstat = (int *) Malloc(sizeof(int)*4*zsize,"Allocating block info");
  input->rpt   = (int *) Malloc(sizeof(int)*(NTHREADS+1),"Allocating first block array");
  input->ngrp  = (int *) Malloc(sizeof(int)*NTHREADS,"Allocating group count array");
  input->gcnt  = (int64 **) Malloc(sizeof(int64 *)*NTHREADS,"Allocating group count header");
  if (input->zrecd == NULL || input->zstat == NULL || input->gcnt == NULL ||
      input->ngrp == NULL || input->rpt == NULL)
    exit (1);

  input->zaux = input->zstat + zsize;
  input->zgrp = input->zaux  + zsize;
  input->zgct = input->zgrp  + zsize;

  { int n;

    for (n = 0; n < NTHREADS; n++)
      input->gcnt[n] = NULL;
  }
}

static void Close_Fastq(Fast_Input *input)
{ int n;

  if (input->recon)
    unlink(input->path);
  free(input->ngrp);
  for (n = 0; n < NTHREADS; n++)
    if (input->gcnt[n] != NULL)
      free(input->gcnt[n]-1);
  free(input->gcnt);
  free(input->rpt);
  free(input->zstat);
  free(input->zrecd);
  free(input->ibufs);
  free(input->zoffs);
  free(input->path);
  free(input->root);
}


/*******************************************************************************************
 *
 *  First pass routines that check syntax and count things
 *
 ********************************************************************************************/

//  In a first pass perform basic .fastq syntax check, count the number of
//    of reads, # of base pairs in all reads, and maximum length of a sequence
//    Count groups and record first/last group that may need to be merged with
//    the group of the previous/next block.

typedef struct
  { int64  nread;       //  # S
    int    maxlen;      //  @ S
    int64  totbp;       //  + S
    int    ngroup;      //  # g
    int64  gmaxc;       //  @ g
    int64  gtotc;       //  + g
    int64  gmaxnr;      //  % g # S
    int64  gmaxbp;      //  % g + S
  } Header_Info;

typedef struct
  { int64  nrs;    //  # of reads seen in group
    int64  bps;    //  # of base pairs seen in group
    char  *name;   //  name of the group
  } Group_Info;

typedef struct
  { int64       beg;    //  Process blocks [beg,end)
    int64       end;
    uint8      *buf;    //  block buffer
    uint8      *zuf;    //  decode buffer
    Fast_Input *inp;    //  file information
    int         genp;   //  account for P line in output estimates
    Header_Info info;   //  statistics accumulate here
    Group_Info  first;  //  statistics for first group processed
    Group_Info  last;   //  statistics for last group processed
    int64      *gcount; //  group counts for groups processed by this thread
    int         etype;  //  error and type if > 0
    int         lineno; //  lineno on which an error occurred
    int         rpt;    //  first block for which this thread processed data
    int         outmax; //  maximum output characters from a block spanned by this thread
    int         outfst; //  output characters from the first block processed by this thread
    int         outlst; //  output characters from the last block processed by this thread
  } Check_Arg;

  //  Automata states (buffer is processed char at a time)

#define UK    0    // Unknown states
#define FK    1
#define FK1   2

#define BEG   3    // Fixed length states
#define QHE   4
#define PLS   5
#define PEN   6

#define SEQ   7   // String states
#define FEQ   8
#define FNL   9
#define QVS  10
#define QHL  11

#if defined(DEBUG_AUTO) || defined(DEBUG_CHECK)

static char *Name[] =
  { "UK", "FK", "FK1", "BEG", "QHE", "PLS", "PEN", 
    "SEQ", "FEQ", "FNL", "QVS", "QHL+0", "QHL+1", "QHL+2" };

#endif

  //  Check blocks [beg,end) starting with the first complete entry and proceeding
  //     into blocks end ... if necessary to complete the last entry started in this
  //     range.

static void *check_thread(void *arg)
{ Check_Arg   *data  = (Check_Arg *) arg;
  int64        beg   = data->beg;
  int64        end   = data->end;
  uint8       *buf   = data->buf;
  uint8       *zuf   = data->zuf;
  Header_Info *info  = &(data->info);
  Fast_Input  *inp   = data->inp;
  int          genp  = data->genp;

  int          fastq = inp->fastq;
  int64       *zoffs = inp->zoffs;
  int64       *brecd = inp->zrecd;
  int         *bstat = inp->zstat;
  int         *baux  = inp->zaux;
  int         *bgrp  = inp->zgrp;
  int         *bgct  = inp->zgct;

  int     fid;
  int     slen;
  int64   lineno, recno;
  int     blen, qlen, glen;
  int     maxlen, ngroup;
  int64   nread, gmaxnr, gnr;
  int64   totbp, gmaxbp, gbp;
  int64   gtotc, gmaxc;
  int     count, lanl, lanmax;
  char   *lanegroup, *lane;
  int64  *gcount, maxgct;
  int     fs, fb;
  int     outmax, olen, outfst;

  int state;
  int nl_1, nl_2, pl_1, alive;
  int b, c;

  fid = open(inp->path,O_RDONLY);

  nread  = -1;
  totbp  = 0;
  maxlen = 0;
  ngroup = -1;
  gmaxnr = 0;    //  max # of records in a group
  gmaxbp = 0;    //  max # of bps in a group
  gtotc  = 0;    //  total chars in group names
  gmaxc  = 0;    //  longest group name
  maxgct = 0;    //  current size of gcount array
  gcount = NULL; //  # of reads in groups seen thus far

  lanmax = 1000;
  lane   = Malloc(1000,"Allocating lane buffer");

  lanegroup = NULL; //  name of current group/lane
  lineno    = 1;
  recno     = 1;
  gnr       = 0;    //  index of 1st read in current group
  gbp       = 0;    //  base pairs prior to 1st read in current group

  if (inp->zipd)
    lseek(fid,zoffs[beg],SEEK_SET);
  else
    lseek(fid,beg*IO_BLOCK,SEEK_SET);
  if (beg == 0)
    { state = BEG;
      data->rpt = 0;
    }
  else if (fastq)
    { state = UK;
      nl_1 = pl_1 = 1;
      alive = 0;
    }
  else
    state = FK;
  fs = -1;
  outfst = -1;
  outmax = 0;

#ifdef DEBUG_CHECK
  printf("\nFrom block %lld - %lld (@ %lld)\n",beg,end,lseek(fid,0,SEEK_CUR));
#endif

  while (1)

    { if (beg == inp->zsize)     //  At the end of the file, no block beyond partition end
        { if (fastq)             //   ==> should be at end of a record
            { if (state != BEG)
                { data->etype = 7;
                  return (NULL);
                }
            }
          else
            { if (state != FNL)
                { data->etype = 7;
                  return (NULL);
                }
              totbp += blen; //  For fasta > terminates internally, not \n, so must push
              if (blen > maxlen)
                maxlen = blen;
              if (fs >= 0)
                { baux[fb] = blen;
                  fs = -1;
                }
            }
          break;
        }

      if (state >= BEG)         //  Processing the start of this block: setup state info for 2nd
        { brecd[beg] = recno;   //    block by block pass
          if (state == QHL)
            bstat[beg] = QHL+count;
          else
            bstat[beg] = state;
          if (ngroup < 0)
            bgrp[beg] = 0;
          else
            bgrp[beg] = ngroup;
          bgct[beg] = nread - gnr;
          if (fs < 0 && state >= SEQ)
            { fs = state;
              fb = beg;
            }
          else
            baux[beg] = 0;
        }

#ifdef DEBUG_CHECK
      printf("  Loading block %lld: @%lld",beg,lseek(fid,0,SEEK_CUR));
#endif
      if (inp->zipd)
        { uint32 dlen, tlen;
          int rez;

          dlen = zoffs[beg+1]-zoffs[beg];
          tlen = IO_BLOCK;
          read(fid,zuf,dlen);
          if ((rez = Gzip_Uncompress(buf,&tlen,zuf,dlen)) != Z_OK)
            { fprintf(stderr,"Decompression not OK!\n");
              exit (1);
            }
          slen = tlen;
#ifdef DEBUG_CHECK
          printf(" %d ->",dlen);
#endif
        }
      else
        slen = read(fid,buf,IO_BLOCK);
#ifdef DEBUG_CHECK
      printf(" %d\n",slen);
#endif

      olen = 0;
      for (b = 0; b < slen; b++)
        { c = buf[b];
#ifdef DEBUG_AUTO
          printf("  %.5s: %c\n",Name[state],c);
#endif
          switch (state)
          { case UK:
              if (c == '@' && alive)
                { if (GROUP)
                    { state = QHL;
                      count = 0;
                      lanl  = 0;
                      olen += 4 + 2*INT_MAXLEN;
                    }
                  else
                    { state = QHE;
                      nread += 1;
                    }
                  data->rpt = beg;
                  break;
                }
              alive = (c == '\n' && ! (nl_2 || pl_1));
              nl_2 = nl_1;
              nl_1 = (c == '\n');
              pl_1 = (c == '+');
#ifdef DEBUG_AUTO
              printf("    n2=%d n1=%d p1=%d a=%d\n",nl_2,nl_1,pl_1,alive);
#endif
              break;
            case BEG:
              if (c != '@')
                { data->etype  = 1;
                  data->lineno = lineno;
                  return (NULL);
                }
              if (GROUP)
                { state = QHL;
                  count = 0;
                  lanl  = 0;
                  olen += 4 + 2*INT_MAXLEN;
                }
              else
                { state = QHE;
                  nread += 1;
                }
              break;
            case QHL:
              if (c == '\n')
                { data->etype  = 2;
                  data->lineno = lineno;
                  return (NULL);
                }
              if (c == ':')
                count += 1;
              if (count == 3)
                { state = QHE;
                  lane[lanl] = '\0';
#ifdef DEBUG_AUTO
                  printf("  lane = %s\n",lane);
#endif
                  nread += 1;
                  if (fs >= 0)
                    { baux[fb] = lanl;
                      fs = -1;
                    }

                  if (ngroup < 0 || strcmp(lanegroup,lane) != 0)
                    { if (ngroup >= 0)
                        { gbp = totbp - gbp;
                          gnr = nread - gnr;
                          if (gbp > gmaxbp)
                            gmaxbp = gbp;
                          if (gnr > gmaxnr)
                            gmaxnr = gnr; 
                          if (ngroup == 0)
                            { data->first.nrs = gnr;
                              data->first.bps = gbp;
                            }

                          if (ngroup >= maxgct)
                            { maxgct = 1.2*ngroup + 10;
                              if (gcount == NULL)
                                gcount = (int64 *) Malloc(sizeof(int64)*(maxgct+1),
                                                          "Allocating group read counts");
                              else
                                gcount = (int64 *) Realloc(gcount-1,sizeof(int64)*(maxgct+1),
                                                           "Allocating group read counts");
                              gcount  += 1;
                            }
                          gcount[ngroup] = gnr;
       
                          glen = strlen(lanegroup);
                          gtotc += glen;
                          if (glen > gmaxc)
                            gmaxc = glen;
#ifdef DEBUG_AUTO
                          printf("   end of %d:%s: %lld %lld\n",ngroup,lanegroup,gbp,gnr);
#endif
                          free(lanegroup);
                        }
                      else
                        data->first.name = strdup(lane);
                      ngroup += 1;
                      lanegroup = strdup(lane);
                      gbp = totbp;
                      gnr = nread;
                    }
                }
              else
                { if (lanl >= lanmax)
                    { lanmax = lanmax*1.2 + 1000;
                      lane = (char *) Realloc(lane,lanmax+1,"Reallocating lane name");
                      if (lane == NULL)
                        exit (1);
                    }
                  lane[lanl++] = c;
                  olen += 1;
                }
              break;
            case QHE:
              if (c == '\n')
                { if (fastq)
                    state = SEQ;
                  else
                    state = FEQ;
                  blen  = 0;
                  lineno += 1;
                  if (genp)
                    olen += 5 + INT_MAXLEN;
                  else
                    olen += 3 + INT_MAXLEN;
                }
              break;
            case SEQ:
              olen += 1;
              if (c == '\n')
                { state = PLS;
                  lineno += 1;
                  if (fs >= 0)
                    { baux[fb] = blen;
                      fs = -1;
                    }
                }
              else if (NotDNA[c])
                { data->etype  = 3;
                  data->lineno = lineno;
                  return (NULL);
                }
              else
                blen += 1;
              break;
            case PLS:
              if (c != '+')
                { data->etype  = 4;
                  data->lineno = lineno;
                  return (NULL);
                }
              state = PEN;
              break;
            case PEN:
              if (c != '\n')
                { data->etype  = 4;
                  data->lineno = lineno;
                  return (NULL);
                }
              state = QVS;
              qlen  = 0;
              lineno += 1;
              if (QVS_OUT)
                olen += 3 + INT_MAXLEN;
              break;
            case QVS:
              if (QVS_OUT)
                olen += 1;
              if (c == '\n')
                { state = BEG;
                  if (qlen != blen)
                    { data->etype  = 5;
                      data->lineno = lineno;
                      return (NULL);
                    }
                  totbp += blen;
                  if (blen > maxlen)
                    maxlen = blen;
                  if (fs >= 0)
                    { baux[fb] = qlen;
                      fs = -1;
                    }
#ifdef DEBUG_AUTO
                  printf("Entry complete: %d (nr = %lld mx = %d bp = %lld)\n",
                         blen,nread,maxlen,totbp);
#endif
                  if (beg > end || (beg == end && b >= 2))
                    goto done;
                  else
                    { lineno += 1;
                      recno  += 1;
                    }
                }
              else if (NotQV[c])
                { data->etype  = 6;
                  data->lineno = lineno;
                  return (NULL);
                }
              else
                qlen += 1;
              break;

            case FK:
              if (c == '\n')
                state = FK1;
              break;
            case FK1:
              if (c == '>')
                { if (GROUP)
                    { state = QHL;
                      count = 0;
                      lanl  = 0;
                      olen += 4 + 2*INT_MAXLEN;
                    }
                  else
                    { state = QHE;
                      nread += 1;
                    }
                  data->rpt = beg;
                }
              else if (c != '\n')
                state = FK;
              break;
            case FEQ:
              if (c == '\n')
                state = FNL;
              else if (NotDNA[c])
                { data->etype  = 3;
                  data->lineno = lineno;
                  return (NULL);
		}
	      else
                { blen += 1;
                  olen += 1;
                }
              break;
            case FNL:
              lineno += 1;
              if (c == '>')
                { if (fs >= 0)
                    { baux[fb] = blen;
                      fs = -1;
                    }
                  recno += 1;
                  if (b == 0)
                    brecd[beg] = recno;
                  olen += 1;
                  totbp += blen;
                  if (blen > maxlen)
                    maxlen = blen;
                  if (beg > end || (beg == end && b >= 1))
                    goto done;
                  else
                    { if (GROUP)
                        { state = QHL;
                          count = 0;
                          lanl  = 0;
                          olen += 4 + 2*INT_MAXLEN;
                        }
                      else
                        { state = QHE;
                          nread += 1;
                        }
                    }
                }
              else if (c != '\n')
                { if (NotDNA[c])
                    { data->etype  = 3;
                      data->lineno = lineno;
                      return (NULL);
                    }
                  blen += 1;
                  olen += 1;
                  state = FEQ;
                }
              break;
          }
        }

      if (olen > outmax)
        outmax = olen;
      if (outfst < 0 && state >= BEG)
        outfst = olen;

      beg += 1;
    }

done:
  nread += 1;
  if (GROUP)
    { gbp = totbp - gbp;
      gnr = nread - gnr;
      if (gbp > gmaxbp)
        gmaxbp = gbp;
      if (gnr > gmaxnr)
        gmaxnr = gnr; 
      if (ngroup == 0)
        { data->first.nrs = gnr;
          data->first.bps = gbp;
        }
      data->last.name = lanegroup;
      data->last.nrs = gnr;
      data->last.bps = gbp;

      if (ngroup >= maxgct)
        { maxgct = 1.2*ngroup + 10;
          if (gcount == NULL)
            gcount = (int64 *) Malloc(sizeof(int64)*(maxgct+1),
                                      "Allocating group read counts");
          else
            gcount = (int64 *) Realloc(gcount-1,sizeof(int64)*(maxgct+1),
                                       "Allocating group read counts");
          gcount  += 1;
        }
      gcount[ngroup++] = gnr;

      glen = strlen(lanegroup);
      gtotc += glen;
      if (glen > gmaxc)
        gmaxc = glen;
    }

  if (olen > outmax)
    outmax = olen;
  if (outfst < 0)
    outfst = olen;

  data->gcount = gcount;
  data->etype  = 0;        //  Non-zero => error return occurred
  data->lineno = lineno;
  data->outmax = outmax;
  data->outfst = outfst;
  data->outlst = olen;

  info->nread  = nread;
  info->totbp  = totbp;
  info->ngroup = ngroup;
  info->gtotc  = gtotc;
  info->gmaxnr = gmaxnr;
  info->maxlen = maxlen;
  info->gmaxc  = gmaxc;
  info->gmaxbp = gmaxbp;

  close(fid);

  return (NULL);
}

  //  Threaded fasta/q file checker

static void Check_FastQ(Fast_Input *input, Header_Info *info, int isforward, int ispair)
{ pthread_t threads[NTHREADS];
  Check_Arg parm[NTHREADS];
  int64     beg;
  int64     lineno;
  int       nthreads, genp;
  int       sgroup, stotc;
  int       n;

  nthreads = NTHREADS;
  if (input->zsize < NTHREADS)
    NTHREADS = input->zsize;

  genp  = (isforward && ispair);

  //  Allocate buffers

  if (input->zipd)
    { input->ibufs = (uint8 *) Malloc(2*IO_BLOCK*NTHREADS,"Allocating IO buffers");
      if (input->ibufs == NULL)
        exit (1);
      input->zbufs = input->ibufs + IO_BLOCK*NTHREADS;
    }
  else
    { input->ibufs = (uint8 *) Malloc(IO_BLOCK*NTHREADS,"Allocating IO buffers");
      if (input->ibufs == NULL)
        exit (1);
    }

  //  Setup input args to each thread including block range [beg,end)
  //    Open a separate file descriptor for each thread !

  beg = 0;
  for (n = 0; n < NTHREADS; n++)
    { parm[n].beg = beg;
      beg = (input->zsize*(n+1))/NTHREADS;
      parm[n].end = beg;
      parm[n].buf = input->ibufs + n*IO_BLOCK;
      if (input->zipd)
        parm[n].zuf = input->zbufs + n*IO_BLOCK;
      else
        parm[n].zuf = NULL;
      parm[n].inp   = input;
      parm[n].genp  = genp;
    }

  for (n = 0; n < NTHREADS; n++)
#ifdef DEBUG_CHECK
    check_thread(parm+n);
#else
    pthread_create(threads+n,NULL,check_thread,parm+n);
#endif

  //  Join threads in order and report the first error if there is one.

  lineno = 0;
  for (n = 0; n < NTHREADS; n++)
    { pthread_join(threads[n],NULL);

      lineno += parm[n].lineno;
      switch (parm[n].etype)
      { case 0:
          break;
        case 1:
          fprintf(stdout,"%s: %s.%lld: Entry header does not start with an @-sign\n",
                         Prog_Name,input->root,lineno);
          exit (1);
        case 2:
          fprintf(stdout,"%s: %s.%lld: Entry header does not appear to be",
                         Prog_Name,input->root,lineno);
          fprintf(stdout," an Illumina header\n");
          exit (1);
        case 3:
          fprintf(stdout,"%s: %s.%lld: Non-DNA symbol in sequence\n",Prog_Name,input->root,lineno);
          exit (1);
        case 4:
          fprintf(stdout,"%s: %s.%lld Divider line does not start with a +-sign\n",
                         Prog_Name,input->root,lineno);
          exit (1);
        case 5:
          fprintf(stdout,"%s: %s.%lld: QV line doesn't have the same length as sequence line\n",
                         Prog_Name,input->root,lineno);
          exit (1);
        case 6:
          fprintf(stdout,"%s: %s.%lld: Non-QV symbol in qv-string\n",Prog_Name,input->root,lineno);
          exit (1);
        case 7:
          fprintf(stdout,"%s: %s: Incomplete last entry\n",Prog_Name,input->root);
          exit (1);
      }
    }

  //  Merge groups across division boundaries as necessary

  if (GROUP)
    { char  *xlast;
      int64 *xgcnt;
      int    xngrp;
      int    x;

#ifdef DEBUG_CHECK
      printf("\nThread: First G :: Last G\n");
      for (n = 0; n < NTHREADS; n++)
        printf("  %d: %s %lld %lld :: %s %lld %lld\n",n,
               parm[n].first.name,parm[n].first.nrs,parm[n].first.bps,
               parm[n].last.name,parm[n].last.nrs,parm[n].last.bps);
#endif

      parm[0].gcount[-1] = 1;
      input->zgrp[0] = -1;
      input->zgct[0] = 0;

      x = 0;
      xlast = parm[x].last.name;
      xgcnt = parm[x].gcount;
      xngrp = parm[x].info.ngroup;

      sgroup = 0;
      stotc  = 0;

      for (n = 1; n < NTHREADS; n++)
        { if (strcmp(xlast,parm[n].first.name) != 0)
            { parm[n].gcount[-1] = 1;
              input->zgrp[parm[n].rpt] = -1;
              input->zgct[parm[n].rpt] = 0;
            }
          else
            { xgcnt[xngrp-1] = xgcnt[xngrp-1] + parm[n].gcount[0];
              input->zgrp[parm[n].rpt] = 0;
	      input->zgct[parm[n].rpt] = -1;
              sgroup += 1;
              stotc  += strlen(parm[n].first.name);
              parm[x].last.nrs += parm[n].first.nrs;
              parm[x].last.bps += parm[n].first.bps;
              if (parm[x].last.nrs > parm[x].info.gmaxnr)
                parm[x].info.gmaxnr = parm[x].last.nrs;
              if (parm[x].last.bps > parm[x].info.gmaxbp)
                parm[x].info.gmaxbp = parm[x].last.bps;
            }
    
          if (strcmp(xlast,parm[n].last.name) != 0)
            { x = n;
              xlast = parm[x].last.name;
              xgcnt = parm[x].gcount;
              xngrp = parm[x].info.ngroup;
            }
        }

#ifdef DEBUG_CHECK
      printf("\nGROUP COUNTS\n");
      for (n = 0; n < NTHREADS; n++)
        { xngrp = parm[n].info.ngroup;
          printf("  Parition %d: rpt = %d\n",n,parm[n].rpt);
          for (x = 0; x < xngrp; x++)
            printf("     %d: %lld\n",x,parm[n].gcount[x]);
        }

      x = 0;
      printf("\nBLOCK TABLE:\n");
      for (n = 0; n < input->zsize; n++)
        { while (x < NTHREADS-1 && n >= parm[x+1].rpt)
            x += 1;
          printf(" %3d: %lld %5s l=%5d g=%d,g%d[%d] = %lld\n",
                 n,input->zrecd[n],Name[input->zstat[n]],input->zaux[n],input->zgct[n],
                   x,input->zgrp[n],parm[x].gcount[input->zgrp[n]]);
        }
      fflush(stdout);
#endif
    }

  for (n = 0; n < NTHREADS; n++)
    { if (parm[n].outmax > MAXOUT)
        MAXOUT = parm[n].outmax;
      input->gcnt[n] = parm[n].gcount;
      input->ngrp[n] = parm[n].info.ngroup;
      input->rpt[n]  = parm[n].rpt;
    }
  input->rpt[NTHREADS] = input->zsize;
  input->zrecd[input->zsize] = 0;

  for (n = 1; n < NTHREADS; n++)
    { if (parm[n-1].outlst + parm[n].outfst > MAXOUT)
        MAXOUT = parm[n-1].outlst + parm[n].outfst;
    }

  //  Accumulate statistics for entire file

  info->nread  = 0;
  info->totbp  = 0;
  info->ngroup = 0;
  info->gtotc  = 0;
  info->gmaxnr = 0;
  info->maxlen = 0;
  info->gmaxc  = 0;
  info->gmaxbp = 0;
  for (n = 0; n < NTHREADS; n++)
    { info->nread  += parm[n].info.nread;
      info->totbp  += parm[n].info.totbp;
      info->ngroup += parm[n].info.ngroup;
      info->gtotc  += parm[n].info.gtotc;
      if (info->gmaxnr < parm[n].info.gmaxnr)
        info->gmaxnr = parm[n].info.gmaxnr;
      if (info->maxlen < parm[n].info.maxlen)
        info->maxlen = parm[n].info.maxlen;
      if (info->gmaxc < parm[n].info.gmaxc)
        info->gmaxc  = parm[n].info.gmaxc;
      if (info->gmaxbp < parm[n].info.gmaxbp)
        info->gmaxbp = parm[n].info.gmaxbp;
    }
  info->ngroup -= sgroup;
  info->gtotc  -= stotc;

  NTHREADS = nthreads;
}


/*******************************************************************************************
 *
 *  Second pass thread routines prepare .seq formated output in buffers for each input block.
 *  Note carefully that there are separate blocks for the forward and reverse files if
 *  processing a pair, the reads in the blocks need to be paired.
 *
 ********************************************************************************************/

  //  Output the contents of a single block using all the information collected during
  //    the first pass to make this possible.  Output just the characters generated from
  //    the data in the block.

typedef struct
  { int64       blk;    //  Process block blk
    int         fid;    //  Input file descriptor (independent for each thread yet same file)
    uint8      *buf;    //  input block buffer
    uint8      *zuf;    //  decode buffer
    char       *out;    //  output buffer
    int         olen;   //  amount of data in output buffer
    char       *seq;    //  string buffer
    Fast_Input *inp;    //  file information
    int         dogrp;  //  build group data structures
    int         genp;   //  account for P line in output estimates
    int64      *gcount; //  gcount array to use for this block
  } Output_Arg;

static void *output_thread(void *arg)
{ Output_Arg  *data   = (Output_Arg *) arg;
  int64        blk    = data->blk;
  int          fid    = data->fid;
  uint8       *buf    = data->buf;
  uint8       *zuf    = data->zuf;
  char        *out    = data->out;
  Fast_Input  *inp    = data->inp;
  char        *seq    = data->seq;
  int64       *gcount = data->gcount;
  int          dogrp  = data->dogrp;
  int          genp   = data->genp;

  int          fastq = inp->fastq;
  int64       *zoffs = inp->zoffs;
  int         *bstat = inp->zstat;
  int         *baux  = inp->zaux;
  int         *bgrp  = inp->zgrp;
  int         *bgct  = inp->zgct;

  int   b, c, slen;
  int   state, count, grpn, grpc;
  int   blen, extn;
  int   olen;

  //  Read the block

  if (inp->zipd)
    { uint32 dlen, tlen;
      int    rez;

      lseek(fid,zoffs[blk],SEEK_SET);
      dlen = zoffs[blk+1]-zoffs[blk];
#ifdef DEBUG_OUT
      printf("Doing block %lld (@ %lld) %d ->",blk,lseek(fid,0,SEEK_CUR),dlen);
#endif
      tlen = IO_BLOCK;
      read(fid,zuf,dlen);
      if ((rez = Gzip_Uncompress(buf,&tlen,zuf,dlen)) != Z_OK)
        { fprintf(stderr,"Decompression not OK!\n");
          exit (1);
        }
      slen = tlen;
    }
  else
    { lseek(fid,blk*IO_BLOCK,SEEK_SET);
#ifdef DEBUG_OUT
      printf("Doing block %lld (@ %lld)",blk,lseek(fid,0,SEEK_CUR));
#endif
      slen = read(fid,buf,IO_BLOCK);
    }
#ifdef DEBUG_OUT
  printf(" %d\n",slen);
#endif

  //  Set up boundary conditions for start saved in the previous pass for each block

  olen  = 0;
  state = bstat[blk];
  if (state >= QHL)
    { count = state-QHL;
      state = QHL;
    }
  grpn  = bgrp[blk];
  grpc  = bgct[blk];
  if (state >= SEQ)
    { blen = 0; 
      extn = 1;
    }
  else
    extn = 0;

  //  Proceed with scan, buffering S, Q, and G strings

  for (b = 0; b < slen; b++)
    { c = buf[b];
#ifdef DEBUG_AUTO
      printf("  %.5s: %c\n",Name[state],c);
#endif
      switch (state)
      { case BEG:
          if (dogrp)
            { state = QHL;
              count = 0;
              blen  = 0;
            }
          else
            state = QHE;
          break;
        case QHL:
          if (c == ':')
            count += 1;
          if (count == 3)
            { seq[blen] = '\0';
              grpc += 1;
              if (grpc == gcount[grpn])
                { grpn += 1;
                  grpc = 0;
                  if (extn)
                    olen += sprintf(out+olen,"%s\n",seq);
                  else
                    olen += sprintf(out+olen,"g %lld %d %s\n",gcount[grpn],blen,seq);
                }
              extn  = 0;
              state = QHE;
            }
          else
            seq[blen++] = c;
          break;
        case QHE:
          if (c == '\n')
            { if (fastq)
                state = SEQ;
              else
                state = FEQ;
              blen = 0;
            }
          break;
        case SEQ:
          if (c == '\n')
            { seq[blen] = '\0';
              if (extn)
                olen += sprintf(out+olen,"%s",seq);
              else if (genp)
                olen += sprintf(out+olen,"P\nS %d %s",blen,seq);
              else
                olen += sprintf(out+olen,"S %d %s",blen,seq);
              if (QVS_OUT)
                out[olen++] = '\n';
              extn  = 0;
              state = PLS;
            }
          else
            seq[blen++] = c;
          break;
        case PLS:
          state = PEN;
          break;
        case PEN:
          state = QVS;
          blen  = 0;
          break;
        case QVS:
          if (c == '\n')
            { if (QVS_OUT)
                { seq[blen] = '\0';
                  if (extn)
                    olen += sprintf(out+olen,"%s",seq);
                  else
                    olen += sprintf(out+olen,"Q %d %s",blen,seq);
                }
              out[olen++] = EOR;
              extn  = 0;
              state = BEG;
            }
          else if (QVS_OUT)
            seq[blen++] = c;
          break;

        case FEQ:
          if (c == '\n')
            state = FNL;
	  else
            seq[blen++] = c;
          break;
        case FNL:
          if (c == '>')
            { seq[blen] = '\0';
              if (extn)
                olen += sprintf(out+olen,"%s%c",seq,EOR);
              else if (genp)
                olen += sprintf(out+olen,"P\nS %d %s%c",blen,seq,EOR);
              else
                olen += sprintf(out+olen,"S %d %s%c",blen,seq,EOR);
              extn = 0;
              if (dogrp)
                { state = QHL;
                  count = 0;
                  blen  = 0;
                }
              else
                state = QHE;
            }
          else if (c != '\n')
            { seq[blen++] = c;
              state = FEQ;
            }
          break;
      }
    }

  //  Flush partial S, Q, or G string if necessary

  if (state >= SEQ)
    { seq[blen] = '\0';
      if (state == QVS)
        { if (QVS_OUT)
            { if (extn)
                olen += sprintf(out+olen,"%s",seq);
              else
                olen += sprintf(out+olen,"Q %d %s",baux[blk+1],seq);
            }
        }
      else if (state == QHL)
        { grpc += 1;
          if (grpc == gcount[grpn])
            { grpn += 1;
              if (extn)
                olen += sprintf(out+olen,"%s",seq);
              else
                olen += sprintf(out+olen,"g %lld %d %s",gcount[grpn],baux[blk+1],seq);
            }
        }
      else
        { seq[blen] = '\0';
          if (extn)
            olen += sprintf(out+olen,"%s",seq);
          else if (genp)
            olen += sprintf(out+olen,"P\nS %d %s",baux[blk+1],seq);
          else
            olen += sprintf(out+olen,"S %d %s",baux[blk+1],seq);
        }
    }

  out[olen] = EOR;
  data->olen = olen;
  return (NULL);
}


/*******************************************************************************************
 *
 *  If there is just a single file, i.e. no read pairing is required, then simply produce
 *    the desired output for the next NTHREAD blocks at a time, and then sequentially output.
 *
 ********************************************************************************************/

static void Output_Sequence(Fast_Input *input, Header_Info *info)
{ Output_Arg parm[NTHREADS];
#ifndef DEBUG_OUT
  pthread_t  thread[NTHREADS];
#endif

  int64   zsize;
  char   *out;
  char   *seq;
  int     sln;
  int     b, c, n;
  int     ng;
  int     idout;

  zsize = input->zsize;
  idout = fileno_unlocked(stdout);   //  Is the unlock helping?

  if (info->gmaxc > info->maxlen)
    sln = info->gmaxc+1;
  else
    sln = info->maxlen+1;
  if (sln > IO_BLOCK)
    sln = IO_BLOCK+1;

  seq  = (char *) Malloc(NTHREADS*(sln+MAXOUT),"Allocating output buffers");
  if (seq == NULL)
    exit (1);
  out = seq + NTHREADS*sln;

  for (n = 0; n < NTHREADS; n++)
    { parm[n].buf   = input->ibufs + n*IO_BLOCK;
      parm[n].zuf   = input->zbufs + n*IO_BLOCK;
      parm[n].seq   = seq + n*sln;
      parm[n].out   = out + n*MAXOUT;
      parm[n].fid   = open(input->path,O_RDWR);
      parm[n].inp   = input;
      parm[n].dogrp = GROUP;
      parm[n].genp  = 0;
    }

  //  Proceed with comnpression and output, accumulating index table

  ng = 0;
  for (b = 0; b < zsize; b += NTHREADS)
    { for (c = b, n = 0; n < NTHREADS; n++, c++)
        if (c < zsize)
          { parm[n].blk = c;
            while (c >= input->rpt[ng+1])
              ng += 1;
            parm[n].gcount = input->gcnt[ng];
#ifdef DEBUG_OUT
            output_thread(parm+n);
            fwrite(parm[n].out,1,parm[n].olen,stdout);
#else
            pthread_create(thread+n,NULL,output_thread,parm+n);
#endif
          }

#ifndef DEBUG_OUT
      for (c = b, n = 0; n < NTHREADS; n++, c++)
        if (c < zsize)
          { pthread_join(thread[n],NULL);
            write(idout,parm[n].out,parm[n].olen);
          }
#endif
    }

  for (n = 0; n < NTHREADS; n++)
    close(parm[n].fid);
  free(seq);
}


/*******************************************************************************************
 *
 *  If a pair of files, then the strategy is to produce the output for the next NTHREADS forward
 *    and reverse blocks based on record order (uses zrecd array) and then in an unthreaded sweep
 *    paired read records are collected in a buffer and streamed to stdout when full.
 *    This is not straight forward when the forward and reverse sequences are not the same length.
 *
 ********************************************************************************************/

static void Output_Pair(Fast_Input *forward, Fast_Input *reverse, Header_Info *info)
{ Output_Arg parm[NTHREADS];
#ifndef DEBUG_OUT
  pthread_t  thread[NTHREADS];
#endif

  int     F_fid[NTHREADS], R_fid[NTHREADS];
  int     idout;
  char   *seq;
  char   *out;
  char   *obuf;

  idout  = fileno_unlocked(stdout);   //  Is the unlock helping?

  EOR = '\01';

  //  Setup buffers and open enough IO units for up to NTHREADS for each input

  { int sln;
    int n;

    if (info->gmaxc > info->maxlen)
      sln = info->gmaxc+1;
    else
      sln = info->maxlen+1;
    if (sln > IO_BLOCK)
      sln = IO_BLOCK+1;

    seq  = (char *) Malloc(NTHREADS*sln,"Allocating output buffers");
    if (seq == NULL)
      exit (1);

    out  = (char *) Malloc((NTHREADS+2)*MAXOUT,"Allocating output buffers");
    if (out == NULL)
      exit (1);

    for (n = 0; n < NTHREADS; n++)
      { parm[n].buf  = forward->ibufs + n*IO_BLOCK;
        parm[n].zuf  = forward->zbufs + n*IO_BLOCK;
        parm[n].seq  = seq + n*sln;
        F_fid[n] = open(forward->path,O_RDWR);
        R_fid[n] = open(reverse->path,O_RDWR);
      }

    obuf = out + (NTHREADS+1)*MAXOUT;
  }

  //  Proceed with converstion to VGP format for each file in properly sequenced blocks
  //    and then interleave the VGP entries to standard out

  { int     nplus;
    int     olen;

    int     F_cblk,  R_cblk;
    int64   F_size,  R_size;
    int64  *F_recd, *R_recd;
    int64   F_off,   R_off;
    int64   F_seg,   R_seg;
    int64  *F_gcnt, *R_gcnt;
    int     F_nblk , R_nblk;
    int     F_thr,   R_thr;
    int     F_rpt,   R_rpt;
    int     F_map[NTHREADS+1], R_map[NTHREADS+1];
    char   *F_loc,  *R_loc;
    char   *F_end,  *R_end;

    F_size = forward->zsize;
    R_size = reverse->zsize;
    F_recd = forward->zrecd;
    R_recd = reverse->zrecd;

    nplus = NTHREADS + 1;
    olen  = 0;

    F_cblk = F_nblk = 0;
    R_cblk = R_nblk = 0;
    R_off  = F_off = 0;
    R_thr  = F_thr = 0;
    R_rpt  = reverse->rpt[1];
    F_rpt  = forward->rpt[1];
    R_seg  = R_recd[R_rpt];
    F_seg  = F_recd[F_rpt];
    R_gcnt = reverse->gcnt[0];
    F_gcnt = forward->gcnt[0];

    while (R_nblk < R_size)   //  While blocks still remaining
      { int ebuf;
        int uthr;
        int forw;
        int n, b;

        //  Determine next NTHREADS blocks of forward & reverse file blocks to convert next

        if (F_cblk < F_nblk)
          ebuf  = F_map[F_cblk % nplus];
        else if (R_cblk < R_nblk)
          ebuf  = R_map[R_cblk % nplus];
        else
          ebuf = NTHREADS;
        b = 0;
        for (n = 0; n < NTHREADS; n++)
          { if (b == ebuf)
              b += 1;
            if (R_nblk >= R_size)
              break;
            parm[n].out = out + b*MAXOUT;
            if (F_nblk >= F_size || R_recd[R_nblk] < F_recd[F_nblk])
              { // Generate rnblk
                parm[n].inp    = reverse;
                parm[n].blk    = R_nblk;
                parm[n].dogrp  = 0;
                parm[n].genp   = 0;
                parm[n].gcount = R_gcnt;
                parm[n].fid    = R_fid[n];
  
                R_map[R_nblk % nplus] = b;
                R_nblk += 1;
                if (R_nblk > R_rpt)
                  { R_thr  += 1;
                    R_rpt  = reverse->rpt[R_thr+1];
                    R_off += R_seg;
                    R_seg  = R_recd[R_rpt];
                    R_gcnt = reverse->gcnt[R_thr];
                  }
                R_recd[R_nblk] += R_off;
              }
            else
              { // Generate fnblk
                parm[n].inp    = forward;
                parm[n].blk    = F_nblk;
                parm[n].dogrp  = GROUP;
                parm[n].genp   = 1;
                parm[n].gcount = F_gcnt;
                parm[n].fid    = F_fid[n];
  
                F_map[F_nblk % nplus] = b;
                F_nblk += 1;
                if (F_nblk > F_rpt)
                  { F_thr += 1;
                    F_rpt  = forward->rpt[F_thr+1];
                    F_off += F_seg;
                    F_seg  = F_recd[F_rpt];
                    F_gcnt = forward->gcnt[F_thr];
                  }
                F_recd[F_nblk] += F_off;
              }
            b += 1;
          }

        //  Compute the VGP format for the next NTHREADS blocks

        uthr = n;
        for (n = 0; n < uthr; n++)
#ifdef DEBUG_OUT
          output_thread(parm+n);
#else
          pthread_create(thread+n,NULL,output_thread,parm+n);
        for (n = 0; n < uthr; n++)
          pthread_join(thread[n],NULL);
#endif

        //  Have blocks F_cblk .. F_nblk-1 and R_cblk .. R_nblk-1
        //    in (F/R)_map assigned output blocks.

#ifdef DEBUG_OUT
        printf("\nNext Batch:\n");
        printf("   Forward: %d\n",ebuf);
        for (b = F_cblk; b < F_nblk; b++)
          { printf("   %5d: %d  %lld",b,F_map[b%nplus],F_recd[b]);
            if (F_map[b%nplus] == ebuf)
              printf(" #");
            printf("\n");
          }
        printf("   Reverse:\n");
        for (b = R_cblk; b < R_nblk; b++)
          { printf("   %5d: %d  %lld",b,R_map[b%nplus],R_recd[b]);
            if (R_map[b%nplus] == ebuf)
              printf(" #");
            printf("\n");
          }
#endif

        //  Interleave VGP sequence records between the forward and reverse outputs
        //    for the blocks.  The F/R_loc/end variables are correctly set to pick up
        //    from where one last left off in the buffered block ebuf.  The interleaving
        //    stops when one or the other of the two traversals runs out of data.

        forw = 1;
        if (F_cblk < F_nblk)
          { b = F_map[F_cblk%nplus];
            if (b != ebuf)
              { F_loc = out + b * MAXOUT;
                if (b < ebuf)
                  F_end = F_loc + parm[b].olen;
                else
                  F_end = F_loc + parm[b-1].olen;
              }
            else
              forw = 0;
          }

        if (R_cblk < R_nblk)
          { b = R_map[R_cblk%nplus];
            if (b != ebuf)
              { R_loc = out + b * MAXOUT;
                if (b < ebuf)
                  R_end = R_loc + parm[b].olen;
                else
                  R_end = R_loc + parm[b-1].olen;
              }
            else
              forw = 1;
          }

        while (1)
          { char *str;
            int   elen;

            if (forw)
              { if (F_cblk >= F_nblk)
                  break;
  
                str = F_loc;
                while (*F_loc != EOR)
                  F_loc += 1;
  
                if (F_loc >= F_end)
                  { elen = F_loc-str;
#ifdef DEBUG_LINES
                    printf("F B %d: '%.*s'\n",F_cblk,elen,str);
#endif
#ifndef DEBUG_OUT
                    if (olen + elen > MAXOUT)
                      { write(idout,obuf,olen);
                        olen = 0;
                      }
                    memcpy(obuf+olen,str,elen);
                    olen += elen;
#endif
                    F_cblk += 1;
                    if (F_cblk >= F_nblk)
                      break;
  
                    b = F_map[F_cblk%nplus];
                    F_loc = out + b * MAXOUT;
                    if (b < ebuf)
                      F_end = F_loc + parm[b].olen;
                    else
                      F_end = F_loc + parm[b-1].olen;
                  }
                else
                  { *F_loc++ = '\n';
                    elen = F_loc-str;
#ifdef DEBUG_LINES
                    printf("F E %d: '%.*s",F_cblk,elen,str);
#endif
#ifndef DEBUG_OUT
                    if (olen + elen > MAXOUT)
                      { write(idout,obuf,olen);
                        olen = 0;
                      }
                    memcpy(obuf+olen,str,elen);
                    olen += elen;
#endif
                    forw = 0;
                  }
              }
            else
              { if (R_cblk >= R_nblk)
                  break;
  
                str = R_loc;
                while (*R_loc != EOR)
                  R_loc += 1;
  
                if (R_loc >= R_end)
                  { elen = R_loc-str;
#ifdef DEBUG_LINES
                    printf("R B %d: '%.*s'\n",R_cblk,elen,str);
#endif
#ifndef DEBUG_OUT
                    if (olen + elen > MAXOUT)
                      { write(idout,obuf,olen);
                        olen = 0;
                      }
                    memcpy(obuf+olen,str,elen);
                    olen += elen;
#endif
                    R_cblk += 1;
                    if (R_cblk >= R_nblk)
                      break;
  
                    b = R_map[R_cblk%nplus];
                    R_loc = out + b * MAXOUT;
                    if (b < ebuf)
                      R_end = R_loc + parm[b].olen;
                    else
                      R_end = R_loc + parm[b-1].olen;
                  }
                else
                  { *R_loc++ = '\n';
                    elen = R_loc-str;
#ifdef DEBUG_LINES
                    printf("R E %d: '%.*s",R_cblk,elen,str);
#endif
#ifndef DEBUG_OUT
                    if (olen + elen > MAXOUT)
                      { write(idout,obuf,olen);
                        olen = 0;
                      }
                    memcpy(obuf+olen,str,elen);
                    olen += elen;
#endif
                    forw = 1;
                  }
              }
          }
      }

    if (olen > 0)
      write(idout,obuf,olen);
  }

  { int n;

    for (n = 0; n < NTHREADS; n++)
      { close(F_fid[n]);
        close(R_fid[n]);
      }
    free(out);
    free(seq);
  }
}


/*******************************************************************************************
 *
 *  Top level.  Process arguments, open files, in a first pass check each file and accumulate
 *    statistics for header, then output header, and in second pass produce .irp output.
 *
 ********************************************************************************************/

int main(int argc, char *argv[])
{ Fast_Input  input1,  input2;
  Header_Info stats1,  stats2;
  int         ispair;

  //  Parse command line options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("VGPseq")

    NTHREADS = 4;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vsg")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE  = flags['v'];
    QVS_OUT  = ! flags['s'];
    GROUP    = flags['g'];

    if (argc != 2 && argc != 3)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, show progress as proceed.\n");
        fprintf(stderr,"      -s: Output sequences only, skip QVs\n");
        fprintf(stderr,"      -g: Output group lanes based on Illumina lanes\n");
        fprintf(stderr,"      -T: Number of threads to use\n");
        exit (1);
      }

    ispair = (argc == 3);
  }

  //  Open input files

  Open_Fastq(argv[1],&input1);

  if (ispair)
    { Open_Fastq(argv[2],&input2);
      if (input1.fastq != input2.fastq)
        { fprintf(stderr,"%s: Pair is not both .fastq or both .fasta\n",Prog_Name);
          exit (1);
        }
      if (input1.zipd != input2.zipd)
        { fprintf(stderr,"%s: Pair should both be zip'd or not\n",Prog_Name);
          exit (1);
        }
    }

  //  Scan 1: check syntax and count

  if (VERBOSE)
    { if (argc == 2)
        fprintf(stderr,"  Checking syntax of file %s\n",argv[1]);
      else
        fprintf(stderr,"  Checking syntax of forward file %s\n",argv[1]);
      fflush(stderr);
    }

  if ( ! input1.fastq)
    QVS_OUT = 0;

  Check_FastQ(&input1,&stats1,1,ispair);

  if (ispair)
    { if (VERBOSE)
        { fprintf(stderr,"  Checking syntax of reverse file %s\n",argv[2]);
          fflush(stderr);
        }

      Check_FastQ(&input2,&stats2,0,ispair);

      if (stats1.nread != stats2.nread)
        { fprintf(stdout,"%s: Number of reads in two files not equal: %lld vs %lld\n",
                         Prog_Name,stats1.nread,stats2.nread);
          exit (1);
        }

      stats1.totbp += stats2.totbp;
      if (stats1.maxlen < stats2.maxlen)
        stats1.maxlen = stats2.maxlen;
      if (stats1.gmaxnr < stats2.gmaxnr)
        stats1.gmaxnr = stats2.gmaxnr;
      stats1.gmaxbp += stats2.gmaxbp;
    }

  //  Output .irp header

  { int    i, clen, optl, len, idout;
    char   buffer[5000];
    char   date[20];
    time_t seconds;

    idout = fileno_unlocked(stdout);   //  Is the unlock helping?

    len = sprintf(buffer,"1 3 seq 1 0\n");
    write(idout,buffer,len);
    if (ispair)
      { len = sprintf(buffer,"2 3 irp\n");
        write(idout,buffer,len);
      }

    optl = VERBOSE + !QVS_OUT + GROUP;
    if (optl == 0)
      clen = -1;
    else
      clen = optl+1;
    for (i = 1; i < argc; i++)
      clen += strlen(argv[i])+1;
    seconds = time(NULL);
    strftime(date,20,"%F_%T",localtime(&seconds));

    len = sprintf(buffer,"! 6 VGPseq 3 1.0 %d",clen);
    if (optl)
      { len += sprintf(buffer+len," -");
        if (VERBOSE)
          len += sprintf(buffer+len,"v");
        if (!QVS_OUT)
          len += sprintf(buffer+len,"s");
        if (GROUP)
          len += sprintf(buffer+len,"g");
      }
    for (i = 1; i < argc; i++)
      len += sprintf(buffer+len," %s",argv[i]);
    len += sprintf(buffer+len," 19 %s\n",date);
    write(idout,buffer,len);

    if (GROUP)
      { len = sprintf(buffer,"# g %d\n",stats1.ngroup);
        write(idout,buffer,len);
      }
    if (ispair)
      { len = sprintf(buffer,"# P %lld\n",stats1.nread);
        write(idout,buffer,len);
      }
    len = sprintf(buffer,"# S %lld\n",(1+ispair)*stats1.nread);
    write(idout,buffer,len);
    if (QVS_OUT)
      { len = sprintf(buffer,"# Q %lld\n",(1+ispair)*stats1.nread);
        write(idout,buffer,len);
      }

    if (GROUP)
      { len = sprintf(buffer,"+ g %lld\n",stats1.gtotc);
        write(idout,buffer,len);
      }
    len = sprintf(buffer,"+ S %lld\n",stats1.totbp);
    write(idout,buffer,len);
    if (QVS_OUT)
      { len = sprintf(buffer,"+ Q %lld\n",stats1.totbp);
        write(idout,buffer,len);
      }

    if (GROUP)
      { len = sprintf(buffer,"@ g %lld\n",stats1.gmaxc);
        write(idout,buffer,len);
      }
    len = sprintf(buffer,"@ S %d\n",stats1.maxlen);
    write(idout,buffer,len);
    if (QVS_OUT)
      { len = sprintf(buffer,"@ Q %d\n",stats1.maxlen);
        write(idout,buffer,len);
      }

    if (GROUP)
      { if (ispair)
          { len = sprintf(buffer,"%% g # P %lld\n",stats1.gmaxnr);
            write(idout,buffer,len);
          }
        len = sprintf(buffer,"%% g # S %lld\n",(1+ispair)*stats1.gmaxnr);
        write(idout,buffer,len);
        if (QVS_OUT)
          { len = sprintf(buffer,"%% g # Q %lld\n",(1+ispair)*stats1.gmaxnr);
            write(idout,buffer,len);
          }
        len = sprintf(buffer,"%% g + S %lld\n",stats1.gmaxbp);
        write(idout,buffer,len);
        if (QVS_OUT)
          { len = sprintf(buffer,"%% g + Q %lld\n",stats1.gmaxbp);
            write(idout,buffer,len);
          }
      }
  }
  
  //  Scan 2: output the data in .irp format

  if (VERBOSE)
    { fprintf(stderr,"  Outputting read pairs in .irp format to stdout\n");
      fflush(stderr);
    }

  if (ispair)
    Output_Pair(&input1,&input2,&stats1);
  else
    Output_Sequence(&input1,&stats1);

  //  Tidy up just for good form

  Close_Fastq(&input1);
  if (ispair)
    Close_Fastq(&input2);
  exit (0);
}
