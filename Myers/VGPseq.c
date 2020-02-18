/*  Last edited: Feb  4 11:52 2020 (rd109) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include <fcntl.h>
#include <sys/stat.h>

#include "LIBDEFLATE/libdeflate.h"
#include "gene_core.h"
#include "../Durbin/VGPlib.h"

#undef   DEBUG_FIND
#undef    DEBUG_OUT
#undef    DEBUG_AUTO

static char *Usage = "[-vsg] [-T<int(4)>] <data:fast[aq][.gz]> ...";

#define IO_BLOCK 10000000ll

typedef struct libdeflate_decompressor DEPRESS;

#define INT_MAXLEN 10

int    VERBOSE;  //  Verbose mode?
int    QVS_OUT;  //  Include QV strings
int    GROUP;    //  Group lanes of reads
int    NTHREADS; //  # of threads to use


/*******************************************************************************************
 *
 *  Routines to open and close a possibly VGPzip'd fasta/q file
 *
 ********************************************************************************************/

  //  File object information

typedef struct
  { char   *path;   //  Full path name
    char   *root;   //  Ascii root name
    int64   fsize;  //  Size of (uncompressed) file in bytes
    int     fastq;  //  Is file a fastq? (else fasta)
    int     zipd;   //  Is file VGPzip'd (has paired .vzi file)
    int     recon;  //  Is file an uncompressed regular zip-file?
    int64   zsize;  //  Size of zip index (in blocks)
    int64  *zoffs;  //  zoffs[i] = offset to compressed block i
  } File_Object;


  //  Open and get info about a possibly VGPzip'd fasta/q file

static void Fetch_Fastq(char *arg, File_Object *input)
{ static char *suffix[5] = { ".fastq.gz", ".fasta.gz", ".fastq", ".fasta" };
  static char *sufidx[5] = { ".fastq.vzi", ".fasta.vzi", "", "" };

  struct stat stats;
  char  *pwd, *root, *path;
  int    fid, i;
  int64  fsize, zsize, *zoffs;
  int    fastq, zipd, recon;

  pwd = PathTo(arg);
  OPEN2(arg,pwd,root,fid,suffix,4)
  if (fid < 0)
    { fprintf(stderr,"%s: Cannot open %s as an .fast[aq] file\n",Prog_Name,arg);
      exit (1);
    }
  path  = Strdup(Catenate(pwd,"/",root,suffix[i]),"Allocating full path name");
  fastq = (i%2 == 0);
  zipd  = (i < 2);

  zoffs = NULL;
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
          zoffs = (int64 *) Malloc(sizeof(int64)*(zsize+1),"Allocating VGPzip index");
          if (zoffs == NULL)
            exit (1);
          read(idx,zoffs+1,sizeof(int64)*zsize);
          zoffs[0] = 0;
          fsize = IO_BLOCK*zsize;   //  An estimate only, upper bound
          close(idx);
        }
    }
  else
    { if (fstat(fid, &stats) == -1)
        { fprintf(stderr,"%s: Cannot get stats for %s\n",Prog_Name,arg);
          exit (1);
        }
      fsize = stats.st_size;
      zsize = (fsize-1)/IO_BLOCK+1;
      free(pwd);
    }

#ifdef DEBUG_FIND
  fprintf(stderr,"\n%s is a %s file, fastq = %d, zipd = %d, fsize = %lld\n",
                 arg,suffix[i],fastq,zipd,fsize);
#endif

  close(fid);

  input->path  = path;
  input->root  = root;
  input->fsize = fsize;
  input->fastq = fastq;
  input->zipd  = zipd;
  input->recon = recon;
  input->zsize = zsize;
  input->zoffs = zoffs;
}

static void Free_Fastq(File_Object *input)
{ if (input->recon)
    unlink(input->path);
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
  { File_Object *fobj;
    uint8       *buf;    //  block buffer
    DEPRESS     *decomp;
    uint8       *zuf;    //  decode buffer
    VgpFile     *vf;
    int          bidx;   //  Scan range is [bidx:beg,eidx:end)
    int64        beg;
    int          eidx;
    int64        end;
    char        *alast;
  } Thread_Arg;

  //  Automata states (buffer is processed char at a time)

#define UK    0
#define FK    1
#define FK1   2
#define QHL   3

#if defined(DEBUG_AUTO)

static char *Name[] =
  { "UK", "FK", "FK1", "QHL" };

#endif

  //  Check blocks [beg,end) starting with the first complete entry and proceeding
  //     into blocks end ... if necessary to complete the last entry started in this
  //     range.

static int64 find_nearest(int fid, Thread_Arg *data)
{ int64        beg    = data->beg;
  uint8       *buf    = data->buf;
  uint8       *zuf    = data->zuf;
  File_Object *inp    = data->fobj + data->bidx;
  DEPRESS     *decomp = data->decomp;

  int          fastq = inp->fastq;
  int64       *zoffs = inp->zoffs;

  int64 blk, off;
  int   slen;
  int   state, count;
  int   nl_1, nl_2, pl_1, alive;
  int   b, c;

  int   gogroup, lanmax, lanl;
  char *lane;

#ifdef DEBUG_FIND
  fprintf(stderr,"\nFind starting at %lld\n",beg);
#endif

  blk = beg / IO_BLOCK;
  off = beg % IO_BLOCK;

  if (inp->zipd)
    lseek(fid,zoffs[blk],SEEK_SET);
  else
    lseek(fid,blk*IO_BLOCK,SEEK_SET);

  if (fastq)
    { state = UK;
      nl_1 = pl_1 = 1;
      alive = 0;
    }
  else
    state = FK;
  count = 0;

  if (GROUP)
    { lanmax  = 1000;
      lane    = Malloc(1000,"Allocating lane buffer");
      lanl    = 0;
      gogroup = 1;
    }
  else
    gogroup = 0;

#ifdef DEBUG_FIND
  fprintf(stderr,"\nFrom block %lld / offset %lld\n",blk,off);
#endif

  while (blk < inp->zsize)
    {
#ifdef DEBUG_FIND
      fprintf(stderr,"  Loading block %lld: @%lld",blk,lseek(fid,0,SEEK_CUR));
#endif
      if (inp->zipd)
        { uint32 dlen, tlen;
          int    rez;
          size_t x;

          dlen = zoffs[blk+1]-zoffs[blk];
          tlen = IO_BLOCK;
          read(fid,zuf,dlen);
          if ((rez = libdeflate_gzip_decompress(decomp,zuf,dlen,buf,tlen,&x)) != 0)
            { fprintf(stderr,"Decompression not OK!\n");
              exit (1);
            }
          slen = (int) x;
#ifdef DEBUG_FIND
          fprintf(stderr," %d ->",dlen);
#endif
        }
      else
        slen = read(fid,buf,IO_BLOCK);
#ifdef DEBUG_FIND
      fprintf(stderr," %d\n",slen);
#endif

      for (b = off; b < slen; b++)
        { c = buf[b];
#ifdef DEBUG_AUTO
          fprintf(stderr,"  %.5s: %c\n",Name[state],c);
#endif
          switch (state)

          { case UK:
              if (c == '@' && alive)
                { if (gogroup)
                    { gogroup = 0;
                      state = QHL;
                      break;
                    }
                  else
                    return (blk*IO_BLOCK + b);
                }
              alive = (c == '\n' && ! (nl_2 || pl_1));
              nl_2 = nl_1;
              nl_1 = (c == '\n');
              pl_1 = (c == '+');
#ifdef DEBUG_AUTO
              fprintf(stderr,"    n2=%d n1=%d p1=%d a=%d\n",nl_2,nl_1,pl_1,alive);
#endif
              break;

            case FK:
              if (c == '\n')
                state = FK1;
              break;

            case FK1:
              if (c == '>')
                { if (gogroup)
                    { gogroup = 0; 
                      state = QHL;
                    }
                  else
                    return (blk*IO_BLOCK + b);
                }
              else if (c != '\n')
                state = FK;
              break;

            case QHL:
              if (c == ':')
                count += 1;
              if (count == 3 || c == '\n')
                { lane[lanl] = '\0';
#ifdef DEBUG_FIND
                  fprintf(stderr,"  lane = %s\n",lane);
#endif
                  data->alast = lane;
                  if (fastq)
                    state = UK;
                  else
                    state = FK;
                }
              else
                { if (lanl >= lanmax)
                    { lanmax = lanmax*1.2 + 1000;
                      lane = (char *) Realloc(lane,lanmax+1,"Reallocating lane name");
                      if (lane == NULL)
                        exit (1);
                    }
                  lane[lanl++] = c;
                }
              break;
           }
        }

      blk += 1;
      off = 0;
    }

  return (-1);
}


/*******************************************************************************************
 *
 *  Second pass thread routines prepare .seq formated output in buffers for each input block.
 *  Note carefully that there are separate blocks for the forward and reverse files if
 *  processing a pair, the reads in the blocks need to be paired.
 *
 ********************************************************************************************/

  //  Automata states (buffer is processed char at a time)

#define QAT     0
#define HEAD    1
#define HSKP    2
#define QSEQ    3
#define QPLS    4
#define QSKP    5
#define QQVS    6
#define ASEQ    7
#define AEOL    8

#if defined(DEBUG_AUTO) || defined(DEBUG_OUTPUT)

static char *Name2[] =
  { "QAT", "HEAD", "HSKP", "QSEQ", "QPLS", "QSKP", "QQVS", "ASEQ", "AEOL" };

#endif

  //  Write fast records in relevant partition

static void *output_thread(void *arg)
{ Thread_Arg  *parm   = (Thread_Arg *) arg;
  File_Object *fobj   = parm->fobj;
  VgpFile     *vf     = parm->vf;
  uint8       *buf    = parm->buf;
  uint8       *zuf    = parm->zuf;
  DEPRESS     *decomp = parm->decomp;
  int          fastq  = fobj->fastq;    //  Is the same for all files (already checked)

  File_Object *inp;
  int          f, fid;
  int64        blk, off;
  int64        epos, eblk, eoff;
  int64       *zoffs;

  int   state, count;

  int   lmax, llen;
  char *line;

  int   omax, olen;
  char *lane, *last;

  //  Do relevant section of each file assigned to this thread in sequence

  if (GROUP)
    { if (strlen(parm->alast) > 1000)
        lmax = strlen(parm->alast);
      else
        lmax = 1000;
      lane = Malloc(lmax+1,"Allocating lane buffer");
      last = Malloc(lmax+1,"Allocating lane buffer");
      llen = 0;
     
      strcpy(last,parm->alast);
      free(parm->alast);
    }

  omax = 100000;
  line = Malloc(omax+1,"Allocating line buffer");

  for (f = parm->bidx; f <= parm->eidx; f++)
    { inp = fobj+f;
      fid = open(inp->path,O_RDONLY);
      if (f < parm->eidx)
        epos = inp->fsize;
      else
        epos = parm->end;
      if (f > parm->bidx)
        parm->beg = 0;

#ifdef DEBUG_OUT
      fprintf(stderr,"Block: %12lld to %12lld --> %8lld\n",parm->beg,epos,epos - parm->beg);
      fflush(stdout);
#endif

      blk   = parm->beg / IO_BLOCK;
      off   = parm->beg % IO_BLOCK;
      eblk  = (epos-1) / IO_BLOCK;
      eoff  = (epos-1) % IO_BLOCK + 1;
      zoffs = inp->zoffs;

      if (inp->zipd)
        lseek(fid,zoffs[blk],SEEK_SET);
      else
        lseek(fid,blk*IO_BLOCK,SEEK_SET);

      state = QAT;
      llen  = 0;
      olen  = 0;
      count = 0;

#ifdef DEBUG_OUT
      fprintf(stderr,"\nFrom block %lld / offset %lld\n",blk,off);
#endif

      while (blk <= eblk)
        { int c, b, slen;

#ifdef DEBUG_OUT
      fprintf(stderr,"  Loading block %lld: @%lld",blk,lseek(fid,0,SEEK_CUR));
#endif
          if (inp->zipd)
            { uint32 dlen, tlen;
              int    rez;
              size_t x;

              dlen = zoffs[blk+1]-zoffs[blk];
              tlen = IO_BLOCK;
              read(fid,zuf,dlen);
              if ((rez = libdeflate_gzip_decompress(decomp,zuf,dlen,buf,tlen,&x)) != 0)
                { fprintf(stderr,"Decompression not OK!\n");
                  exit (1);
                }
              slen = (int) x;
#ifdef DEBUG_OUT
              fprintf(stderr," %d ->",dlen);
#endif
            }
          else
            slen = read(fid,buf,IO_BLOCK);
#ifdef DEBUG_OUT
          fprintf(stderr," %d\n",slen);
#endif

          if (blk == eblk && eoff < slen)
            slen = eoff;

          for (b = off; b < slen; b++)
            { c = buf[b];
#ifdef DEBUG_AUTO
              fprintf(stderr,"  %.5s: %c\n",Name2[state],c);
#endif
              switch (state)

              { case QAT:
                  if (GROUP)
                    state = HEAD;
                  else
                    state = HSKP;
                  break;

                case HSKP:
                  if (c == '\n')
                    { if (fastq)
                        state = QSEQ;
                      else
                        state = ASEQ;
                    }
                  break;
                  
                case HEAD:
                  if (c == ':')
                    count += 1;
                  if (count == 3 || c == '\n')
                    { lane[llen] = '\0';
#ifdef DEBUG_AUTO
                      fprintf(stderr,"  lane = %s\n",lane);
#endif
                      if (strcmp(lane,last) != 0)
                        { vgpInt(vf,0) = 0;
                          vgpInt(vf,1) = llen;
                          vgpWriteLine(vf,'g',lane);

                          { char *x; x = lane; lane = last; last = x; }
                        }
                      count = 0;
                      llen  = 0;
                      if (c == '\n')
                        { if (fastq)
                            state = QSEQ;
                          else
                            state = ASEQ;
                        }
                      else
                        state = HSKP;
                    }
                  else
                    { if (llen >= lmax)
                        { lmax = lmax*1.2 + 1000;
                          lane = (char *) Realloc(lane,lmax+1,"Reallocating lane name");
                          last = (char *) Realloc(last,lmax+1,"Reallocating lane name");
                          if (lane == NULL || last == NULL)
                            exit (1);
                        }
                      if (llen > 0 || c != '@' || !fastq)
                        lane[llen++] = c;
                    }
                  break;

                case QSEQ:
                  if (c != '\n')
                    { if (olen >= omax)
                        { omax = omax*1.2 + 1000;
                          line = (char *) Realloc(line,omax+1,"Reallocating line buffer");
                          if (line == NULL)
                            exit (1);
                        }
                      line[olen++] = c;
                    }
                  else
                    { vgpInt(vf,0) = olen;
                      vgpWriteLine(vf,'S',line);
                      olen  = 0;
                      state = QPLS;
                    }
                  break;

                case QPLS:
                  if (c == '\n')
                    { if (QVS_OUT)
                        state = QQVS;
                      else
                        state = QSKP;
                    }
                  break;

                case QSKP:
                  if (c == '\n')
                    state = QAT;
                  break;

                case QQVS:
                  if (c != '\n')
                    { if (olen >= omax)
                        fprintf(stderr,"%s: Fastq QV line longer than sequence line?\n",Prog_Name);
                      line[olen++] = c;
                    }
                  else
                    { vgpInt(vf,0) = olen;
                      vgpWriteLine(vf,'Q',line);
                      olen  = 0;
                      state = QAT;
                    }
                  break;

                case ASEQ:
                  if (c == '\n')
                    state = AEOL;
                  else
                    { if (olen >= omax)
                        { omax = omax*1.2 + 1000;
                          line = (char *) Realloc(line,omax+1,"Reallocating line buffer");
                          if (line == NULL)
                            exit (1);
                        }
                      line[olen++] = c;
                    }
                  break;

                case AEOL:
                  if (c == '>')
                    { vgpInt(vf,0) = olen;
                      vgpWriteLine(vf,'S',line);
                      olen  = 0;
                      if (GROUP)
                        state = HEAD;
                      else
                        state = HSKP;
                    }
                  else
                    state = ASEQ;
                  break;
              }
            }
          blk += 1;
          off = 0;
        }
      if (state == AEOL)
        { vgpInt(vf,0) = olen;
          vgpWriteLine(vf,'S',line);
          olen  = 0;
        }
      close(fid);
    }

  free(line);
  free(last);
  free(lane);

  return (NULL);
}


/****************************************************************************************
 *
 *  The top-level program
 *
 ****************************************************************************************/

int main(int argc, char *argv[])
{ char       *command;
  int         nfiles;

  //  Capture command line for provenance

  { int   n, i;
    char *c;

    n = -1;
    for (i = 1; i < argc; i++)
      n += strlen(argv[i])+1;

    command = Malloc(n+1,"Allocating command string");
    if (command == NULL)
      exit (1);

    c = command;
    if (argc >= 1)
      { c += sprintf(c,"%s",argv[1]);
        for (i = 2; i < argc; i++)
          c += sprintf(c," %s",argv[i]);
      }
    *c = '\0';
  }

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

    if (argc < 2)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: verbose mode, output progress as proceed\n");
        fprintf(stderr,"      -s: Output sequences only, skip QVs\n");
        fprintf(stderr,"      -g: Output group lanes based on Illumina lanes\n");
        fprintf(stderr,"      -T: Number of threads to use\n");
        exit (1);
      }
  }


  //  Determine number of .las files, nfiles

  nfiles = argc-1;
  //  Find partition points dividing data in all files into NTHREADS roughly equal parts
  //    and then in parallel threads produce the output for each part.

  { File_Object fobj[nfiles];
    Thread_Arg  parm[NTHREADS];

#ifndef DEBUG_OUT
    pthread_t   threads[NTHREADS];
#endif

    { uint8        *bf;
      int           f, i, c;
      int64         b, work, wper;

      //  Allocate IO buffer space for threads

      bf = Malloc(2*NTHREADS*IO_BLOCK,"Allocating IO_Buffer\n");
      if (bf == NULL)
        exit (1);

      //  Get name and size of each file in 'fobj[]', determine type, etc.

      work   = 0;
      for (c = 1; c < argc; c++)
        { Fetch_Fastq(argv[c],fobj+(c-1));

          if (c > 1 && fobj[1].fastq != fobj[c-1].fastq)
            { fprintf(stderr,"%s: All files must either be .fastq or .fastq but not a mix\n",
                             Prog_Name);
              exit (1);
            }

          work += fobj[c-1].fsize;
        }

      if (VERBOSE)
        { fprintf(stderr,"  Partitioning %d %sfast%c files into %d parts\n",
                         nfiles,fobj->zipd?"compressed ":"",fobj->fastq?'q':'a',NTHREADS);
          fflush(stderr);
        }

      //  Allocate work evenly amongst threads, setting up search start
      //    point for each thread.  Also find the beginning of data in
      //    each file that a thread will start in (place in end.fpos)

      wper = work / NTHREADS;

      f = 0;
      b = 0;
      work = fobj[f].fsize;
      for (i = 0; i < NTHREADS; i++)
        { parm[i].fobj   = fobj;
          parm[i].buf    = bf + 2*i*IO_BLOCK;
          parm[i].zuf    = bf + (2*i+1)*IO_BLOCK;
          parm[i].decomp = libdeflate_alloc_decompressor();

          parm[i].end = 0;   //  Could special case, but keep general template

#ifdef DEBUG_FIND
          fprintf(stderr," %2d: %1d %10lld (%10lld)\n",i,f,b,parm[i].end);
          fflush(stdout);
#endif

          parm[i].beg  = b;
          parm[i].bidx = f;

          work -= wper;
          while (work < .01*IO_BLOCK)
            { if (f == nfiles)
                { NTHREADS = i+1;
                  break;
                }
              work += fobj[++f].fsize;
            }
          b = fobj[f].fsize - work;
          if (b < 0)
            { work += b;
              b = 0;
            }
        }
    }

    { int i, f, fid;

      //  For each non-zero start point find synchronization point in
      //    bam/sam file.  If can't find it then start at beginning of
      //    next file, and if at or before first data line then signal
      //    start at beginning by zero'ing the synch point.

      for (i = 0; i < NTHREADS; i++)
        { if (parm[i].beg != 0)
            { fid = open(fobj[parm[i].bidx].path,O_RDONLY);
              parm[i].beg = find_nearest(fid,parm+i);
              if (parm[i].beg < 0 || parm[i].beg >= fobj[parm[i].bidx].fsize)
                { parm[i].beg   = 0;
                  parm[i].bidx += 1;
                }
              else if (parm[i].beg <= parm[i].end)
                parm[i].beg = 0;
              close(fid);
            }
          if (parm[i].beg == 0)
            parm[i].alast = Strdup("","Allocating empty string");;
        }

      //  Paranoid: if one thread's synch point overtakes the next one (will almost
      //    certainly never happen unless files very small and threads very large),
      //    remove the redundant threads.

      f = 0;
      for (i = 1; i < NTHREADS; i++)
        if (parm[i].bidx > parm[f].bidx || parm[i].beg > parm[f].beg)
          parm[++f] = parm[i];
        else
          { free(parm[i].alast);
            libdeflate_free_decompressor(parm[i].decomp);
          }
      NTHREADS = f+1;

      //  Develop end points of each threads work using the start point of the next thread

      for (i = 1; i < NTHREADS; i++)
        if (parm[i].beg == 0)
          { parm[i-1].end  = fobj[parm[i].bidx-1].fsize;
            parm[i-1].eidx = parm[i].bidx-1;
          }
        else
          { parm[i-1].end  = parm[i].beg;
            parm[i-1].eidx = parm[i].bidx;
          }
      parm[NTHREADS-1].end  = fobj[nfiles-1].fsize;
      parm[NTHREADS-1].eidx = nfiles-1;

#if defined(DEBUG_FIND) || defined(DEBUG_OUT)
      fprintf(stderr,"\nPartition:\n");
      for (i = 0; i < NTHREADS; i++)
        { fprintf(stderr," %2d: %2d / %12lld",i,parm[i].bidx,parm[i].beg);
          fprintf(stderr,"  -  %2d / %12lld\n",parm[i].eidx,parm[i].end);
        }
      fflush(stdout);
#endif
    }

    //  Produce output in parallel threads based on partition

    { VgpFile *vf;
      int      i;

      vf = vgpFileOpenWriteNew("-",SEQ,0,TRUE,NTHREADS);

      vgpAddProvenance(vf,Prog_Name,"1.0",command,NULL);

      vgpWriteHeader(vf);

#ifdef DEBUG_OUT
      fprintf(stderr,"Opened\n");
      fflush(stdout);
#endif

      if (VERBOSE)
        { fprintf(stderr,"  Producing .seq segments in parallel\n");
          fflush(stderr);
        }

      //  Generate the data lines in parallel threads

#ifdef DEBUG_OUT
      for (i = 1; i < NTHREADS; i++)
        { parm[i].vf = vf+i;
          output_thread(parm+i);
        }
#else
      for (i = 0; i < NTHREADS; i++)
        { parm[i].vf = vf+i;
          pthread_create(threads+i,NULL,output_thread,parm+i);
        }

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);
#endif

      if (VERBOSE)
        { fprintf(stderr,"  Cat'ing .seq segments\n");
          fflush(stderr);
        }

      vgpFileClose(vf);
    }

    //  Free everything as a matter of good form

    { int i, f;

      for (i = 0; i < NTHREADS; i++)
        libdeflate_free_decompressor(parm[i].decomp);
      free(parm[0].buf);
      for (f = 0; f < nfiles; f++)
        Free_Fastq(fobj+f);
      free(command);
    }
  }

  if (VERBOSE)
    { fprintf(stderr,"  Done\n");
      fflush(stderr);
    }

  exit (0);
}
