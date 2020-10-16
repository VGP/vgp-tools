/*******************************************************************************************
 *
 *  VGPseq can read fasta, fastq, sam, bam, and cram files and outputs to the standard output
 *    either a .seq (default) or a .irp (if -p set) file.  If paired reads are requested
 *    then the paired reads must be adjacent in the given file type.  VGPzip'd fasta or
 *    fastq files are processed directly.  Gzip'd fasta and fastq files are supported but
 *    be advised that they are first decompressed into a temporary file.
 *
 *  Author:    Gene Myers
 *  Creation:  April 2020
 *
 *******************************************************************************************/

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
#include <stdbool.h>

#include "VGPschema.h"

#include "gene_core.h"
#include "../Core/ONElib.h"
#include "LIBDEFLATE/libdeflate.h"
#include "HTSLIB/htslib/hts.h"
#include "HTSLIB/htslib/hfile.h"
#include "HTSLIB/cram/cram.h"

#undef    DEBUG_FIND
#undef    DEBUG_OUT

#undef    DEBUG_AUTO

#undef    DEBUG_BAM_IO
#undef    DEBUG_BAM_RECORD

static char *Usage = "[-viqp] [-g#x] [-T<int(4)>] <data:cram|[bs]am|f{ast}[aq][.gz]> ...";

#define IO_BLOCK 10000000ll

typedef struct libdeflate_decompressor DEPRESS;

#define INT_MAXLEN 10

int    NTHREADS; //  # of threads to use
int    IS_BIG_ENDIAN;   //  Is machine big-endian?

int    VERBOSE;  //  Verbose mode?
int    QVS_OUT;  //  Include QV strings
int    QNAME;    //  Include identifier
int    GROUP;    //  Group lanes of reads
int      GROUP_REP;
int      GROUP_CHAR;
int    PAIRING;  //  Pair reads if possible

#define CRAM   0
#define BAM    1
#define SAM    2
#define FASTQ  3
#define FASTA  4

static char *Tstring[5] = { "cram", "bam", "sam", "fastq", "fasta" };

typedef struct
  { char   *path;   //  Full path name
    char   *root;   //  Ascii root name
    int64   fsize;  //  Size of (uncompressed) file in bytes
    int     ftype;  //  Type of file
                    //  Fasta/q and Cram specific:
    int64   zsize;  //    Size of zip index (in blocks)
    int64  *zoffs;  //    zoffs[i] = offset to compressed block i
                    //  Fasta/q specific:
    int     zipd;   //    Is file VGPzip'd (has paired .vzi file)
    int     recon;  //    Is file an uncompressed regular zip-file?
  } File_Object;

typedef struct
  { int64  fpos;   //  offset in source file of desired unit (gzip block, cram contaianer)
    uint32 boff;   //  offset within a unit block of desired record
  } Location;

typedef struct
  { File_Object *fobj;   //  array of file objects
    uint8       *buf;    //  block buffer
    int          fid;    //  fid of fobj[bidx].path
    OneFile     *vf;     //  One-file for output
    int          bidx;   //  Scan range is [bidx:beg,eidx:end)
    Location     beg; 
    int          eidx;
    Location     end;
    int          error;  //    pairing not making sense
                         //  fasta/q & bam specific:
    DEPRESS     *decomp; //    decompressor
                         //  fasta/q specific:
    char        *alast;  //    last group name
    uint8       *zuf;    //    decode buffer (for zip'd files)
  } Thread_Arg;

static void do_nothing(Thread_Arg *parm)
{ (void) parm; }



/*******************************************************************************************
 *
 *  Routines to open and close the sequence file types supported
 *
 ********************************************************************************************/

  //  Open and get info about each input file

static void Fetch_File(char *arg, File_Object *input)
{ static char *suffix[11] = { ".cram", ".bam", ".sam",
                              ".fastq.gz",   ".fasta.gz", ".fastq", ".fasta",
                              ".fq.gz",  ".fa.gz", ".fq", ".fa" };
  static char *sufidx[11] = { "", "", "", ".fastq.vzi", ".fasta.vzi", "", "", ".fq.vzi", ".fa.vzi", "", "" };

  int64 *genes_cram_index(char *path, int64 fsize, int64 *zsize);

  struct stat stats;
  char  *pwd, *root, *path;
  int    fid, i;
  int64  fsize, zsize, *zoffs;
  int    ftype, zipd, recon;

  pwd = PathTo(arg);
  OPEN2(arg,pwd,root,fid,suffix,11)
  if (fid < 0)
    { fprintf(stderr,"%s: Cannot open %s as a .cram|[bs]am|f{ast}[aq][.gz] file\n",Prog_Name,arg);
      exit (1);
    }
  path  = Strdup(Catenate(pwd,"/",root,suffix[i]),"Allocating full path name");
  ftype = i;
  if (i >= 3)
    if (i >= 7)
      { ftype = 4 - (i%2);
        zipd  = (i < 9);
      }
    else
      { ftype = 4 - (i%2);
        zipd  = (i < 5);
      }
  else
    { ftype = i;
      zipd  = 0;
    }

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
      if (ftype == CRAM)
        zoffs = genes_cram_index(path,fsize,&zsize);
      else
        zsize = (fsize-1)/IO_BLOCK+1;
      free(pwd);
    }

#ifdef DEBUG_FIND
  fprintf(stderr,"\n%s is a %s file, ftype = %d, zipd = %d, fsize = %lld\n",
                 arg,suffix[i],ftype,zipd,fsize);
#endif

  close(fid);

  input->path  = path;
  input->root  = root;
  input->fsize = fsize;
  input->ftype = ftype;
  input->zipd  = zipd;
  input->recon = recon;
  input->zsize = zsize;
  input->zoffs = zoffs;
}

static void Free_File(File_Object *input)
{ if (input->recon)
    unlink(input->path);
  free(input->zoffs);
  free(input->path);
  free(input->root);
}



/*******************************************************************************************
 *
 *  FASTA / FASTQ SPECIFIC CODE
 *
 ********************************************************************************************/

/*******************************************************************************************
 *
 *  fast_nearest: Routine to find next entry given arbitray start point.
 *
 ********************************************************************************************/

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

static void fast_nearest(Thread_Arg *data)
{ int          fid    = data->fid;
  int64        beg    = data->beg.fpos;
  uint8       *buf    = data->buf;
  uint8       *zuf    = data->zuf;
  File_Object *inp    = data->fobj + data->bidx;
  DEPRESS     *decomp = data->decomp;

  int          fastq = (inp->ftype == FASTQ);
  int64       *zoffs = inp->zoffs;

  int64 blk, off;
  int   slen;
  int   state;
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

  gogroup = GROUP;
  if (GROUP)
    { lanmax  = 1000;
      lane    = Malloc(1000,"Allocating lane buffer");
      lanl    = 0;
    }

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
                      state   = QHL;
                      break;
                    }
                  else
                    { data->beg.fpos = blk*IO_BLOCK + b;
                      data->beg.boff = 0;
                      return;
                    }
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
                      state   = QHL;
                      break;
                    }
                  else
                    { data->beg.fpos = blk*IO_BLOCK + b;
                      data->beg.boff = 0;
                      return;
                    }
                }
              else if (c != '\n')
                state = FK;
              break;

            case QHL:
              if (c == '\n' || isspace(c))
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

  data->beg.fpos = -1;
  data->beg.boff = 0;
  return;
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

#if defined(DEBUG_AUTO) || defined(DEBUG_OUT)

static char *Name2[] =
  { "QAT", "HEAD", "HSKP", "QSEQ", "QPLS", "QSKP", "QQVS", "ASEQ", "AEOL" };

#endif

  //  Write fast records in relevant partition

static void *fast_output_thread(void *arg)
{ Thread_Arg  *parm   = (Thread_Arg *) arg;
  File_Object *fobj   = parm->fobj;
  OneFile     *vf     = parm->vf;
  uint8       *buf    = parm->buf;
  uint8       *zuf    = parm->zuf;
  DEPRESS     *decomp = parm->decomp;
  int          fastq  = (fobj->ftype == FASTQ);    //  Is the same for all files (already checked)

  File_Object *inp;
  int          f, fid;
  int64        blk, off;
  int64        epos, eblk, eoff;
  int64       *zoffs;

  int   state, count;

  int   lmax, llen, glen, ilen;
  char *lane, *last;
  int   first;

  int   omax, olen;
  char *line;

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

  first = 2;
  for (f = parm->bidx; f <= parm->eidx; f++)
    { inp = fobj+f;
      fid = open(inp->path,O_RDONLY);
      if (f < parm->eidx)
        epos = inp->fsize;
      else
        epos = parm->end.fpos;
      if (f > parm->bidx)
        parm->beg.fpos = 0;

#ifdef DEBUG_OUT
      fprintf(stderr,"Block: %12lld to %12lld --> %8lld\n",parm->beg,epos,epos - parm->beg);
      fflush(stdout);
#endif

      blk   = parm->beg.fpos / IO_BLOCK;
      off   = parm->beg.fpos % IO_BLOCK;
      eblk  = (epos-1) / IO_BLOCK;
      eoff  = (epos-1) % IO_BLOCK + 1;
      zoffs = inp->zoffs;

      if (inp->zipd)
        lseek(fid,zoffs[blk],SEEK_SET);
      else
        lseek(fid,blk*IO_BLOCK,SEEK_SET);

      state = QAT;
      olen  = 0;
      llen  = 0;
      glen  = 0;
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
                  if (GROUP || QNAME)
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
                  if (c == GROUP_CHAR)
                    { count += 1;
                      if (count == GROUP_REP)
                        glen = llen;
                    }
                  if (c == '\n' || isspace(c))
                    { char *x;

                      lane[llen] = '\0';
                      if (glen <= 0)
                        glen = llen;
#ifdef DEBUG_AUTO
                      fprintf(stderr,"  header = %s\n",lane);
#endif

                      if (GROUP)
                        { if (strncmp(lane,last,glen) != 0)
                            { oneInt(vf,0) = 0;
                              oneInt(vf,1) = glen;
                              oneWriteLine(vf,'g',glen,lane);
                            }
                        }

                      if (PAIRING)
                        { if (strncmp(lane,last,llen) != 0)
                            { oneWriteLine(vf,'P',0,NULL);
                              if (first == 0)
                                { parm->error = 1;
                                  return (NULL);
                                }
                              first = 0;
                            }
                          else
                            { if (first == 1)
                                { parm->error = 1;
                                  return (NULL);
                                }
                              first = 1;
                            }
                        }

                      x = lane;
                      lane = last;
                      last = x;
                      ilen = llen;

                      count = 0;
                      glen  = 0;
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
                    { oneInt(vf,0) = olen;
                      oneWriteLine(vf,'S',olen,line);
                      if (QNAME)
                        { oneInt(vf,0) = ilen;
                          oneWriteLine(vf,'I',ilen,last);
                        }
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
                    { if (QVS_OUT)
                        { oneInt(vf,0) = olen;
                          oneWriteLine(vf,'Q',olen,line);
                        }
                      olen  = 0;
                      state = QAT;
                    }
                  break;

                case AEOL:
                  if (c == '>')
                    { oneInt(vf,0) = olen;
                      oneWriteLine(vf,'S',olen,line);
                      if (QNAME)
                        { oneInt(vf,0) = ilen;
                          oneWriteLine(vf,'I',ilen,last);
                        }
                      olen  = 0;
                      if (GROUP || QNAME)
                        state = HEAD;
                      else
                        state = HSKP;
                      break;
                    }

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
              }
            }
          blk += 1;
          off = 0;
        }
      if (state == AEOL)
        { oneInt(vf,0) = olen;
          oneWriteLine(vf,'S',olen,line);
          olen  = 0;
        }
      close(fid);
    }

  free(line);
  free(last);
  free(lane);

  parm->error = 0;
  return (NULL);
}


/*******************************************************************************************
 *
 *  BAM/SAM SPECIFIC CODE
 *
 ********************************************************************************************/

#define BAM_BLOCK  0x10000
#define HEADER_LEN      36
#define SEQ_RUN         40

static int DNA[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
  };

 //  Get value of little endian integer of n-bytes

static inline uint32 getint(uint8 *buf, int n)
{ uint32 val;
  int    k;

  val = 0;
  for (k = n-1; k >= 0; k--)
    val = (val << 8) | buf[k];
  return (val);
}

 //  Next len chars are printable and last is zero?

static inline int valid_name(char *name, int len)
{ int i;

  if (len < 1)
    return (0);
  len -= 1;
  for (i = 0; i < len; i++)
    if ( ! isprint(name[i]))
      return (0);
  return (name[len] == 0);
}

 //  Next len ints are plausible cigar codes

static inline int valid_cigar(int32 *cigar, int len, int range)
{ int i, c;

  for (i = 0; i < len; i++)
    { c = cigar[i];
      if ((c & 0xf) > 8)
        return (0);
      c >>= 4;
      if (c < 0 || c > range)
        return (0);
    }
  return (1);
}


/*******************************************************************************************
 *
 *  Routines to manage a BAM_FILE stream object
 *
 ********************************************************************************************/

  //   Open BAM stream where compressed block start is known.  Compressed BAM blocks are buffered
  //     in a very big IO buffer and the current uncompressed block is in a 64Kbp array.

typedef struct
  { int       fid;              //  file descriptor
    int       last;             //  last block of data in file has been read
    uint8    *buf;              //  IO buffer (of IO_BLOCK bytes, supplied by caller)
    int       blen;             //  # of bytes currently in IO buffer
    int       bptr;             //  start of next BAM block in IO buffer
    uint8     bam[BAM_BLOCK+1]; //  uncompressed bam block
    uint32    bsize;            //  length of compressed bam block
    uint32    ssize;            //  length of uncompressed bam block
    Location  loc;              //  current location in bam file
    DEPRESS  *decomp;
  } BAM_FILE;

  //  Load next len bytes of uncompressed BAM data into array data

static void bam_get(BAM_FILE *file, uint8 *data, int len)
{ int    chk, off;
  uint8 *bam  = file->bam;
  int    boff = file->loc.boff;

  off = 0;
  chk = file->ssize - boff;
  while (len >= chk)
    { int    bptr, blen, bsize;
      uint8 *block, *buf;
      uint32 ssize;
      size_t tsize;

#ifdef DEBUG_BAM_IO
      printf("Move %d bytes to data+%d from bam+%d\n",chk,off,boff);
#endif
      if (data != NULL)
        memcpy(data+off,bam+boff,chk);
      off += chk;
      len -= chk;

      file->loc.fpos += file->bsize;
#ifdef DEBUG_BAM_IO
      printf("File pos %lld\n",file->fpos);
#endif

      if (chk == 0 && len == 0)
        { file->loc.boff = 0;
          return;
        }

      buf   = file->buf;
      bptr  = file->bptr;
      blen  = file->blen;
      block = buf+bptr;
#ifdef DEBUG_BAM_IO
      printf("Block at buffer %d\n",bptr);
#endif
      while (bptr + 18 > blen || bptr + (bsize = getint(block+16,2) + 1) > blen)
        { chk = blen-bptr;
          if (file->last)
            { fprintf(stderr, "%s: Corrupted BAM file\n",Prog_Name);
              exit (1);
            }
          memmove(buf,block,chk);
          blen = chk + read(file->fid,buf+chk,IO_BLOCK-chk);
#ifdef DEBUG_BAM_IO
          printf("Loaded %d to buf+%d for a total of %d\n",IO_BLOCK-chk,chk,blen);
#endif
          if (blen < IO_BLOCK)
            file->last = 1;
          file->blen = blen;
          bptr = 0;
          block = buf;
        }

      //  Fetch and uncompress next Bam block

      if (libdeflate_gzip_decompress(file->decomp,block,bsize,bam,BAM_BLOCK,&tsize) != 0)
        { fprintf(stderr,"%s: Bad gzip block\n",Prog_Name);
          exit (1);
        }
      ssize = tsize;
      boff = 0;
#ifdef DEBUG_BAM_IO
      printf("Loaded gzip block of size %d into %d\n",bsize,ssize);
#endif

      file->bsize = bsize;
      file->ssize = ssize;
      file->bptr  = bptr + bsize;
      chk = ssize;
    }

#ifdef DEBUG_BAM_IO
  printf("Xfer %d bytes to data+%d from bam+%d\n",len,off,boff);
#endif
  if (data != NULL)
    memcpy(data+off,bam+boff,len);
  file->loc.boff = boff+len;
}

  //  Startup a bam stream, the location must be valid.

static void bam_start(BAM_FILE *file, int fid, uint8 *buf, Location *loc)
{ file->fid   = fid;
  file->ssize = 0;
  file->bsize = 0;
  file->buf   = buf;
  file->bptr  = 0;
  file->blen  = 0;
  file->last  = 0;
  lseek(fid,loc->fpos,SEEK_SET);
  file->loc.fpos = loc->fpos;
  file->loc.boff = 0;
  bam_get(file,NULL,loc->boff);
}

static int bam_eof(BAM_FILE *file)
{ return (file->loc.boff == file->ssize && file->bptr == file->blen && file->last); }


/*******************************************************************************************
 *
 *  Routines to manage a SAM stream, but as a BAM_FILE (only select fields are used)
 *
 ********************************************************************************************/

  //  Get next line of SAM input if possible

static uint8 *sam_getline(BAM_FILE *file)
{ int    rem;
  int    blen = file->blen;
  int    bptr = file->bptr;
  uint8 *buf  = file->buf;
  uint8  *b, *d;

  b = buf + bptr;
  rem = blen-bptr;
  if (rem == 0)
    d = NULL;
  else
    d = memchr(b,'\n',rem);
  if (d == NULL)
    { if (file->last)
        { fprintf(stderr,"%s: Corrupted SAM file",Prog_Name);
          exit (1);
        }
      memmove(buf,buf+bptr,rem);
      blen = rem + read(file->fid,buf+rem,IO_BLOCK-rem);
      if (blen < IO_BLOCK)
        file->last = 1;
      file->blen = blen;
      bptr = 0;
      b = buf;
      d = memchr(b,'\n',blen);
      if (d == NULL)
        { if (blen < IO_BLOCK)
            fprintf(stderr,"%s: Corrupted SAM file",Prog_Name);
          else
            fprintf(stderr,"%s: SAM-line is longer than max %lld\n",Prog_Name,IO_BLOCK);
          exit (1);
        }
    }
  d += 1;
  file->bptr = d-buf;
  file->loc.fpos += d-b;
  return (b);
}

  //  Startup a sam stream

static void sam_start(BAM_FILE *file, int fid, uint8 *buf, Location *loc)
{ file->fid   = fid;
  file->buf   = buf;
  file->bptr  = 0;
  file->blen  = 0;
  file->last  = 0;
  lseek(fid,loc->fpos,SEEK_SET);
  file->loc.fpos = loc->fpos;
  file->loc.boff = 0;
}


/*******************************************************************************************
 *
 *  Routines to find bam blocks and valid locations to start scan threads for first pass
 *
 ********************************************************************************************/

  //  Find first record location (skip over header) in parm->fid
  //    Return value is in parm->beg

static void skip_bam_header(Thread_Arg *parm)
{ uint8    *buf = parm->buf;
  int       fid = parm->fid;

  BAM_FILE  _bam, *bam = &_bam;
  Location  zero = { 0ll, 0 };
  uint8     data[4];
  int       i, ntxt, ncnt, nlen;

  //  At start of file so can use BAM stream

  bam->decomp = parm->decomp;
  bam_start(bam,fid,buf,&zero);

#ifdef DEBUG_FIND
  fprintf(stderr,"Header seek\n");
  fflush(stderr);
#endif

  bam_get(bam,data,4);
  if (memcmp(data,"BAM\1",4) != 0)
    { fprintf(stderr, "%s: Corrupted BAM header %.4s\n",Prog_Name,data);
      exit (1);
    }

  bam_get(bam,data,4);
  ntxt = getint(data,4);
  bam_get(bam,NULL,ntxt);

  bam_get(bam,data,4);
  ncnt = getint(data,4);
  for (i = 0; i < ncnt; i++)
    { bam_get(bam,data,4);
      nlen = getint(data,4);
      bam_get(bam,NULL,nlen+4);
    }

  parm->beg = bam->loc;

#ifdef DEBUG_FIND
  fprintf(stderr,"  Begin @ %lld/%d\n",parm->beg.fpos,parm->beg.boff);
  fflush(stderr);
#endif
}

  //  Find next identifiable entry location forward of parm->fpos in parm->fid
  //    Return value is in parm->beg

static void bam_nearest(Thread_Arg *parm)
{ uint8       *buf  = parm->buf;
  int          fid  = parm->fid;
  DEPRESS     *decomp = parm->decomp;
  int64        fpos = parm->beg.fpos;

  BAM_FILE  _bam, *bam = &_bam;

  uint32 bptr, blen;
  int    last, notfound;

  uint8 *block;
  uint32 bsize, ssize;
  size_t tsize;

#ifdef DEBUG_FIND
  fprintf(stderr,"Searching from %lld\n",fpos);
  fflush(stderr);
#endif

  lseek(fid,fpos,SEEK_SET);
  blen = 0;
  bptr = 0;
  last = 0;

  //  Search until find a gzip block header

  notfound = 1;
  while (notfound)
    { int    j;
      uint32 isize, crc;

      fpos += bptr;      //   Get more data at level of IO blocks
      if (last)
        { fprintf(stderr,"%s: Could not find bam block structure!\n",Prog_Name);
          exit (1);
        }
      else
        { uint32 x = blen-bptr;
          memmove(buf,buf+bptr,x);
          blen = x + read(fid,buf+x,IO_BLOCK-x);
          if (blen < IO_BLOCK)
            last = 1;
#ifdef DEBUG_FIND
          fprintf(stderr,"Loading %d(last=%d)\n",blen,last);
          fflush(stderr);
#endif
          bptr = 0;
        }

      while (bptr < blen)          //  Search IO block for Gzip block start
        { if (buf[bptr++] != 31)
            continue;
          if ( buf[bptr] != 139)
            continue;
          bptr += 1;
          if ( buf[bptr] != 8)
            continue;
          bptr += 1;
          if ( buf[bptr] != 4)
            continue;
  
#ifdef DEBUG_FIND
          fprintf(stderr,"  Putative header @ %d\n",bptr-3);
          fflush(stderr);
#endif

          if (bptr + 12 > blen)
            { if (last)
                continue;
              bptr -= 3;
              break;
            }
  
          j = bptr+9;
          if (buf[j] != 66)
		  continue;
          j += 1;
          if (buf[j] != 67)
            continue;
          j += 1;
          if (buf[j] != 2)
            continue;
          j += 1;
          if (buf[j] != 0)
            continue;
          j += 1;
    
          bsize = getint(buf+j,2)+1;
          block = buf+(bptr-3);

          if ((bptr-3) + bsize > blen)
            { if (last)
                continue;
              bptr -= 3;
              break;
            }

#ifdef DEBUG_FIND
          fprintf(stderr,"    Putative Extra %d\n",bsize);
          fflush(stderr);
#endif
  
          isize = getint(block+(bsize-4),4);
          crc   = getint(block+(bsize-8),4);
  
          if (libdeflate_gzip_decompress(decomp,block,bsize,bam->bam,BAM_BLOCK,&tsize) != 0)
            continue;
          ssize = tsize;

          if (ssize == isize && crc == libdeflate_crc32(0,bam->bam,ssize))
            { bptr -= 3;
              fpos  += bptr;
              notfound = 0;

#ifdef DEBUG_FIND
              fprintf(stderr,"    First block at %lld (%d)\n",fpos,ssize);
              fflush(stderr);
#endif
	      break;
            }
        }
    }

  //  Have found a gzip/bam block start, now scan blocks until can identify the start
  //    of a sequence entry

  bam->fid      = fid;      //  Kick-start BAM stream object
  bam->last     = last;
  bam->buf      = buf;
  bam->blen     = blen;
  bam->bptr     = bptr+bsize;
  bam->bsize    = bsize;
  bam->ssize    = ssize;
  bam->loc.fpos = fpos;
  bam->loc.boff = 0;
  bam->decomp   = decomp;

  while ( ! bam_eof(bam))
    { int    j, k;
      int    run, out, del;
      int    lname, lcigar, lseq, ldata;

      block = bam->bam;
      ssize = bam->ssize;

      run = HEADER_LEN-1;
      out = 1;
      for (j = HEADER_LEN; j < 10000; j++)
        if (DNA[block[j]])
          { if (out && j >= run+SEQ_RUN)
              {
#ifdef DEBUG_FIND
                fprintf(stderr,"      Possible seq @ %d\n",run+1);
                fflush(stderr);
#endif
                for (k = run-(HEADER_LEN-1); k >= 0; k--)
                  { ldata  = getint(block+k,4);
                    lname  = block[k+12];
                    lcigar = getint(block+(k+16),2);
                    lseq   = getint(block+(k+20),4);
                    if (lname > 0 && lcigar >= 0 && lseq > 0 &&
                        (lseq+1)/2+lseq+lname+(lcigar<<2) < ldata)
                      { del = (k+35+lname+(lcigar<<2)) - run;
                        if (del >= 0 && del < SEQ_RUN/2)
                          { if (valid_name((char *) (block+(k+36)),lname) &&
                                valid_cigar((int32 *) (block+(k+36+lname)),lcigar,lseq))
                              { parm->beg.fpos = bam->loc.fpos;
                                parm->beg.boff = k;
#ifdef DEBUG_FIND
                                fprintf(stderr,"      Found @ %d: '%s':%d\n",k,block+(k+36),lseq);
                                fflush(stderr);
#endif

                                close(fid);

                                return;
                              }
                          }
                      }
                  }
                out = 0;
              }
          }
        else
          { run = j;
            out = 1;
          }

      bam_get(bam,NULL,ssize);
    }

  parm->beg.fpos = -1;
}

  //  Find next identifiable sam entry location forward of parm->fpos in parm->fid
  //    Return location is in parm->beg.  NB: works to skip sam header as well

static void sam_nearest(Thread_Arg *parm)
{ uint8       *buf  = parm->buf;
  int          fid  = parm->fid;

  BAM_FILE  _bam, *bam = &_bam;

  bam->decomp = parm->decomp;
  sam_start(bam,fid,buf,&(parm->beg));

  sam_getline(bam);
  if (parm->beg.fpos == 0)
    { while (buf[bam->bptr] == '@')
        sam_getline(bam);
    }
  parm->beg = bam->loc;
}


/*******************************************************************************************
 *
 *  Routines to scan and parse bam and sam entries
 *
 ********************************************************************************************/

typedef struct
  { int    hlen;
    char  *header;
    uint32 flags;
    int    len;
    char  *seq;
    char  *arr;
    char  *qvs;
    int    lmax;     //  current max size for seq, arr, and qvs
    int    dmax;     //  current max size for data
    uint8 *data;     //  data buffer
  } samRecord;

static char INT_2_IUPAC[16] = "=acmgrsvtwyhkdbn";

  //  Scan next bam entry and load PacBio info in record 'theR'
 
static int bam_record_scan(BAM_FILE *sf, samRecord *theR)
{ int ldata, lname, lseq, lcigar, aux;

  { uint8  x[36];     //  Process 36 byte header

    bam_get(sf,x,36);

    ldata  = getint(x,4) - 32;
    lname  = getint(x+12,1);
    lcigar = getint(x+16,2);
    lseq   = getint(x+20,4);

    if (ldata < 0 || lseq < 0 || lname < 1)
      { fprintf(stderr,"%s: Non-sensical BAM record, file corrupted?\n",Prog_Name);
        exit (1);
      }

    aux = lname + ((lseq + 1)>>1) + lseq + (lcigar<<2);
    if (aux > ldata)
      { fprintf(stderr,"%s: Non-sensical BAM record, file corrupted?\n",Prog_Name);
        exit (1);
      }

    if (lseq > theR->lmax)
      { theR->lmax = 1.2*lseq + 1000;
        theR->seq  = (char *) Realloc(theR->seq,3*theR->lmax,"Reallocating sequence buffer");
        if (theR->seq == NULL)
          exit (1);
        theR->arr  = theR->seq + theR->lmax;
        theR->qvs  = theR->arr + theR->lmax;
      }

    if (ldata > theR->dmax)
      { theR->dmax = 1.2*ldata + 1000;
        theR->data = (uint8 *) Realloc(theR->data,theR->dmax,"Reallocating data buffer");
        if (theR->data == NULL)
          exit (1);
      }

    bam_get(sf,theR->data,ldata);

    if ((getint(x+18,2) & 0x900) != 0)
      { theR->len = 0;
        return (0);
      }

    if (lseq <= 0)
      fprintf(stderr,"%s: WARNING: no sequence for subread !?\n",Prog_Name);

    theR->header = (char *) theR->data;
    theR->hlen   = lname;
    theR->len    = lseq;

    { uint8 *t;
      char  *s;
      int    i, e;

      t = theR->data + (lname + (lcigar<<2));
      s = theR->seq;
      lseq -= 1;
      for (e = i = 0; i < lseq; e++)
        { s[i++] = INT_2_IUPAC[t[e] >> 4];
          s[i++] = INT_2_IUPAC[t[e] & 0xf];
        }
      if (i <= lseq)
        s[i] = INT_2_IUPAC[t[e++] >> 4];
      lseq += 1;

      t += e;
      s  = theR->qvs;
      if (t[0] == 0xff)
        return (0);
      for (i = 0; i < lseq; i++)
        s[i] = t[i] + 33;
    }
  }

  return (1);
}

  //  Scan next bam entry and load PacBio info in record 'theR'

static char  IUPAC_2_DNA[256] =
  { 'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',   'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
    'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',   'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
    'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',   'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',
    'a', 'c', 'g', 't', 'a', 'a', 'a', 'a',   'a', 'a', 'a', 'a', 'a', 'a', 'a', 'a',

    'a', 'a', 'c', 'c', 'a', 'a', 'a', 'g',   'a', 'a', 'a', 'g', 'a', 'a', 'a', 'a',
    'a', 'a', 'a', 'c', 't', 'a', 'a', 'a',   'a', 'c', 'a', 'a', 'a', 'a', 'a', 'a',
    'a', 'a', 'c', 'c', 'a', 'a', 'a', 'g',   'a', 'a', 'a', 'g', 'a', 'a', 'a', 'a',
    'a', 'a', 'a', 'c', 't', 'a', 'a', 'a',   'a', 'c', 'a', 'a', 'a', 'a', 'a', 'a',
  };

#define CHECK(cond, msg)                                \
{ if ((cond))                                           \
    { fprintf(stderr, "%s: %s\n", Prog_Name, msg);      \
       exit (1); 	                                \
    }                                                   \
}

#define NEXT_ITEM(b,e)                                  \
{ b = e;                                                \
  while (*e != '\0' && *e != '\t')                      \
    e++;                                                \
  CHECK( *e == '\0', "Missing one or more fields")      \
  *e = 0;                                               \
}

static int sam_record_scan(BAM_FILE *sf, samRecord *theR)
{ char      *p;
  int        qlen, flags;

  //  read next line

  theR->data = sam_getline(sf);

  p = (char *) theR->data;

  { char *q, *seq;     //  Load header and sequence from required fields
    int   i;

    NEXT_ITEM(q,p)
    qlen = p-q;
    CHECK( qlen <= 1, "Empty header name")
    CHECK( qlen > 255, "Header is too long")

    theR->header = q;

    p = index(p+1,'\t');
    flags = strtol(q=p+1,&p,0);
    CHECK( p == q, "Cannot parse flags")

    for (i = 0; i < 7; i++)   // Skip next 8 required fields
      { p = index(p+1,'\t');
        CHECK( p == NULL, "Too few required fields in SAM record, file corrupted?")
      }
    p += 1;

    NEXT_ITEM(q,p)
    qlen = p-q;
    CHECK (*q == '*', "No sequence for read?");

    if (qlen > theR->lmax)
      { theR->lmax = 1.2*qlen + 1000;
        theR->seq  = (char *) Realloc(theR->seq,3*theR->lmax,"Reallocating sequence buffer");
        if (theR->seq == NULL)
          exit (1);
        theR->arr  = theR->seq + theR->lmax;
        theR->qvs  = theR->arr + theR->lmax;
      }

    if ((flags & 0x900) != 0)
      { theR->len = 0;
        return (0);
      }

    if (qlen <= 0)
      fprintf(stderr,"%s: WARNING: no sequence for subread !?\n",Prog_Name);

    theR->len = qlen;
    seq = theR->seq;
    for (i = 0; i < qlen; i++)
      seq[i] = IUPAC_2_DNA[(int) (*q++)];

    q = ++p;
    p = index(p,'\t');
    CHECK( p == NULL, "No auxilliary tags in SAM record, file corrupted?")

    if (*q == '*')
      return (0);
    qlen = p-q;
    seq = theR->qvs;
    for (i = 0; i < qlen; i++)
      seq[i] = *q++;
  }

  return (1);
}

/*******************************************************************************************
 *
 *  Parallel:  Each thread processes a contiguous stripe across the input files
 *               sending the compressed binary data lines to their assigned OneFile.
 *
 ********************************************************************************************/

  //  Write subread data in samRecord rec to non-NULL file types

static void *bam_output_thread(void *arg)
{ Thread_Arg  *parm  = (Thread_Arg *) arg;
  File_Object *fobj  = parm->fobj;
  uint8       *buf   = parm->buf;
  OneFile     *vf    = parm->vf;

  samRecord    _theR, *theR = &_theR;
  BAM_FILE     _bam, *bam = &_bam;

  int64        epos;
  uint32       eoff;
  int          isbam;
  int          f, fid;
  int          hasQV;
  char        *last;
  int          lmax, extra, skip1;
  int          first;

  //  Know the max size of sequence and data from pass 1, so set up accordingly

  if (fobj->ftype == BAM)
    { theR->dmax = 50000;
      theR->data = Malloc(theR->dmax,"Allocating sequence array");
      if (theR->data == NULL)
        exit (1);
    }
  theR->lmax = 75000;
  theR->seq  = Malloc(3*theR->lmax,"Allocating sequence array");
  if (theR->seq == NULL)
    exit (1);
  theR->arr = theR->seq + theR->lmax;
  theR->qvs = theR->arr + theR->lmax;

  bam->decomp = parm->decomp;

  if (GROUP)
    { lmax = 1000;
      last = Malloc(lmax+1,"Allocating lane buffer");
      last[0] = '\0';
    }

  first = 2;
  for (f = parm->bidx; f <= parm->eidx; f++)
    { fid   = open(fobj[f].path,O_RDONLY);
      isbam = (fobj[f].ftype == BAM);
      if (f < parm->eidx || parm->end.fpos >= fobj[f].fsize)
        { epos  = fobj[f].fsize;
          eoff  = 0;
          extra = 0;
        }
      else 
        { epos  = parm->end.fpos;
          eoff  = parm->end.boff;
          extra = GROUP;
        }
      if (f > parm->bidx || parm->beg.fpos == 0)
        { parm->beg.fpos = 0;
          parm->fid      = fid;
          if (isbam)
            skip_bam_header(parm);
          else
            sam_nearest(parm);
          skip1 = 0;
        }
      else
        skip1 = GROUP;

      if (isbam)
        bam_start(bam,fid,buf,&(parm->beg));
      else
        sam_start(bam,fid,buf,&(parm->beg));

#ifdef DEBUG_OUT
      printf("Block: %12lld / %5d to %12lld / %5d --> %8lld\n",bam->loc.fpos,bam->loc.boff,
                                                               epos,eoff,epos - bam->loc.fpos);
      fflush(stdout);
#endif

      while (1)
        { if (bam->loc.fpos >= epos && bam->loc.boff >= eoff)
            { if (extra)
                extra = 0;
              else
                break;
            }

          if (isbam)
            hasQV = bam_record_scan(bam,theR);
          else
            hasQV = sam_record_scan(bam,theR);

          if (theR->len <= 0)
            continue;

#ifdef DEBUG_BAM_RECORD
          fprintf(stderr,"S = '%s'\n",theR->seq);
          if (hasQV)
            fprintf(stderr,"Q = '%.*s'\n",theR->len,theR->qvs);
#endif

          if (GROUP)
            { char *s;
              int   i, glen;

              s = theR->header-1;
              for (i = 0; i < GROUP_REP; i++)
                { s = index(s+1,GROUP_CHAR);
                  if (s == NULL)
                    { s = theR->header + strlen(theR->header);
                      break;
                    }
                }
#ifdef DEBUG_AUTO
                  fprintf(stderr,"  group = %s\n",theR->header);
#endif
              glen = s-theR->header;
              if (strncmp(theR->header,last,glen) != 0)
                { if (glen >= lmax)
                    { lmax = lmax*1.2 + 1000;
                      last = (char *) Realloc(last,lmax+1,"Reallocating lane name");
                      if (last == NULL)
                        exit (1);
                    }
                  strncpy(last,theR->header,glen);
                  last[glen] = '\0';

                  if (skip1)
                    { skip1 = 0;
                      continue;
                    }

                  oneInt(vf,0) = 0;
                  oneInt(vf,1) = glen;
                  oneWriteLine(vf,'g',glen,theR->header);
                }
            }

          if (PAIRING)
            { if ((theR->flags & 0xc0) == 0x40)
                { oneWriteLine(vf,'P',0,NULL);
                  if (first == 0)
                    { parm->error = 1;
                      return (NULL);
                    }
                  first = 0;
                }
              else
                { if (first == 1)
                    { parm->error = 1;
                      return (NULL);
                    }
                  first = 1;
                }
            }

          oneInt(vf,0) = theR->len;
          oneWriteLine(vf,'S',theR->len,theR->seq);
    
          if (hasQV && QVS_OUT)
            { oneInt(vf,0) = theR->len;
              oneWriteLine(vf,'Q',theR->len,theR->qvs);
            }

          if (QNAME)
            { oneInt(vf,0) = theR->hlen;
              oneWriteLine(vf,'I',theR->hlen,theR->header);
            }
        }

      close(fid);
    }

  if (GROUP)
    free(last);
  free(theR->seq);
  if (fobj->ftype == BAM)
    free(theR->data);

  parm->error = 0;
  return (NULL);
}


/*******************************************************************************************
 *
 *  CRAM SPECIFIC CODE
 *
 ********************************************************************************************/

/*******************************************************************************************
 *
 *  cram_nearest: Routine to find next entry given arbitray start point.
 *
 ********************************************************************************************/

#include "htslib/hts.h"
#include "htslib/hfile.h"
#include "cram/cram.h"

int ITF_LEN[16] = { 0, 0, 0, 0,
                    0, 0, 0, 0,
                    1, 1, 1, 1,
                    2, 2, 3, 4 };

void scan_itf8(hFILE *fp)
{ char buf[4];
  int  nb = ITF_LEN[hgetc(fp)>>4];
  if (nb > 0)
    hread(fp,buf,nb);
}

const int LTF_LEN[256] = {
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,

    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
    0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,

    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
    1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,

    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,  2, 2, 2, 2,
    3, 3, 3, 3,  3, 3, 3, 3,  3, 3, 3, 3,  3, 3, 3, 3,
    4, 4, 4, 4,  4, 4, 4, 4,  5, 5, 5, 5,  6, 6, 7, 8 };

void scan_ltf8(hFILE *fp)
{ char buf[8];
  int  nb = LTF_LEN[hgetc(fp)];
  if (nb > 0)
    hread(fp,buf,nb);
}

int int32_decode(hFILE *fp, int32 *val)
{ char *v = (char *) val;

  if (hread(fp,v,4) != 4)
    return (1);
#if __ORDER_LITTLE_ENDIAN__ != __BYTE_ORDER__
  { char  x;
    x = v[0];
    v[0] = v[3];
    v[3] = x;
    x = v[1];
    v[1] = v[2];
    v[2] = x;
  }
#endif
  return (0);
}

// reading cram block, header is a block so wrapped.

int64 scan_container(cram_fd *fd)
{ int     i, len;
  int32   nslice;
  hFILE  *fp = fd->fp;

  if (int32_decode(fp,&len))                    //  total length
    return (-1);
  scan_itf8(fp);                                //  ref seq id
  scan_itf8(fp);                                //  start pos on ref
  scan_itf8(fp);                                //  align span
  scan_itf8(fp);                                //  # of records
  if (CRAM_MAJOR_VERS(fd->version) > 1)
    { if (CRAM_MAJOR_VERS(fd->version) >= 3)    //  record counter
        scan_ltf8(fp);
      else
        scan_itf8(fp);
      scan_ltf8(fp);                            //  # bps
    }
  scan_itf8(fp);                                //  # of blocks
  itf8_decode(fd,&nslice);                      //  # of slices
  for (i = 0; i < nslice; i++)                  //  landmarks
    scan_itf8(fp);
  hread(fp,&nslice,4);                          //  crc code
  hseek(fp,len,SEEK_CUR);                       //  skip contents of container
  return (htell(fp));
}

int64 *genes_cram_index(char *path, int64 fsize, int64 *zsize)
{ cram_fd *fd;
  int64    cash[10], *zoff;
  int64    s, e;
  int      i, j;
 
  fd = cram_open(path,"r");

  cash[0] = htell(fd->fp);
  for (i = 1; i < 10; i++)
    { cash[i] = s = scan_container(fd);
      if (s == fsize)
        { zoff = Malloc(sizeof(int64)*i,"Allocating cram index"); 
          if (zoff == NULL)
            exit (1);
          for (j = 0; j < i; i++)
            zoff[j] = cash[i];
          *zsize = i-1;
          return (zoff);
        }
    }

  e = ((fsize-(cash[0]+38))/(cash[9]-cash[0])+1)*9 + 100;
  zoff = Malloc(sizeof(int64)*e,"Allocating cram index"); 
  if (zoff == NULL)
    exit (1);
  for (j = 0; j < i; j++)
    zoff[j] = cash[j];

  while (s != fsize)
    { zoff[i++] = s = scan_container(fd);
      if (i >= e)
        { e = ((fsize-(zoff[0]+38.))/(zoff[i-1]-zoff[0]))*(i-1) + 100;
          zoff = Realloc(zoff,sizeof(int64)*e,"Allocating cram index"); 
          if (zoff == NULL)
            exit (1);
        }
    }

  cram_close(fd);

  *zsize = i-2;
  return (zoff);
}       

void cram_nearest(Thread_Arg *data)
{ File_Object *inp   = data->fobj + data->bidx;
  int64       *zoffs = inp->zoffs;
  int64       zsize  = inp->zsize;
  int i;

  for (i = 0; i < zsize; i++)
    if (zoffs[i] >= data->beg.fpos)
      break;
  if (i == zsize)
    data->beg.fpos = -1;
  else
    data->beg.fpos = zoffs[i];
}

static void *cram_output_thread(void *arg)
{ Thread_Arg  *parm  = (Thread_Arg *) arg;
  File_Object *fobj  = parm->fobj;
  OneFile     *vf    = parm->vf;

  File_Object *inp;
  int          f;
  cram_fd     *fid;
  int64        bpos, epos;

  char        *last;
  int          lmax;
  int          extra, skip1;
  int          first;

  if (GROUP)
    { lmax = 1000;
      last = Malloc(lmax+1,"Allocating lane buffer");
      last[0] = '\0';
    }

  first = 2;
  for (f = parm->bidx; f <= parm->eidx; f++)
    { inp = fobj+f;
      fid = cram_open(inp->path,"r");
      if (f < parm->eidx || parm->end.fpos >= inp->zoffs[inp->zsize])
        { epos  = inp->zoffs[inp->zsize];
          extra = 0;
        }
      else
        { epos  = parm->end.fpos;
          extra = GROUP;
        }
      if (f > parm->bidx || parm->beg.fpos < inp->zoffs[0])
        { bpos  = inp->zoffs[0];
          skip1 = 0;
        }
      else
        { bpos  = parm->beg.fpos;
          skip1 = GROUP;
        }
      hseek(fid->fp,bpos,SEEK_SET);

#ifdef DEBUG_OUT
      fprintf(stderr,"Block: %12lld to %12lld --> %8lld\n",bpos,epos,epos - bpos);
      fflush(stderr);
#endif

      while (1)
        { cram_record *rec;
          uint8       *qual;

          rec = cram_get_seq(fid);
          if (rec == NULL)
            break;

          if (htell(fid->fp) > epos)
            { if (extra)
                extra = 0;
              else
                break;
            }

          if (GROUP)
            { char *s, *h;
              int   i, glen;

              h = (char *) rec->s->name_blk->data + rec->name;
              s = h-1;
              for (i = 0; i < GROUP_REP; i++)
                { s = index(s+1,GROUP_CHAR);
                  if (s == NULL)
                    { s = h + rec->name_len;
                      break;
                    }
                }
#ifdef DEBUG_AUTO
              fprintf(stderr,"  group = %s\n",h);
#endif
              glen = s-h;
              if (strncmp(h,last,s-h) != 0)
                { if (glen >= lmax)
                    { lmax = glen*1.2 + 1000;
                      last = (char *) Realloc(last,lmax+1,"Reallocating lane name");
                      if (last == NULL)
                        exit (1);
                    }
                  strncpy(last,h,glen);
                  last[glen] = '\0';

                  if (skip1)
                    { skip1 = 0;
                      continue;
                    }

                  oneInt(vf,0) = 0;
                  oneInt(vf,1) = glen;
                  oneWriteLine(vf,'g',glen,h);
                }
            }

          if (PAIRING)
            { if ((rec->flags & 0xc0) == 0x40)
                { oneWriteLine(vf,'P',0,NULL);
                  if (first == 0)
                    { parm->error = 1;
                      return (NULL);
                    }
                  first = 0;
                }
              else
                { if (first == 1)
                    { parm->error = 1;
                      return (NULL);
                    }
                  first = 1;
                }
            }

          oneInt(vf,0) = rec->len;
          oneWriteLine(vf,'S',rec->len,rec->s->seqs_blk->data+rec->seq);

          qual = rec->s->qual_blk->data+rec->qual;
          if (QVS_OUT && qual[0] != 0xff)
            { int i;
              for (i = 0; i < rec->len; i++)
                qual[i] += 33;
              oneInt(vf,0) = rec->len;
              oneWriteLine(vf,'Q',rec->len,qual);
            }

          if (QNAME)
            { oneInt(vf,0) = rec->name_len;
              oneWriteLine(vf,'I',rec->name_len,rec->s->name_blk->data+rec->name);
            }
        }
    }

  if (GROUP)
    free(last);

  parm->error = 0;
  return (NULL);
}


/****************************************************************************************
 *
 *  The top-level program
 *
 ****************************************************************************************/

int main(int argc, char *argv[])
{ OneSchema  *schema;
  char       *command;
  int         nfiles;
  int         ftype;
  int         need_decon;
  int         need_zuf;
  void *    (*output_thread)(void *);
  void      (*scan_header)(Thread_Arg *);
  void      (*find_nearest)(Thread_Arg *);

  //  Capture command line for provenance

  { int   n, i;
    char *c;

    n = 0;
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

    schema = oneSchemaCreateFromText(vgpSchemaText);
  }

  //  Parse command line options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;
    int    one = 1;

    ARG_INIT("VGPseq")

    NTHREADS = 4;
    GROUP    = 0;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("viqp")
            break;
          case 'g':
            GROUP_REP = strtol(argv[i]+2,&eptr,10);
            if (eptr <= argv[i]+2)
              { fprintf(stderr,"%s: -g group repetition count not present\n",Prog_Name);
                exit (1);
              }
            if (GROUP_REP <= 0)
              { fprintf(stderr,"%s: -g group repetition count must be positive\n",Prog_Name);
                exit (1);
              }
            k = *eptr++;
            if (k == '\\')
              { k = *eptr++;
                switch (k)
                { case 'f':
                    GROUP_CHAR = 0x0c;
                    break; 
                  case 'n':
                    GROUP_CHAR = 0x0a;
                    break; 
                  case 'r':
                    GROUP_CHAR = 0x0d;
                    break; 
                  case 't':
                    GROUP_CHAR = 0x09;
                    break; 
                  default:
                    GROUP_CHAR = k;
                    break; 
                }
              }
            else
              GROUP_CHAR = k;
            GROUP = 1;
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE  = flags['v'];
    QNAME    = flags['i'];
    QVS_OUT  = flags['q'];
    PAIRING  = flags['p'];

    if (argc < 2)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: verbose mode, output progress as proceed\n");
        fprintf(stderr,"      -i: Output identify\n");
        fprintf(stderr,"      -q: Output QV string\n");
        fprintf(stderr,"      -g: Output group where names = identifier prefix\n");
        fprintf(stderr,"                   to #'th instance of character x\n");
        fprintf(stderr,"      -p: If pairing information, then produce .irp\n");
        fprintf(stderr,"      -T: Number of threads to use\n");
        exit (1);
      }

    IS_BIG_ENDIAN = ( *((char *) (&one)) == 0);

    nfiles = argc-1;
  }

  //  Find partition points dividing data in all files into NTHREADS roughly equal parts
  //    and then in parallel threads produce the output for each part.

  { File_Object fobj[nfiles];
    Thread_Arg  parm[NTHREADS];

#ifndef DEBUG_OUT
    pthread_t   threads[NTHREADS];
#endif

    { uint8        *bf;
      int           f, i;
      int64         b, work, wper;

      //  Get name and size of each file in 'fobj[]', determine type, etc.

      need_decon = 0;
      need_zuf   = 0;

      work = 0;
      for (f = 0; f < nfiles; f++)
        { Fetch_File(argv[f+1],fobj+f);

          if (f > 0)
            { if (fobj[f].ftype != ftype)
                { fprintf(stderr,"%s: All files must be of the same type\n",Prog_Name);
                  exit (1);
                }
            }
          else
            { ftype = fobj[0].ftype;
              if (ftype == FASTA || ftype == FASTQ)
                { output_thread = fast_output_thread;
                  scan_header   = do_nothing;
                  find_nearest  = fast_nearest;
                }
              else if (ftype == BAM)
                { output_thread = bam_output_thread;
                  scan_header   = skip_bam_header;
                  find_nearest  = bam_nearest;
                }
              else if (ftype == SAM)
                { output_thread = bam_output_thread;
                  scan_header   = sam_nearest;
                  find_nearest  = sam_nearest;
                }
              else // ftype == CRAM
                { output_thread = cram_output_thread;
                  scan_header   = do_nothing;
                  find_nearest  = cram_nearest;
                }
            }
          if (ftype == BAM)
            need_decon = 1;
          else if (fobj[f].zipd)
            { need_decon = 1;
              need_zuf   = 1;
            }

          work += fobj[f].fsize;
        }

      //  Allocate IO buffer space for threads

      if (need_zuf)
        bf = Malloc(2*NTHREADS*IO_BLOCK,"Allocating IO_Buffer\n");
      else
        bf = Malloc(NTHREADS*IO_BLOCK,"Allocating IO_Buffer\n");
      if (bf == NULL)
        exit (1);

      if (VERBOSE)
        { if (nfiles > 1)
            fprintf(stderr,"  Partitioning %d %s%s files into %d parts\n",
                           nfiles,fobj->zipd?"compressed ":"",Tstring[fobj->ftype],NTHREADS);
          else
            fprintf(stderr,"  Partitioning %d %s%s file into %d parts\n",
                           nfiles,fobj->zipd?"compressed ":"",Tstring[fobj->ftype],NTHREADS);
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
          parm[i].buf    = bf + i*IO_BLOCK;
          if (need_decon)
            parm[i].decomp = libdeflate_alloc_decompressor();
          if (need_zuf)
            parm[i].zuf = bf + (NTHREADS+i)*IO_BLOCK;

          if (b != 0)
            { parm[i].fid = open(fobj[f].path,O_RDONLY);
              parm[i].beg.fpos = 0;
              parm[i].beg.boff = 0;
              scan_header(parm+i);
              parm[i].end = parm[i].beg;
            }
          else
            parm[i].end.fpos = parm[i].end.boff = 0;

#ifdef DEBUG_FIND
          fprintf(stderr," %2d: %1d %10lld (%10lld / %5d)\n",
                         i,f,b,parm[i].end.fpos,parm[i].end.boff);
          fflush(stderr);
#endif

          parm[i].beg.fpos = b;
          parm[i].bidx = f;

          work -= wper;
          while (work < .01*IO_BLOCK)
            { if (f == nfiles)
                { if (VERBOSE && NTHREADS != i+1)
                    { if (i > 0)
                        fprintf(stderr,"  File is so small dividing into only %d parts\n",i+1);
                      else
                        fprintf(stderr,"  File is so small will not divide it\n");
                    }
                  NTHREADS = i+1;
                  break;
                }
              work += fobj[++f].fsize;
            }
          b = fobj[f].fsize - work;
          if (b < 0)
            { work += b;
              b = 0;
              if (f == nfiles)
                { if (VERBOSE && NTHREADS != i+1)
                    { if (i > 0)
                        fprintf(stderr,"  File is so small dividing into only %d parts\n",i+1);
                      else
                        fprintf(stderr,"  File is so small will not divide it\n");
                    }
                  NTHREADS = i+1;
                  break;
                }
            }
        }
    }

    { int i, f;

      //  For each non-zero start point find synchronization point in
      //    sequence file.  If can't find it then start at beginning of
      //    next file, and if at or before first data line then signal
      //    start at beginning by zero'ing the synch point.

      for (i = 0; i < NTHREADS; i++)
        { if (parm[i].beg.fpos != 0)
           { find_nearest(parm+i);
              if (parm[i].beg.fpos < 0)
                { parm[i].beg.fpos = 0;
                  parm[i].bidx += 1;
                }
              else if (parm[i].beg.fpos <= parm[i].end.fpos)
                parm[i].beg.fpos = 0;
              close(parm[i].fid);
            }
          if (parm[i].beg.fpos == 0 && ftype >= FASTQ)
            parm[i].alast = Strdup("","Allocating empty string");;
        }

      //  Paranoid: if one thread's synch point overtakes the next one (will almost
      //    certainly never happen unless files very small and threads very large),
      //    remove the redundant threads.

      f = 0;
      for (i = 1; i < NTHREADS; i++)
        if (parm[i].bidx > parm[f].bidx || parm[i].beg.fpos > parm[f].beg.fpos)
          parm[++f] = parm[i];
        else
          { if (ftype >= FASTQ)
              free(parm[i].alast);
            if (need_decon)
              libdeflate_free_decompressor(parm[i].decomp);
          }
      if (VERBOSE && NTHREADS != f+1)
        { if (f > 0)
            fprintf(stderr,"  File is so small will divide into only %d parts\n",f+1);
          else
            fprintf(stderr,"  File is so small will not divide it\n");
        }
      NTHREADS = f+1;

      //  Develop end points of each threads work using the start point of the next thread

      for (i = 1; i < NTHREADS; i++)
        if (parm[i].beg.fpos == 0)
          { parm[i-1].end.fpos = fobj[parm[i].bidx-1].fsize;
            parm[i-1].end.boff = 0;
            parm[i-1].eidx     = parm[i].bidx-1;
          }
        else
          { parm[i-1].end.fpos = parm[i].beg.fpos;
            parm[i-1].end.boff = parm[i].beg.boff;
            parm[i-1].eidx     = parm[i].bidx;
          }
      parm[NTHREADS-1].end.fpos = fobj[argc-2].fsize;
      parm[NTHREADS-1].end.boff = 0;
      parm[NTHREADS-1].eidx     = argc-2;

#if defined(DEBUG_FIND) || defined(DEBUG_OUT)
      fprintf(stderr,"\nPartition:\n");
      for (i = 0; i < NTHREADS; i++)
        { fprintf(stderr," %2d: %2d / %12lld / %5d",
                         i,parm[i].bidx,parm[i].beg.fpos,parm[i].beg.boff);
          fprintf(stderr,"  -  %2d / %12lld / %5d\n",
                           parm[i].eidx,parm[i].end.fpos,parm[i].end.boff);
        }
      fflush(stderr);
#endif
    }

    //  Produce output in parallel threads based on partition

    { OneFile *vf;
      int      i, error;

      if (PAIRING)
        vf = oneFileOpenWriteNew("-",schema,"irp",true,NTHREADS);
      else
        vf = oneFileOpenWriteNew("-",schema,"seq",true,NTHREADS);
      oneAddProvenance(vf,Prog_Name,"1.0",command,NULL);
      oneWriteHeader(vf);

#ifdef DEBUG_OUT
      fprintf(stderr,"Opened\n");
      fflush(stderr);
#endif

      if (VERBOSE)
        { if (NTHREADS > 1)
            fprintf(stderr,"  Producing .seq segments in parallel\n");
          else
            fprintf(stderr,"  Producing .seq segments\n");
          fflush(stderr);
        }

      //  Generate the data lines in parallel threads

#ifdef DEBUG_OUT
      for (i = 0; i < NTHREADS; i++)
        { parm[i].vf = vf+i;
          fprintf(stderr,"Thread %d\n",i);
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

      //  If asked for paired reads make the files support it

      error = 0;
      for (i = 0; i < NTHREADS; i++)
        error |= parm[i].error;

      if (error)
        { fprintf(stderr,"%s: Input file(s) are not properly sorted for pairing\n",Prog_Name);
          exit (1);
        }

      if (VERBOSE)
        { fprintf(stderr,"  Cat'ing .seq segments\n");
          fflush(stderr);
        }

      oneFileClose(vf);
    }

    //  Free everything as a matter of good form

    { int i, f;

      if (need_decon)
        for (i = 0; i < NTHREADS; i++)
          libdeflate_free_decompressor(parm[i].decomp);
      free(parm[0].buf);
      for (f = 0; f < nfiles; f++)
        Free_File(fobj+f);
      free(command);
    }
  }

  oneSchemaDestroy(schema);

  if (VERBOSE)
    { fprintf(stderr,"  Done\n");
      fflush(stderr);
    }

  exit (0);
}
