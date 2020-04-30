#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include <fcntl.h>
#include <dirent.h>
#include <sys/mman.h>
#include <pthread.h>

#include "gene_core.h"
#include "lsd.sort.h"
#include "exsort.h"
#include "../Core/VGPlib.h"

#undef DEBUG

int    VERBOSE;    //  Verbose mode?
char  *SORT_PATH;  //  Directory to do external sort in
int    NTHREADS;   //  # of threads to use for parallized sorts

#define VALID_THRESH 100

static char *Usage = "[-v] [-P<dir(/tmp)>] [-T<int(4)>] <clouds:irp>";

static uint8 bit[8] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };


/****************************************************************************************
 *
 *  Fixed length compression / decompression for Q-strings & decompression for DNA S-strings
 *
 *****************************************************************************************/

static void Compress_QV(int len, uint8 *qvec, int qbits, uint8 *map)
{ uint8 *buf, x;
  int    i, rem;

  if (qbits == 8)
    return;

  buf = qvec;
  *buf = map[qvec[0]];
  rem = qbits;
  for (i = 1; i < len; i++)
    { x = map[qvec[i]];
      *buf |= (x << rem);
      rem += qbits;
      if (rem >= 8)
        { rem -= 8;
          if (rem > 0)
            *++buf = (x >> (qbits-rem));
          else
            *++buf = 0;
        }
    }
}

static void Uncompress_QV(int len, uint8 *qvec, int qbits, uint8 qmask, uint8 *inv)
{ uint8 *buf, x, v, bmask;
  int    i, rem;

  bmask = (1 << qbits) - 1;
  buf = qvec + (len*qbits-1)/8;
  rem = (len*qbits-1) % 8 + 1;
  x   = *buf--;
  for (i = len-1; i >= 0; i--)
    { rem -= qbits;
      if (rem < 0)
        { v = (x << (-rem)) & qmask;
          rem += 8;
          x = *buf--;
          v |= (x >> rem);
          qvec[i] = inv[v];
        }
      else
        { v = (x >> rem) & qmask;
          qvec[i] = inv[v];
        }
    }
}

  //  Uncompress read from 2-bits per base into [0-3] per byte representation

static char Base[4] = { 'a', 'c', 'g', 't' };

static void Uncompress_SEQ(char *s, int len)
{ int   i, byte;
  char *t0, *t1, *t2, *t3;

  t0 = s;
  t1 = t0+1;
  t2 = t1+1;
  t3 = t2+1;

  i  = len/4;
  s += i;
  i = 4*i;
  switch (len-i)
  { case 3:
      byte = *s--;
      t0[i] = Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 4) & 0x3];
      t2[i] = Base[(byte >> 2) & 0x3];
      break;
    case 2:
      byte = *s--;
      t0[i] = Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 4) & 0x3];
      break;
    case 1:
      byte = *s--;
      t0[i] = Base[(byte >> 6) & 0x3];
      break;
    default:
      s -= 1;
      break;
  }

  for ( ; i > 0; i -= 4)
    { byte = *s--;
      t0[i] = Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 4) & 0x3];
      t2[i] = Base[(byte >> 2) & 0x3];
      t3[i] = Base[byte & 0x3];
    }
}


/****************************************************************************************
 *
 *  Thread for building unsorted bar code list and collecting Q-vector symbol usage
 *
 *****************************************************************************************/

//  Optimized ReadLine optimized to read binary S-line without decompression

static char read_raw_seq(VgpFile *vf)
{ U8        x;
  char      t;
  LineInfo *li;

  x = getc(vf->f);
  if (feof(vf->f) || x == '\n')
    return (0);

  if (x & 0x80)
    { t = vf->binaryTypeUnpack[x];
      if (t != 'S')
        return (0);
    }
  else
    return (0);
  vf->lineType = t;

  li = vf->lineInfo['S'];

  { int64 nBits;

    if (x & 0x1)
      { nBits = (U8) getc(vf->f);
        if (fread(vf->codecBuf,((nBits+7) >> 3),1,vf->f) != 1)
          die("fail to read compressed fields");
      }
    else
      { if (fread(vf->field,sizeof(Field),1,vf->f) != (unsigned long) 1)
          die("fail to read fields");
      }

    if (fread(&nBits,sizeof(I64),1,vf->f) != 1)
      die("fail to read list nBits");
    if (fread(vf->codecBuf,((nBits+7) >> 3),1,vf->f) != 1)
      die("fail to read compressed list");

    vf->field[0].i = (nBits >> 1);
  }

  { U8 peek = getc(vf->f);    // check if next line is a comment - if so then read it
    ungetc(peek,vf->f);
    if (peek & 0x80)
      peek = vf->binaryTypeUnpack[peek];
    if (peek == '/')
      { Field keepField0 = vf->field[0];
        vgpReadLine(vf);
        vf->lineType = t;
        vf->field[0] = keepField0;
      }
  }


  return (t);
}


typedef struct
  { VgpFile     *vf;       //  VgpFile for input
    int64        beg;      //  Range of reads to process
    int64        end;
    uint32      *count;
    int64        usedqvs[256];
    int          error;
    int          flen;
    int          rlen;
  } BarCode_Arg;

#define ERROR_PLINE 1
#define ERROR_SLINE 2
#define ERROR_SLEN  3
#define ERROR_QLINE 4
#define ERROR_QLEN  5

static void error_report(int error, int64 line)
{ printf("Error Called %d\n",error);
  switch (error)
  { case ERROR_PLINE:
      fprintf(stderr,"%s: Expecting P-line, line %lld\n",Prog_Name,line);
      break;
    case ERROR_SLINE:
      fprintf(stderr,"%s: Expecting S-line, line %lld\n",Prog_Name,line);
      break;
    case ERROR_SLEN:
      fprintf(stderr,"%s: S-strings are not all the same size, line %lld\n",Prog_Name,line);
      break;
    case ERROR_QLINE:
      fprintf(stderr,"%s: Expecting Q-line, line %lld\n",Prog_Name,line);
      break;
    case ERROR_QLEN:
      fprintf(stderr,"%s: Q-string is not the same length as S-string, line %lld\n",Prog_Name,line);
      break;
  }
}

static void *barcodes_thread(void *arg)
{ BarCode_Arg *parm  = (BarCode_Arg *) arg;
  VgpFile     *vf    = parm->vf;
  int64       *used  = parm->usedqvs;
  int64        beg   = 2*parm->beg;
  int64        end   = 2*parm->end;
  uint32      *count = parm->count;

  uint32 *fcmp;
  uint8  *fqvs;
  int64   o, coll;
  int     t, n;
  int     flen, rlen;

  fqvs = vf->lineInfo['Q']->buffer;
  fcmp = (uint32 *) (vf->codecBuf);

  bzero(used,sizeof(int64)*256);

  vgpGotoObject(vf,beg);

  coll  = beg + 10000000/NTHREADS;
  rlen  = 0;
  flen  = 0;
  for (o = beg; o < end; o += 2)
    { t = read_raw_seq(vf);
      if (t != 'S')
        { parm->error = ERROR_SLINE;
          return (NULL);
	}

      *count++ = *fcmp;

      if (o > coll)
        { vgpGotoObject(vf,o+2);
          continue;
        }

      n = vgpInt(vf,0);
      if (flen == 0)
        flen = n;
      else if (flen != n)
        { parm->error = ERROR_SLEN;
          return (NULL);
        }

      t = vgpReadLine(vf);
      if (t != 'Q')
        { parm->error = ERROR_QLINE;
          return (NULL);
        }
      if (vgpInt(vf,0) != flen)
        { parm->error = ERROR_QLEN;
          return (NULL);
        }
      for (n = 0; n < flen; n++)
        used[(int) fqvs[n]] += 1;

      t = read_raw_seq(vf);
      if (t != 'S')
        { parm->error = ERROR_SLINE;
          return (NULL);
        }
      n = vgpInt(vf,0);
      if (rlen == 0)
        rlen = n;
      else if (rlen != n)
        { parm->error = ERROR_SLEN;
          return (NULL);
        }

      t = vgpReadLine(vf);
      if (t != 'Q')
        { parm->error = ERROR_QLINE;
          return (NULL);
        }
      if (vgpInt(vf,0) != rlen)
        { parm->error = ERROR_QLEN;
          return (NULL);
        }
      for (n = 0; n < rlen; n++)
        used[(int) fqvs[n]] += 1;

      t = vgpReadLine(vf);
      if (t == 'g')
        t = vgpReadLine(vf);

      if (t != 'P')
        { parm->error = ERROR_PLINE;
          return (NULL);
        }
    }

  parm->error = 0;
  parm->flen  = flen;
  parm->rlen  = rlen;
  return (NULL);
}


/****************************************************************************************
 *
 *  4 thread passes for building list of 1-corrections
 *
 *****************************************************************************************/

typedef struct
  { int64        beg;      //  Range of reads to process
    int64        end;
    uint32      *count;
    uint8       *vector;
    uint8       *valid;
    int64        ngood;
    int64        ndist;
  } Correct_Arg;

static void *zero_thread(void *arg)
{ Correct_Arg *parm    = (Correct_Arg *) arg;
  int64         beg    = parm->beg;
  int64         end    = parm->end;
  uint8        *valid  = parm->valid;
  uint8        *vector = parm->vector;

  bzero(vector,0x20000000ll);
  bzero(valid+beg,end-beg);

  return (NULL);
}

static void *good_thread(void *arg)
{ Correct_Arg *parm  = (Correct_Arg *) arg;
  int64        beg   = parm->beg;
  int64        end   = parm->end;
  uint8       *valid = parm->valid;
  uint32      *count = parm->count;

  int64  g, f, i;
  uint32 v;
  int64  ngood, ndist;

  ngood = 0;
  ndist = 0;
  f = beg;
  g = beg;
  v = count[g];
  for (i = beg+1; i <= end; i++)
    if (count[i] != v)
      { g = i-g;
#ifdef DEBUG_GOOD
        printf(" %8lld: %4lld %08x",i,g,v);
#endif
        if (g > VALID_THRESH)
          { valid[v] = 0x1;
            ngood += g;
            ndist += 1;
#ifdef DEBUG_GOOD
            printf(" ++\n");
#endif
          }
#ifdef DEBUG_GOOD
        else
          printf("\n");
#endif
        count[f++] = v;
        g = i;
        v = count[g];
      }

  parm->end   = f; 
  parm->ngood = ngood;
  parm->ndist = ndist;

  return (NULL);
}

static void *mark_thread(void *arg)
{ Correct_Arg *parm  = (Correct_Arg *) arg;
  int64        beg   = parm->beg;
  int64        end   = parm->end;
  uint8       *valid = parm->valid;
  uint8       *vector = parm->vector;
  uint32      *count = parm->count;

  int64  i;
  int    j, m;
  uint32 v;
  uint32 var, v3, vb;
  uint32 chg, msk, u, t;

  // Precondition: valid[v] & 0x3 == 0 if not valid, 1 if valid, bit_t[v] = 0 for all t, v

  for (i = beg; i < end; i++)
    { v = count[i];
      if ((valid[v] & 0x3) == 0x1)
        { chg = 0x1;
          msk = 0x3;
#ifdef DEBUG_GOOD
          printf("  %08x\n",v);
#endif
          for (j = 0; j < 16; j++)
            { t = (v & msk); 
              v -= t;
              u = 0;
              for (m = 0; m < 4; m++)
                { if (t != u)
                    { var = v + u;
#ifdef DEBUG_MARK
                      printf("   %08x (%d/%d)",var,j,m);
#endif
                      v3 = (var >> 3);
                      vb = bit[var & 0x7];
                      if (vector[v3] & vb)
                        { if ((valid[var] & 0x1) == 0)
                            { valid[var] |= 0x3;
#ifdef DEBUG_MARK
                              printf(" V");
#endif
                            }
                        }
                      else
                        { vector[v3] |= vb;
#ifdef DEBUG_MARK
                          printf(" B");
#endif
                        }
#ifdef DEBUG_MARK
                      printf("\n");
#endif
                    }
                  u += chg;
                }
              chg <<= 2;
              msk <<= 2;
              v += t;
            }
        }
    }

  // valid[v] & 0x3 == 1 if valid, 3 if not valid & touched >= 2x by a thread, 0 otherwise
  // bit_t[v] = 1 iff thread t touched v at least once

  return (NULL);
}

static void *fix_thread(void *arg)
{ Correct_Arg *parm   = (Correct_Arg *) arg;
  int          beg    = parm->beg;
  int          end    = parm->end;
  uint8       *valid  = parm->valid;
  uint8       *vector = parm->vector;
  uint32      *count  = parm->count;

  uint8 *vec[NTHREADS];

  int64  i, ngood;
  int    j, m, n, cnt;
  uint32 v, var, v3, vb;
  uint32 chg, msk, u, t, c;

  for (n = 0; n < NTHREADS; n++)
    vec[n] = vector + n*0x20000000ll;

  // valid[v] & 0x3 == 1 if valid, 3 if not valid & touched >=2x by a thread, 0 otherwise
  // bit_t[v] = 1 iff thread t touched v at least once

  ngood = 0;
  for (i = beg; i < end; i++)
    { v = count[i];
      if ((valid[v] & 0x3) == 0x1)
        { chg = 0x1;
          msk = 0x3;
#ifdef DEBUG_FIXES
          printf("%08x\n",v);
#endif
          for (j = 0; j < 16; j++)
            { t = (v & msk); 
              v -= t;
              c = (t >> (2*j));
              u = 0;
              for (m = 0; m < 4; m++)
                { if (t != u)
                    { var = v + u;
#ifdef DEBUG_FIXES
                      printf("   %08x (%d/%d)\n",var,j,m);
#endif
                      v3 = (var >> 3);
                      vb = bit[var & 0x7];
                      if ((valid[var] & 0x3) == 0)
                        { cnt = 0;
                          for (n = 0; n < NTHREADS; n++)
                            if (vec[n][v3] & vb)
                              cnt += 1;
                          if (cnt == 1)
                            { ngood += 1;
                              valid[var] = (0x2 | (c << 2) | (j << 4));
#ifdef DEBUG_FIXES
                              printf("      OK %x  %d %d\n",valid[var],c,j);
#endif
                            }
                          else
                            { valid[var] |= 0x3;
#ifdef DEBUG_FIXES
                              printf("      C-%d\n",cnt);
#endif
                            }
                        }
#ifdef DEBUG_FIXES
                      else
                        printf("      Z %x\n",valid[var]&0x3);
#endif
                    }
                  u += chg;
                }
              chg <<= 2;
              msk <<= 2;
              v += t;
            }
        }
    }

  // valid[v] & 0x3 == 0  not validd & never touched
  //                   1  valid
  //                   2  not valid & touched once => uniquely correctable (with correction)
  //                   3  not valid & touched 2 or more times

  parm->ngood = ngood;

  return (NULL);
}


/****************************************************************************************
 *
 *  4 thread passes for building list of 1-corrections
 *
 *****************************************************************************************/

typedef struct
  { int          tidx;     //  Thread index [0,NTHREADS)
    VgpFile     *vf;       //  VgpFile for input
    int64        beg;      //  Range of reads to process
    int64        end;
    int          qbits;
    uint8       *map;
    int          flen, fclen, fqlen;
    int          rlen, rclen, rqlen;
    uint8       *valid;
    uint8       *array;
    int64        nused;
    int64        size[256];
  } Load_Arg;

static pthread_mutex_t Allocator = PTHREAD_MUTEX_INITIALIZER;
static int64           Alloc_Offset;
static int64           Block_Size;

static void *load_thread(void *arg)
{ Load_Arg *parm  = (Load_Arg *) arg;
  VgpFile  *vf    = parm->vf;
  int       tidx  = parm->tidx;
  int64     beg   = 2*parm->beg;
  int64     end   = 2*parm->end;
  uint8    *valid = parm->valid;
  int       qbits = parm->qbits;
  uint8    *map   = parm->map;
  int       flen  = parm->flen;
  int       rlen  = parm->rlen;
  int       fclen = parm->fclen;
  int       fqlen = parm->fqlen;
  int       rclen = parm->rclen;
  int       rqlen = parm->rqlen;
  uint8    *array = parm->array;
  int64    *size  = parm->size;

  uint8 *bptr[256];
  uint8 *bend[256];

  uint32 *dcode;
  uint8  *bcode, *fqvs;
  int64   o;
  uint32  v, c, bar, byte;
  int64   nused;
  uint8  *aptr, *aend;
  int     t;

  fqvs  = vf->lineInfo['Q']->buffer;
  bcode = (uint8 *)  (vf->codecBuf);
  dcode = (uint32 *) (vf->codecBuf);

  aptr = array + tidx*256*Block_Size;
  for (t = 0; t < 256; t++)
    { bptr[t] = aptr;
      aptr   += Block_Size;
      bend[t] = aptr - sizeof(uint8 *);
      size[t] = 0;
    }

  vgpGotoObject(vf,beg);

  nused = 0;
  for (o = beg; o < end; o += 2)
    { read_raw_seq(vf);

      bar = *dcode;

      v   = valid[bar];
      c   = v&0x3;
      if (c == 2)
        { int i = (v>>3) & 0x1e;
          int x = (v>>2) & 0x3;
          bar -= (bar & (0x3 << i));
          bar |= (x << i);
          *dcode = bar;
          c = 1;
        }
      if (c != 1)
        { vgpGotoObject(vf,o+2);
          continue;
        }

      byte = *bcode;
      aptr = bptr[byte];
      aend = bend[byte];
      if (aptr >= aend)
        { pthread_mutex_lock(&Allocator);
            aptr          = array + Alloc_Offset;
            Alloc_Offset += Block_Size;
          pthread_mutex_unlock(&Allocator);
          *((uint8 **) aend) = aptr; 
          bend[byte] = aptr + (Block_Size-sizeof(uint8 *));
        }

      memcpy(aptr,dcode,fclen);
      aptr += fclen;

      vgpReadLine(vf);
      Compress_QV(flen,fqvs,qbits,map);
      memcpy(aptr,fqvs,fqlen);
      aptr += fqlen;

      read_raw_seq(vf);
      memcpy(aptr,dcode,rclen);
      aptr += rclen;

      vgpReadLine(vf);
      Compress_QV(rlen,fqvs,qbits,map);
      memcpy(aptr,fqvs,rqlen);
      aptr += rclen;

      t = vgpReadLine(vf);
      if (t == 'g')
        t = vgpReadLine(vf);

      nused += 1;
      size[byte] += 1;
      bptr[byte]  = aptr;
    }

  parm->nused = nused;
  return (NULL);
}

/****************************************************************************************
 *
 *  MAIN
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ char    *command;
  char    *fname;
  VgpFile *vf;

  uint32 *count;        // Barcode list and eventually valid barcode list & their counts
  int64   npairs;       // # of reads in each file
  int64   nreads;       // # of reads in each file

  uint8   map[256];     // Map from QV char to bit code
  uint8   inv[256];     // Map from bit code to QV char
  int     qbits;        // # of bits to encode QV values
  uint8   qmask;        // mask of lowest qbits ina a byte

  uint8  *valid;        // 4^16 bit vector of good codes
  int64   ndist;
  int64   ngood;
  int64   nused;

  int     reclen;       // Compressed record length
  int     flen, rlen;
  int     fclen, rclen; // Length of compressed fields
  int     fqlen, rqlen;

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
    DIR   *dirp;

    ARG_INIT("VGPcloud")

    SORT_PATH = "/tmp";
    NTHREADS  = 4;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'P':
            SORT_PATH = argv[i]+2;
            if ((dirp = opendir(SORT_PATH)) == NULL)
              { fprintf(stderr,"%s: -P option: cannot open directory %s\n",Prog_Name,SORT_PATH);
                exit (1);
              }
            closedir(dirp);
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc != 2)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, show progress as proceed.\n");
        fprintf(stderr,"      -P: Do external sorts in this directory.\n");
        fprintf(stderr,"      -T: Use -T threads.\n");
 
        exit (1);
      }
  }

  //  Open .irp VGPfile file with NTHREADS

  { char *pwd, *root;
    int   i;

    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".irp");
    fname = Strdup(Catenate(pwd,"/",root,".irp"),"Allocating full path name");
    vf    = vgpFileOpenRead(fname,SEQ,NTHREADS);
    if (vf == NULL)
      { fprintf(stderr,"%s: Cannot open %s as an .irp file\n",Prog_Name,argv[1]);
        exit (1);
      }
    free(root);
    free(pwd);

    if ( ! vf->isBinary)
      { fprintf(stderr,"%s: Input is not a binary file\n",Prog_Name);
        exit (1);
      }
    if (vf->subType != IRP)
      { fprintf(stderr,"%s: Input is not a .irp file\n",Prog_Name);
        exit (1);
      }
    if (vf->lineInfo['Q']->given.count <= 0)
      { fprintf(stderr,"%s: Input file does not have QV vectors\n",Prog_Name);
        exit (1);
      }
    for (i = 'A'; i < 'Z'; i++)
      if (vf->lineInfo[i] != NULL && vf->lineInfo[i]->given.count > 0)
        if ( ! (i == 'S' || i == 'Q' || i == 'P'))
          { fprintf(stderr,"%s: Input contains data lines other than S, Q, P, or g\n",Prog_Name);
            exit (1);
          }
    if (vf->lineInfo['S']->given.count != vf->lineInfo['Q']->given.count)
      { fprintf(stderr,"%s: The number of sequences and QV's are not equal\n",Prog_Name);
        exit (1);
      }
  }

  //  Pass 1: Build list of barcodes & build fixed-bit QV code on first 10 million symbols

  if (VERBOSE)
    { fprintf(stderr,"  Scanning barcodes in file %s\n",argv[1]);
      fflush(stderr);
    }

  nreads = vf->lineInfo['S']->given.count;
  npairs = vf->lineInfo['P']->given.count;

  count = (uint32 *) Malloc(sizeof(uint32)*(2*npairs+1),"Allocating barcode array");

  { BarCode_Arg parm[NTHREADS];
    pthread_t   threads[NTHREADS];
    
    int64     usedqvs[256];
    int64     line;
    int       i, n, c, p;

    for (i = 0; i < NTHREADS; i++)
      { parm[i].beg   = (npairs * i) / NTHREADS;
        parm[i].end   = (npairs * (i+1)) / NTHREADS;
        parm[i].vf    = vf+i;
        parm[i].count = count + parm[i].beg;
#ifdef DEBUG
        barcodes_thread(parm+i);
#else
        pthread_create(threads+i,NULL,barcodes_thread,parm+i);
#endif
      }

    bzero(usedqvs,sizeof(int64)*256);

    line = 0;
    for (i = 0; i < NTHREADS; i++)
      { 
#ifndef DEBUG
        pthread_join(threads[i],NULL);
#endif
        if (parm[i].error > 0)
          { error_report(parm[i].error,line + vf[i].line);
            exit (1);
          }
        if (i == 0)
          flen = parm[0].flen;
        else
          { if (flen != parm[i].flen)
              { error_report(ERROR_SLEN,line);
                exit (1);
              }
          }
        if (i == 0)
          rlen = parm[0].rlen;
        else
          { if (rlen != parm[i].rlen)
              { error_report(ERROR_SLEN,line);
                exit (1);
              }
          }
        line += vf[i].line;
        for (n = 0; n < 256; n++)
          usedqvs[n] += parm[i].usedqvs[n];
      }

    n = 0;                         //  Map non-empy QV histogram entries to a bit code, and
    for (i = 0; i < 256; i++)      //    all empty entries to the nearest code,
      if (usedqvs[i])              //    and determine # of bits to encode largest
        { map[i] = n;
          inv[n] = i;
          if (n == 0)
            { for (c = 0; c < i; c++)
                map[c] = n;
            }
          else
            { p = (inv[n-1] + i)/2;
              for (c = inv[n-1]+1; c < p; c++)
                map[c] = n-1;
              for (c = p; c < i; c++)
                map[c] = n;
            }
          n += 1;
        }
    n -= 1;
    for (c = inv[n]+1; c < 256; c++)
      map[c] = n;
    for (qbits = 0; n > 0; qbits++)
      n >>= 1;
    qmask  = (1 << qbits) - 1;

    if (VERBOSE)
      { fprintf(stderr,"    Forward / reverse lengths = %d / %d\n",flen,rlen);
        fprintf(stderr,"    QVs will take %d bits per base\n",qbits);
        fflush(stderr);
      }
  }

  //  Sort barcode list,then compress to list of valid codes, and finally
  //    build bit vector 'valid' of good codes and those uniquely 1-away from a good code.

  if (VERBOSE)
    { fprintf(stderr,"  About to sort ");
      Print_Number((int64) npairs,0,stderr);
      fprintf(stderr," barcodes\n");
      fflush(stderr);
    }

  { Correct_Arg parm[NTHREADS];
    pthread_t   threads[NTHREADS];

    uint8 *vector;
    int    i, barsort[5];

#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
    for (i = 0; i < 4; i++)
      barsort[i] = i;
#else
    for (i = 0; i < 4; i++)
      barsort[i] = 3-i;
#endif
    barsort[4] = -1;

    Set_LSD_Params(NTHREADS,VERBOSE);

    LSD_Sort(npairs,count,count+npairs,sizeof(uint32),sizeof(uint32),barsort);

    count[npairs] = count[npairs-1] + 1;

    // Generate 1-neighborhood vector

    count  = Realloc(count,sizeof(uint32)*(npairs+1),"Resizing valid codes vector");
    valid  = (uint8 *) Malloc(0x100000000ll,"Bit vectors");
    vector = (uint8 *) Malloc(0x20000000ll*NTHREADS,"Bit vectors");
    if (count == NULL || valid == NULL || vector == NULL)
      exit (1);

    if (VERBOSE)
      { fprintf(stderr,"  Sorted, now anlyzing ...\n");
        fflush(stderr);
      }

    parm[i].beg = 0;
    for (i = 0; i < NTHREADS; i++)
      { int64 e;

        if (i == 0)
          parm[i].beg = 0;
        else
          parm[i].beg = parm[i-1].end;

        parm[i].valid  = valid;
        parm[i].vector = vector + i * 0x20000000ll;
        parm[i].count  = count;

        e = (npairs*(i+1)) / NTHREADS;
        if (e > 0)
          while (count[e-1] == count[e])
            e += 1;
        parm[i].end = e;
      }

    // Zero all bit vectors

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,zero_thread,parm+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

    //  Find valid bar codes

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,good_thread,parm+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

    ndist = 0;
    ngood = 0;
    for (i = 0; i < NTHREADS; i++)
      { ndist += parm[i].ndist;
        ngood += parm[i].ngood;
      }

    if (VERBOSE)
      { fprintf(stderr,"  There are ");
        Print_Number(ngood,0,stderr);
        fprintf(stderr," (%.1f%%) valid codes with ",(100.*ngood)/npairs);
        Print_Number(ndist,0,stderr);
        fprintf(stderr," distinct values.\n");
        fflush(stderr);
      }

    //  Set neighborhood vectors

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,mark_thread,parm+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

     //  Set 1-correction where possible

    for (i = 0; i < NTHREADS; i++)
      { parm[i].vector = vector;
        pthread_create(threads+i,NULL,fix_thread,parm+i);
      }

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

exit (1);
  }


  //  Pass 2: Correct barcodes when possible, output compressed pairs for sorting

  { int64  Array_Size;
    uint8 *array;
    int    fd;

    if (VERBOSE)
      { fprintf(stderr,"  Scan to compressing data for external sort\n");
        fprintf(stderr,"         and repair barcodes where possible.\n");
        fflush(stderr);
      }

    fd = open(Catenate(SORT_PATH,"/",fname,".sort"),O_RDWR);
    if (fd == -1)
      { fprintf(stderr,"%s: Cannot open %s/%s.sort\n",Prog_Name,SORT_PATH,fname);
        exit (1);
      }

    Array_Size   = reclen*nreads;
    Block_Size   = (reclen << 12);
    Array_Size  += ((Array_Size+(Block_Size-1))/Block_Size)*sizeof(uint8 *);
    Block_Size  += sizeof(uint8 *);
    Alloc_Offset = NTHREADS*256*Block_Size;
    Array_Size  += Alloc_Offset;

    array = (uint8 *) mmap(NULL,Array_Size,PROT_READ | PROT_WRITE,MAP_PRIVATE,fd,0);
    if (array == NULL)
      { close(fd);
        fprintf(stderr,"%s: Cannot memory map %s/%s.sort\n",Prog_Name,SORT_PATH,fname);
        exit (1);
      }
  
    madvise(array, Array_Size, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);

    { Load_Arg  parm[NTHREADS];
      pthread_t threads[NTHREADS];

      int i;

      for (i = 0; i < NTHREADS; i++)
        { parm[i].tidx   = i;
          parm[i].vf     = vf+i;
          parm[i].beg    = (npairs * i) / NTHREADS;
          parm[i].end    = (npairs * (i+1)) / NTHREADS;
          parm[i].array  = valid;
          parm[i].qbits  = qbits;
          parm[i].map    = map;
          parm[i].flen   = flen;
          parm[i].rlen   = rlen;
          parm[i].fclen  = (flen*2-1)/8 + 1; 
          parm[i].rclen  = (rlen*2-1)/8 + 1; 
          parm[i].fqlen  = (flen*qbits-1)/8 + 1; 
          parm[i].rqlen  = (rlen*qbits-1)/8 + 1; 
          parm[i].array  = array;
#ifdef DEBUG
          load_thread(parm+i);
#else
          pthread_create(threads+i,NULL,load_thread,parm+i);
#endif
        }

      for (i = 0; i < NTHREADS; i++)
        { 
#ifndef DEBUG
          pthread_join(threads[i],NULL);
#endif
        }

      if (VERBOSE)
        { fprintf(stderr,"  ");
          Print_Number((int64) ngood,0,stderr);
          fprintf(stderr," (%.1f%%) pairs have good codes, ",(100.*ngood)/nreads);
          Print_Number((int64) (nused-ngood),0,stderr);
          fprintf(stderr," (%.1f%%) have repairable codes, and ",(100.*(nused-ngood))/nreads);
          Print_Number((int64) (nreads-nused),0,stderr);
          fprintf(stderr," (%.1f%%) were dropped.\n",(100.*(nreads-nused))/nreads);
          fflush(stderr);
        }
    }
  }

  //  External sort on barcode (1st 4 bytes of each record)

  if (VERBOSE)
    { fprintf(stderr,"  Performing external sort.\n");
      fflush(stderr);
    }

  Ex_sort(Catenate(SORT_PATH,"/",fname,".sort"),reclen,4,NTHREADS);

  //  Read in sorted array and output to standard out

  if (VERBOSE)
    { fprintf(stderr,"  Final scan to produce cloud grouped pairs on stdout.\n");
      fflush(stderr);
    }

  //  Output clouds

  { FILE    *sfile;
    VgpFile *vg;
    int      gc;
    uint8   *bp, *fqvs;
    uint32   val;
    char    *forw, *fcode;
    char    *pwd, *root, *gname;

    sfile = fopen(Catenate(SORT_PATH,"/",fname,".sort"),"r");
    if (sfile == NULL)
      { fprintf(stderr,"%s: Cannot create %s.sort in directory %s\n",Prog_Name,fname,SORT_PATH);
        exit (1);
      }

    pwd   = PathTo(argv[1]);
    root  = Root(argv[1],".irp");
    gname = Strdup(Catenate(pwd,"/",root,".10x"),"Allocating full path name");
    vg    = vgpFileOpenWriteNew(gname,SEQ,X10,TRUE,NTHREADS);
    if (vg == NULL)
      { fprintf(stderr,"%s: Cannot open %s.10x for writing\n",Prog_Name,root);
        exit (1);
      }
    free(gname);
    free(root);
    free(pwd);

    vgpInheritProvenance(vg,vf);
    vgpAddProvenance(vg,"VGPcloud","1.0",command,NULL);

    vgpGotoObject(vf,0);

    forw  = vf->lineInfo['S']->buffer;
    fqvs  = vf->lineInfo['Q']->buffer;
    fcode = (char *) Malloc(fclen,"Allocating compressed DNA buffer");

    gc  = 0;
    bp  = (uint8 *) count;
    val = (bp[0] << 24 | bp[1] << 16 | bp[2] << 8 | bp[3]);
    while (fread(forw,sizeof(uint8),fclen,sfile) >= (uint32) fclen)
      { if ( *((uint32 *) forw) != val)
          { printf("g %d",count[gc+1]);
            gc += 2;
            bp += 8;
            val = (bp[0] << 24 | bp[1] << 16 | bp[2] << 8 | bp[3]);

            Uncompress_SEQ(forw,flen);

            vgpInt(vf,0) = 0;
            vgpInt(vf,1) = 16;
            vgpWriteLine(vf,'g',16,forw);
          }
        else
          Uncompress_SEQ(forw,flen);

        vgpWriteLine(vf,'P',0,NULL);

        vgpInt(vf,0) = flen-23;
        vgpWriteLine(vf,'S',flen-23,forw+23);

        fread(fqvs,sizeof(uint8),fqlen,sfile);
        Uncompress_QV(flen,fqvs,qbits,qmask,inv);

        vgpWriteLine(vf,'Q',flen-23,fqvs+23);

        fread(forw,sizeof(uint8),rclen,sfile);
        Uncompress_SEQ(forw,flen);

        vgpInt(vf,0) = rlen;
        vgpWriteLine(vf,'S',rlen,forw);

        fread(fqvs,sizeof(uint8),rqlen,sfile);
        Uncompress_QV(rlen,fqvs,qbits,qmask,inv);

        vgpWriteLine(vf,'Q',rlen,fqvs);
      }

    vgpFileClose(vg);
    fclose(sfile);
    free(fcode);
  }

  //  Tidy up just for good form

  vgpFileClose(vf);
  free(fname);
  exit (0);
}
