/*******************************************************************************************
 *
 *  VGPpacbio: Converte Pacbio subreads.bam file to a VGP .pbr file
 *
 *  Author:  Gene Myers
 *  Date  :  Aug. 9, 2019
 *
 ********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include <fcntl.h>
#include <sys/stat.h>

#include "gene_core.h"
#include "pb_expr.h"

#include "LIBDEFLATE/libdeflate.h"
typedef  struct libdeflate_decompressor Depress;

#undef   DEBUG_CORE
#undef   DEBUG_FIND
#undef   DEBUG_CHECK
#undef   DEBUG_RECORDS
#undef   DEBUG_PART
#undef   DEBUG_OUT

#define LOWER_OFFSET 32
#define INT_MAXLEN   10
#define FLOAT_MAXLEN 11

#define IO_BLOCK  10000000
#define BAM_BLOCK  0x10000
#define MAX_PREFIX     292
#define HEADER_LEN      36
#define SEQ_RUN         40

static int     VERBOSE;
static int     NTHREADS;
static int     ARROW;           //  Output arrow/A lines?
static int     UPPER;           //  Output in upper case?
static Filter *EXPR;            //  Filter expression
static int     IS_BIG_ENDIAN;   //  Is machine big-endian?

static char *Usage = "[-vaU] [-e<expr(ln>=500 && rq>=750)> [-T<int(4)>] <input:pacbio> ...";

static int DNA[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
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


/*******************************************************************************************
 *
 *  Routines to manage a BAM_FILE stream object
 *
 ********************************************************************************************/

  //   Open BAM stream where compressed block start is known.  Compressed BAM blocks are buffered
  //     in a very big IO buffer and the current uncompressed block is in a 64Kbp array.

typedef struct
  { int64  fpos;   //  offset in compressed file of desired bam block
    uint32 boff;   //  offset within bam block of desired record
  } Location;

typedef struct
  { int      fid;              //  file descriptor
    int      last;             //  last block of data in file has been read
    uint8   *buf;              //  IO buffer (of IO_BLOCK bytes, supplied by caller)
    int      blen;             //  # of bytes currently in IO buffer
    int      bptr;             //  start of next BAM block in IO buffer
    uint8    bam[BAM_BLOCK+1]; //  uncompressed bam block
    uint32   bsize;            //  length of compressed bam block
    uint32   ssize;            //  length of uncompressed bam block
    Location loc;              //  current location in bam file
    Depress *decomp;
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

#ifdef DEBUG_CORE
      printf("Move %d bytes to data+%d from bam+%d\n",chk,off,boff);
#endif
      if (data != NULL)
        memcpy(data+off,bam+boff,chk);
      off += chk;
      len -= chk;

      file->loc.fpos += file->bsize;
#ifdef DEBUG_CORE
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
#ifdef DEBUG_CORE
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
#ifdef DEBUG_CORE
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
#ifdef DEBUG_CORE
      printf("Loaded gzip block of size %d into %d\n",bsize,ssize);
#endif

      file->bsize = bsize;
      file->ssize = ssize;
      file->bptr  = bptr + bsize;
      chk = ssize;
    }

#ifdef DEBUG_CORE
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
            fprintf(stderr,"%s: SAM-line is longer than max %d\n",Prog_Name,IO_BLOCK);
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

typedef struct
  { int       fid;      //  File id (distinct for each thread)
    uint8    *buf;      //  IO buffer for this thread
    int       isbam;    //  Is a bam file (vs. sam file)
    Location  beg;      //  Scan range is [beg,end)
    Location  end;
                      //  Outputs from first pass
    int       hlen;     //  length of SMRT cell name
    int       nr;       //  # of reads
    int64     bp;       //  # of base pairs
    int       maxbp;    //  length of longest read
    int64     maxout;   //  Largest possible VGP output for a second pass thread interval
    int       maxdata;  //  Maximum bam entry data block
    int       nparts;   //  Number of ~10Mb intervals for second scan
    Location *parts;    //  Start point of each interval for second scan
                      //  Outputs from second pass
    char     *out;      //  Output buffer
    int       olen;     //  Current length of output buffer
    int       head;     //  Should the g line start the output?

    Depress  *decomp;
  } Thread_Arg;

  //  Find first record location (skip over header)

static void *header_thread(void *arg)
{ Thread_Arg *parm = (Thread_Arg *) arg;
  int         fid  = parm->fid;
  uint8      *buf  = parm->buf;

  Location zero = { 0ll, 0 };
  BAM_FILE _bam, *bam = &_bam;
  uint8    data[4];
  int      i, ntxt, ncnt, nlen;

  //  At start of file so can use BAM stream

  bam->decomp = parm->decomp;
  bam_start(bam,fid,buf,&zero);

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
      bam_get(bam,NULL,nlen);
    }

  parm->beg = bam->loc;

  return (NULL);
}

  //  Find next identifiable entry location forward of data->fpos in data->fid

static void *find_thread(void *arg)
{ Thread_Arg *parm = (Thread_Arg *) arg;
  int         fid  = parm->fid;
  uint8      *buf  = parm->buf;
  int64       fpos = parm->beg.fpos;

  Depress *decomp = parm->decomp;

  uint32 bptr, blen;
  int    last, notfound;

  uint8 *block;
  uint32 bsize, ssize;
  size_t tsize;

  BAM_FILE _bam, *bam = &_bam;

#ifdef DEBUG_FIND
  printf("Searching from %lld\n",fpos);
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
          printf("Loading %d(last=%d)\n",blen,last);
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
          printf("  Putative header @ %d\n",bptr-3);
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
          printf("    Putative Extra %d\n",bsize);
#endif
  
          isize = getint(block+(bsize-4),4);
          crc   = getint(block+(bsize-8),4);
  
          if (libdeflate_gzip_decompress(decomp,block,bsize,bam,BAM_BLOCK,&tsize) != 0)
            continue;
          ssize = tsize;

          if (ssize == isize && crc == libdeflate_crc32(0,bam->bam,ssize))
            { bptr -= 3;
              fpos  += bptr;
              notfound = 0;

#ifdef DEBUG_FIND
              printf("    First block at %lld (%d)\n",fpos,ssize);
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
    { int    j, k, beg;
      int    run, out;
      int    lname, lcigar, lseq, ldata;

      block = bam->bam;
      ssize = bam->ssize;

      run = HEADER_LEN-1;
      out = 1;
      for (j = HEADER_LEN; j < (int) 10000; j++)
        if (DNA[block[j]])
          { if (out && j >= run+SEQ_RUN)
              {
#ifdef DEBUG_FIND
                printf("Possible start @ %d (%d)\n",run+1,run-MAX_PREFIX);
#endif
                if (run < MAX_PREFIX)
                  beg = 0;
                else
                  beg = run-MAX_PREFIX;
                for (k = run-(HEADER_LEN-1); k >= beg; k--)
                  { ldata  = getint(block+k,4);
                    lname  = block[k+12];
                    lcigar = getint(block+(k+16),2);
                    lseq   = getint(block+(k+20),4);
                    if (lcigar == 0 && lseq+lname < ldata && k + 35 + lname == run)
                      if (lseq > 0 && ldata > 0 && valid_name((char *) (block+(k+36)),lname))
                        { parm->beg.fpos = bam->loc.fpos;
                          parm->beg.boff = k;
#ifdef DEBUG_FIND
                          printf("  Found '%s':%d at %d\n",block+(k+36),lseq,k);
#endif
                          return (NULL);
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

  fprintf(stderr,"%s: Could not find bam record structure!\n",Prog_Name);
  exit (1);
}

  //  Find next identifiable sam entry location forward of data->fpos in data->fid

static void *sam_thread(void *arg)
{ Thread_Arg *parm = (Thread_Arg *) arg;
  int         fid  = parm->fid;
  uint8      *buf  = parm->buf;

  BAM_FILE _bam, *bam = &_bam;

  bam->decomp = parm->decomp;
  sam_start(bam,fid,buf,&(parm->beg));

  sam_getline(bam);
  if (parm->beg.fpos == 0)
    { while (buf[bam->bptr] == '@')
        sam_getline(bam);
    }
  parm->beg = bam->loc;

  return (NULL);
}


/*******************************************************************************************
 *
 *  Routines to scan and parse a bam entry
 *
 ********************************************************************************************/

static int is_integer[128] =
  { 0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,

    0, 0,  0, 1, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 1, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 1, 0, 0, 0, 0,   0, 1, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 1, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
  };

static int bam_tag_size[128] =
  { 0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,

    0, 1, 10, 1, 0, 0, 0, 0,   9, 4, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 2, 0, 0, 0, 0,   0, 0, 9, 0, 0, 0, 0, 0,
    0, 0,  0, 1, 8, 0, 4, 0,   0, 4, 0, 0, 0, 0, 0, 0,
    0, 0,  0, 2, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0,
  };

static char INT_2_IUPAC[16] = "=acmgrsvtwyhkdbn";

  // flip_auxilliary converts big endian auxilliary data to little endian

#define SWAP(i,j)  (x = v[i], v[i] = v[j], v[j] = x)

static void flip_short(void *w)
{ uint8 x, *v = (uint8 *) w;
  SWAP(0,1);
}

static void flip_int(void *w)
{ uint8 x, *v = (uint8 *) w;
  SWAP(0,3);
  SWAP(1,2);
}

static void flip_double(void *w)
{ uint8 x, *v = (uint8 *) w;
  SWAP(0,7);
  SWAP(1,6);
  SWAP(2,5);
  SWAP(3,4);
}

static void flip_auxilliary(uint8 *s, uint8 *e)
{ int  i, n, x;

  while (s < e)
    { s += 2;
      x = bam_tag_size[*s++];
      switch (x)
      { case 1:
          s += 1;
          break;
        case 2:
          flip_short(s);
          s += 2;
          break;
        case 4:
          if (is_integer[s[-1]])
            flip_int(s);
          s += 4;
          break;
        case 8:
          flip_double(s);
          s += 8;
          break;
        case 9:   //  Z or H
          while (*s != 0)
            s += 1;
          s += 1;
          break;
        case 10:  //  B
          x = bam_tag_size[*s++];
          flip_int(s);
          n = *((int *) s);
          s += 4;
          switch (x)
          { case 1:
              s += n;
              break;
            case 2:
              for (i = 0; i < n; i++, s += 2)
                flip_short(s);
              break;
            case 4:
              for (i = 0; i < n; i++, s += 4)
                flip_int(s);
              break;
            case 8:
              for (i = 0; i < n; i++, s += 8)
                flip_double(s);
              break;
          }
          break;
       }
    }
}

  //  Scan next bam entry and load PacBio info in record 'theR'
 
static void bam_record_scan(BAM_FILE *sf, samRecord *theR, int pass2)
{ int ldata, lname, lseq, aux;

  { uint8  x[36];     //  Process 36 byte header, if pass2 get sequence
    char  *eoh;

    bam_get(sf,x,36);

    ldata  = getint(x,4) - 32;
    lname  = getint(x+12,1);
    lseq   = getint(x+20,4);

    if (ldata < 0 || lseq < 0 || lname < 1)
      { fprintf(stderr,"%s: Non-sensical BAM record, file corrupted?\n",Prog_Name);
        exit (1);
      }

    aux = lname + ((lseq + 1)>>1) + lseq;
    if (aux > ldata)
      { fprintf(stderr,"%s: Non-sensical BAM record, file corrupted?\n",Prog_Name);
        exit (1);
      }

    if (ldata > theR->dmax)
      { theR->dmax = 1.2*ldata + 1000;
        theR->data = (uint8 *) Realloc(theR->data,theR->dmax,"Reallocating data buffer");
        if (theR->data == NULL)
          exit (1);
      }

    bam_get(sf,theR->data,ldata);

    if (lseq <= 0)
      { fprintf(stderr,"%s: no sequence for subread !?\n",Prog_Name);
        exit (1);
      }

    theR->header = (char *) theR->data;
    theR->len = lseq;
    eoh = index(theR->header,'/');
    if (eoh != NULL)
      *eoh = 0;

    if (pass2)
      { uint8 *t;
        char  *s;
        int    i, e;

        e = getint(x+16,2);
        t = theR->data + (lname + (e<<2));
        s = theR->seq;
        lseq -= 1;
        for (e = i = 0; i < lseq; e++)
          { s[i++] = INT_2_IUPAC[t[e] >> 4];
            s[i++] = INT_2_IUPAC[t[e] & 0xf];
          }
        if (i <= lseq)
          s[i] = INT_2_IUPAC[t[e] >> 4];
        lseq += 1;
      }
  }

#define PULSES(type)                    \
for (i = 0; i < len; i++, p += size)    \
  { uint32 x = *((type *) p);           \
    if (x >= 5)                         \
      x = 4;                            \
    arr[i] = x+'0';  		        \
  }

#define BARCODES(type)                  \
for (i = 0; i < len; i++, p += size)    \
  { uint32 x = *((type *) p);           \
    theR->bc[i] = x;                    \
  }

#define GET_UINT(name,target)                                                   \
{ switch (p[2])                                                                 \
  { case 'C':                                                                   \
    case 'c':                                                                   \
      target = p[3];								\
      p += 4;									\
      break;                                                                    \
    case 'S':                                                                   \
    case 's':                                                                   \
      target = *((short *) (p+3));						\
      p += 5;									\
      break;                                                                    \
    case 'I':                                                                   \
    case 'i':                                                                   \
      target = *((int *) (p+3));						\
      p += 7;									\
      break;                                                                    \
    default:                                                                    \
      fprintf(stderr,"%s: %s-tag is not of integer type\n",Prog_Name,name);     \
      exit (1);                                                                 \
  }                                                                             \
}

  { int      size, len;    //  Get sn, pw, bc, bq, zm, qs, qe, rq, and np from auxilliary tags
    uint8   *p, *e;
    char    *arr;
    int      i;

    arr = theR->arr;
    p = theR->data + aux;
    e = theR->data + ldata;

    if (IS_BIG_ENDIAN)
      flip_auxilliary(p,e);

    while (p < e)
      { switch (*p)
        { case 's':
            if (ARROW && memcmp(p+1,"nBf",3) == 0)
              { len = *((int *) (p+4));
                if (len != 4)
                  { fprintf(stderr,"%s: sn-tag does not have 4 floats\n",Prog_Name);
                    exit (1);
                  }
                if (pass2)
                  { theR->snr[0] = *((float *) (p+8));
                    theR->snr[1] = *((float *) (p+12));
                    theR->snr[2] = *((float *) (p+16));
                    theR->snr[3] = *((float *) (p+20));
                  }
                p += 24;
                continue;
              }
            break;
          case 'p':
            if (ARROW && memcmp(p+1,"wB",2) == 0)
              { if (!is_integer[p[3]])
                  { fprintf(stderr,"%s: pw-tag is not of integer type\n",Prog_Name);
                    exit (1);
                  }
                len = *((int *) (p+4));
                if (len != lseq)
                  { fprintf(stderr,"%s: pw-tag is not the same length as sequence\n",Prog_Name);
                    exit (1);
                  }
                size = bam_tag_size[p[3]];
                if (pass2)
                  { p += 8;
                    switch (size)
                    { case 1:
                        PULSES(uint8)
                        break;
                      case 2:
                        PULSES(uint16)
                        break;
                      case 4:
                        PULSES(uint32)
                        break;
                    }
                  }
                else
                  p += 8 + size*len;
                continue;
              }
            break;
          case 'b':
            if (memcmp(p+1,"cB",2) == 0)
              { if (! is_integer[p[3]])
                  { fprintf(stderr,"%s: bc-tag is not of integer type\n",Prog_Name);
                    exit (1);
                  }
                size = bam_tag_size[p[3]];
                len  = *((int *) (p+4));
                if (len > 2)
                  { fprintf(stderr,"%s: More than two barcode values\n",Prog_Name);
                    exit (1);
                  }
                p += 8;
                switch (size)
                { case 1:
                    BARCODES(uint8)
                    break;
                  case 2:
                    BARCODES(uint16)
                    break;
                  case 4:
                    BARCODES(uint32)
                    break;
                }
                continue;
              }
            else if (p[1] == 'q')
              { GET_UINT("bq",theR->bqual)
                continue;
              }
            break;
          case 'z':
            if (p[1] == 'm')
              { GET_UINT("zm",theR->well)
                continue;
              }
            break;
          case 'q':
            if (p[1] == 's')
              { GET_UINT("qs",theR->beg)
                continue;
              }
            else if (p[1] == 'e')
              { GET_UINT("qe",theR->end)
                continue;
              }
            break;
          case 'r':
            if (memcmp(p+1,"qf",2) == 0)
              { theR->qual = *((float *) (p+3));
                p += 7;
                continue;
              }
            break;
          case 'n':
            if (p[1] == 'p')
              { GET_UINT("np",theR->nump)
                continue;
              }
          default:
            break;
        }

        size = bam_tag_size[p[2]];
        if (size <= 8)
          p += 3+size;
        else if (size == 9)
          p += strlen(((char *) p)+3) + 4;
        else
          { size = bam_tag_size[p[3]];
            len  = *((int *) (p+4));
            p += 8 + size*len;
          }
      }
  }
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

static void sam_record_scan(BAM_FILE *sf, samRecord *theR, int pass2)
{ char      *p;
  int        qlen;

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
    q = index(q,'/');         // Truncate pacbio well & pulse numbers
    if (q != NULL && q < p)
      *q = 0;

    for (i = 0; i < 8; i++)   // Skip next 8 required fields
      { p = index(p+1,'\t');
        CHECK( p == NULL, "Too few required fields in SAM record, file corrupted?")
      }
    p += 1;

    NEXT_ITEM(q,p)
    qlen = p-q;
    CHECK (*q == '*', "No sequence for read?");

    theR->len = qlen;
    if (pass2)
      { seq = theR->seq;
        for (i = 0; i < qlen; i++)
          seq[i] = IUPAC_2_DNA[(int) (*q++)];
      }

    p = index(p+1,'\t');  // Skip qual
    CHECK( p == NULL, "No auxilliary tags in SAM record, file corrupted?")
  }

  { char *q, *arr;       //  Get zm, qs, qe, rq, sn, and pw from auxilliary tags
    int   x, cnt;

    arr = theR->arr;
    while (*p++ == '\t')
      { switch (*p)
        { case 's':
            if (ARROW && strncmp(p+1,"n:B:f,",6) == 0)
              { p += 6;
                for (cnt = 0; *p == ',' && cnt < 4; cnt++)
                  { theR->snr[cnt] = strtod(q=p+1,&p);
                    CHECK( p == q, "Cannot parse snr value")
                  }
                CHECK ( *p == ',' || cnt < 4, "Expecting 4 snr values")
                CHECK ( *p != '\t' && *p != '\n', "Cannot parse snr values")
                continue;
              }
            break;
          case 'p':
            if (ARROW && strncmp(p+1,"w:B:",4) == 0)
              { if (p[5] == 'f' || p[5] == 'd' || p[5] == 'A' || p[5] == 'a')
                  { fprintf(stderr, "Type of pulse width tag is not integer\n");
                    exit (1);
                  }
                p += 6;
                if (pass2)
                  for (cnt = 0; *p == ',' && cnt < qlen; cnt++)
                    { x = strtol(q=p+1,&p,0);
                      CHECK( p == q, "Cannot parse pulse width value")
                      if (x >= 5)
                        x = 4;
                      arr[cnt] = x + '0';
                    }
                else
                  for (cnt = 0; *p == ',' && cnt < qlen; cnt++)
                    { x = strtol(q=p+1,&p,0);
                      CHECK( p == q, "Cannot parse pulse width value")
                    }
                CHECK ( *p == ',' || cnt < qlen, "pulse width arraw has different length than read")
                CHECK ( *p != '\t' && *p != '\n', "Cannot parse pulse width values")
                continue;
              }
            break;
          case 'b':
            if (strncmp(p+1,"c:B:",4) == 0)
              { if (p[5] == 'f' || p[5] == 'd' || p[5] == 'A' || p[5] == 'a')
                  { fprintf(stderr, "Type of barcode tag is not integer\n");
                    exit (1);
                  }
                p += 6;
                for (cnt = 0; *p == ',' && cnt < 2; cnt++)
                  { x = strtol(q=p+1,&p,0);
                    CHECK( p == q, "Cannot parse barcode value")
                    theR->bc[cnt] = x;
                  }
                CHECK ( *p == ',' || cnt < 2, "more than 2 barcode values?")
                CHECK ( *p != '\t' && *p != '\n', "Cannot parse barcode values")
                continue;
              }
            else if (strncmp(p+1,"q:",2) == 0 && is_integer[(int) p[3]] && p[4] == ':')
              { theR->bqual = strtol(p+5,&q,0);
                CHECK (p+5 == q, "Could not parse integer barcode quality")
                p = q;
                continue;
              }
            break;
          case 'z':
            if (strncmp(p+1,"m:",2) == 0 && is_integer[(int) p[3]] && p[4] == ':')
              { theR->well = strtol(p+5,&q,0);
                CHECK (p+5 == q, "Could not parse integer well number")
                p = q;
                continue;
              }
            break;
          case 'q':
            if (strncmp(p+1,"s:",2) == 0 && is_integer[(int) p[3]] && p[4] == ':')
              { theR->beg = strtol(p+5,&q,0);
                CHECK (p+5 == q, "Could not parse integer start pulse")
                p = q;
                continue;
              }
            else if (strncmp(p+1,"e:",2) == 0 && is_integer[(int) p[3]] && p[4] == ':')
              { theR->end = strtol(p+5,&q,0);
                CHECK (p+5 == q, "Could not parse integer end pulse")
                p = q;
                continue;
              }
            break;
          case 'r':
            if (strncmp(p+1,"q:f:",4) == 0)
              { theR->qual = strtod(p+5,&q);
                CHECK (p+5 == q, "Could not parse floating point quality value")
                p = q;
                continue;
              }
            break;
          case 'n':
            if (strncmp(p+1,"np:",3) == 0 && is_integer[(int) p[3]] && p[4] == ':')
              { theR->nump = strtol(p+5,&q,0);
                CHECK (p+5 == q, "Could not parse integer number of passes")
                p = q;
                continue;
              }
            break;
          default:
            break;
        }

        while (*p != '\t' && *p != '\n')
          p += 1;
      }
  }
}


/*******************************************************************************************
 *
 *  Pass 1: In NTHREAD partition of bam file, check record and collect information of second pass
 *
 ********************************************************************************************/

static void *check_thread(void *arg)
{ Thread_Arg  *parm  = (Thread_Arg *) arg;
  int          fid   = parm->fid;
  uint8       *buf   = parm->buf;
  int64        epos  = parm->end.fpos;
  uint32       eoff  = parm->end.boff;
  int          isbam = parm->isbam;

  BAM_FILE     _bam, *bam = &_bam;
  samRecord    _theR, *theR = &_theR;

  Location lpos, ldiv, *parts;
  int64    divpt, lbp, maxout;
  int      lnr, nparts, npmax;

  int   nr, maxbp, hlen, first; 
  int64 bp;

  theR->dmax = 0;
  theR->data = NULL;

  nr    = 0;
  bp    = 0;
  first = 1;
  maxbp = 0;

  maxout = 0;
  nparts = 0;
  npmax  = 1.2 * ((epos - parm->beg.fpos) / (IO_BLOCK - BAM_BLOCK)) + 100;
  parts  = (Location *) Malloc(sizeof(Location)*npmax,"Allocating location array");
  if (parts == NULL)
    exit (1);

  bam->decomp = parm->decomp;
  if (isbam)
    bam_start(bam,fid,buf,&(parm->beg));
  else
    sam_start(bam,fid,buf,&(parm->beg));

#ifdef DEBUG_CHECK
  printf("\nStart: %lld %d to %lld %d (%d)\n",bam->loc.fpos,bam->loc.boff,epos,eoff,npmax);
#endif

  lnr  = 0;
  lbp  = 0;
  ldiv = parm->beg;
#ifdef DEBUG_PART
  printf("  Divpt @  %12lld\n",ldiv.fpos);
#endif
  divpt = ldiv.fpos + IO_BLOCK;
  parts[nparts++] = ldiv;
  lpos = parm->beg;

  while (bam->loc.fpos != epos || bam->loc.boff != eoff)
    { if (isbam)
        bam_record_scan(bam,theR,0);
      else
        sam_record_scan(bam,theR,0);

      if (bam->loc.fpos > divpt)
        { int64 dnr, dbp, omax;

          if (nr == lnr)
            { dnr = 1;
              dbp = (bp-lbp)+theR->len;
#ifdef DEBUG_PART
              printf("  Divpt+ @ %12lld   s=%lld r=%lld bp=%lld\n",
                     bam->loc.fpos+bam->loc.boff,bam->loc.fpos-ldiv.fpos,dnr,dbp);
#endif
              lnr  = nr+1;
              lbp  = bp+theR->len;
              ldiv = bam->loc;
            }
          else
            { dnr = nr-lnr;
              dbp = bp-lbp;
#ifdef DEBUG_PART
              printf("  Divpt  @ %12lld   s=%lld r=%lld bp=%lld\n",
                     bam->loc.fpos+bam->loc.boff,bam->loc.fpos-ldiv.fpos,dnr,dbp);
#endif
              lnr = nr;
              lbp = bp;
              ldiv = lpos;
            }
          divpt = ldiv.fpos + IO_BLOCK;

          if (nparts >= npmax)
            { npmax = 1.2*nparts + 100;
              parts  = (Location *)
                          Realloc(parts,sizeof(Location)*npmax,"Allocating location array");
              if (parts == NULL)
                exit (1);
            }
          parts[nparts++] = ldiv;

          omax = dnr*(15 + 4*INT_MAXLEN) + dbp;
          if (ARROW)
            omax += dnr*(10 + INT_MAXLEN + 4*FLOAT_MAXLEN) + dbp;
          if (omax > maxout)
            maxout = omax;
        }
      lpos = bam->loc;

      if ( ! evaluate_bam_filter(EXPR,theR))
        continue;

      if (first)
        { first = 0;
          hlen  = strlen(theR->header);
        }

      nr += 1;
      bp += theR->len;
      if (theR->len > maxbp)
        maxbp = theR->len;
    }

  if (nr > lnr)
    { int64 dnr, dbp, omax;

      dnr = nr-lnr;
      dbp = bp-lbp;
      omax = dnr*(15 + 4*INT_MAXLEN) + dbp;
      if (ARROW)
        omax += dnr*(10 + INT_MAXLEN + 4*FLOAT_MAXLEN) + dbp;
      if (omax > maxout)
        maxout = omax;
    }
  else
    nparts -= 1;

  if (isbam)
    free(theR->data);

#ifdef DEBUG_CHECK
  printf("  nr = %d bp = %lld maxbp = %d\n",nr,bp,maxbp);
  printf("  parts = %d   max = %lld\n",nparts,maxout);
#endif

  parm->nr      = nr;
  parm->bp      = bp;
  parm->hlen    = hlen;
  parm->maxbp   = maxbp;
  parm->maxout  = maxout;
  parm->maxdata = theR->dmax;
  parm->nparts  = nparts;
  parm->parts   = parts;

  return (NULL);
}


/*******************************************************************************************
 *
 *  Pass 2: Produe the data NTHREAD "blocks" at a time in buffer in parallel and then
 *             stream the buffer to standard out serially.
 *
 ********************************************************************************************/

  //  Write subread data in samRecord rec to non-NULL file types

static void *output_thread(void *arg)
{ Thread_Arg  *parm  = (Thread_Arg *) arg;
  int          fid   = parm->fid;
  uint8       *buf   = parm->buf;
  int          isbam = parm->isbam;
  int64        epos  = parm->end.fpos;
  uint32       eoff  = parm->end.boff;
  char        *out   = parm->out;
  int          head  = parm->head;

  BAM_FILE     _bam, *bam = &_bam;
  samRecord    _theR, *theR = &_theR;

  int          i, olen;

  //  Know the max size of sequence and data from pass 1, so set up accordingly

  theR->dmax = parm->maxdata;
  theR->data = Malloc(theR->dmax,"Allocating sequence array");
  theR->seq  = Malloc((1+ARROW)*parm->maxbp,"Allocating sequence array");
  if (theR->seq == NULL || theR->data == NULL)
    exit (1);
  if (ARROW)
    theR->arr = theR->seq + parm->maxbp;

  bam->decomp = parm->decomp;
  if (isbam)
    bam_start(bam,fid,buf,&(parm->beg));
  else
    sam_start(bam,fid,buf,&(parm->beg));

#ifdef DEBUG_OUT
  printf("Block: %12lld / %5d to %12lld / %5d --> %8lld\n",bam->loc.fpos,bam->loc.boff,epos,eoff,
                                                           epos - bam->loc.fpos);
  fflush(stdout);
#endif


  olen = 0;
  while (bam->loc.fpos != epos || bam->loc.boff != eoff)
    { if (isbam)
        bam_record_scan(bam,theR,1);
      else
        sam_record_scan(bam,theR,1);

#ifdef DEBUG_RECORDS
      printf("S = '%s'\n",theR->seq);
      printf("A = '%s'\n",theR->arr);
      printf("zm = %d\n",theR->well);
      printf("ln = %d\n",theR->len);
      printf("rq = %g\n",theR->qual);
      printf("bc1 = %d\n",theR->bc[0]);
      printf("bc2 = %d\n",theR->bc[1]);
      printf("bq = %d\n",theR->bqual);
      printf("np = %d\n",theR->nump);
      printf("qs = %d\n",theR->beg);
      printf("qe = %d\n",theR->end);
#endif

      if (head > 0)
        { int len = strlen(theR->header);
          olen += sprintf(out+olen,"g %d %d %.*s\n",head,len,len,theR->header);
          head = 0;
        }

      if ( ! evaluate_bam_filter(EXPR,theR))
        continue;

      if (UPPER)
        { if (islower(theR->seq[0]))
            for (i = 0; i < theR->len; i++)
              theR->seq[i] -= LOWER_OFFSET;
        }
      else
        { if (isupper(theR->seq[0]))
            for (i = 0; i < theR->len; i++)
              theR->seq[i] += LOWER_OFFSET;
        }
    
      olen += sprintf(out+olen,"S %d %.*s\n",theR->len,theR->len,theR->seq);
      olen += sprintf(out+olen,"W %d %d %d 0.%0d\n",
                               theR->well,theR->beg,theR->end,(int) (theR->qual*1000.));
    
      if (ARROW)
        { olen += sprintf(out+olen,"N %.2f %.2f %.2f %.2f\n",
                                   theR->snr[0],theR->snr[1],theR->snr[2],theR->snr[3]);
          olen += sprintf(out+olen,"A %d %.*s\n",theR->len,theR->len,theR->arr);
        }
    }

  parm->olen = olen;
  return (NULL);
}


/*******************************************************************************************
 *
 *  Main routine: parameter, 2 passes to stdout. 
 *
 *******************************************************************************************/

  //  Main

typedef struct
  { char      *fname;   //  Full path name of file
    int64      fsize;   //  size of file
    int        isbam;   //  bam or sam?
    int        nread;   //  number of reads in this file
    int       *nparts;  //  partition and location point computed in pass 1 for pass 2
    Location **parts;
  } File_Object;

int main(int argc, char* argv[])
{ int    Oargc;
  char **Oargv; 

  Oargc = argc;
  Oargv = argv;

  //  Process command line arguments

  { int   i, j, k;
    int   flags[128];
    char *eptr;
    int   one = 1;

    ARG_INIT("VGPpacbio")

    EXPR     = NULL;
    NTHREADS = 4;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vaU")
            break;
          case 'e':
            EXPR = parse_filter(argv[i]+2);
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    ARROW   = flags['a'];
    UPPER   = flags['U'];

    if (EXPR == NULL)
      EXPR = parse_filter("ln>=500 && rq>=750");

    if (argc == 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: verbose mode, output progress as proceed\n");
        fprintf(stderr,"      -a: extract Arrow information on N- and A-lines.\n");
        fprintf(stderr,"      -U: use upper-case for DNA, default is lower-case\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -e: subread selection expression.  Possible variables are:\n");
        fprintf(stderr,"           zm  - well number\n");
        fprintf(stderr,"           ln  - length of subread\n");
        fprintf(stderr,"           rq  - quality value of subread (normalized to [0,1000])\n");
        fprintf(stderr,"           bc1 - # of first barcode\n");
        fprintf(stderr,"           bc2 - # of second barcode\n");
        fprintf(stderr,"           bq  - quality of barcode detection (normalized to [0,100])\n");
        fprintf(stderr,"           np  - number of passes producing subread\n");
        fprintf(stderr,"           qs  - start pulse of subread\n");
        fprintf(stderr,"           qe  - last pulse of subread\n");
        exit (1);
      }

    IS_BIG_ENDIAN = ( *((char *) (&one)) == 0);
  }

  //  Begin main code block

  { File_Object fobj[argc-1];
    int         idout;

#if ! defined(DEBUG_CORE) || ! defined(DEBUG_FIND) || ! defined(DEBUG_CHECK)
    pthread_t  threads[NTHREADS];
#endif
    Thread_Arg parm[NTHREADS];

    //  Setup the initial file objects for each file and thread IO buffers

    { int       *np;
      Location **pt;
      uint8     *bf;
      int        f, i;

      np = (int *) Malloc((argc-1)*NTHREADS*sizeof(int),"Allocating part vector");
      pt = (Location **) Malloc((argc-1)*NTHREADS*sizeof(Location *),"Allocating part vector");
      bf = Malloc(NTHREADS*IO_BLOCK,"Allocating IO_Buffer\n");
      if (np == NULL || pt == NULL || bf == NULL)
        exit (1);

      for (f = 0; f < argc-1; f++)
        { struct stat stats;
          char *suffix[4] = { ".subreads.bam", ".bam", ".subreads.sam", ".sam" };
          char *pwd, *root;
          int   fid;

          pwd = PathTo(argv[f+1]);
          OPEN2(argv[f+1],pwd,root,fid,suffix,4);
          if (fid < 0)
            { fprintf(stderr,"%s: Cannot open %s as an .subreads.bam/sam file\n",
                             Prog_Name,argv[f+1]);
              exit (1);
            }

          if (fstat(fid, &stats) == -1)
            { fprintf(stderr,"%s: Cannot get stats for %s\n",Prog_Name,argv[f+1]);
              exit (1);
            }

          fobj[f].fname  = Strdup(Catenate(pwd,"/",root,suffix[i]),"Allocating full path name");
          fobj[f].fsize  = stats.st_size;
          fobj[f].nparts = np + f*NTHREADS;
          fobj[f].parts  = pt + f*NTHREADS;
          fobj[f].isbam  = (i < 2);

          free(pwd);
          free(root);
          close(fid);
        }

      for (i = 0; i < NTHREADS; i++)
        { parm[i].buf = bf + i*IO_BLOCK;
          parm[i].decomp = libdeflate_alloc_decompressor();
        }

      idout = fileno_unlocked(stdout);   //  Is the unlock helping?
    }

    //  First pass, find valid locations for first pass NTHREAD partition, and then perform
    //    scan checking syntax, accumulating statistics, and recording valid locations for
    //    finer grained second pass.

    { int      nread, rg_nread, gmaxc;
      int64    totbp, rg_totbp, gtotc;
      int      maxbp, maxdata, maxout;

      int      f, i;

      totbp    = 0;
      nread    = 0;
      rg_nread = 0;
      rg_totbp = 0;
      gtotc    = 0;
      gmaxc    = 0;
      maxbp    = 0;
      maxout   = 0;
      maxdata  = 0;
      for (f = 0; f < argc-1; f++)
        { if (VERBOSE)
            { fprintf(stderr,"  Checking syntax of file %s\n",fobj[f].fname);
              fflush(stderr);
            }

          for (i = 0; i < NTHREADS; i++)
            { parm[i].beg.fpos = (fobj[f].fsize*i)/NTHREADS;
              parm[i].fid      = open(fobj[f].fname,O_RDONLY);
              parm[i].isbam    = fobj[f].isbam;
            }

          //  Find valid locations for each thread

          if (fobj[f].isbam)
            { for (i = 1; i < NTHREADS; i++)
#ifdef DEBUG_FIND
                find_thread(parm+i);
#else
                pthread_create(threads+i,NULL,find_thread,parm+i);
#endif
              header_thread(parm);
            }
          else
            { for (i = 1; i < NTHREADS; i++)
#ifdef DEBUG_FIND
                sam_thread(parm+i);
#else
                pthread_create(threads+i,NULL,sam_thread,parm+i);
#endif
              sam_thread(parm);
            }
#ifndef DEBUG_FIND
          for (i = 1; i < NTHREADS; i++)
            pthread_join(threads[i],NULL);
#endif

          for (i = 1; i < NTHREADS; i++)
            { parm[i-1].end.fpos = parm[i].beg.fpos;
              parm[i-1].end.boff = parm[i].beg.boff;
            }
          parm[NTHREADS-1].end.fpos = fobj[f].fsize;
          parm[NTHREADS-1].end.boff = 0;

#ifdef DEBUG_FIND
          for (i = 0; i < NTHREADS; i++)
            printf(" %2d: %12lld / %5d\n",i,parm[i].beg.fpos,parm[i].beg.boff);
          printf(" %2d: %12lld / %5d\n",i,parm[NTHREADS-1].end.fpos,parm[NTHREADS-1].end.boff);
#endif

          //  Perform first pass scans between valid locations just found

          for (i = 1; i < NTHREADS; i++)
#ifdef DEBUG_CHECK
            check_thread(parm+i);
#else
            pthread_create(threads+i,NULL,check_thread,parm+i);
#endif
          check_thread(parm);
#ifndef DEBUG_CHECK
          for (i = 1; i < NTHREADS; i++)
            pthread_join(threads[i],NULL);
#endif

          //  Record valid block locations for second pass in per file object

          for (i = 0; i < NTHREADS; i++)
            { fobj[f].nparts[i] = parm[i].nparts;
              fobj[f].parts[i]  = parm[i].parts;
              close(parm[i].fid);
            }

          //  Accumulate stats for VGP header

          { int   fnr;
            int64 fbp; 

            fnr = 0;
            fbp = 0;
            for (i = 0; i < NTHREADS; i++)
              { fnr += parm[i].nr;
                fbp += parm[i].bp;
                if (parm[i].maxbp > maxbp)
                  maxbp = parm[i].maxbp;
                if (parm[i].maxout > maxout)
                  maxout = parm[i].maxout;
                if (parm[i].maxdata > maxdata)
                  maxdata = parm[i].maxdata;
              }

            fobj[f].nread = fnr;
   
            nread += fnr;
            totbp += fbp;
            gtotc += parm[0].hlen;
            if (parm[0].hlen > gmaxc)
              gmaxc = parm[0].hlen;
            if (fnr > rg_nread)
              rg_nread = fnr;
            if (fbp > rg_totbp)
              rg_totbp = fbp;
          }
        }

      { int   n, i;
        char *out;
        int    clen, olen;
        char   date[20];
        time_t seconds;

        //  Use cumulative metrics to inform 2nd pass

        if (maxout*NTHREADS < 5000)
          out = Malloc(5000,"Allocating IO_Buffer\n");
        else
          out = Malloc(NTHREADS*maxout,"Allocating IO_Buffer\n");
        if (out == NULL)
          exit (1);

        for (n = 0; n < NTHREADS; n++)
          { parm[n].maxbp   = maxbp;
            parm[n].maxdata = maxdata;
            parm[n].out = out + n*maxout;
          }

        //  Output size headers

        olen = sprintf(out,"1 3 seq 1 0\n");
        olen += sprintf(out+olen,"2 3 pbr\n");

        clen = -1;
        for (i = 1; i < Oargc; i++)
          clen += strlen(Oargv[i])+1;
  
        olen += sprintf(out+olen,"! 9 VGPpacbio 3 1.0 %d",clen);
        for (i = 1; i < Oargc; i++)
          olen += sprintf(out+olen," %s",Oargv[i]);
        seconds = time(NULL);
        strftime(date,20,"%F_%T",localtime(&seconds));
        olen += sprintf(out+olen," 19 %s\n",date);

        olen += sprintf(out+olen,"# g %d\n",argc-1);
        olen += sprintf(out+olen,"# S %d\n",nread);
        olen += sprintf(out+olen,"# W %d\n",nread);
        if (ARROW)
          { olen += sprintf(out+olen,"# N %d\n",nread);
            olen += sprintf(out+olen,"# A %d\n",nread);
          }
        olen += sprintf(out+olen,"+ g %lld\n",gtotc);
        olen += sprintf(out+olen,"+ S %lld\n",totbp);
        if (ARROW)
          olen += sprintf(out+olen,"+ A %lld\n",totbp);

        olen += sprintf(out+olen,"@ g %d\n",gmaxc);
        olen += sprintf(out+olen,"@ S %d\n",maxbp);
        if (ARROW)
          olen += sprintf(out+olen,"@ A %d\n",maxbp);
  
        olen += sprintf(out+olen,"%% g # S %d\n",rg_nread);
        olen += sprintf(out+olen,"%% g # W %d\n",rg_nread);
        if (ARROW)
          { olen += sprintf(out+olen,"%% g # N %d\n",rg_nread);
            olen += sprintf(out+olen,"%% g # A %d\n",rg_nread);
          }
  
        olen += sprintf(out+olen,"%% g + S %lld\n",rg_totbp);
        if (ARROW)
          olen += sprintf(out+olen,"%% g + A %lld\n",rg_totbp);

        write(idout,out,olen);
      }
    }
  
    //  Second pass: Scan files and output .pbr

    { int   f, i, p, n;
      int   uthreads;

      for (f = 0; f < argc-1; f++)
        { for (n = 0; n < NTHREADS; n++)
            { parm[n].fid   = open(fobj[f].fname,O_RDONLY);
              parm[n].head  = 0;
              parm[n].isbam = fobj[f].isbam;
            }
          uthreads = n;

          if (VERBOSE)
            { fprintf(stderr,"  Computing output for file %s\n",fobj[f].fname);
              fflush(stderr);
            }
 
          parm[0].head = fobj[f].nread;
          i = 0;
          p = 0;
          while (i < NTHREADS)
            { for (n = 0; n < NTHREADS; n++)
                { parm[n].beg = fobj[f].parts[i][p];
                  p += 1;
                  if (p >= fobj[f].nparts[i])
                    { p = 0;
                      i += 1;
                      if (i >= NTHREADS)
                        { parm[n].end.fpos = fobj[f].fsize;
                          parm[n].end.boff = 0;
                          uthreads = n+1;
                          break;
                        }
                    }
                  parm[n].end = fobj[f].parts[i][p];
                }

#ifdef DEBUG_OUT
              for (n = 0; n < uthreads; n++)
                output_thread(parm+n);
#else
              for (n = 1; n < uthreads; n++)
                pthread_create(threads+n,NULL,output_thread,parm+n);
              output_thread(parm);

              for (n = 0; n < uthreads; n++)
                { pthread_join(threads[n],NULL);
                  write(idout,parm[n].out,parm[n].olen);
                }
#endif
              parm[0].head = 0;
            }

          for (n = 0; n < NTHREADS; n++)
            close(parm[n].fid);
        }
    }

    //  Free everything as a matter of good form

    { int f, n;

      for (f = 0; f < argc-1; f++)
        { free(fobj[f].fname);
          free(fobj[f].nparts);
          for (n = 0; n < NTHREADS; n++)
            free(fobj[f].parts[n]);
          free(fobj[f].parts);
        }
  
      for (n = 0; n < NTHREADS; n++)
        libdeflate_free_decompressor(parm[n].decomp);
      free(parm[0].out);
      free(parm[0].buf);
    }
  }

  exit (0);
}
