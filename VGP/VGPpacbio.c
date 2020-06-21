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
#include <stdbool.h>

#include "gene_core.h"
#include "pb_expr.h"
#include "../Core/ONElib.h"

#include "LIBDEFLATE/libdeflate.h"

#include "VGPschema.h"

typedef  struct libdeflate_decompressor Deflator;

#undef   DEBUG_CORE
#undef   DEBUG_FIND
#undef   DEBUG_RECORDS
#undef   DEBUG_OUT

#define IO_BLOCK  10000000
#define BAM_BLOCK  0x10000
#define HEADER_LEN      36
#define SEQ_RUN         40

static int     VERBOSE;
static int     NTHREADS;
static int     QUALITY;         //  Output quality vector
static int     ARROW;           //  Output arrow/A lines?
static Filter *EXPR;            //  Filter expression
static int     EFLAGS;          //  Fields needed in filter expression
static int     IS_BIG_ENDIAN;   //  Is machine big-endian?

static char *Usage = "[-vaq] [-e<expr(ln>=500 && rq>=750)> [-T<int(4)>] <input:pacbio> ...";

typedef struct
  { char      *fname;   //  Full path name of file
    int64      fsize;   //  size of file
    int        isbam;   //  bam or sam?
  } File_Object;

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
  { int64  fpos;   //  offset in compressed file of desired bam block
    uint32 boff;   //  offset within uncompressed bam block of desired record
  } Location;

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
    Deflator *decomp;
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
  { File_Object *fobj;     //  Reference to array of file objects
    uint8       *buf;      //  IO buffer for this thread
    Deflator    *decomp;   //  Decompression codec (LIBDEFLATE)
    OneFile     *vf;       //  OneFile for output
    BAM_FILE     bam;
    int          fid;      //  fid of fobj[bidx].fname
    int          bidx;     //  Scan range is [bidx:beg,eidx:end)
    Location     beg;
    int          eidx;
    Location     end;
    int          error;
  } Thread_Arg;

  //  Find first record location (skip over header) in parm->fid
  //    Return value is in parm->beg

static void skip_header(Thread_Arg *parm)
{ uint8    *buf = parm->buf;
  BAM_FILE *bam = &(parm->bam);
  int       fid = parm->fid;

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

static void find_nearest(Thread_Arg *parm)
{ uint8       *buf  = parm->buf;
  BAM_FILE    *bam  = &(parm->bam);
  int          fid  = parm->fid;
  Deflator    *decomp = parm->decomp;
  int64        fpos = parm->beg.fpos;

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
  BAM_FILE    *bam  = &(parm->bam);
  int          fid  = parm->fid;

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
 
static int bam_record_scan(BAM_FILE *sf, samRecord *theR)
{ int ldata, lname, lseq, lcigar, aux;

  { uint8  x[36];     //  Process 36 byte header
    char  *eoh;

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
    theR->len = lseq;
    eoh = index(theR->header,'/');
    if (eoh != NULL)
      *eoh = 0;

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

      if (QUALITY)
        { t += e;
          s  = theR->qvs;
          if (t[0] == 0xff)
            return (2);
          for (i = 0; i < lseq; i++)
            s[i] = t[i] + 33;
        }
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
    int      i, alines, contents;

    contents = 0;

    arr = theR->arr;
    p = theR->data + aux;
    e = theR->data + ldata;

    if (IS_BIG_ENDIAN)
      flip_auxilliary(p,e);

    alines = 0;
    while (p < e)
      { switch (*p)
        { case 's':
            if (ARROW && memcmp(p+1,"nBf",3) == 0)
              { len = *((int *) (p+4));
                if (len != 4)
                  { fprintf(stderr,"%s: sn-tag does not have 4 floats\n",Prog_Name);
                    exit (1);
                  }
                theR->snr[0] = *((float *) (p+8));
                theR->snr[1] = *((float *) (p+12));
                theR->snr[2] = *((float *) (p+16));
                theR->snr[3] = *((float *) (p+20));
                p += 24;
                alines += 1;
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
                p += 8;
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
                alines += 1;
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
                contents |= HAS_BC;
                continue;
              }
            else if (p[1] == 'q')
              { GET_UINT("bq",theR->bqual)
                contents |= HAS_BQ;
                continue;
              }
            break;
          case 'z':
            if (p[1] == 'm')
              { GET_UINT("zm",theR->well)
                contents |= HAS_ZM;
                continue;
              }
            break;
          case 'q':
            if (p[1] == 's')
              { GET_UINT("qs",theR->beg)
                contents |= HAS_QS;
                continue;
              }
            else if (p[1] == 'e')
              { GET_UINT("qe",theR->end)
                contents |= HAS_QE;
                continue;
              }
            break;
          case 'r':
            if (memcmp(p+1,"qf",2) == 0)
              { theR->qual = *((float *) (p+3));
                p += 7;
                contents |= HAS_RQ;
                continue;
              }
            break;
          case 'n':
            if (p[1] == 'p')
              { GET_UINT("np",theR->nump)
                contents |= HAS_NP;
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
    if (ARROW && alines < 2)
      return (1);
    theR->defined = contents;
  }

  return (0);
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
    q = index(q,'/');         // Truncate pacbio well & pulse numbers
    if (q != NULL && q < p)
      *q = 0;

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

    if (QUALITY)
      { if (*q == '*')
          return (2);
        qlen = p-q;
        seq = theR->qvs;
        for (i = 0; i < qlen; i++)
          seq[i] = *q++;
      }
  }

  { char *q, *arr;       //  Get zm, qs, qe, rq, sn, and pw from auxilliary tags
    int   x, cnt, aline;
    int   contents;

    contents = 0;
    aline = 0;
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
                aline += 1;
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
                for (cnt = 0; *p == ',' && cnt < qlen; cnt++)
                  { x = strtol(q=p+1,&p,0);
                    CHECK( p == q, "Cannot parse pulse width value")
                    if (x >= 5)
                      x = 4;
                    arr[cnt] = x + '0';
                  }
                CHECK ( *p == ',' || cnt < qlen, "pulse width arraw has different length than read")
                CHECK ( *p != '\t' && *p != '\n', "Cannot parse pulse width values")
                aline += 1;
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
                contents |= HAS_BC;
                continue;
              }
            else if (strncmp(p+1,"q:",2) == 0 && is_integer[(int) p[3]] && p[4] == ':')
              { theR->bqual = strtol(p+5,&q,0);
                CHECK (p+5 == q, "Could not parse integer barcode quality")
                p = q;
                contents |= HAS_BQ;
                continue;
              }
            break;
          case 'z':
            if (strncmp(p+1,"m:",2) == 0 && is_integer[(int) p[3]] && p[4] == ':')
              { theR->well = strtol(p+5,&q,0);
                CHECK (p+5 == q, "Could not parse integer well number")
                p = q;
                contents |= HAS_ZM;
                continue;
              }
            break;
          case 'q':
            if (strncmp(p+1,"s:",2) == 0 && is_integer[(int) p[3]] && p[4] == ':')
              { theR->beg = strtol(p+5,&q,0);
                CHECK (p+5 == q, "Could not parse integer start pulse")
                p = q;
                contents |= HAS_QS;
                continue;
              }
            else if (strncmp(p+1,"e:",2) == 0 && is_integer[(int) p[3]] && p[4] == ':')
              { theR->end = strtol(p+5,&q,0);
                CHECK (p+5 == q, "Could not parse integer end pulse")
                p = q;
                contents |= HAS_QE;
                continue;
              }
            break;
          case 'r':
            if (strncmp(p+1,"q:f:",4) == 0)
              { theR->qual = strtod(p+5,&q);
                CHECK (p+5 == q, "Could not parse floating point quality value")
                p = q;
                contents |= HAS_RQ;
                continue;
              }
            break;
          case 'n':
            if (strncmp(p+1,"np:",3) == 0 && is_integer[(int) p[3]] && p[4] == ':')
              { theR->nump = strtol(p+5,&q,0);
                CHECK (p+5 == q, "Could not parse integer number of passes")
                p = q;
                contents |= HAS_NP;
                continue;
              }
            break;
          default:
            break;
        }

        while (*p != '\t' && *p != '\n')
          p += 1;
      }

    if (ARROW && aline < 2)
      return (1);
    theR->defined = contents;
  }

  return (0);
}

/*******************************************************************************************
 *
 *  Parallel:  Each thread processes a contiguous stripe across the input files
 *               sending the compressed binary data lines to their assigned OneFile.
 *
 ********************************************************************************************/

  //  Write subread data in samRecord rec to non-NULL file types

static void *output_thread(void *arg)
{ Thread_Arg  *parm  = (Thread_Arg *) arg;
  File_Object *fobj  = parm->fobj;
  uint8       *buf   = parm->buf;
  BAM_FILE    *bam   = &(parm->bam);
  OneFile     *vf    = parm->vf;

  samRecord    _theR, *theR = &_theR;

  int64        epos;
  uint32       eoff;
  int          head, isbam;
  int          f, fid;
  int          error;

  //  Know the max size of sequence and data from pass 1, so set up accordingly

  theR->dmax = 50000;
  theR->data = Malloc(theR->dmax,"Allocating sequence array");
  theR->lmax = 75000;
  theR->seq  = Malloc(3*theR->lmax,"Allocating sequence array");
  if (theR->seq == NULL || theR->data == NULL)
    exit (1);
  theR->arr = theR->seq + theR->lmax;
  theR->qvs = theR->arr + theR->lmax;

  bam->decomp = parm->decomp;

  parm->error = 0;

  for (f = parm->bidx; f <= parm->eidx; f++)
    { fid   = open(fobj[f].fname,O_RDONLY);
      isbam = fobj[f].isbam;
      if (f < parm->eidx)
        { epos = fobj[f].fsize;
          eoff = 0;
        }
      else 
        { epos = parm->end.fpos;
          eoff = parm->end.boff;
        }
      if (f > parm->bidx || parm->beg.fpos == 0)
        { head = 1;
          parm->beg.fpos = 0;
          parm->fid = fid;
          if (isbam)
            skip_header(parm);
          else
            sam_nearest(parm);
        }
      else
        head = 0;
      if (isbam)
        bam_start(bam,fid,buf,&(parm->beg));
      else
        sam_start(bam,fid,buf,&(parm->beg));

#ifdef DEBUG_OUT
      fprintf(stderr,"Block: %12lld / %5d to %12lld / %5d --> %8lld\n",bam->loc.fpos,bam->loc.boff,
                                                               epos,eoff,epos - bam->loc.fpos);
      fflush(stderr);
#endif

      while (bam->loc.fpos != epos || bam->loc.boff != eoff)
        { if (isbam)
            { error = bam_record_scan(bam,theR);
              if (error)
                { parm->error = error;
                  return (NULL);
                }
            }
          else
            { error = sam_record_scan(bam,theR);
              if (error)
                { parm->error = error;
                  return (NULL);
                }
            }
          if (EFLAGS & ~ theR->defined)
            { parm->error = 4;
              return (NULL);
            }

#ifdef DEBUG_RECORDS
          fprintf(stderr,"S = '%s'\n",theR->seq);
          if (ARROW)
            fprintf(stderr,"A = '%.*s'\n",theR->len,theR->arr);
          if (QUALITY)
            fprintf(stderr,"Q = '%.*s'\n",theR->len,theR->qvs);
          fprintf(stderr,"zm = %d\n",theR->well);
          fprintf(stderr,"ln = %d\n",theR->len);
          fprintf(stderr,"rq = %g\n",theR->qual);
          fprintf(stderr,"bc1 = %d\n",theR->bc[0]);
          fprintf(stderr,"bc2 = %d\n",theR->bc[1]);
          fprintf(stderr,"bq = %d\n",theR->bqual);
          fprintf(stderr,"np = %d\n",theR->nump);
          fprintf(stderr,"qs = %d\n",theR->beg);
          fprintf(stderr,"qe = %d\n",theR->end);
#endif

          if (head)
            { int len = strlen(theR->header);
    
              oneInt(vf,0) = 0;     //  don't actually know, OK for binary
              oneInt(vf,1) = len;
              oneWriteLine(vf,'g',len,theR->header);
              head = 0;
            }

          if (theR->len <= 0 || ! evaluate_bam_filter(EXPR,theR))
            continue;
    
          oneInt(vf,0) = theR->len;
          oneWriteLine(vf,'S',theR->len,theR->seq);
    
          oneInt(vf,0) = theR->well;
          oneInt(vf,1) = theR->beg;
          oneInt(vf,2) = theR->end;
          oneReal(vf,3) = theR->qual;
          oneWriteLine(vf,'W',0,NULL);

          if (QUALITY)
            { oneInt(vf,0) = theR->len;
              oneWriteLine(vf,'Q',theR->len,theR->qvs);
            }
    
          if (ARROW)
            { oneInt(vf,0) = theR->len;
              oneWriteLine(vf,'A',theR->len,theR->arr);

              oneReal(vf,0) = theR->snr[0];
              oneReal(vf,1) = theR->snr[1];
              oneReal(vf,2) = theR->snr[2];
              oneReal(vf,3) = theR->snr[3];
              oneWriteLine(vf,'N',0,NULL);
            }
        }

      close(fid);
    }

  free(theR->data);
  free(theR->seq);

  return (NULL);
}


/*******************************************************************************************
 *
 *  Main routine: parameter, 2 passes to stdout. 
 *
 *******************************************************************************************/

  //  Main

int main(int argc, char* argv[])
{ OneSchema *schema;
  char      *command; 

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
            ARG_FLAGS("vaq")
            break;
          case 'e':
            EXPR = parse_filter(argv[i]+2,&EFLAGS);
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
    QUALITY = flags['q'];

    if (EXPR == NULL)
      EXPR = parse_filter("ln>=500 && rq>=750",&EFLAGS);
    EFLAGS |= HAS_ZM | HAS_QS | HAS_QE | HAS_RQ;

    if (argc == 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: verbose mode, output progress as proceed\n");
        fprintf(stderr,"      -a: extract Arrow information on N- and A-lines.\n");
        fprintf(stderr,"      -q: extract QV information on a Q-line.\n");
        fprintf(stderr,"      -T: Number of threads to use\n");
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

  { File_Object fobj[10];
    Thread_Arg  parm[NTHREADS];

#if ! defined(DEBUG_OUT)
    pthread_t   threads[NTHREADS];
#endif

    if (VERBOSE)
      { fprintf(stderr,"  Partitioning %d .subreads files into %d parts\n",argc-1,NTHREADS);
        fflush(stderr);
      }

    { uint8 *bf;
      int    f, i;
      int64  b, work, wper;

      //  Allocate IO buffer space for threads

      bf = Malloc(NTHREADS*IO_BLOCK,"Allocating IO_Buffer\n");
      if (bf == NULL)
        exit (1);

      //  Get name, size, and type (sam or bam) of each file in 'fobj[]'

      work = 0;
      for (f = 0; f < argc-1; f++)
        { struct stat stats;
          char *suffix[6] = { ".ccs.bam", ".ccs.sam", ".subreads.bam", ".subreads.sam",
                              ".bam", ".sam" };
          char *pwd, *root;
          int   fid;

          pwd = PathTo(argv[f+1]);
          OPEN2(argv[f+1],pwd,root,fid,suffix,6);
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
          fobj[f].isbam  = (i % 2 == 0);

          work += stats.st_size;

          free(pwd);
          free(root);
          close(fid);
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
          parm[i].decomp = libdeflate_alloc_decompressor();

          if (b != 0)
            { parm[i].fid      = open(fobj[f].fname,O_RDONLY);
              parm[i].beg.fpos = 0;
              if (fobj[f].isbam)
                skip_header(parm+i);
              else
                sam_nearest(parm+i);
              parm[i].end = parm[i].beg;
            }
          else
            parm[i].end.fpos = 0;

#ifdef DEBUG_FIND
          fprintf(stderr," %2d: %1d %10lld (%10lld,%d)\n",i,f,b,parm[i].end.fpos,parm[i].end.boff);
          fflush(stderr);
#endif

          parm[i].beg.fpos = b;
          parm[i].bidx     = f;

          work -= wper;
          while (work < 2*BAM_BLOCK)
            { if (f == argc-2)
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

    { int i, f;

      //  For each non-zero start point find synchronization point in
      //    bam/sam file.  If can't find it then start at beginning of
      //    next file, and if at or before first data line then signal
      //    start at beginning by zero'ing the synch point.

      for (i = 0; i < NTHREADS; i++)
        if (parm[i].beg.fpos != 0)
          { if (fobj[parm[i].bidx].isbam)
              find_nearest(parm+i);
            else
              sam_nearest(parm+i);
            if (parm[i].beg.fpos < 0)
              { parm[i].beg.fpos = 0;
                parm[i].bidx    += 1;
              }
            else if (parm[i].beg.fpos <= parm[i].end.fpos)
              parm[i].beg.fpos = 0;
            close(parm[i].fid);
          }

      //  Paranoid: if one thread's synch point overtakes the next one (will almost
      //    certainly never happen unless files very small and threads very large),
      //    remove the redundant threads.

      f = 0;
      for (i = 1; i < NTHREADS; i++)
        if (parm[i].bidx > parm[f].bidx || parm[i].beg.fpos > parm[f].beg.fpos
            || parm[i].beg.boff > parm[f].beg.boff)
          parm[++f] = parm[i];
        else
          libdeflate_free_decompressor(parm[i].decomp);
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
      for (i = 0; i < NTHREADS; i++)
        { fprintf(stderr," %2d: %2d / %12lld / %5d",
                         i,parm[i].bidx,parm[i].beg.fpos,parm[i].beg.boff);
          fprintf(stderr,"  -  %2d / %12lld / %5d\n",
                         parm[i].eidx,parm[i].end.fpos,parm[i].end.boff);
        }
      fflush(stderr);
#endif
    }

    //  Setup a OneFile for each thread, put the header in the first one

    { OneFile *vf;
      int      i, error;

      vf = oneFileOpenWriteNew("-",schema,"pbr",true,NTHREADS);
      oneAddProvenance(vf,Prog_Name,"1.0",command,NULL);
      oneWriteHeader(vf);
#ifdef DEBUG_OUT
      fprintf(stderr,"Opened\n");
      fflush(stderr);
#endif

      if (VERBOSE)
        { fprintf(stderr,"  Producing .pbr segements in parallel\n");
          fflush(stderr);
        }

      //  Generate the data lines in parallel threads

#ifdef DEBUG_OUT
      for (i = 0; i < NTHREADS; i++)
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
        { fprintf(stderr,"  Cat'ing .pbr segments\n");
          fflush(stderr);
        }

      oneFileClose(vf);

      //  If asked for arrow / qv vectors but not in input, then error will be 1 / 2.

      error = 0;
      for (i = 0; i < NTHREADS; i++)
        error |= parm[i].error;

      if (error)
        { if (error >= 4)
            fprintf(stderr,"%s: Bam file does not have auxiliary info of a PacBio file\n",
                           Prog_Name);
          else if (error == 1)
            fprintf(stderr,"%s: Bam file does not contain pulse information for -a option\n",
                           Prog_Name);
          else if (error == 2)
            fprintf(stderr,"%s: Bam file does not contain qv information for -q option\n",
                           Prog_Name);
          else
            fprintf(stderr,"%s: Bam file does not contain the information for -a & -q options\n",
                           Prog_Name);
          exit (1);
        }
    }

    //  Free everything as a matter of good form

    { int f, n;

      for (f = 0; f < argc-1; f++)
        free(fobj[f].fname);

      for (n = 0; n < NTHREADS; n++)
        libdeflate_free_decompressor(parm[n].decomp);
      free(parm[0].buf);

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
