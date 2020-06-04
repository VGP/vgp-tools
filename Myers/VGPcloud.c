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
#include <stdbool.h>

#include "gene_core.h"
#include "msd.sort.h"
#include "../Core/ONElib.h"

#include "VGPschema.h"

#undef DEBUG

int    VERBOSE;        //  Verbose mode?
int    NTHREADS;       //  # of threads to use for parallized sorts
int    HISTOGRAM;      //  Histogram bar-code counts
int    VALID;          //  Threshold above which bar-code is considered valid

static char *Usage = "[-vH] [-t<int(100)>] [-T<int(4)>] <clouds:irp>";

static uint8 bit[8] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };


/****************************************************************************************
 *
 *  Fixed length compression / decompression for Q-strings & decompression for DNA S-strings
 *
 *****************************************************************************************/

static uint8 *Compress_QV(uint8 *trg, int len, uint8 *qvec, int qbits, uint8 *map)
{ uint8  x;
  int    i, rem;

  if (qbits == 8)
    { memcpy(trg,qvec,len);
      return (trg+len);
    }

  *trg = map[qvec[0]];
  rem = qbits;
  for (i = 1; i < len; i++)
    { x = map[qvec[i]];
      *trg |= (x << rem);
      rem += qbits;
      if (rem >= 8)
        { rem -= 8;
          if (rem > 0)
            *++trg = (x >> (qbits-rem));
          else
            *++trg = 0;
        }
    }
  if (rem == 0)
    return (trg);
  else
    return (trg+1);
}

static uint8 *Uncompress_QV(uint8 *src, int len, char *qvec, int qbits, uint8 qmask, uint8 *inv)
{ uint8  x, v;
  int    i, rem;

  if (qbits == 8)
    { memcpy(qvec,src,len);
      return (src+len);
    }

  x = *src;
  rem = 0;
  for (i = 0; i < len; i++)
    { v = x & qmask;
      x >>= qbits;
      rem += qbits;
      if (rem >= 8)
        { rem -= 8;
          x = *++src;
          if (rem > 0)
            { v |= (x << (qbits-rem)) & qmask;
              x >>= rem;
            }
        }
      qvec[i] = inv[v];
    }
  if (rem == 0)
    return (src);
  else
    return (src+1);
}

  //  Uncompress read from 2-bits per base into [0-3] per byte representation

static char Base[4] = { 'a', 'c', 'g', 't' };

static uint8 *Uncompress_SEQ(uint8 *src, char *s, int len)
{ int   i, tlen, byte;
  char *t0, *t1, *t2, *t3;

  t0 = s;
  t1 = t0+1;
  t2 = t1+1;
  t3 = t2+1;

  tlen = len-3;
  for (i = 0; i < tlen; i += 4)
    { byte = *src++;
      t0[i] = Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 4) & 0x3];
      t2[i] = Base[(byte >> 2) & 0x3];
      t3[i] = Base[byte & 0x3];
    }

  switch (i-tlen)
  { case 0:
      byte = *src++;
      t0[i] = Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 4) & 0x3];
      t2[i] = Base[(byte >> 2) & 0x3];
      break;
    case 1:
      byte = *src++;
      t0[i] = Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 4) & 0x3];
      break;
    case 2:
      byte = *src++;
      t0[i] = Base[(byte >> 6) & 0x3];
      break;
    default:
      break;
  }

  return (src);
}


/****************************************************************************************
 *
 *  Thread for building unsorted bar code list and collecting Q-vector symbol usage
 *
 *****************************************************************************************/

//  Optimized ReadLine optimized to read binary S-line without decompression

static char read_raw_seq(OneFile *vf)
{ U8        x;
  char      t;
  OneInfo  *li;

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

  li = vf->info['S'];

  { int64 nBits;

    if (x & 0x1)
      { nBits = (U8) getc(vf->f);
        if (fread(vf->codecBuf,((nBits+7) >> 3),1,vf->f) != 1)
          { fprintf(stderr,"%s: *FATAL* fail to read compressed fields\n",Prog_Name);
            exit (1);
          }
      }
    else
      { if (fread(vf->field,sizeof(OneField),1,vf->f) != (unsigned long) 1)
          { fprintf(stderr,"%s: *FATAL* fail to read fields\n",Prog_Name);
            exit (1);
          }
      }

    if (fread(&nBits,sizeof(I64),1,vf->f) != 1)
      { fprintf(stderr,"%s: *FATAL* fail to read list nBits\n",Prog_Name);
        exit (1);
      }
    if (fread(vf->codecBuf,((nBits+7) >> 3),1,vf->f) != 1)
      { fprintf(stderr,"%s: *FATAL* fail to read compressed list\n",Prog_Name);
        exit (1);
      }

    vf->field[0].i = (nBits >> 1);
  }

  { U8 peek = getc(vf->f);    // check if next line is a comment - if so then read it
    ungetc(peek,vf->f);
    if (peek & 0x80)
      peek = vf->binaryTypeUnpack[peek];
    if (peek == '/')
      { OneField keepField0 = vf->field[0];
        oneReadLine(vf);
        vf->lineType = t;
        vf->field[0] = keepField0;
      }
  }

  return (t);
}

typedef struct
  { OneFile     *vf;       //  OneFile for input
    int64        beg;      //  Range of reads to process
    int64        end;
    uint32      *codes;
    int64        usedqvs[256];
    int          error;
    int64        where;
    int          flen;
    int          rlen;
  } BarCode_Arg;

#define ERROR_PLINE 1
#define ERROR_SLINE 2
#define ERROR_SLEN  3
#define ERROR_QLINE 4
#define ERROR_QLEN  5

#define ERROR_SLINE2 6

static void error_report(int error, int64 line)
{ switch (error)
  { case ERROR_PLINE:
      fprintf(stderr,"%s: Expecting P-line, at or after S-line %lld\n",Prog_Name,line);
      break;
    case ERROR_SLINE:
    case ERROR_SLINE2:
      fprintf(stderr,"%s: Expecting S-line, at or after S-line %lld (%d)\n",Prog_Name,line,error);
      break;
    case ERROR_SLEN:
      fprintf(stderr,"%s: S-strings are not all the same size, at or after S-line %lld\n",
                     Prog_Name,line);
      break;
    case ERROR_QLINE:
      fprintf(stderr,"%s: Expecting Q-line, at or after S-line %lld\n",
                     Prog_Name,line);
      break;
    case ERROR_QLEN:
      fprintf(stderr,"%s: Q-string is not the same length as S-string, at or after S-line %lld\n",
                     Prog_Name,line);
      break;
  }
}

static void *barcodes_thread(void *arg)
{ BarCode_Arg *parm  = (BarCode_Arg *) arg;
  OneFile     *vf    = parm->vf;
  int64       *used  = parm->usedqvs;
  int64        beg   = 2*parm->beg;
  int64        end   = 2*parm->end;
  uint32      *codes = parm->codes;

  uint32 *fcmp;
  uint8  *fqvs;
  int64   o, coll;
  int     t, n;
  int     flen, rlen;

  fqvs = vf->info['Q']->buffer;
  fcmp = (uint32 *) (vf->codecBuf);

  bzero(used,sizeof(int64)*256);

  oneGotoObject(vf,beg);

  coll  = beg + 10000000/NTHREADS;
  rlen  = 0;
  flen  = 0;
  for (o = beg; o < end; o += 2)
    { t = read_raw_seq(vf);
      if (t != 'S')
        { parm->error = ERROR_SLINE;
          parm->where = o;
          return (NULL);
	}

      *codes++ = *fcmp;

      if (o > coll)
        { oneGotoObject(vf,o+2);
          continue;
        }

      n = oneInt(vf,0);
      if (flen == 0)
        flen = n;
      else if (flen != n)
        { parm->error = ERROR_SLEN;
          parm->where = o;
          return (NULL);
        }

      t = oneReadLine(vf);
      if (t != 'Q')
        { parm->error = ERROR_QLINE;
          parm->where = o;
          return (NULL);
        }
      if (oneInt(vf,0) != flen)
        { parm->error = ERROR_QLEN;
          parm->where = o;
          return (NULL);
        }

      for (n = 0; n < flen; n++)
        used[(int) fqvs[n]] += 1;

      t = read_raw_seq(vf);
      if (t != 'S')
        { parm->error = ERROR_SLINE;
          parm->where = o+1;
          return (NULL);
        }
      n = oneInt(vf,0);
      if (rlen == 0)
        rlen = n;
      else if (rlen != n)
        { parm->error = ERROR_SLEN;
          parm->where = o+1;
          return (NULL);
        }

      t = oneReadLine(vf);
      if (t != 'Q')
        { parm->error = ERROR_QLINE;
          parm->where = o+1;
          return (NULL);
        }
      if (oneInt(vf,0) != rlen)
        { parm->error = ERROR_QLEN;
          parm->where = o+1;
          return (NULL);
        }
      for (n = 0; n < rlen; n++)
        used[(int) fqvs[n]] += 1;

      t = oneReadLine(vf);
      if (t == 'g')
        t = oneReadLine(vf);

      if (t != 'P')
        { parm->error = ERROR_PLINE;
          parm->where = o+1;
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
 *  Threaded LSD radix sort and unique+count compression:
 *    Each thread sorts codes[beg..end) and then compacts to unique codes adjusting end.
 *    If count > 255 then multiple entries summing to count so that # of bar-codes with
 *    a given code is preserved.
 *
 *****************************************************************************************/

typedef struct
  { int          tid;      //  # of the thread (in [0,NTHREAD) )

    uint32      *codes;
    uint32      *vlist;
    uint8       *count;
    uint8       *ctype;
    uint8       *first;

    int64        beg;      //  Range of reads to process on input
    int64        end;      //  Indexes end of compacted list upon return

    int64       *finger;   // For the merge_thread
    int64        vbeg;
    int64        vlen;
    int64        vend;
    int64        nrepair;

    int64        hist10[101];
    int64        hist1[50];

    int64       *bucks;   //  To seed top level sort
  } Sort_Arg;

  //  Simple LSD radix sort of uint32's

#undef TEST_LSORT

void LSD_sort(int64 nelem, uint32 *src, uint32 *trg)
{ int64  ptrs[512];

  int64 *Cptr = ptrs;
  int64 *Nptr = ptrs + 256;

  int64    i;
  int      j, b;

  //  If first pass, then explicitly sweep to get Cptr counts

  { uint32  v;
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
    uint8  *p = ((uint8 *) &v);
#else
    uint8  *p = ((uint8 *) &v) + 3;
#endif

    for (j = 0; j < 256; j++)
      Cptr[j] = 0;

    for (i = 0; i < nelem; i++)
      { v = src[i];
        Cptr[*p] += 1;
      }
  }

  //  For each requested byte b in order, radix sort

  for (b = 0; b < 4; b++)
    { int Cbyte, Nbyte;

#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
      Cbyte  = b;
      Nbyte  = b+1;
#else
      Cbyte  = 3-b;
      Nbyte  = 2-b;
#endif

#ifdef TEST_LSORT
      printf("\nSorting byte %d (next %d)\n",Cbyte,Nbyte);
      fflush(stdout);
      printf("\nBUCKETS %d\n",Cbyte);
      for (j = 0; j < 256; j++)
        printf(" %3d: %10lld\n",j,Cptr[j]);
      fflush(stdout);
#endif

      //  Convert Cptr from counts to fingers

      { int64 x, y;

        x = 0;
        for (j = 0; j < 256; j++)
          { y = Cptr[j];
            Cptr[j] = x;
            x += y;
          }
      }

      //  Radix sort from src to trg

      { uint32  v;
        uint8  *p = ((uint8 *) &v) + Cbyte;
        uint8  *n = ((uint8 *) &v) + Nbyte;
        uint8   d;
        int64   x;

        if (b == 3)
          for (i = 0; i < nelem; i++)
            { v = src[i];
              d = *p;
              x = Cptr[d]++;
              trg[x] = v;
            }
        else
          { for (j = 0; j < 256; j++)
              Nptr[j] = 0;
            for (i = 0; i < nelem; i++)
              { v = src[i]; 
                d = *p;
                x = Cptr[d]++;
                trg[x] = v;
                Nptr[*n] += 1;
              }
          }
      }
    
      //  Flip roles of data and ptr vectors

      { int64  *p;
        uint32 *d;

        d   = src;
        src = trg;
        trg = d;

        p    = Cptr;
        Cptr = Nptr;
        Nptr = p;
      }

#ifdef TEST_LSORT
      { int64  c;

        printf("\nLSORT %d\n",Cbyte);
        for (c = 0; c < 300; c++)
          { for (j = 0; j < 4; j++)
              printf(" %02x",((uint8 *) (src+c))[j]);
            printf("\n");
          }
        printf("\n");
        for (c = nelem/2-150; c < nelem/2+150; c++)
          { for (j = 0; j < 4; j++)
              printf(" %02x",((uint8 *) (src+c))[j]);
            printf("\n");
          }
        printf("\n");
        for (c = nelem-300; c < nelem; c++)
          { for (j = 0; j < 4; j++)
              printf(" %02x",((uint8 *) (src+c))[j]);
            printf("\n");
          }
        printf("\n");

        for (c = 1; c < nelem; c++)
          { for (j = Cbyte; j >= 0; j--)
              if (((uint8 *) (src+c))[j] > ((uint8 *) (src+(c-1)))[j])
                break;
              else if (((uint8 *) (src+c))[j] < ((uint8 *) (src+(c-1)))[j])
                { printf("  Order: %lld",c);
                  for (j = Cbyte; j >= 0; j--)
                    printf(" %02x",((uint8 *) (src+(c+1)))[j]);
                  printf(" vs");
                  for (j = Cbyte; j >= 0; j--)
                    printf(" %02x",((uint8 *) (src+c))[j]);
                  printf("\n");
                  break;
                }
          }
      }
#endif
    }
}

static void *sort_thread(void *arg)
{ Sort_Arg *parm  = (Sort_Arg *) arg;
  int64     beg   = parm->beg;
  int64     end   = parm->end;
  uint32   *codes = parm->codes;
  uint32   *vlist = parm->vlist;
  uint8    *count = parm->count;

  int64  g, f, i;
  uint32 v;

  LSD_sort(end-beg,codes+beg,vlist+beg);

#ifdef CHECK_SORT
  for (i = beg+1; i < end; i++)
    if (codes[i-1] > codes[i])
      printf("Out of order %lld: %u vs %u\n",i,codes[i-1],codes[i]);
#endif

  codes[end] = codes[end-1]+1;  // != codes[end-1]

  f = beg;
  g = beg;
  v = codes[g];
  for (i = beg+1; i <= end; i++)
    if (codes[i] != v)
      { g = i-g;
#ifdef DEBUG_COMPRESS
        printf(" %8lld: %4lld %08x",i,g,v);
#endif
        while (g > 255)
          { count[f] = 255;
            codes[f] = v;
            f += 1;
            g -= 255;
          }
        count[f] = g;
        codes[f] = v;
        f += 1;

        g = i;
        v = codes[g];
      }

#ifdef CHECK_SORT
  { int64 sum;

    sum = count[beg];
    for (i = beg+1; i < f; i++)
      { if (codes[i-1] > codes[i] || (codes[i] == codes[i-1] && count[i-1] != 255))
          printf("Out of order %lld: %u vs %u\n",i,codes[i-1],codes[i]);
        sum += count[i];
      }
    if (sum != end-beg)
      printf("Lost counts\n");
  }
#endif

  count[f] = 0;
  parm->end = f;
  return (NULL);
}


/****************************************************************************************
 *
 *  Merge value range slabs of each of the NTHREAD segments into order into vlist[vbeg,vend)
 *
 *****************************************************************************************/

static int64 find(uint32 *codes, int64 n, uint32 x)
{ int64 l, r, m;

  // smallest k s.t. codes[k] >= x (or n if does not exist)

  l = 0;
  r = n;
  while (l < r)
    { m = ((l+r) >> 1);
      if (codes[m] < x)
        l = m+1;
      else
        r = m;
    }
#ifdef CHECK_SORT
  if (l > 0 && codes[l-1] >= x)
    printf("Find failed\n");
  if (l < n && codes[l] < x) 
    printf("Find failed\n");
#endif
  return (l);
}

  //  Heap sort of codes

static void reheap(int s, int64 *heap, uint32 *code, int hsize)
{ int    c, l, r;
  int64  hs, hl, hr;
  uint32 vs, vl, vr;

  vs = code[hs = heap[s]];
  c  = s;
  while ((l = 2*c) <= hsize)
    { r  = l+1;
      vl = code[hl = heap[l]];
      if (r > hsize)
        vr = 0xffffffffu;
      else
        vr = code[hr = heap[r]];
      if (vr >= vl)
        { if (vs > vl)
            { heap[c] = hl;
              c = l;
            }
          else
            break;
        }
      else
        { if (vs > vr)
            { heap[c] = hr;
              c = r;
            }
          else
            break;
        }
    }
  if (c != s)
    heap[c] = hs;
}

static void *merge_thread(void *arg)
{ Sort_Arg *parm   = (Sort_Arg *) arg;
  uint32   *codes  = parm->codes;
  uint32   *vlist  = parm->vlist;
  uint8    *count  = parm->count;
  uint8    *ctype  = parm->ctype;
  uint8    *first  = parm->first;
  int64    *heap   = parm->finger-1;
  int64     vbeg   = parm->vbeg; 
  int64     vlen   = parm->vlen;
  int64    *hist10 = parm->hist10;
  int64    *hist1  = parm->hist1;

  int    v, cum, hsize;
  int64  f, i, p;
  uint32 bar, x;
#ifdef CHECK_SORT
  uint32 a;
#endif
  int64 s, e;

  s = (0x100000000ll * parm->tid) / NTHREADS;
  e = (0x100000000ll * (parm->tid+1)) / NTHREADS;
  bzero(ctype+s,e-s);

  if (HISTOGRAM)
    { bzero(hist10,101*sizeof(int64));
      bzero(hist1,50*sizeof(int64));
    }

  hsize = NTHREADS;
  for (i = hsize/2; i >= 1; i--)
    reheap(i,heap,codes,hsize);

  f = vbeg;
  bar = codes[heap[1]] + 1;
  cum = 0;
#ifdef CHECK_SORT
  a = codes[heap[1]];
  if (a < s)
    printf("Merge start not in range\n");
#endif
  for (i = 0; i < vlen; i++)
    { p = heap[1];
      v = count[p];
      while (v == 0)
        { heap[1] = heap[hsize];
          hsize -= 1;
          reheap(1,heap,codes,hsize);
          p = heap[1];
          v = count[p];
        }
      x = codes[p];
#ifdef CHECK_SORT
      if (x < bar && i > 0)
        printf("Out of order heap\n");
#endif
      if (x != bar)
        { if (cum >= VALID)
            { vlist[f++] = bar;
              ctype[bar] = 0x40;
            }
          if (HISTOGRAM)
            { if (cum >= 1000)
                hist10[100] += 1;
              else if (cum >= 50)
                hist10[cum/10] += 1;
              else
                hist1[cum] += 1;
            }
          bar = x;
          cum = v;
        }
      else
        cum += v;
      heap[1] = p+1;
      reheap(1,heap,codes,hsize);
    }
  if (cum >= VALID)
    { vlist[f++] = bar;
      ctype[bar] = 0x40;
      if (HISTOGRAM)
        { if (cum >= 1000)
            hist10[100] += 1;
          else if (cum >= 50)
            hist10[cum/10] += 1;
          else
            hist1[cum] += 1;
        }
    }

#ifdef CHECK_SORT
  if (bar >= e)
    printf("Merge end beyond range\n");
  for (i = 1; i <= hsize; i++)
    if (count[heap[i]] > 0 && codes[heap[i]] < e)
      printf("Heap is not exhausted\n");
  printf("Code Range: %08x - %08x [%08llx - %08llx] (%lld / %lld / %d)\n",
          a,bar,s,e,vlen,f-vbeg,hsize);
#endif

  bzero(first,0x20000000ll);

  parm->vend  = f;

  return (NULL);
}


/****************************************************************************************
 *
 *  Threads to mark 1-neighborhood of every valid code in ctype & thread bit vector first
 *
 *****************************************************************************************/

static void *mark_thread(void *arg)
{ Sort_Arg *parm  = (Sort_Arg *) arg;
  int64     vbeg  = parm->vbeg;
  int64     vend  = parm->vend;
  uint8    *ctype = parm->ctype;
  uint8    *first = parm->first;
  uint32   *vlist = parm->vlist;

  int64  i;
  int    j, m;
  uint32 v;
  uint32 var, v3, vb;
  uint32 chg, msk, u, t;

  // Precondition: ctype[v] & 0xc0 == 0x00 if not valid, 0x40 if valid, first_t[v] = 0 for all t, v

  for (i = vbeg; i < vend; i++)
    { v = vlist[i];
      chg = 0x1;
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
                { var = v | u;
#ifdef DEBUG_MARK
                  printf("   %08x (%d/%d)",var,j,m);
#endif
                  v3 = (var >> 3);
                  vb = bit[var & 0x7];
                  if (first[v3] & vb)
                    { if ((ctype[var] & 0xc0) == 0x00)
                        { ctype[var] |= 0xc0;
#ifdef DEBUG_MARK
                          printf(" V");
#endif
                        }
                    }
                  else
                    { first[v3] |= vb;
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
          v |= t;
        }
    }
    
  // ctype[v] & 0xc0 == 0x40 if valid, 0xc0 if not valid & touched >= 2x by a thread, 0x00 otherwise
  // first_t[v] = 1 iff thread t touched v at least once

  return (NULL);
}


/****************************************************************************************
 *
 *  Threads to examine marks of the previous pass and identify uniquely correctable codes
 *    and their correction (position,symbol)
 *
 *****************************************************************************************/

static void *fix_thread(void *arg)
{ Sort_Arg *parm  = (Sort_Arg *) arg;
  int       vbeg  = parm->vbeg;
  int       vend  = parm->vend;
  uint8    *ctype = parm->ctype;
  uint8    *first = parm->first;
  uint32   *vlist = parm->vlist;

  uint8 *vec[NTHREADS];

  int64  i;
  int    m, n, cnt;
  uint32 v, var, v3, vb;
  uint32 chg, msk, u, j, t, c;

  for (n = 0; n < NTHREADS; n++)
    vec[n] = first + n*0x20000000ll;

  // ctype[v] & 0xc0 == 0x40 if valid, 0xc0 if not valid & touched >=2x by a thread, 0x00 otherwise
  // bit_t[v] = 1 iff thread t touched v at least once

  for (i = vbeg; i < vend; i++)
    { v = vlist[i];
      chg = 0x1;
      msk = 0x3;
#ifdef DEBUG_FIXES
      printf("%08x\n",v);
#endif
      for (j = 0; j < 32; j += 2)
        { t = (v & msk); 
          v -= t;
          c = (t >> j);
          u = 0;
          for (m = 0; m < 4; m++)
            { if (t != u)
                { var = v | u;
#ifdef DEBUG_FIXES
                  printf("   %08x (%d/%d)\n",var,j,m);
#endif
                  v3 = (var >> 3);
                  vb = bit[var & 0x7];
                  if ((ctype[var] & 0xc0) == 0x00)
                    { cnt = 0;
                      for (n = 0; n < NTHREADS; n++)
                        if (vec[n][v3] & vb)
                          cnt += 1;
                      if (cnt == 1)
                        { ctype[var] = (0x80 | (j << 1) | c);
#ifdef DEBUG_FIXES
                          printf("      OK %x  %d %d\n",ctype[var],c,j);
#endif
                        }
                      else
                        { ctype[var] |= 0xc0;
#ifdef DEBUG_FIXES
                          printf("      C-%d\n",cnt);
#endif
                        }
                    }
#ifdef DEBUG_FIXES
                  else
                    printf("      Z %d\n",ctype[var]>>6);
#endif
                }
              u += chg;
            }
          chg <<= 2;
          msk <<= 2;
          v |= t;
        }
    }

  // ctype[v] & 0xc0 = 0x00  not validd & never touched
  //                   0x40  valid
  //                   0x80  not valid & touched once => uniquely correctable (with correction)
  //                   0xc0  not valid & touched 2 or more times

  return (NULL);
}


/****************************************************************************************
 *
 *  Threads to count first byte frequencies of all valid & correctable codes
 *
 *****************************************************************************************/

static uint32 off[64] = { 0xfffffffc, 0xfffffffc, 0xfffffffc, 0xfffffffc,
                          0xfffffff3, 0xfffffff3, 0xfffffff3, 0xfffffff3,
                          0xffffffcf, 0xffffffcf, 0xffffffcf, 0xffffffcf, 
                          0xffffff3f, 0xffffff3f, 0xffffff3f, 0xffffff3f, 

                          0xfffffcff, 0xfffffcff, 0xfffffcff, 0xfffffcff, 
                          0xfffff3ff, 0xfffff3ff, 0xfffff3ff, 0xfffff3ff, 
                          0xffffcfff, 0xffffcfff, 0xffffcfff, 0xffffcfff, 
                          0xffff3fff, 0xffff3fff, 0xffff3fff, 0xffff3fff, 

                          0xfffcffff, 0xfffcffff, 0xfffcffff, 0xfffcffff,
                          0xfff3ffff, 0xfff3ffff, 0xfff3ffff, 0xfff3ffff,
                          0xffcfffff, 0xffcfffff, 0xffcfffff, 0xffcfffff, 
                          0xff3fffff, 0xff3fffff, 0xff3fffff, 0xff3fffff, 

                          0xfcffffff, 0xfcffffff, 0xfcffffff, 0xfcffffff, 
                          0xf3ffffff, 0xf3ffffff, 0xf3ffffff, 0xf3ffffff, 
                          0xcfffffff, 0xcfffffff, 0xcfffffff, 0xcfffffff, 
			  0x3fffffff, 0x3fffffff, 0x3fffffff, 0x3fffffff
                        };

static uint32 sym[64] = { 0x00000000, 0x00000001, 0x00000002, 0x00000003,
                          0x00000000, 0x00000004, 0x00000008, 0x0000000c,
                          0x00000000, 0x00000010, 0x00000020, 0x00000030,
                          0x00000000, 0x00000040, 0x00000080, 0x000000c0,

                          0x00000000, 0x00000100, 0x00000200, 0x00000300,
                          0x00000000, 0x00000400, 0x00000800, 0x00000c00,
                          0x00000000, 0x00001000, 0x00002000, 0x00003000,
                          0x00000000, 0x00004000, 0x00008000, 0x0000c000,

                          0x00000000, 0x00010000, 0x00020000, 0x00030000,
                          0x00000000, 0x00040000, 0x00080000, 0x000c0000,
                          0x00000000, 0x00100000, 0x00200000, 0x00300000,
                          0x00000000, 0x00400000, 0x00800000, 0x00c00000,

                          0x00000000, 0x01000000, 0x02000000, 0x03000000,
                          0x00000000, 0x04000000, 0x08000000, 0x0c000000,
                          0x00000000, 0x10000000, 0x20000000, 0x30000000,
			  0x00000000, 0x40000000, 0x80000000, 0xc0000000
                        };

static void *bucket_thread(void *arg)
{ Sort_Arg *parm   = (Sort_Arg *) arg;
  int       beg    = parm->beg;
  int       end    = parm->end;
  uint8    *ctype  = parm->ctype;
  uint8    *count  = parm->count;
  uint32   *codes  = parm->codes;
  int64    *bucks  = parm->bucks;

  int64  i, nrepair;
  int    p;
  uint8  x, t;
  uint32 v;
 
  nrepair = 0;
  for (p = 0; p < 256; p++)
    bucks[p] = 0;
  for (i = beg; i < end; i++)
    { v = codes[i];
      x = ctype[v];
      t = (x & 0xc0);
      if (t == 0x80)
        { p = (x & 0x3f);
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
          if (p < 16)
            v = (v & off[p]) | sym[p];
          bucks[v&0xff] += count[i];
#else
          if (p >= 48)
            v = (v & off[p]) | sym[p];
          bucks[v>>24] += count[i];
#endif
          nrepair += count[i];
        }
      else if (t == 0x40)
#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
        bucks[v&0xff] += count[i];
#else
        bucks[v>>24] += count[i];
#endif
    }

  parm->nrepair = nrepair;

  return (NULL);
}


/****************************************************************************************
 *
 *  Thread to load pairs into array sorted on highest order byte
 *
 *****************************************************************************************/

typedef struct
  { OneFile     *vf;       //  OneFile for input
    int64        beg;      //  Range of reads to process
    int64        end;
    int          qbits;
    uint8       *map;
    int          flen, fclen, fqlen;
    int          rlen, rclen, rqlen;
    uint8       *ctype;
    uint8       *array;
    int64       *bptr;
  } Load_Arg;

static void *load_thread(void *arg)
{ Load_Arg *parm  = (Load_Arg *) arg;
  OneFile  *vf    = parm->vf;
  int64     beg   = 2*parm->beg;
  int64     end   = 2*parm->end;
  uint8    *ctype = parm->ctype;
  int       qbits = parm->qbits;
  uint8    *map   = parm->map;
  int       flen  = parm->flen;
  int       rlen  = parm->rlen;
  int       fclen = parm->fclen;
  int       rclen = parm->rclen;
  uint8    *array = parm->array;
  int64    *bptr  = parm->bptr;

  uint32 *dcode;
  uint8  *bcode, *fqvs;
  int64   o;
  uint32  v, c, p, bar, byte;
  uint8  *aptr;
  int     t;

  fqvs  = vf->info['Q']->buffer;
  bcode = (uint8 *)  (vf->codecBuf);
  dcode = (uint32 *) (vf->codecBuf);

  oneGotoObject(vf,beg);

  o = beg;
  while (1)
    { read_raw_seq(vf);

      bar = *dcode;

      v   = ctype[bar];
      c   = (v & 0xc0);
      if (c == 0x80)
        { p = (v & 0x3f);
          *dcode = (bar & off[p]) | sym[p];
          c = 0x40;
        }
      if (c != 0x40)
        { o += 2;
          if (o >= end)
            break;
          oneGotoObject(vf,o);
          continue;
        }

      byte = *bcode;
      aptr = array + bptr[byte];

      memcpy(aptr,bcode,fclen);
      aptr += fclen;

      oneReadLine(vf);
      aptr = Compress_QV(aptr,flen,fqvs,qbits,map);

      read_raw_seq(vf);
      memcpy(aptr,bcode,rclen);
      aptr += rclen;

      oneReadLine(vf);
      aptr = Compress_QV(aptr,rlen,fqvs,qbits,map);

      bptr[byte] = aptr - array;

      o += 2;
      if (o >= end)
        break;

      t = oneReadLine(vf);
      if (t == 'g')
        oneReadLine(vf);
    }

  return (NULL);
}


/****************************************************************************************
 *
 *  Thread to output sorted array of compressed pairs into a .10x file
 *
 *****************************************************************************************/

typedef struct
  { OneFile     *vg;       //  OneFile for output
    int64        beg;      //  Range of read pairs to process
    int64        end;
    int          qbits;
    uint8       *inv;
    int          flen, rlen;
    int          reclen;
    uint8       *array;
  } Out_Arg;

static void *output_thread(void *arg)
{ Out_Arg  *parm   = (Out_Arg *) arg;
  OneFile  *vg     = parm->vg;
  int64     beg    = parm->beg;
  int64     end    = parm->end;
  int       qbits  = parm->qbits;
  uint8    *inv    = parm->inv;
  int       flen   = parm->flen;
  int       rlen   = parm->rlen;
  int       reclen = parm->reclen;
  uint8    *array  = parm->array;

  uint8  qmask = (1<<qbits)-1;
  char  *data;
  char  *forw, *revr;
  char  *fqvs, *rqvs;
  uint32 code, last;
  uint8 *aptr;
  OneInfo *ls, *lx;
  int    i;

  data = Malloc(2*(flen+rlen),"Allocating decompression buffers");
  forw = data;
  fqvs = forw + flen;
  revr = fqvs + flen;
  rqvs = revr + rlen;

  lx = vg->info['&'];
  ls = vg->info['S'];

  if (lx->buffer != NULL)
    { fprintf(stderr,"Surprised\n"); fflush(stdout);
      free(lx->buffer);
    }
  lx->buffer  = Malloc((end-beg)*sizeof(I64),"Sequence index setup");
  lx->bufSize = end-beg; 

  if ( ! vg->isLastLineBinary)      // terminate previous ascii line
    fputc ('\n', vg->f);
  vg->byte = ftello (vg->f);
  vg->isLastLineBinary = true;

  aptr = array + beg*reclen;
  if (beg == 0)
    last = *((uint32 *) array) - 1;
  else
    last = *((uint32 *) (aptr-reclen));

  for (i = beg; i < end; i++)
    { code = *((uint32 *) aptr);

      aptr = Uncompress_SEQ(aptr,forw,flen);
      aptr = Uncompress_QV(aptr,flen,fqvs,qbits,qmask,inv);
      aptr = Uncompress_SEQ(aptr,revr,rlen);
      aptr = Uncompress_QV(aptr,rlen,rqvs,qbits,qmask,inv);

      if (code != last)
        { oneInt(vg,0) = 0;
          oneInt(vg,1) = 16;
          oneWriteLine(vg,'g',16,forw);
          last = code;
        }

      oneWriteLine(vg,'P',0,NULL);

      oneInt(vg,0) = flen-23;
      oneWriteLine(vg,'S',flen-23,forw+23);
      oneWriteLine(vg,'Q',flen-23,fqvs+23);

      oneInt(vg,0) = rlen;
      oneWriteLine(vg,'S',rlen,revr);
      oneWriteLine(vg,'Q',rlen,rqvs);
    }

  free(data);

  return (NULL);
}

/****************************************************************************************
 *
 *  MAIN
 *
 *****************************************************************************************/

int main(int argc, char *argv[])
{ OneSchema *schema;
  char      *command;
  char      *fname;
  OneFile   *vf;

  uint32 *codes;        // Barcode list
  int64   npairs;       // # of reads in each file
  int64   nreads;       // # of reads in each file

  uint8   map[256];     // Map from QV char to bit code
  uint8   inv[256];     // Map from bit code to QV char
  int     qbits;        // # of bits to encode QV values
  uint8   qmask;        // mask of lowest qbits ina a byte

  uint8  *ctype;        // 4^16 bit vector of good codes and corrections
  int64   ndist;
  int64   ngood;
  int64   nrepair;
  int64  *bucks;        // count buckets/fingers for top level radix sort of data

  int flen, rlen;       // forward and reverse read lengths

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

    ARG_INIT("VGPcloud")

    NTHREADS  = 4;
    VALID     = 100;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vH")
            break;
          case 't':
            ARG_POSITIVE(VALID,"Valid code threshhold")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE   = flags['v'];
    HISTOGRAM = flags['H'];

    if (argc != 2)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, show progress as proceed.\n");
        fprintf(stderr,"      -H: Display histogram of all bar code counts.\n");
        fprintf(stderr,"      -t: Threshold for valid bar-codes.\n");
        fprintf(stderr,"      -T: Use -T threads.\n");
 
        exit (1);
      }
  }

  //  Open .irp VGPfile file with NTHREADS

  { char *pwd;
    int   i;

    pwd   = PathTo(argv[1]);
    fname = Root(argv[1],".irp");
    vf    = oneFileOpenRead(Catenate(pwd,"/",fname,".irp"),schema,"seq",NTHREADS);
    if (vf == NULL)
      { fprintf(stderr,"%s: Cannot open %s as an .irp file\n",Prog_Name,argv[1]);
        exit (1);
      }
    free(pwd);

    if ( ! vf->isBinary)
      { fprintf(stderr,"%s: Input is not a binary file\n",Prog_Name);
        exit (1);
      }
    if (vf->info['Q']->given.count <= 0)
      { fprintf(stderr,"%s: Input file does not have QV vectors\n",Prog_Name);
        exit (1);
      }
    for (i = 'A'; i < 'Z'; i++)
      if (vf->info[i] != NULL && vf->info[i]->given.count > 0)
        if ( ! (i == 'S' || i == 'Q' || i == 'P'))
          { fprintf(stderr,"%s: Input contains data lines other than S, Q, P, or g\n",Prog_Name);
            exit (1);
          }
    if (vf->info['S']->given.count != vf->info['Q']->given.count)
      { fprintf(stderr,"%s: The number of sequences and QV's are not equal\n",Prog_Name);
        exit (1);
      }
    if (2*vf->info['P']->given.count != vf->info['S']->given.count)
      { fprintf(stderr,"%s: The sequences are not all paired\n",Prog_Name);
        exit (1);
      }
  }


  //  Allocate vectors for code sorting (codes), code correction (ctype), and first radix
  //    sort bucket counts (bucks)

  nreads = vf->info['S']->given.count;
  npairs = vf->info['P']->given.count;

  bucks = (int64 *) Malloc(sizeof(int64)*256*NTHREADS,"Allocating bucket counters");
  ctype = (uint8 *) Malloc(0x100000000ll,"Allocating correction vector");
  codes = (uint32 *) Malloc(sizeof(uint32)*2*(npairs+NTHREADS),"Allocating bar-code array");
  if (bucks == NULL || ctype == NULL || codes == NULL)
    exit (1);

  if (VERBOSE)
    { fprintf(stderr,"  %.1fGB memory needed for bar-code sorting & correction\n",
                     (npairs+NTHREADS) * (2*sizeof(uint32) + sizeof(uint8)) / 1073741824.
                     + 4. + .5*NTHREADS);
      fprintf(stderr,"  Extracting bar-codes in scan of file %s ...\n",argv[1]);
      fflush(stderr);
    }

  //  Pass 1: Build list of bar-codes & build fixed-bit QV code on first 10 million symbols

  { BarCode_Arg parm[NTHREADS];
    pthread_t   threads[NTHREADS];
    
    int64     usedqvs[256];
    int       i, n, c, p;

    for (i = 0; i < NTHREADS; i++)
      { parm[i].beg   = (npairs * i) / NTHREADS;
        parm[i].end   = (npairs * (i+1)) / NTHREADS;
        parm[i].vf    = vf+i;
        parm[i].codes = codes + (parm[i].beg + i);  // leave room for a terminal entry
#ifdef DEBUG
        barcodes_thread(parm+i);
#else
        pthread_create(threads+i,NULL,barcodes_thread,parm+i);
#endif
      }

#ifndef DEBUG
    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);
#endif

    //  Merge QV histograms and check for errors

    bzero(usedqvs,sizeof(int64)*256);
    
    flen = parm[0].flen;
    rlen = parm[0].rlen;
    for (i = 0; i < NTHREADS; i++)
      { if (parm[i].error > 0)
          { error_report(parm[i].error,parm[i].where);
            exit (1);
          }
        if (flen != parm[i].flen)
          { error_report(ERROR_SLEN,parm[i].beg);
            exit (1);
          }
        if (rlen != parm[i].rlen)
          { error_report(ERROR_SLEN,parm[i].beg);
            exit (1);
          }
        for (n = 0; n < 256; n++)
          usedqvs[n] += parm[i].usedqvs[n];
      }

    n = 0;                         //  Map non-emty QV histogram entries to consecutive integers
    for (i = 0; i < 256; i++)      //    and all empty entries to the nearest code,
      if (usedqvs[i])              //    and determine # of bits to encode the # of ints used
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

  //  5 passes to produce correction vector and 1st byte sort bucket counts

  if (VERBOSE)
    { fprintf(stderr,"  About to sort & analyzea ");
      Print_Number((int64) npairs,0,stderr);
      fprintf(stderr," bar-codes ...\n");
      fflush(stderr);
    }

  { pthread_t  threads[NTHREADS];
    Sort_Arg   sparm[NTHREADS];

    int64     parts[NTHREADS+1][NTHREADS];
    int64     slice[NTHREADS];

    uint8 *first;   //  bit-vector of all codes for each thread
    uint8 *count;   //  small counts of comnpressed unique codes
    int    i;

    count = (uint8 *) Malloc(sizeof(uint8)*(npairs+NTHREADS),"Allocating count array");
    first = (uint8 *) Malloc(sizeof(uint8)*NTHREADS*0x20000000ll,"Allocating first array");
    if (count == NULL || first == NULL)
      exit (1);

    //  1. Sort and then build unique/count lists in codes/count[beg..end) for each thread panel

    for (i = 0; i < NTHREADS; i++)
      { sparm[i].tid   = i;
        sparm[i].beg   = (npairs * i) / NTHREADS + i;
        sparm[i].end   = (npairs * (i+1)) / NTHREADS + i;
        sparm[i].codes = codes;
        sparm[i].vlist = codes + (npairs + NTHREADS);
        sparm[i].count = count;
	sparm[i].ctype = ctype;
        sparm[i].first = first + i * 0x20000000ll;
        pthread_create(threads+i,NULL,sort_thread,sparm+i);
      }

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

    //  2. Merge NTHREAD unique+count lists and make list of valid codes in vlist

    { int64 s, e, n;
      int64 hist10[101], hist1[50];
      int   t, j;

      //  parts[i][t] = index of 1st code of input panel t >= (2^32)*(i/NTHREADS)

      for (t = 0; t < NTHREADS; t++)
        { s = sparm[t].beg;
          e = sparm[t].end;
          parts[0][t] = s;
          for (i = 1; i < NTHREADS; i++)
            parts[i][t] = s + find(codes+s,e-s,(uint32) ((0x100000000llu * i) / NTHREADS));
          parts[NTHREADS][t] = e;
        }

      //  slice[i] = # of codes in the i (out of NTHREADS) slices of [0,2^32]

      for (i = 0; i < NTHREADS; i++)
        { n = 0;
          for (t = 0; t < NTHREADS; t++) 
            n += parts[i+1][t] - parts[i][t];
          slice[i] = n;
        }

      n = 0;
      for (i = 0; i < NTHREADS; i++)
        { sparm[i].finger = parts[i]; 
          sparm[i].vbeg   = n;
          sparm[i].vlen   = slice[i];
	  n += slice[i];
          pthread_create(threads+i,NULL,merge_thread,sparm+i);
        }

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);

      if (HISTOGRAM)
        { for (j = 0; j <= 100; j++)
            hist10[j] = sparm[0].hist10[j];
          for (j = 0; j < 50; j++)
            hist1[j] = sparm[0].hist1[j];
          for (i = 1; i < NTHREADS; i++)
            { for (j = 0; j <= 100; j++)
                hist10[j] += sparm[i].hist10[j];
              for (j = 0; j < 50; j++)
                hist1[j] += sparm[i].hist1[j];
            }
     
          fprintf(stderr,"  Histogram:\n");
          fprintf(stderr,"    >999: %7lld\n",hist10[100]);
          for (t = 99; t >= 5; t--)
            if (hist10[t] > 0)
              fprintf(stderr,"    %4d: %7lld\n",t*10,hist10[t]);
          for (t = 49; t >= 1; t--)
            fprintf(stderr,"    %4d: %7lld\n",t,hist1[t]);
        }
    }

#ifdef DEBUG_COUNTS

  { int64  nvalid, ncorrect, nwaste;
    int    th, i, j, n, m;
    int    ex;
    uint32 v, t, u;
    uint32 chg, msk, var;

    nvalid   = 0;
    ncorrect = 0;
    nwaste   = 0;
    for (th = 0; th < NTHREADS; th++)
      for (i = sparm[th].beg; i < sparm[th].end; i++)
        { v = codes[i];
          n = count[i];
          if (ctype[v] > 0)
            nvalid += n;
          else
            { ex  = 0;
              chg = 0x1;
              msk = 0x3;
#ifdef DEBUG_MARK
              printf("  %08x\n",v);
#endif
              for (j = 0; j < 16; j++)
                { t = (v & msk);
                  v -= t;
                  u = 0;
                  for (m = 0; m < 4; m++)
                    { if (t != u)
                        { var = v | u;
#ifdef DEBUG_MARK
                          printf("   %08x (%d/%d)",var,j,m);
#endif
                          if (ctype[var] > 0)
                            { ex += 1;
#ifdef DEBUG_MARK
                              printf(" V");
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
                  v |= t;
                }

              if (ex == 1)
                ncorrect += n;
              else
                nwaste += n;
            }
        }
      printf(" %lld  %lld  %lld\n",nvalid,ncorrect,nwaste);
    }
#endif

    ndist = 0;
    for (i = 0; i < NTHREADS; i++)
      ndist += sparm[i].vend - sparm[i].vbeg;

    // 3. Generate 1-neighborhood and mark

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,mark_thread,sparm+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

    // 4. Set 1-correctable entries and the correction in ctype

    for (i = 0; i < NTHREADS; i++)
      { sparm[i].first = first;
        pthread_create(threads+i,NULL,fix_thread,sparm+i);
      }

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

    free(first);

    // 5. Setup bucket vectors for simultaneous loading and byte-1 sort by each thread

    for (i = 0; i < NTHREADS; i++)
      { sparm[i].bucks = bucks + i*256;
        pthread_create(threads+i,NULL,bucket_thread,sparm+i);
      }

    for (i = 0; i < NTHREADS; i++)
       pthread_join(threads[i],NULL);
	
    nrepair = 0;
    for (i = 0; i < NTHREADS; i++)
      nrepair += sparm[i].nrepair;

    { int64 x, y;
      int   j;

      x = 0;
      for (j = 0; j < 256; j++)
        for (i = 0; i < NTHREADS; i++)
          { y = sparm[i].bucks[j];
            sparm[i].bucks[j] = x;
            x += y;
          }
      ngood = x;
    }

    if (VERBOSE)
      { fprintf(stderr,"    There are ");
        Print_Number(ndist,0,stderr);
        fprintf(stderr," distinct valid bar-codes\n");
        fprintf(stderr,"      ");
        Print_Number((int64) (ngood-nrepair),12,stderr);
        fprintf(stderr," (%.1f%%) bar-codes are valid\n",(100.*(ngood-nrepair))/npairs);
        fprintf(stderr,"      ");
        Print_Number((int64) nrepair,12,stderr);
        fprintf(stderr," (%.1f%%) will be repaired\n",(100.*nrepair)/npairs);
        fprintf(stderr,"      ");
        Print_Number((int64) (npairs-ngood),12,stderr);
        fprintf(stderr," (%.1f%%) will be dropped.\n",(100.*(npairs-ngood))/npairs);
        fflush(stderr);
      }

    free(count);
    free(codes);
  }


  //  Pass 2: Correct bar-codes when possible, build compressed pairs in a mem-mapped array in
  //             sorted order of the 1st byte

  { int64  asize;
    uint8 *array;
    int    fclen, rclen;
    int    fqlen, rqlen;
    int    reclen;

    fclen  = (flen*2+7)/8; 
    rclen  = (rlen*2+7)/8; 
    fqlen  = (flen*qbits+7)/8; 
    rqlen  = (rlen*qbits+7)/8; 
    reclen = fclen + fqlen + rclen + rqlen;
    asize  = reclen*ngood;

    if (VERBOSE)
      { fprintf(stderr,"  Building compressed array of read pairs with valid bar-codes\n");
        fprintf(stderr,"    sorted on their 1st byte in a second scan of file\n");
        fprintf(stderr,"    Memory-mapped array size is %.1fGB (%d bytes/pair) ...\n",
                       asize/1073741824.,reclen);
        fflush(stderr);
      }

    array = (uint8 *) mmap(NULL,asize,PROT_READ | PROT_WRITE,MAP_ANON | MAP_PRIVATE,-1,0);
    if (array == MAP_FAILED)
      { fprintf(stderr,"%s: Cannot create memory map\n",Prog_Name);
        exit (1);
      }
  
    madvise(array, asize, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);

    { Load_Arg  parm[NTHREADS];
      pthread_t threads[NTHREADS];

      int i;

      for (i = NTHREADS*256-1; i >= 0; i--)
        bucks[i] *= reclen;

      for (i = 0; i < NTHREADS; i++)
        { parm[i].vf     = vf+i;
          parm[i].beg    = (npairs * i) / NTHREADS;
          parm[i].end    = (npairs * (i+1)) / NTHREADS;
          parm[i].bptr   = bucks + i*256;
          parm[i].array  = array;
          parm[i].ctype  = ctype;
          parm[i].qbits  = qbits;
          parm[i].map    = map;
          parm[i].flen   = flen;
          parm[i].rlen   = rlen;
          parm[i].fclen  = fclen;
          parm[i].rclen  = rclen;
          parm[i].fqlen  = fqlen;
          parm[i].rqlen  = rqlen;
          parm[i].array  = array;
#ifdef DEBUG
          load_thread(parm+i);
#else
          pthread_create(threads+i,NULL,load_thread,parm+i);
#endif
        }

#ifndef DEBUG
      for (i = 0; i < NTHREADS; i++)
          pthread_join(threads[i],NULL);
#endif

      free(ctype);

      if (VERBOSE)
        { fprintf(stderr,"    Fill & 1st byte sort complete.\n");
          fprintf(stderr,"  Commencing MSD sort of remaining 3 bytes ...\n");
          fflush(stderr);
        }
    }
    //  External sort on remaining 3 bytes of bar-code (1st 4 bytes of each record)

    { int64 *nucks = bucks + (NTHREADS-1)*256;
      int64  x, y;
      int    i;

      x = 0;
      for (i = 0; i < 256; i++)
        { y = nucks[i];
          bucks[i] = y-x;
          x = y;
        }

      MSD_Sort(array,ngood,reclen,reclen,4,bucks,NTHREADS);

      free(bucks);
    }

    //  Output clouds to <input>.10x trimming 1st 23bp of forward reads.

    { Out_Arg   parm[NTHREADS];
      pthread_t threads[NTHREADS];
      OneFile  *vg;
      char     *pwd, *root, *gname;
      int       i;

      pwd   = PathTo(argv[1]);
      root  = Root(argv[1],".irp");
      gname = Strdup(Catenate(pwd,"/",root,".10x"),"Allocating full path name");

      if (VERBOSE)
        { fprintf(stderr,"  Final scan to produce cloud grouped pairs in %s.10x ...\n",root);
          fflush(stderr);
        }

      vg = oneFileOpenWriteNew(gname,schema,"x10",true,NTHREADS);
      if (vg == NULL)
        { fprintf(stderr,"%s: Cannot open %s.10x for writing\n",Prog_Name,root);
          exit (1);
        }

      free(gname);
      free(root);
      free(pwd);

      oneInheritProvenance(vg,vf);
      oneAddProvenance(vg,"VGPcloud","1.0",command,NULL);

      oneWriteHeader(vg);

      oneFileClose(vf);

      if (VERBOSE)
        { if (NTHREADS > 1)
            fprintf(stderr,"  Producing .10x segments in parallel\n");
          else
            fprintf(stderr,"  Producing .10x segments\n");
          fflush(stderr);
        }

        //  Generate the data lines in parallel threads

      for (i = 0; i < NTHREADS; i++)
        { parm[i].vg     = vg+i;
          parm[i].beg    = (ngood * i) / NTHREADS;
          parm[i].end    = (ngood * (i+1)) / NTHREADS;
          parm[i].qbits  = qbits;
          parm[i].inv    = inv;
          parm[i].flen   = flen;
          parm[i].rlen   = rlen;
          parm[i].reclen = reclen;
          parm[i].array  = array;
        }

#ifdef DEBUG_OUT
      for (i = 0; i < NTHREADS; i++)
        { fprintf(stderr,"Thread %d\n",i); fflush(stderr);
          output_thread(parm+i);
        }
#else
      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,output_thread,parm+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);
#endif

      if (VERBOSE)
        { fprintf(stderr,"  Cat'ing .10x segments\n");
          fflush(stderr);
        }
  
      oneFileClose(vg);
    }
  }

  oneSchemaDestroy(schema);

  free(command);

  exit (0);
}
