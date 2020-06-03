/*  Last edited: Apr 16 20:16 2020 (rd109) */
/*******************************************************************************************
 *
 *  Length limited Huffman Compressor/decompressor with special 2-bit compressor for DNA.
 *
 *  Author:   Gene Myers
 *  Date:     June 27, 2019
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

typedef void OneCodec;  //  Compressor details are supressed

#include "compression.h"

#undef  DEBUG
#undef  TEST

typedef unsigned long long uint64;
typedef unsigned int       uint32;
typedef unsigned short     uint16;
typedef unsigned char      uint8;

#define HUFF_CUTOFF  12     //  This cannot be larger than 16 !

  //  Endian flipping macros

#define FLIP64(p)	\
{ uint8 x = p[0];	\
  p[0] = p[7];		\
  p[7] = x;		\
  x     = p[1];		\
  p[1] = p[6];		\
  p[6] = x;		\
  x     = p[2];		\
  p[2] = p[5];		\
  p[5] = x;		\
  x     = p[3];		\
  p[3] = p[4];		\
  p[4] = x;		\
}

#define FLIP32(p)	\
{ uint8 x = p[0];	\
  p[0] = p[3];		\
  p[3] = x;		\
  x     = p[1];		\
  p[1] = p[2];		\
  p[2] = x;		\
}

#define FLIP16(p)	\
{ uint8 x = p[0];	\
  p[0] = p[1];		\
  p[1] = x;		\
}

/*******************************************************************************************
 *
 *  Routines for computing a length-limited Huffman Encoding Scheme
 *
 ********************************************************************************************/

#define EMPTY        0      //  Compressor just created, histogram zero'd
#define FILLED       1      //  Compressor histogram being filled, no codec
#define CODED_WITH   2      //  Compressor has a codec (can no longer accumulate histogram)
#define CODED_READ   3      //  Compressor has codec but no histogram as was created by read

typedef struct
  { int    state;            //  1 of the 4 states immediately above
    int    isbig;            //  endian of the current machine
    uint16 codebits[256];    //  Code esc_code is the special code for
    uint8  codelens[256];    //    non-Huffman exceptions
    char   lookup[0x10000];  //  Lookup table (just for decoding)
    int    esc_code;         //  The special escape code (-1 if not partial)
    int    esc_len;          //  The length in bits of the special code (if present)
    uint64 hist[256];        //  Byte distribution for codec
  } _OneCodec;

  //  The special "predefined" DNA compressor

static _OneCodec _DNAcodec = { .state = CODED_READ };
OneCodec  *DNAcodec = (OneCodec *) &_DNAcodec;

  //  Create an EMPTY compressor object with zero'd histogram and determine machine endian

OneCodec *vcCreate()
{ _OneCodec *v;
  int i;

  v = (_OneCodec *) malloc(sizeof(_OneCodec));
  if (v == NULL)
    { fprintf(stderr,"vcCreate: Could not allocate compressor\n");
      exit (1);
    }

  v->state = EMPTY;
  for (i = 0; i < 256; i++)
    v->hist[i] = 0;

  { uint32 t;
    uint8 *b;

    t = 1;
    b = (uint8 *) (&t);
    v->isbig = (b[0] == 0);
  }

  return ((OneCodec *) v);
}

  //  Free a compressor object

void vcDestroy(OneCodec *vc)
{ _OneCodec *v = (_OneCodec *) vc;
  if (vc != DNAcodec)
    free(v);
}

  //  Add the frequencies of bytes in bytes[0..len) to vc's histogram
  //    State becomes FILLED

void vcAddToTable(OneCodec *vc, int len, char *bytes)
{ _OneCodec *v = (_OneCodec *) vc;
  uint8 *data = (uint8 *) bytes;
  int i;

  for (i = 0; i < len; i++)
    v->hist[(int) data[i]] += 1;
  if (v->state < FILLED)
    v->state = FILLED;
}

  //  Add the frequencies of bytes in bytes[0..len) to vc's histogram
  //    State becomes FILLED

void vcAddHistogram(OneCodec *vc, OneCodec *vh)
{ _OneCodec *v = (_OneCodec *) vc;
  _OneCodec *h = (_OneCodec *) vh;
  int i;

  if (v->state >= CODED_WITH)
    { fprintf(stderr,"vcAddHistogram: Compressor already has a codec\n");
      exit (1);
    }
  if (h->state == CODED_READ)
    { fprintf(stderr,"vcAddHistogram: Source compressor doesn't have a histogram\n");
      exit (1);
    }

  for (i = 0; i < 256; i++)
    v->hist[i] += h->hist[i];
  v->state = FILLED;
}

  //  Check vc has a non-empty distribution histogram and if so then build
  //    length-limited Huffman tables for the bytes that occur in the histogram,
  //    plus a special escape code if partial is set and there is at least one byte
  //    with a zero count in the histogram.  The algorithm is by Larmore & Hirschberg,
  //    JACM 73, 3 (1990).

uint64 *HIST;

int HSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);
  return (HIST[x] - HIST[y]);
}

void vcCreateCodec(OneCodec *vc, int partial)
{ _OneCodec *v = (_OneCodec *) vc;

  uint64  *hist;
  char    *look;
  uint8   *lens;
  uint16  *bitv;

  int      code[256];
  int      leng[256];
  uint16   bits[256];
  int      ncode, dcode, ecode;

  int      i;

  if (v->state >= CODED_WITH)
    { fprintf(stderr,"vcCreateCoder: Compressor already has a codec\n");
      exit (1);
    }
  if (v->state == EMPTY)
    { fprintf(stderr,"vcCreateCoder: Compressor has no byte distribution data\n");
      exit (1);
    }

  hist  = v->hist;
  look  = v->lookup;
  lens  = v->codelens;
  bitv  = v->codebits;

  ecode = -partial;
  ncode = 0;
  for (i = 0; i < 256; i++)
    if (hist[i] > 0)
      code[ncode++] = i;
    else if (ecode < 0)
      { ecode = i;
        code[ncode++] = i;
      }
  dcode = 2*ncode;

  if (ecode < 0)
    partial = 0;

  HIST = hist;
  qsort(code,ncode,sizeof(int),HSORT);

#ifdef DEBUG
  fprintf(stderr,"\nSorted Codes %d:\n",ncode);
  for (i = 0; i < ncode; i++)
    fprintf(stderr," %3d: %3d %10llu\n",i,code[i],hist[code[i]]);
#endif

  { uint8   matrix[HUFF_CUTOFF][dcode];
    uint64  count1[dcode], count2[dcode], countb[ncode];
    uint64 *lcnt, *ccnt, *swp;
    int     llen, span;
    int     j, k, n, L;

    for (n = 0; n < ncode; n++)
      { count1[n] = countb[n] = hist[code[n]];
        leng[n] = 0;
      }

#ifdef DEBUG
    fprintf(stderr,"\nCoin Filter:\n");
    fprintf(stderr,"  Row %2d:",HUFF_CUTOFF);
    for (n = 0; n < ncode; n++)
      fprintf(stderr," %lld*",countb[n]);
    fprintf(stderr,"\n");
#endif

    lcnt = count1;
    ccnt = count2;
    llen = ncode-1;
    for (L = HUFF_CUTOFF-1; L > 0; L--)
      { j = 0;
        k = 0;
        for (n = 0; j < ncode || k < llen; n++)
          { if (k >= llen || (j < ncode && countb[j] <= lcnt[k] + lcnt[k+1]))
              { ccnt[n] = countb[j];
                matrix[L][n] = 1;
                j += 1;
              }
            else
              { ccnt[n] = lcnt[k] + lcnt[k+1];
                matrix[L][n] = 0;
                k += 2;
              }
          }
        llen = n-1;
        swp  = lcnt;
        lcnt = ccnt;
        ccnt = swp;

#ifdef DEBUG
        fprintf(stderr,"  Row %2d:",L);
        for (n = 0; n <= llen; n++)
          fprintf(stderr," %lld%c",lcnt[n],matrix[L][n]?'*':'+');
        fprintf(stderr,"\n");
#endif
      }

    span = 2*(ncode-1);
    for (L = 1; L < HUFF_CUTOFF; L++)
      { j = 0;
        for (n = 0; n < span; n++)
          { if (matrix[L][n])
              leng[j++] += 1;
          }
        span = 2*(span-j);
      }
    for (n = 0; n < span; n++)
      leng[n] += 1;

#ifdef DEBUG
    fprintf(stderr,"\nBack Trace:\n");
    span = 2*(ncode-1);
    for (L = 1; L < HUFF_CUTOFF; L++)
      { j = 0;
        fprintf(stderr,"  Row %2d:",L);
        for (n = 0; n < span; n++)
          { if (matrix[L][n])
              j += 1;
            fprintf(stderr," %c",matrix[L][n]?'*':'+');
          }
        fprintf(stderr,"\n");
        span = 2*(span-j);
      }
    fprintf(stderr,"  Length:");
    for (n = 0; n < ncode; n++)
      fprintf(stderr," %d",leng[n]);
    fprintf(stderr,"\n");
#endif
  }

  { int    n, llen;
    uint16 lbits;

    llen  = leng[0];
    lbits = bits[0] = (1 << llen) - 1;
    for (n = 1; n < ncode; n++)
      { while ((lbits & 0x1) == 0)
          { lbits >>= 1;
            llen -= 1;
          }
        lbits -= 1;
        while (llen < leng[n])
          { lbits = (lbits << 1) | 0x1;
            llen += 1;
          }
        bits[n] = lbits;
      }

#ifdef DEBUG
    { int j;

      fprintf(stderr,"\nCodes:\n");
      for (n = 0; n < ncode; n++)
        { fprintf(stderr,"   %3d: %2d ",code[n],leng[n]);
          for (j = leng[n]-1; j >= 0; j--)
            fprintf(stderr,"%x",(bits[n]>>j)&0x1);
          fprintf(stderr,"\n");
        }
    }
#endif
  }

  for (i = 0; i < 256; i++)
    { lens[i] = 0;
      bitv[i] = 0;
    }

  for (i = 0; i < ncode; i++)
    { lens[code[i]] = leng[i];
      bitv[code[i]] = bits[i];
    }

  { int    j, powr;    //  Fill in a decoder table giving the next Huffman code
    uint16 base;       //    that is a prefix of the next 16 bits

    for (i = 0; i < 256; i++)
      { if (lens[i] > 0)
          { base = (bitv[i] << (16-lens[i]));
            powr = (1 << (16-lens[i]));
            for (j = 0; j < powr; j++)
              look[base+j] = i;
          }
      }
  }

  if (partial)
    { v->esc_code = ecode;
      v->esc_len  = lens[ecode];
      lens[ecode] = 0;
    }
  else
    v->esc_code = -1;
  v->state = CODED_WITH;
}

  //  For debug, give a nice print out of the distribution histogram (if present)
  //     and the Huffman codec

void vcPrint(OneCodec *vc, FILE *to)
{ _OneCodec *v = (_OneCodec *) vc;

  uint64  total_bits, ucomp_bits, count;
  uint16  mask, code, *bits;
  uint64 *hist;
  uint8  *lens;
  int     clen;
  int     hashist;
  int     i, k;

  if (vc == DNAcodec)
    { fprintf(to,"    DNAcompressor\n");
      return;
    }

  if (v->state < CODED_WITH)
    { fprintf(stderr,"vcPrint: Compressor has no codec\n");
      exit (1);
    }
  hashist = (v->state == CODED_WITH);

  bits = v->codebits;
  lens = v->codelens;

  if (hashist)
    { hist = v->hist;
      total_bits = 0;
      ucomp_bits = 0;

      count = 0;
      for (i = 0; i < 256; i++)
        count += hist[i];

      fprintf(to,"\nHistogram:\n");
      for (i = 0; i < 256; i++)
        if (hist[i] > 0)
          { if (isprint(i))
              fprintf(to,"      %c: %12llu %5.1f%%\n",i,hist[i],(hist[i]*100.)/count);
            else
              fprintf(to,"    %3d: %12llu %5.1f%%\n",i,hist[i],(hist[i]*100.)/count);
          }
    }

  fprintf(to,"\nCode Table:\n");
  for (i = 0; i < 256; i++)
    { clen = lens[i];
      if (i == v->esc_code)
        clen = v->esc_len;
      if (clen > 0)
        { mask = (1 << clen);
          code = bits[i];
          if (isprint(i))
            fprintf(to,"   %c: %2d ",i,clen);
          else
            fprintf(to," %3d: %2d ",i,clen);
          for (k = 0; k < clen; k++)
            { mask >>= 1;
              if (code & mask)
                fprintf(to,"1");
              else
                fprintf(to,"0");
            }
          if (i == v->esc_code)
            fprintf(to," ***\n");
          else
            { fprintf(to,"\n");
              if (hashist)
                { total_bits += clen*hist[i];
                  ucomp_bits += (hist[i]<<3);
                }
            }
        }
    }
  if (hashist)
    fprintf(to,"\nTotal Bytes = %lld (%.2f%%)\n",(total_bits-1)/8+1,(100.*total_bits)/ucomp_bits);
}


/*******************************************************************************************
 *
 *  Read and Write Huffman Schemes (actually just (de)serialize)
 *
 ********************************************************************************************/

  //  Maximum # of bytes in a serialized compressor code

int vcMaxSerialSize()
{ return (257 + 2*sizeof(int) + 256*sizeof(uint16)); }

  //  Code the compressor into blob 'out' and return number of bytes in the code

int vcSerialize(OneCodec *vc, void *out)
{ _OneCodec *v = (_OneCodec *) vc;
  
  int     i;
  uint16 *bits;
  uint8  *lens, *o;

  if (vc == DNAcodec)
    return (0);

  if (v->state < CODED_WITH)
    { fprintf(stderr,"vcWrite: Compressor does not have a codec\n");
      exit (1);
    }

  lens = v->codelens;
  bits = v->codebits;
  o    = (uint8 *) out;

  //  Only need to record endian, escape code, code lengths, and codes for those
  //    with non-zero length

  *o++ = v->isbig;
  memcpy(o,&(v->esc_code),sizeof(int));
  o += sizeof(int);
  memcpy(o,&(v->esc_len),sizeof(int));
  o += sizeof(int);
  for (i = 0; i < 256; i++)
    { *o++ = lens[i];
      if (lens[i] > 0 || i == v->esc_code)
        { memcpy(o,bits+i,sizeof(uint16));
          o += sizeof(uint16);
        }
    }
  return (o - (uint8 *) out);
}

  //  Create a compressor object from the serialized code in blob 'in'.
  //    The compressor does not have the original histogram from which
  //    its codec was created.  If the endian of the current machine and
  //    the one that serialized the compressor don't match, then all relevant
  //    items are byte-flipped.

OneCodec *vcDeserialize(void *in)
{ _OneCodec *v;

  char    *look;
  uint8   *lens, *ip;
  uint16  *bits, base;
  int      i, j, powr;

  v = (_OneCodec *) malloc(sizeof(_OneCodec));
  if (v == NULL)
    { fprintf(stderr,"vcRead: Could not allocate compressor\n");
      exit (1);
    }

  v->state = CODED_READ;
  lens = v->codelens;
  bits = v->codebits;
  look = v->lookup;
  ip   = (uint8 *) in;

  { uint32 t;
    uint8 *b;

    t = 1;
    b = (uint8 *) (&t);
    v->isbig = (b[0] == 0);
  }

  if (v->isbig != *ip++)  // If endians out and in don't match then flip item bytes as needed
    { FLIP32(ip)
      memcpy(&(v->esc_code),ip,sizeof(int));
      ip += sizeof(int);
      FLIP32(ip)
      memcpy(&(v->esc_len),ip,sizeof(int));
      ip += sizeof(int);
      for (i = 0; i < 256; i++)
        { lens[i] = *ip++;
          if (lens[i] > 0 || i == v->esc_code)
            { FLIP16(ip)
              memcpy(bits+i,ip,sizeof(uint16));
              ip += sizeof(uint16);
            }
          else
            bits[i] = 0;
        }
    }
  else
    { memcpy(&(v->esc_code),ip,sizeof(int));
      ip += sizeof(int);
      memcpy(&(v->esc_len),ip,sizeof(int));
      ip += sizeof(int);
      for (i = 0; i < 256; i++)
        { lens[i] = *ip++;
          if (lens[i] > 0 || i == v->esc_code)
            { memcpy(bits+i,ip,sizeof(uint16));
              ip += sizeof(uint16);
            }
          else
            bits[i] = 0;
        }
    }

  if (v->esc_code >= 0)
    lens[v->esc_code] = v->esc_len;
  for (i = 0; i < 256; i++)
    { if (lens[i] > 0)
        { base = (bits[i] << (16-lens[i]));
          powr = (1 << (16-lens[i]));
          for (j = 0; j < powr; j++)
            look[base+j] = i;
        }
    }
  if (v->esc_code >= 0)
    lens[v->esc_code] = 0;

  return ((OneCodec *) v);
}


/*******************************************************************************************
 *
 *  Encoders and Decoders
 *
 ********************************************************************************************/

static uint8 Number[128] =
    { 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
    };

  //  Compress DNA into 2-bits per base

int Compress_DNA(int len, char *s, char *t)
{ int    i, j;
  uint8 *s0, *s1, *s2, *s3;

  s0 = (uint8 *) s;
  s1 = s0+1;
  s2 = s1+1;
  s3 = s2+1;

  len -= 3;
  for (i = j = 0; i < len; i += 4)
    t[j++] = (Number[s0[i]] << 6) | (Number[s1[i]] << 4) | (Number[s2[i]] << 2) | Number[s3[i]];
  switch (i-len)
  { case 0:
      t[j++] = (Number[s0[i]] << 6) | (Number[s1[i]] << 4) | (Number[s2[i]] << 2);
      break;
    case 1:
      t[j++] = (Number[s0[i]] << 6) | (Number[s1[i]] << 4);
      break;
    case 2:
      t[j++] = (Number[s0[i]] << 6);
      break;
    default:
      break;
  }

  return ((len+3)<<1);
}

  //  Encode ibytes[0..ilen) according to compressor vc and place in obytes
  //  Return the # of bits used.

int vcEncode(OneCodec *vc, int ilen, char *ibytes, char *obytes)
{ _OneCodec *v = (_OneCodec *) vc;

  uint64  c, ocode, *ob;
  int     n, k, rem, tbits, ibits, esc, elen;
  uint8  *clens, x, *bcode, *bb;
  uint16 *cbits;

  if (vc == DNAcodec)
    return (Compress_DNA(ilen,ibytes,obytes));

  if (v->state < CODED_WITH)
    { fprintf(stderr,"vcEncode: Compressor does not have a codec\n");
      exit (1);
    }

  esc   = v->esc_code;
  elen  = v->esc_len;
  clens = v->codelens;
  cbits = v->codebits;
  ibits = (ilen << 3);
  bcode = (uint8 *) &ocode;

#define OCODE(L,C)				\
{ rem -= L;					\
  if (rem <= 0)					\
    { ocode |= (C >> (-rem));			\
      *ob++ = ocode;				\
      if (rem < 0)				\
        { rem   += 64;				\
          ocode = (C << rem);			\
        }					\
      else					\
        { rem   = 64;				\
          ocode = 0;				\
        }					\
    } 						\
  else						\
    ocode |= (C << rem);			\
}

  ob    = (uint64 *) obytes;
  tbits = 2;
  rem   = 62;
  if (v->isbig)
    ocode = 0x4000000000000000llu;
  else
    ocode = 0;
  for (k = 0; k < ilen; k++)
    { x = ibytes[k];
      n = clens[x];
      if (n == 0)
        { if (esc < 0)
            { fprintf(stderr,"Compression lib: No code for %c(%x) and no escape code\n",x,x);
              exit (1);
            }
          c = cbits[esc];
          tbits += 8+elen;
          if (tbits > ibits)
            break;
          OCODE(elen,c);
          c = x;
          OCODE(8,c);
        }
      else
        { tbits += n;
          if (tbits > ibits)
            break;
          c = cbits[x];
          OCODE(n,c);
        }
    }
  
  if (k < ilen)
    { *obytes = 0xff;
      memcpy(obytes+1,ibytes,ilen);
      return (ibits+8);
    }

  bb = (uint8 *) ob;
  if (v->isbig)
    { rem = ((71-rem)>>3);
      for (k = 0; k < rem; k++)
        *bb++ = bcode[k];
    }
  else
    { rem = 7 - ((63-rem)>>3);
      for (k = 7; k >= rem; k--)
        *bb++ = bcode[k];
    }

  if (tbits >= 64 && !v->isbig)
    { x = obytes[7];
      obytes[7] = obytes[0];
      obytes[0] = x;
    }

  return (tbits);
}

  //  Uncompress read from 2-bits per base into [0-3] per byte representation

static char Base[4] = { 'a', 'c', 'g', 't' };

int Uncompress_DNA(char *s, int len, char *t)
{ int   i, tlen, byte;
  char *t0, *t1, *t2, *t3;

  t0 = t;
  t1 = t0+1;
  t2 = t1+1;
  t3 = t2+1;

  tlen = len-3;
  for (i = 0; i < tlen; i += 4)
    { byte = *s++;
      t0[i] = Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 4) & 0x3];
      t2[i] = Base[(byte >> 2) & 0x3];
      t3[i] = Base[byte & 0x3];
    }

  switch (i-tlen)
  { case 0:
      byte = *s++;
      t0[i] = Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 4) & 0x3];
      t2[i] = Base[(byte >> 2) & 0x3];
      break;
    case 1:
      byte = *s++;
      t0[i] = Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 4) & 0x3];
      break;
    case 2:
      byte = *s++;
      t0[i] = Base[(byte >> 6) & 0x3];
      break;
    default:
      break;
  }

  return (len);
}

  //  Decode ilen bits in ibytes, into obytes according to vc's codec
  //  Return the number of bytes decoded.

int vcDecode(OneCodec *vc, int ilen, char *ibytes, char *obytes)
{ _OneCodec *v = (_OneCodec *) vc;

  char   *look;
  uint8  *lens, *q;
  uint64  icode, ncode, *p;
  int     rem, nem;
  uint8   c, *o;
  int     n, k, elen, inbig, esc;

  if (vc == DNAcodec)
    return (Uncompress_DNA(ibytes,ilen>>1,obytes));

  if (v->state < CODED_WITH)
    { fprintf(stderr,"vcDecode: Compressor does not have a codec\n");
      exit (1);
    }

  if (*((uint8 *) ibytes) == 0xff)
    { int olen = (ilen>>3)-1;
      memcpy(obytes,ibytes+1,olen);
      return (olen);
    }

  p = (uint64 *) ibytes;

  inbig = (*ibytes & 0x40);
  if (!inbig && ilen >= 64)
    { uint8 x = ibytes[7];
      ibytes[7] = ibytes[0];
      ibytes[0] = x;
    }

  if (inbig != v->isbig)
    { q = (uint8 *) ibytes;
      for (k = 64; k <= ilen; k += 64)
        { FLIP64(q)
          q += 8;
        }
    }

  lens = v->codelens;
  look = v->lookup;
  esc  = v->esc_code;
  elen = v->esc_len;

#define GET(n)						\
  ilen  -= n;						\
  icode <<= n;						\
  rem   -= n;						\
  while (rem < 16)					\
    { int z = 64-rem;					\
      icode |= (ncode >> rem);				\
      if (nem > z)					\
        { nem -= z;					\
          ncode <<= z;					\
          rem = 64;					\
          break;					\
        }						\
      else						\
        { rem += nem; 					\
          if (rem >= ilen)				\
            break;					\
          else if (ilen-rem < 64)			\
            { nem = ilen-rem;				\
              q = (uint8 *) p;				\
              ncode = 0;				\
              for (k = 0; k < nem; k += 8)		\
                ncode |= (((uint64) (*q++)) << (56-k));	\
            }						\
          else						\
            { ncode = *p++;				\
              nem   = 64;				\
            }						\
	}						\
    }
 
  if (ilen < 64)
    { q = (uint8 *) ibytes;
      icode = 0;
      for (k = 0; k < ilen; k += 8)
        icode |= (((uint64) (*q++)) << (56-k));
    }
  else
    icode = *p++;
  o = (uint8 *) obytes;
  icode <<= 2;
  ilen -= 2;
  rem   = 62;
  if (rem > ilen)
    rem = ilen;
  ncode = 0;
  nem   = 0;
  while (ilen > 0)
    { c = look[icode >> 48];
      if (c == esc)
        { GET(elen)
          c = (icode >> 56);
          GET(8);
        }
      else
        { n = lens[(int) c];
          GET(n)
        }
      *o++ = c;
    }

  return (o - (uint8 *) obytes);
}

#ifdef TEST

static char *bits[256] =
  { "00000000", "00000001", "00000010", "00000011", "00000100", "00000101", "00000110", "00000111",
    "00001000", "00001001", "00001010", "00001011", "00001100", "00001101", "00001110", "00001111",
    "00010000", "00010001", "00010010", "00010011", "00010100", "00010101", "00010110", "00010111",
    "00011000", "00011001", "00011010", "00011011", "00011100", "00011101", "00011110", "00011111",
    "00100000", "00100001", "00100010", "00100011", "00100100", "00100101", "00100110", "00100111",
    "00101000", "00101001", "00101010", "00101011", "00101100", "00101101", "00101110", "00101111",
    "00110000", "00110001", "00110010", "00110011", "00110100", "00110101", "00110110", "00110111",
    "00111000", "00111001", "00111010", "00111011", "00111100", "00111101", "00111110", "00111111",
    "01000000", "01000001", "01000010", "01000011", "01000100", "01000101", "01000110", "01000111",
    "01001000", "01001001", "01001010", "01001011", "01001100", "01001101", "01001110", "01001111",
    "01010000", "01010001", "01010010", "01010011", "01010100", "01010101", "01010110", "01010111",
    "01011000", "01011001", "01011010", "01011011", "01011100", "01011101", "01011110", "01011111",
    "01100000", "01100001", "01100010", "01100011", "01100100", "01100101", "01100110", "01100111",
    "01101000", "01101001", "01101010", "01101011", "01101100", "01101101", "01101110", "01101111",
    "01110000", "01110001", "01110010", "01110011", "01110100", "01110101", "01110110", "01110111",
    "01111000", "01111001", "01111010", "01111011", "01111100", "01111101", "01111110", "01111111",
    "10000000", "10000001", "10000010", "10000011", "10000100", "10000101", "10000110", "10000111",
    "10001000", "10001001", "10001010", "10001011", "10001100", "10001101", "10001110", "10001111",
    "10010000", "10010001", "10010010", "10010011", "10010100", "10010101", "10010110", "10010111",
    "10011000", "10011001", "10011010", "10011011", "10011100", "10011101", "10011110", "10011111",
    "10100000", "10100001", "10100010", "10100011", "10100100", "10100101", "10100110", "10100111",
    "10101000", "10101001", "10101010", "10101011", "10101100", "10101101", "10101110", "10101111",
    "10110000", "10110001", "10110010", "10110011", "10110100", "10110101", "10110110", "10110111",
    "10111000", "10111001", "10111010", "10111011", "10111100", "10111101", "10111110", "10111111",
    "11000000", "11000001", "11000010", "11000011", "11000100", "11000101", "11000110", "11000111",
    "11001000", "11001001", "11001010", "11001011", "11001100", "11001101", "11001110", "11001111",
    "11010000", "11010001", "11010010", "11010011", "11010100", "11010101", "11010110", "11010111",
    "11011000", "11011001", "11011010", "11011011", "11011100", "11011101", "11011110", "11011111",
    "11100000", "11100001", "11100010", "11100011", "11100100", "11100101", "11100110", "11100111",
    "11101000", "11101001", "11101010", "11101011", "11101100", "11101101", "11101110", "11101111",
    "11110000", "11110001", "11110010", "11110011", "11110100", "11110101", "11110110", "11110111",
    "11111000", "11111001", "11111010", "11111011", "11111100", "11111101", "11111110", "11111111",
  };

static char *test[4] =
  { "llkllkjllkllkjithlhlkl",
    "llkllkjllkllkjithlhlkk",
    "llkllkjllkllkjithlhlkkl",
    "mnopq",
  };

int main(int argc, char *argv[])
{ OneCodec *scheme;
  uint8     junk[10];
  FILE     *io;
  void     *blob;
  int       size, olen, ilen, tlen;
  uint8     obuf[300];
  char      ibuf[300];
  long long xbuf[50] =
       { 0, 0, 0, 0, 1, 1, 0, 0, 0, 1,
         0, 1, 1, 0, 0, 0, 0, 1, 0, 2,
         1, 1, 0, 0, 2, 0, 0, 1, 0, 0,
         0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 1, 1, 0, 0, 0, 13, 42 };

  int i, j, t;

  scheme = vcCreate();

  vcAddToTable(scheme,12,"abcdefghijkl");
  vcAddToTable(scheme,10,"cdefghijkl");
  for (i = 0; i < 2; i++)
    vcAddToTable(scheme,9,"defghijkl");
  for (i = 0; i < 4; i++)
    vcAddToTable(scheme,8,"efghijkl");
  for (i = 0; i < 8; i++)
    vcAddToTable(scheme,7,"fghijkl");
  for (i = 0; i < 16; i++)
    vcAddToTable(scheme,6,"ghijkl");
  for (i = 0; i < 32; i++)
    vcAddToTable(scheme,5,"hijkl");
  for (i = 0; i < 64; i++)
    vcAddToTable(scheme,4,"ijkl");
  for (i = 0; i < 128; i++)
    vcAddToTable(scheme,3,"jkl");
  for (i = 0; i < 256; i++)
    vcAddToTable(scheme,2,"kl");
  for (i = 0; i < 512; i++)
    vcAddToTable(scheme,1,"l");

  vcCreateCodec(scheme,1);

  vcPrint(scheme);

  blob = malloc(vcMaxSerialSize());

  size = vcSerialize(scheme,blob);
  printf("\nWriting %d bytes to .xxtest\n",size);

  io = fopen(".xxtest","w");
  fwrite(junk,1,10,io);
  fwrite(blob,1,size,io);
  fwrite(junk,1,10,io);
  fclose(io);

  vcDestroy(scheme);

  io = fopen(".xxtest","r");
  fread(junk,1,10,io);
  fread(blob,1,size,io);
  fclose(io);

  scheme = vcDeserialize(blob);

  vcPrint(scheme);

  for (t = 0; t < 4; t++)
    { strcpy(ibuf,test[t]);

      olen = vcEncode(scheme,strlen(ibuf),ibuf,(char *) obuf);

      printf("\nIn: %s\n",ibuf);
      printf("\nEncode: %d",olen);
      for (i = 0; i < ((olen+7)>>3); i++)
        printf(" %s",bits[obuf[i]]);
      printf("\n");
    
      ilen = vcDecode(scheme,olen,(char *) obuf,ibuf);
      ibuf[ilen] = '\0';
      printf("Decode: %s\n",ibuf);
    }
    
  vcDestroy(scheme);

  io = fopen("X.vcx","r");
  fread(blob,1,500,io);
  fclose(io);

  scheme = vcDeserialize(blob);

  vcPrint(scheme);

  olen = vcEncode(scheme,400,(char *) xbuf,(char *) obuf);

  printf("\nIn:\n");
  for (i = 0; i < 400; i++)
    printf(" %s",bits[((char *) xbuf)[i]]);
  printf("\nEncode: %d\n",olen);
  tlen = (olen+7)>>3;
  for (i = 0; i < tlen; i += 8)
    { if (i+8 <= tlen)
        for (j = i+7; j >= i; j--)
          printf(" %s",bits[obuf[j]]);
      else
        for (j = i; j < tlen; j++)
          printf(" %s",bits[obuf[j]]);
      if (i%8 == 0)
        printf("\n");
    }
  printf("\n");
 
  ilen = vcDecode(scheme,olen,(char *) obuf,ibuf);
  printf("Decode:\n");
  for (i = 0; i < ilen; i++)
    printf(" %s",bits[ibuf[i]]);

  vcDestroy(scheme);

  exit (0);
}

#endif
