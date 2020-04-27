/*******************************************************************************************
 *
 *  Simple single-thread radix sort
 *
 *  Author :  Gene Myers
 *  First  :  Feb 2020
 *
 ********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "lsd.sort.h"

typedef unsigned char uint8;
typedef long long     int64;

#undef TEST_LSORT

//  Simple LSD radix sort of uint32's

void LSD_Sort(int64 nelem, uint32 *src, uint32 *trg)
{ int64  ptrs[512];
  int    bytes[5];

  int64    i;
  int      j, b;

  Cptr  = ptrs;
  Nptr  = ptrs + 256;

#if __ORDER_LITTLE_ENDIAN__ == __BYTE_ORDER__
  for (i = 0; i < 4; i++)
    bytes[i] = i;
#else
  for (i = 0; i < 4; i++)
    bytes[i] = 3-i;
#endif
  bytes[4] = -1;

  //  For each requested byte b in order, radix sort

  for (b = 0; b < 4; b++)
    { int Cbyte, Nbyte;

      Cbyte  = bytes[b];
      Nbyte  = bytes[b+1];

#ifdef TEST_LSORT
      printf("\nSorting byte %d\n",Cbyte);
      fflush(stdout);
#endif

      //  Zero Cptr counters

      for (j = 0; j < 256; j++)
        Cptr[j] = 0;

      //  If first pass, then explicitly sweep to get Cptr counts

      if (b == 0)
        { uint8 *dig = srcD + Cbyte;

          for (i = 0; i < asize; i += rsize)
            Cptr[dig[i]] += 1;
        }

#ifdef TEST_LSORT
      printf("\nBUCKETS %d\n",Cbyte);
      for (j = 0; j < 255; j++)
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

      //  Radix sort from srcD to trgD

      { uint32  v;
        uint8  *p = ((uint8 *) &v) + Cbyte;
        uint8  *n = ((uint8 *) &v) + Nbyte;

        if (Nbyte < 0)
          for (i = 0; i < nelem; i++)
            { v = src[i];
              d = *p;
              x = Cptr[d]++;
              trg[x] = v;
            }
        else
          for (i = 0; i < nelem; i++)
            { v = src[i]; 
              d = *p;
              x = Cptr[d]++;
              trg[x] = v;
              Nptr[*n] += 1;
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
        uint8 *psort = srcD-rsize;

        printf("\nLSORT %d\n",Cbyte);
        for (c = 10000*rsize; c < 11000*rsize; c += rsize)
          { printf(" %4lld: ",c/rsize);
            for (j = 0; j < dsize; j++)
              printf(" %02x",srcD[c+j]);
            printf("\n");
          }

        for (c = rsize; c < asize; c += rsize)
          { for (j = Cbyte; j >= 2; j--)
              if (srcD[c+j] > psort[c+j])
                break;
              else if (srcD[c+j] < psort[c+j])
                { printf("  Order: %lld",c/rsize);
                  for (x = 2; x <= Cbyte; x++)
                    printf(" %02x",psort[c+x]);
                  printf(" vs");
                  for (x = 2; x <= Cbyte; x++)
                    printf(" %02x",srcD[c+x]);
                  printf("\n");
                  break;
                }
          }
      }
#endif
    }
}
