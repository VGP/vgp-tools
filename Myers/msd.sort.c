#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>

#include "gene_core.h"
#include "msd.sort.h"

#undef  DEBUG

#define SMAX  6

#define THR0 15
#define THR1 15
#define THR2  8
#define GAP1  7
#define GAP2  3

static int S_thr0, S_thr1, S_thr2;
static int S_gap1, S_gap2;

static int    RSIZE;
static int    DSIZE;
static int    KSIZE;
static int    PSIZE;

static int    COFF;
static int64  CMAX;

#ifdef SHOW

static void print_kmer(uint8 *x)
{ int i;
  for (i = 0; i < KSIZE; i++)
    printf(" %02x",x[i]);
}

#endif

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n--)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}

static inline void mycpy(uint8 *a, uint8 *b, int n)
{ while (n--)
    *a++ = *b++;
}

#ifdef DEBUG

static inline void sorted(uint8 *array, int64 asize, int level)
{ int64 p, i;

  for (p = RSIZE; p < asize; p += RSIZE)
    if (mycmp(array + (p-RSIZE),array + p,KSIZE) > 0)
      { printf("Not sorted %12lld: ",p);
        for (i = level; i < KSIZE; i++)
          printf(" %02x",array[(p-RSIZE)+i]);
        printf(" vs ");
        for (i = level; i < KSIZE; i++)
          printf(" %02x",array[p+i]);
        printf("\n");
      }
}

#endif


static inline void gap_sort(uint8 *array, int asize, int gap, int cmp)
{ int    i, j, rem;
  uint8  temp[RSIZE];
  uint8 *garray;

  rem    = cmp + PSIZE;
  garray = array + gap;
  for (i = gap; i < asize; i += RSIZE)
    { j = i-gap;
      if (mycmp(array+j,array+i,cmp) <= 0)
        continue;
      mycpy(temp,array+i,rem);
      mycpy(array+i,array+j,rem);
      for(j -= gap; j >= 0; j -= gap)
        { if (mycmp(array+j,temp,cmp) <= 0)
            break;
          mycpy(garray+j,array+j,rem);
        }
      mycpy(garray+j,temp,rem);
    }
}

static inline void shell_sort(uint8 *array, int asize, int digit)
{ int    cmp;

  cmp    = KSIZE-digit;
  array += digit;

  if (asize > S_thr1)
    gap_sort(array,asize,S_gap1,cmp);
  if (asize > S_thr2)
    gap_sort(array,asize,S_gap2,cmp);
  gap_sort(array,asize,RSIZE,cmp);
}

#ifdef PERMUTE_SORT_IDEA

static inline void gap_perm(uint8 *array, int *elem, int n, int gap, int cmp)
{ int i, j, rem;

  for (i = gap; i < n; i += RSIZE)
    { j = i-gap;
      if (mycmp(el[j],el[i],cmp) <= 0)
        continue;
      t = el[i];
      el[i] = el[j];
      for(j -= gap; j >= 0; j -= gap)
        { if (mycmp(el[j],t,cmp) <= 0)
            break;
          el[j+gap] = el[j];
        }
      el[j+gap] = t;
    }
}

static inline void shell_sort(uint8 *array, int asize, int digit)
{ int    cmp;
  int    i, j;
  uint16 cnt;

  cmp    = KSIZE-digit;
  array += digit;

  n = 0;
  for (i = 0; i < asize; i += RSIZE)
    elem[n++] = array + i;

  if (n > THR1)
    gap_sort(array,elem,n,S_GAP1,cmp);
  if (n > THR2)
    gap_sort(array,elem,n,S_GAP2,cmp);
  gap_sort(array,elem,n,1,cmp);

  rem = cmp + PSIZE;
  for (i = 0; i < n; i++)
    if (elem[i] != NULL)
      { mycpy(temp,elem[i]),
        for (j = i; elem[j] != NULL; j = elem[j]) 
          { mycpy(array+j*RSIZE,elem[j],rem);
            elem[j] = NULL;
          }
      }
}

#endif

static void radix_sort(uint8 *array, int64 asize, int digit)
{ int64  n, len[256];
  int    x;

  { uint8 *end[256];
    uint8 *u, *arrow = array + digit;
    int64  o;
    int    rems;
    int    e;

    uint8 *off[256];
    uint8  temp[DSIZE];
    uint8 *stack[SMAX];

    while (1)
      { e = arrow[0];
        for (o = RSIZE; o < asize; o += RSIZE)
          { if (arrow[o] != e)
              break;
          }
        if (o < asize)
          break;

        digit += 1;
        if (digit >= KSIZE)
          return;
        arrow += 1;
      }

    for (x = 0; x < 256; x++)
      len[x] = 0;

    len[e] = o;
    for (; o < asize; o += RSIZE)
      len[arrow[o]] += RSIZE;

    u = arrow;
    for (x = 0; x < 256; x++)
      { off[x] = u;
        end[x] = u += len[x];
      }

    rems = DSIZE-digit;
    for (x = 0; x < 256; x++)
      { uint8   *p;
        int      t, s;
        int      z;

        while (off[x] < end[x])
          { t = *off[x];

            if (t == x)
              off[x] += RSIZE;
	    else
              { s = 0;
                stack[s++] = off[x];
                while (s < SMAX)
                  if (t == x)
                    { off[x] += RSIZE;
                      break;
                    }
                  else
                    { u = off[t];
                      while ((z = *u) == t)
                        u += RSIZE;
                      off[t] = u+RSIZE;
                      stack[s++] = u;
                      t = z;
                    }

                u = stack[--s];
                mycpy(temp,u,rems);
	        while (s > 0)
                  { p = stack[--s];
                    mycpy(u,p,rems);
                    u = p;
                  }
                mycpy(u,temp,rems);
              }
          }
      }
  }
  
  digit += 1;

  if (digit < KSIZE)
    for (x = 0; x < 256; x++)
      { n = len[x];
        if (n > S_thr0)
          radix_sort(array, n, digit);
        else if (n > RSIZE)
          shell_sort(array, n, digit);
        array += n;
      }
}

typedef struct
  { uint8 *array;
    int    npart;
    int64 *parts;
    int64  offset;
  } Arg;

static int64 *Parts;

static void *sort_thread(void *arg) 
{ Arg *param = (Arg *) arg;

  uint8   *array = param->array;
  int      npart = param->npart;
  int64   *parts = param->parts;
  int64    off;
  int      x;

  if (KSIZE <= 0)
    return (NULL);

  off = param->offset;
  for (x = 0; x < npart; x++)
    { if (parts[x] == 0)
        continue;

#ifdef DEBUG
      printf("Thread %3ld: %12lld - %12lld\n",x+(parts-Parts),off,off+parts[x]);
#endif

      radix_sort(array + off, parts[x], 1);

      off += parts[x];
    }

  return (NULL);
}

void MSD_Sort(uint8 *array, int64 nelem, int rsize, int dsize, int ksize, int64 *part, int nthreads)
{ pthread_t   threads[nthreads];
  Arg         parms[nthreads];

  int   x, n, beg;
  int64 sum, thr, off;
  int64 asize;

  asize    = nelem*rsize;

  Parts    = part;
  RSIZE    = rsize;
  DSIZE    = dsize;
  KSIZE    = ksize;
  PSIZE    = rsize - ksize;

  COFF     = KSIZE-2;
  CMAX     = RSIZE*0x7fff;

  S_thr0 = THR0*RSIZE;
  S_thr1 = THR1*RSIZE;
  S_thr2 = THR2*RSIZE;
  S_gap1 = GAP1*RSIZE;
  S_gap2 = GAP2*RSIZE;

  n   = 0;
  thr = asize / nthreads;
  off = 0;
  sum = 0;
  beg = 0;
  for (x = 0; x < 256; x++)
    { sum += part[x];
      if (sum >= thr)
        { parms[n].array  = array;
          parms[n].npart  = (x+1) - beg;
          parms[n].parts  = Parts + beg;
          parms[n].offset = off;
          n  += 1;
          thr = (asize * (n+1))/nthreads;
          beg = x+1;
          off = sum;
        }
    }

#ifdef DEBUG
  for (x = 0; x < nthreads; x++)
    sort_thread(parms+x);
#else
  for (x = 1; x < nthreads; x++)
    pthread_create(threads+x,NULL,sort_thread,parms+x);

  sort_thread(parms);

  for (x = 1; x < nthreads; x++)
    pthread_join(threads[x],NULL);
#endif

#ifdef DEBUG
  sum = 0;
  for (x = 0; x < 256; x++)
    { sorted(array+sum,part[x],0);
      sum += part[x];
    }
#endif
}
