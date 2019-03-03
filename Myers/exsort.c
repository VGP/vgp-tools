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
#include "exsort.h"

#define SHELL 24
#define SMAX  20
#define MOVE(A,B) memcpy(A,B,rsize)

static int64 Shell;

static inline void gap_sort(uint8 *array, int asize, int rsize, int ksize, int gap)
{ int    i, j, step;
  uint8  temp[rsize];
  uint8 *garray;

  step   = gap*rsize;
  garray = array + step;
  for (i = step; i < asize; i += rsize)
    { j = i-step;
      if (memcmp(array+j,array+i,ksize) <= 0)
        continue;
      MOVE(temp,array+i);
      MOVE(array+i,array+j);
      for(j -= step; j >= 0; j -= step)
        { if (memcmp(array+j,temp,ksize) <= 0)
            break;
          MOVE(garray+j,array+j);
        }
      MOVE(garray+j,temp);
    }
}

static inline void shell_sort(uint8 *array, int asize, int rsize, int ksize)
{ gap_sort(array,asize,rsize,ksize,10);
  gap_sort(array,asize,rsize,ksize,4);
  gap_sort(array,asize,rsize,ksize,1);
}

typedef struct
  { uint8 *array;
    int    rsize;
    int    ksize;
    int    npart;
    int64 *parts;
  } Arg;

static int        Nthreads;
static pthread_t *Threads;
static Arg       *Parms;
static int64      Parts[256];
static void       radix_sort(uint8 *, int64, int, int, int);

static void *sort_thread(void *arg) 
{ Arg *param = (Arg *) arg;

  uint8 *array = param->array;
  int    rsize = param->rsize;
  int    ksize = param->ksize;
  int    npart = param->npart;
  int64 *parts = param->parts;
  int64  off;
  int    x;

  off = 0;
  for (x = 0; x < npart; x++)
    { if (parts[x] > Shell)
        radix_sort(array + off, parts[x], rsize, ksize, 1);
      else if (parts[x] > rsize)
        shell_sort(array + off, parts[x], rsize, ksize);
      off += parts[x];
    }

  return (NULL);
}

static void radix_sort(uint8 *array, int64 asize, int rsize, int ksize, int digit)
{ int64  len[256];
  int64  beg[256];
  int64  end[256];
  uint8 *darray = array + digit;

  int64  off[256];
  uint8  temp[rsize];
  int64  stack[SMAX];

  int x;

  if (digit == 0)
    fprintf(stderr, "radixify(asize=%lldd, digit=%d, rsize=%d, ksize=%d)\n",
                   asize,digit,rsize,ksize);

  { int64 o;

    for (x = 0; x < 256; x++)
      len[x] = 0;
    for (o = 0; o < asize; o += rsize)
      len[darray[o]] += rsize;

    o = 0;
    for (x = 0; x < 256; x++)
      { beg[x] = off[x] = o;
        end[x] = o += len[x];
      }
  }

  for (x = 0; x < 256; x++)
    { uint8 *o, *p;
      int    t, s;
      int64  u;

      while (off[x] < end[x])
        { t = darray[off[x]];
          if (t == x)
            off[x] += rsize;
          else
            { s = 0;
              stack[s++] = off[x];
              while (s < SMAX)
                { u = off[t];
                  off[t] = u + rsize;
                  if (t == x)
                    break;
                  stack[s++] = u;
                  t = darray[u];
                }

              o = array + stack[--s];
              MOVE(temp,o);
	      while (s > 0)
                { p = array + stack[--s];
                  MOVE(o,p);
                  o = p;
                }
              MOVE(o,temp);
            }
        }
    }
  
  digit += 1;
  if (digit >= ksize)
    return;

  if (digit == 1)
    { int   n, beg;
      int64 sum, thr;

      n   = 0;
      thr = asize / Nthreads;
      sum = 0;
      beg = 0;
      for (x = 0; x < 256; x++)
        { Parts[x] = len[x];
          sum += len[x];
          if (sum >= thr)
            { Parms[n].array = array;
              Parms[n].rsize = rsize;
              Parms[n].ksize = ksize;
              Parms[n].npart = (x+1) - beg;
              Parms[n].parts = Parts + beg;
              n  += 1;
              thr = (asize * (n+1))/Nthreads;
              beg = x+1;
            }
        }

      for (x = 0; x < Nthreads; x++)
        pthread_create(Threads+x,NULL,sort_thread,Parms+x);

      for (x = 0; x < Nthreads; x++)
        pthread_join(Threads[x],NULL);
    }
  else
    { for (x = 0; x < 256; x++)
        if (len[x] > Shell)
          radix_sort(array + beg[x], len[x], rsize, ksize, digit);
        else if (len[x] > rsize)
          shell_sort(array + beg[x], len[x], rsize, ksize);
    }
}


void Ex_sort(char *path, int rsize, int ksize, int nthreads)
{ pthread_t   threads[nthreads];
  Arg         parms[nthreads];

  void       *array = NULL;
  struct stat stats;
  int64       asize;
  int         fd;

  fd = open(path, O_RDWR);
  if (fd == -1)
    { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,path);
      exit (1);
    }

  if (fstat(fd, &stats) == -1)
    { close(fd);
      fprintf(stderr,"%s: Cannot get stats for %s\n",Prog_Name,path);
      exit (1);
    }
  asize = stats.st_size;

  array = mmap(NULL,asize,PROT_READ | PROT_WRITE,MAP_SHARED,fd,0);
  if (array == NULL)
    { close(fd);
      fprintf(stderr,"%s: Cannot memory map %s\n",Prog_Name,path);
      exit (1);
    }

  madvise(array, asize, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);

  Shell    = SHELL * rsize;
  Nthreads = nthreads;
  Threads  = threads;
  Parms    = parms;

  radix_sort(array, asize, rsize, ksize, 0);

  munmap(array,asize);
  close(fd);
}
