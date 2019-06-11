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

#undef DEBUG

#define SHELL 24
#define SMAX   6
#define MOVE(A,B) memcpy(A,B,rsize)

static int64 Shell;

#ifdef DEBUG

static inline void sorted(uint8 *array, int asize, int rsize, int ksize)
{ int p;
  for (p = rsize; p < asize; p += rsize)
    if (memcmp(array+(p-rsize),array+p,ksize) > 0)
      printf("Not sorted\n");
}

#endif

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
{ if (memcmp(array,array+(asize-rsize),ksize) == 0)
    return;
  gap_sort(array,asize,rsize,ksize,10);
  gap_sort(array,asize,rsize,ksize,4);
  gap_sort(array,asize,rsize,ksize,1);
}

typedef struct
  { uint8 *array;
    int64  asize;
    int    rsize;
    int    ksize;
    int64  offset;
    int    npart;
    int64 *parts;
    int64  len[256];
  } Arg;

static pthread_t *Threads;
static Arg       *Parms;
static int        Nthreads;
static int64      Parts[256];

static void *split_thread(void *arg) 
{ Arg *param = (Arg *) arg;

  uint8 *array = param->array;
  int64  asize = param->asize;
  int    rsize = param->rsize;
  int64 *len   = param->len;

  int   x;
  int64 o;

  for (x = 0; x < 256; x++)
    len[x] = 0;
  for (o = 0; o < asize; o += rsize)
    len[array[o]] += rsize;

  return (NULL);
}

static void radix_sort(uint8 *array, int64 asize, int rsize, int ksize, int digit)
{ int64  len[256];
  int64  beg[256];
  int    x;

  { int64  end[256];
    uint8 *darray = array + digit;
    int64  o;

    int64  off[256];
    uint8  temp[rsize];
    int64  stack[SMAX];

    for (x = 0; x < 256; x++)
       len[x] = 0;

    if (digit == 0)
      { int64 s;
        int   t;

        s = ((asize/rsize)/Nthreads)*rsize;
        o = 0;
        for (t = 0; t < Nthreads; t++)
          { Parms[t].array = array + o;
            Parms[t].asize = s;
            Parms[t].rsize = rsize;
            o += s;
          }
        Parms[Nthreads-1].asize = asize - (o-s);

        for (t = 0; t < Nthreads; t++)
          pthread_create(Threads+t,NULL,split_thread,Parms+t);

        for (t = 0; t < Nthreads; t++)
          pthread_join(Threads[t],NULL);

        for (t = 0; t < Nthreads; t++)
          for (x = 0; x < 256; x++)
            len[x] += Parms[t].len[x];
      }
    else
      for (o = 0; o < asize; o += rsize)
        len[darray[o]] += rsize;

    o = 0;
    for (x = 0; x < 256; x++)
      { beg[x] = off[x] = o;
        end[x] = o += len[x];
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
  }

  digit += 1;
  if (digit >= ksize)
    return;

  if (digit == 1)
    { { int   t, beg;
        int64 sum, thr, off;

        t   = 0;
        thr = asize / Nthreads;
        off = 0;
        sum = 0;
        beg = 0;
        for (x = 0; x < 256; x++)
          { Parts[x] = len[x];
            sum += len[x];
            if (sum >= thr)
              { Parms[t].offset = off;
                Parms[t].npart  = (x+1) - beg;
                Parms[t].parts  = Parts + beg;
                Parms[t].array  = array;
                Parms[t].asize  = asize;
                Parms[t].rsize  = rsize;
                Parms[t].ksize  = ksize;
                t  += 1;
                thr = (asize * (t+1))/Nthreads;
                beg = x+1;
                off = sum;
              }
          }
      }
    }
  else
    { for (x = 0; x < 256; x++)
        if (len[x] > Shell)
          radix_sort(array + beg[x], len[x], rsize, ksize, digit);
        else if (len[x] > rsize)
          shell_sort(array + beg[x], len[x], rsize, ksize);
    }
}

static void *sort_thread(void *arg) 
{ Arg *param = (Arg *) arg;

  uint8 *array = param->array;
  int    rsize = param->rsize;
  int    ksize = param->ksize;
  int    npart = param->npart;
  int64 *parts = param->parts;

  int64  off;
  int    x;

  off = param->offset;
  for (x = 0; x < npart; x++)
    {
#ifdef DEBUG
      if (parts[x] > 0)
        printf("Thread %3ld: %12lld - %12lld\n",x+(parts-Parts),off,off+parts[x]);
#endif
      if (parts[x] > Shell)
        radix_sort(array + off, parts[x], rsize, ksize, 1);
      else if (parts[x] > rsize)
        shell_sort(array + off, parts[x], rsize, ksize);
      off += parts[x];
    }

  return (NULL);
}

void Ex_sort(char *path, int rsize, int ksize, int nthreads)
{ pthread_t threads[nthreads];
  Arg       parms[nthreads];

  int         fd;
  uint8      *array;
  int64       asize;

  fd = open(path, O_RDWR);
  if (fd == -1)
    { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,path);
      exit (1);
    }

  { struct stat stats;
    if (fstat(fd, &stats) == -1)
      { close(fd);
        fprintf(stderr,"%s: Cannot get stats for %s\n",Prog_Name,path);
        exit (1);
      }
    asize = stats.st_size;
  }

  array = (uint8 *) mmap(NULL,asize,PROT_READ | PROT_WRITE,MAP_PRIVATE,fd,0);
  if (array == NULL)
    { close(fd);
      fprintf(stderr,"%s: Cannot memory map %s\n",Prog_Name,path);
      exit (1);
    }

  madvise(array, asize, POSIX_MADV_WILLNEED | POSIX_MADV_SEQUENTIAL);

  Threads  = threads;
  Parms    = parms;
  Shell    = SHELL * rsize;
  Nthreads = nthreads;

  radix_sort(array, asize, rsize, ksize, 0);

  if (ksize > 1)
    { int x;

      for (x = 0; x < Nthreads; x++)
        pthread_create(Threads+x,NULL,sort_thread,Parms+x);

      for (x = 0; x < Nthreads; x++)
        pthread_join(Threads[x],NULL);
    }

#ifdef DEBUG
  sorted(array,asize,rsize,ksize);
#endif

  { int   fe;
    int64 t;

    unlink(path);

    fe = open(path, O_WRONLY | O_CREAT, S_IRWXU);
    if (fe == -1)
      { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,path);
        exit (1);
      }

    for (t = 0; t < asize; t += 0x40000000ll)
      { if (asize-t < 0x40000000ll)
          write(fe,array+t,asize-t);
        else
          write(fe,array+t,0x40000000ll);
      }

    close(fe);

    munmap(array,asize);
  }
}
