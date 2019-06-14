#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <time.h>
#include <fcntl.h>
#include <zlib.h>
#include <pthread.h>
#include <sys/stat.h>

#include "gene_core.h"

int    VERBOSE;  //  Verbose mode?
int    NTHREADS; //  Do not include QV strings
int    CLEVEL;   //  Compression level (in [1,9]);

#define IN_BLOCK  1000000

static int64 MAXOUT;

static char *Usage = "[-v] [-T<int(4)>] [-C<int(6)>] <input>";

typedef struct
  { int   inp;
    uint8 *in;
    uint8 *out;
    int64  seek;
    uLong  dlen;
    int    eof;
  } Deflate_Arg;

static void *deflate_thread(void *arg)
{ Deflate_Arg *data  = (Deflate_Arg *) arg;
  int          input = data->inp;
  uint8       *in    = data->in;
  uint8       *out   = data->out;
  uLong        dlen;
  uint32       rlen;

  lseek(input,data->seek,SEEK_SET);
  rlen = read(input,in,IN_BLOCK);
// printf("Read %d\n",rlen);
  if (rlen > 0)
    { dlen = MAXOUT;
      if (rlen < IN_BLOCK)
        { compress2(out,&dlen,in,rlen,CLEVEL);
// printf("  Comp %d\n",rlen);
          data->eof = 1;
        }
      else
        { compress2(out,&dlen,in,IN_BLOCK,CLEVEL);
// printf("  Comp %d\n",IN_BLOCK);
        }
      data->dlen = dlen;
    }
  else
    { data->dlen = 0;
      data->eof = 1;
    }
  return (NULL);
}

int main(int argc, char *argv[])
{ int output, table;

  //  Parse command line options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("VGPzip")

    NTHREADS = 4;
    CLEVEL   = Z_DEFAULT_COMPRESSION;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("v")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
          case 'C':
            ARG_NON_NEGATIVE(CLEVEL,"Number of threads")
            if (CLEVEL < 0 || CLEVEL > 9)
              { fprintf(stderr,"%s: Compression level must be in [0,9]\n",Prog_Name);
                exit (1);
              }
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE  = flags['v'];

    if (argc != 2)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, show progress as proceed.\n");
        fprintf(stderr,"      -T: Number of threads to use\n");
        exit (1);
      }
  }

  //  Open input files

  { char *fname;

    fname = Malloc(strlen(argv[1])+5,"");
    sprintf(fname,"%s.gz",argv[1]);
    output = open(fname, O_WRONLY | O_CREAT, S_IRWXU);

    sprintf(fname,"%s.vzi",argv[1]);
    table = open(fname, O_WRONLY | O_CREAT, S_IRWXU);

    free(fname);
  }

  { pthread_t   thread[NTHREADS];
    Deflate_Arg parm[NTHREADS];
    struct stat stats;

    int64   fsize, isize, *idx;
    uint8  *in, *out;
    int64   b, d;
    int     n;

    MAXOUT = compressBound(IN_BLOCK);

    in  = Malloc(IN_BLOCK*NTHREADS,"Allocating input buffer");
    out = Malloc(MAXOUT*NTHREADS,"Allocating output buffer");
    if (in == NULL || out == NULL)
      exit (1);

    for (n = 0; n < NTHREADS; n++)
      { parm[n].in   = in + n*IN_BLOCK;
        parm[n].out  = out + n*MAXOUT;
        parm[n].seek = n*IN_BLOCK;
        parm[n].inp  = open(argv[1],O_RDWR);
        parm[n].eof  = 0;
      }

    if (fstat(parm[0].inp, &stats) == -1)
      { fprintf(stderr,"%s: Cannot get stats for %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    fsize = stats.st_size;

    isize = (fsize-1)/IN_BLOCK + 1;
    idx = Malloc(sizeof(int64)*isize,"Allocating table");

    b = 0;
    d = 0;
    while (1)
      { for (n = 0; n < NTHREADS; n++)
          pthread_create(thread+n,NULL,deflate_thread,parm+n);

        for (n = 0; n < NTHREADS; n++)
          { pthread_join(thread[n],NULL);
            write(output,out,parm[n].dlen);
// printf("    Writ %ld\n",dlen);
            parm[n].seek += NTHREADS*IN_BLOCK;
            if (parm[n].seek > fsize)
              parm[n].seek = fsize;
            idx[b++] = d;
            d += parm[n].dlen;
          }

        if (parm[NTHREADS-1].eof)
          break;
      }

    for (n = 0; n < NTHREADS; n++)
      close(parm[n].inp);

    write(table,&isize,sizeof(int64));
    write(table,idx,sizeof(int64)*isize);

    free(idx);
    free(out);
    free(in);
  }

  close(table);
  close(output);
  exit (0);
}
