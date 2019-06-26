/*******************************************************************************************
 *
 *  VGPzip [-v] [-T<int(4)>] [-C<int(6)>] <input>
 *
 *  Block compression of a file with index
 *
 *  Author:  Gene Myers
 *  Date  :  June 2019
 *
 ********************************************************************************************/

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

#define IN_BLOCK  10000000

#define GZIP_HEAD 10
#define GZIP_TAIL  8

static int64 MAXOUT;
static int64 OUT_BLOCK;
static int64 SEEK_STEP;

static char *Usage = "[-v] [-T<int(4)>] [-C<int(6)>] <input>";

typedef struct
  { int   inp;     //  Input file descriptor (independent for each thread even though same file)
    uint8 *in;     //  Input buffer (IO_BLOCK bytes)
    uint8 *out;    //  Output buffer (MAXOUT bytes)
    int64  seek;   //  Location in file to start next read+compress
    uint32 dlen;   //  Size of compressed block
    int    eof;    //  Hit end of file on last go
  } Deflate_Arg;

static uint8 gzip_header[GZIP_HEAD] = "\037\213\010\0\0\0\0\0\0\377";

static int compress3(uint8 *dest, uint32 *destLen, uint8 *source, uint32 sourceLen, int level)
{ z_stream stream;
  int      err;

  stream.zalloc = (alloc_func) NULL;
  stream.zfree  = (free_func) NULL;
  stream.opaque = (voidpf) NULL;

  err = deflateInit2(&stream,level,Z_DEFLATED,-15,9,Z_DEFAULT_STRATEGY);
  if (err != Z_OK)
    return (err);

  stream.next_out  = dest;
  stream.avail_out = *destLen;
  stream.next_in   = source;
  stream.avail_in  = sourceLen;

  err = deflate(&stream,Z_FINISH);
  if (err != Z_STREAM_END)
    return (err);

  *destLen = stream.total_out;
  deflateEnd(&stream);
  return (Z_OK);
}

static inline void pack32(uint32 x, uint8 *m)
{ m[0] = x;
  m[1] = x >> 8;
  m[2] = x >> 16;
  m[3] = x >> 24;
}

static void *deflate_thread(void *arg)
{ Deflate_Arg *data  = (Deflate_Arg *) arg;
  int          input = data->inp;
  uint8       *in    = data->in;
  uint8       *out   = data->out;
  uint32       crc;
  uint32       dlen, rlen;

  lseek(input,data->seek,SEEK_SET);
  rlen = read(input,in,IN_BLOCK);
  if (rlen > 0)
    { dlen = MAXOUT;
      if (compress3(out+GZIP_HEAD,&dlen,in,rlen,CLEVEL) != Z_OK)
        { fprintf(stderr,"Compression not OK\n");
          exit (1);
        }
      dlen += GZIP_HEAD;
      crc = crc32(0L,in,rlen);
      pack32(crc,out+dlen);
      pack32(rlen,out+dlen+4);
      data->dlen = dlen + GZIP_TAIL;
      if (rlen < IN_BLOCK)
        data->eof = 1;
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

  //  Open output files

  { char *fname;
    int   input;

    input = open(argv[1], O_RDONLY);
    if (input < 0)
      { fprintf(stderr,"%s: Cannot open %s for reading\n",Prog_Name,argv[1]);
        exit (1);
      }
    close(input);

    fname = Malloc(strlen(argv[1])+5,"");
    sprintf(fname,"%s.gz",argv[1]);
    output = open(fname, O_WRONLY | O_CREAT, S_IRWXU);

    sprintf(fname,"%s.vzi",argv[1]);
    table = open(fname, O_WRONLY | O_CREAT, S_IRWXU);

    free(fname);
  }

  //  Create NTHREADS independent reader/compressors and have them produce
  //    compressed outputs of the next NTHREAD IO_BLOCK sized chunks at a time
  //    and output the compressed blocks sequentially.

  { pthread_t   thread[NTHREADS];
    Deflate_Arg parm[NTHREADS];
    struct stat stats;

    int64   fsize, isize, *idx;
    uint8  *in, *out;
    int64   b, d;
    int     n;

    MAXOUT    = compressBound(IN_BLOCK);
    OUT_BLOCK = MAXOUT + GZIP_HEAD + GZIP_TAIL;
    SEEK_STEP = NTHREADS * IN_BLOCK;

    in  = Malloc(IN_BLOCK*NTHREADS,"Allocating input buffer");
    out = Malloc(MAXOUT*NTHREADS,"Allocating output buffer");
    if (in == NULL || out == NULL)
      exit (1);

    for (n = 0; n < NTHREADS; n++)
      { parm[n].in   = in + n*IN_BLOCK;
        parm[n].out  = out + n*OUT_BLOCK;
        parm[n].seek = n*IN_BLOCK;
        parm[n].inp  = open(argv[1],O_RDWR);
        parm[n].eof  = 0;
        memcpy(parm[n].out,gzip_header,GZIP_HEAD);
      }

    //  Allocate index

    if (fstat(parm[0].inp, &stats) == -1)
      { fprintf(stderr,"%s: Cannot get stats for %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    fsize = stats.st_size;

    isize = (fsize-1)/IN_BLOCK + 1;
    idx = Malloc(sizeof(int64)*isize,"Allocating table");

    //  Proceed with comnpression and output, accumulating index table

    b = 0;
    d = 0;
    while (1)
      { for (n = 0; n < NTHREADS; n++)
          pthread_create(thread+n,NULL,deflate_thread,parm+n);

        for (n = 0; n < NTHREADS; n++)
          { pthread_join(thread[n],NULL);
            write(output,parm[n].out,parm[n].dlen);
            parm[n].seek += SEEK_STEP;
            if (parm[n].seek > fsize)
              parm[n].seek = fsize;
            d += parm[n].dlen;
            idx[b++] = d;
          }

        if (parm[NTHREADS-1].eof)
          break;
      }

    for (n = 0; n < NTHREADS; n++)
      close(parm[n].inp);

    //  Output index

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
