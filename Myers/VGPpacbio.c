/*******************************************************************************************
 *
 *  Dextract: pullls requested info out of subreads.[bs]am and .bax.h5 files produced by
 *               Pacbio sequencing instruments and software
 *
 *  Author:  Gene Myers
 *  Date  :  Oct. 9, 2016
 *
 ********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include "gene_core.h"
#include "pb_sam.h"
#include "pb_expr.h"

#define LOWER_OFFSET 32

static int VERBOSE;
static int ARROW;
static int UPPER;

static char *Usage = "[-vaU] [-e<expr(ln>=500 && rq>=750)> <input:pacbio> ...";

  //  Write subread data in samRecord rec to non-NULL file types

static void writeSamRecord(samRecord *rec)
{ int i;

  printf("L %d %d %d 0.%0d\n",rec->well,rec->beg,rec->end,(int) (rec->qual*1000.));

  if (UPPER)
    { if (islower(rec->seq[0]))
      for (i = 0; i < rec->len; i++)
        rec->seq[i] -= LOWER_OFFSET;
    }
  else
    { if (isupper(rec->seq[0]))
      for (i = 0; i < rec->len; i++)
        rec->seq[i] += LOWER_OFFSET;
    }

  printf("S %d %.*s\n",rec->len,rec->len,rec->seq);

  if (ARROW)
    { printf("N %.2f %.2f %.2f %.2f\n",rec->snr[0],rec->snr[1],rec->snr[2],rec->snr[3]);
      printf("A %d %.*s\n",rec->len,rec->len,rec->arr);
    }
}

  //  Main

int main(int argc, char* argv[])
{ char   *path, *core;
  Filter *EXPR;

  //  Process command line arguments

  { int   i, j, k;
    int   flags[128];

    ARG_INIT("dextract")

    path   = NULL;
    core   = NULL;
    EXPR   = NULL;

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
  }

  //  Scan each file for size statistics

  { int      i;
    samFile *in;
    int      nread, rg_nread, gmaxc;
    int64    totbp, rg_totbp, gtotc;
    int      maxbp;
    int      gsize[argc];

    totbp    = 0;
    nread    = 0;
    maxbp    = 0;
    rg_nread = 0;
    rg_totbp = 0;
    gtotc    = 0;
    gmaxc    = 0;
    for (i = 1; i < argc; i++)
      { FILE *file;
        int   status, intype;
        int   nr, bp;

        //  Determine file type

#define IS_BAM 1
#define IS_SAM 2

        path  = PathTo(argv[i]);
        core  = Root(argv[i],".subreads.bam");
        if ((file = fopen(Catenate(path,"/",core,".subreads.bam"),"r")) == NULL)
          { core  = Root(argv[i],".subreads.sam");
            if ((file = fopen(Catenate(path,"/",core,".subreads.sam"),"r")) == NULL)
              { fprintf(stderr,"%s: Cannot find %s/%s with a Pacbio extension\n",
                               Prog_Name,path,core);
                goto error;
              }
            else
              intype = IS_SAM;
          }
        else
          intype = IS_BAM;
        fclose(file);

        //  Extract from a .bam or .sam

        if (VERBOSE)
          { fprintf(stderr, "Scanning file : %s ...\n", core); fflush(stderr); }

        if (intype == IS_BAM)
          { if ((in = sam_open(Catenate(path,"/",core,".subreads.bam"))) == NULL)
              { fprintf(stderr, "%s: can't open %s as a Bam file\n", Prog_Name, argv[i]);
                goto error;
              }
          }
        else
          { if ((in = sam_open(Catenate(path,"/",core,".subreads.sam"))) == NULL)
              { fprintf(stderr, "%s: can't open %s as a Sam file\n", Prog_Name, argv[i]);
                goto error;
              }
          }

        status = sam_header_process(in,0);
        if (status < 0)
          goto error;
        else if ((status & HASPW) == 0 && ARROW)
          { fprintf(stderr, "%s: %s does not have Arrow information\n", Prog_Name, argv[i]);
            goto error;
          }

        { samRecord *rec;
          int        first, hlen;

          first = 1;
          nr = 0;
          bp = 0;
          while (1)
            { rec = sam_scan(in);
              if (rec == NULL)
                goto error;
              if (rec == SAM_EOF)
                break;

              if ( ! evaluate_bam_filter(EXPR,rec))
                continue;

              if (first)
                { first = 0;
                  hlen  = strlen(rec->header);
                }
 
              nr += 1;
              bp += rec->len;
              if (rec->len > maxbp)
                maxbp = rec->len;
            }

          if (bp > rg_totbp)
            rg_totbp = bp;
          if (nr > rg_nread)
            rg_nread = nr;
          totbp += bp;
          nread += nr;
          gsize[i] = nr;
          if (!first)
            { gtotc += hlen;
              if (hlen > gmaxc)
                gmaxc = hlen;
            }
        }

        if (sam_close(in))
          { fprintf(stderr, "%s: Error closing file %s\n", Prog_Name, core);
            goto error;
          }

        free(path);
        free(core);
      }

    //  Output size headers

    { int    i, clen, optl;
      char   date[26];
      time_t seconds;

      optl = ARROW+UPPER+VERBOSE;
      if (optl > 0)
        clen = optl+1;
      else
        clen = -1;
      for (i = 1; i < argc; i++)
        clen += strlen(argv[i])+1;

      printf("1 3 seq 1 0\n");
      printf("2 3 pbr\n");
      printf("# ! 1\n");
      printf("# g %d\n",argc-1);
      printf("# L %d\n",nread);
      printf("# S %d\n",nread);
      if (ARROW)
        { printf("# N %d\n",nread);
          printf("# A %d\n",nread);
        }
      printf("+ ! %d\n",clen+36);
      printf("+ g %lld\n",gtotc);
      printf("+ S %lld\n",totbp);
      if (ARROW)
        printf("+ A %lld\n",totbp);

      if (clen > 24)
        printf("@ ! %d\n",clen);
      else
        printf("@ ! 24\n");
      printf("@ g %d\n",gmaxc);
      printf("@ S %d\n",maxbp);
      if (ARROW)
        printf("@ A %d\n",maxbp);

      printf("%% g # L %d\n",rg_nread);
      printf("%% g # S %d\n",rg_nread);
      if (ARROW)
        { printf("%% g # N %d\n",rg_nread);
          printf("%% g # A %d\n",rg_nread);
        }

      printf("%% g + S %lld\n",rg_totbp);
      if (ARROW)
        printf("%% g + A %lld\n",rg_totbp);

      printf("\n! 9 VGPpacbio 3 1.0 %d",clen);
      for (i = 1; i < argc; i++)
        printf(" %s",argv[i]);
      seconds = time(NULL);
      ctime_r(&seconds,date);
      date[24] = '\0';
      printf(" 24 %s\n",date);
    }

    //  Scan files and output .pbr

    for (i = 1; i < argc; i++)
      { FILE *file;
        int   status, intype;

        //  Determine file type

        path  = PathTo(argv[i]);
        core  = Root(argv[i],".subreads.bam");
        if ((file = fopen(Catenate(path,"/",core,".subreads.bam"),"r")) == NULL)
          { core  = Root(argv[i],".subreads.sam");
            intype = IS_SAM;
          }
        else
          intype = IS_BAM;
        fclose(file);

        //  Extract from a .bam or .sam

        if (VERBOSE)
          { fprintf(stderr, "Processing file : %s ...\n", core); fflush(stderr); }

        if (intype == IS_BAM)
          in = sam_open(Catenate(path,"/",core,".subreads.bam"));
        else
          in = sam_open(Catenate(path,"/",core,".subreads.sam"));

        status = sam_header_process(in,0);

        { samRecord *rec;
          int        first;

          first = 1;
          while (1)
            { rec = sam_record_extract(in, status);
              if (rec == SAM_EOF)
                break;

              if ( ! evaluate_bam_filter(EXPR,rec))
                continue;

              if (first)
                { first = 0;
                  printf("\ng %d %ld %s\n",gsize[i],strlen(rec->header),rec->header);
                }

              writeSamRecord(rec);
            }
        }

        sam_close(in);

        free(path);
        free(core);

        if (VERBOSE)
          { fprintf(stderr, "Done\n"); fflush(stdout); }
      }
  }

  exit (0);

  //  An error occured, carefully undo any files in progress

error:
  free(path);
  free(core);
  exit (1);
}
