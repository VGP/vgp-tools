#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <time.h>
#include <zlib.h>

#include "gene_core.h"

int    VERBOSE;  //  Verbose mode?
int    SEQ_ONLY; //  Do not include QV strings
int    GROUP;    //  Group lanes of reads

static char *Usage = "[-vsg] <forward:fastq> <reverse:fastq>";

//  Read next line into a buffer and return a pointer to the buffer and set *plen
//    the length of the line.  NB: replaces '\n' with '\0'.

static char *read_line(gzFile input, int *plen)
{ static char *buffer;
  static int   bmax = 0;
  int len;

  if (bmax == 0)
    { bmax = 500;
      buffer = (char *) Malloc(bmax,"Allocating read buffer");
      if (buffer == NULL)
        exit (1);
    }

  if (gzgets(input,buffer,bmax) == NULL)
    return (NULL);

  len = strlen(buffer);
  while (buffer[len-1] != '\n')
    { bmax = ((int) (1.4*bmax)) + 100;
      buffer = (char *) Realloc(buffer,bmax,"Reallocating read buffer");
      if (buffer == NULL)
        exit (1);
      if (gzgets(input,buffer+len,bmax-len) == NULL)
        { fprintf(stdout,"%s: Last line of file does not end with new-line\n",Prog_Name);
          exit (1);
        }
      len += strlen(buffer+len);
    }
  buffer[--len] = '\0';

  if (plen != NULL)
    *plen = len;
  return (buffer);
}

  //  In a first pass perform basic .fastq syntax check, count the number of
  //    of reads, # of base pairs in all reads, and maximum length of a sequence

typedef struct
  { int64  nread;
    int    maxlen;
    int64  totbp;
    int    ngroup;
    int64  gmaxnr;
    int64  gmaxbp;
    int64 *gsize;
    int64  gmaxc;
    int64  gtotc;
  } Header_Info;

static void Check_FastQ(gzFile input, char *file_name, Header_Info *info)
{ char   *line;
  int64   lineno;
  int     len, qlen;
  int     maxlen, ngroup, maxgroup;
  int64   nr, nread, gmaxnr, *gcount;
  int64   bp, totbp, gmaxbp, gtotc, gmaxc;
  char   *readgroup, *p;
  int     i;

  readgroup = NULL;
  maxgroup  = 0;
  gcount    = NULL;

  lineno  = 0;
  nread   = 0;
  totbp   = 0;
  maxlen  = 0;
  ngroup  = 0;
  gmaxnr  = 0;
  gmaxbp  = 0;
  gtotc   = 0;
  gmaxc   = 0;

  nr = 0;
  bp = 0;
  while ((line = read_line(input,NULL)) != NULL)
    { lineno += 1;
      if (line[0] != '@')
        { fprintf(stdout,"%s: %s.%lld: Entry header does not start with an @-sign\n",
                         Prog_Name,file_name,lineno);
          exit (1);
        }

      if (GROUP)
        { p = line+1;
          for (i = 0; i < 4; i++)
            { p = index(p+1,':');
              if (p == NULL)
                { fprintf(stdout,"%s: %s.%lld: Entry header does not appear to be",
                                 Prog_Name,file_name,lineno);
                  fprintf(stdout," an Illumina header\n");
                  exit (1);
                }
            }
          *p = '\0';

          if (readgroup == NULL || strcmp(readgroup,line) != 0)
            { free(readgroup);
    
              if (readgroup != NULL)
                { totbp += bp;
                  nread += nr;
                  if (bp > gmaxbp)
                    gmaxbp = bp;
                  if (nr > gmaxnr)
                    gmaxnr = nr; 
                  if (ngroup >= maxgroup)
                    { maxgroup = 1.2*ngroup + 10;
                      gcount   = (int64 *) Realloc(gcount,sizeof(int64)*maxgroup,
                                                   "Allocating group read counts");
                    }
                  gcount[ngroup++] = nr;

                  len = strlen(readgroup+1);
                  gtotc += len;
                  if (len > gmaxc)
                    gmaxc = len;
                }
    
              readgroup = strdup(line);
              bp = 0;
              nr = 0;
            } 
        }

      line = read_line(input,&len);
      lineno += 1;
      if (line == NULL)
        { fprintf(stdout,"%s: %s.%lld: Incomplete fastq entry\n",Prog_Name,file_name,lineno);
          exit (1);
        }

      nr  += 1;
      bp  += len;
      if (len > maxlen)
        maxlen = len;

      line = read_line(input,NULL);
      lineno += 1;
      if (line == NULL)
        { fprintf(stdout,"%s: %s.%lld Incomplete fastq entry\n",Prog_Name,file_name,lineno);
          exit (1);
        }
      if (line[0] != '+')
        { fprintf(stdout,"%s: %s.%lld Divider line does not start with a +-sign\n",
                         Prog_Name,file_name,lineno);
          exit (1);
        }
      line = read_line(input,&qlen);
      lineno += 1;
      if (line == NULL)
        { fprintf(stdout,"%s: %s.%lld Incomplete fastq entry\n",Prog_Name,file_name,lineno);
          exit (1);
        }
      if (len != qlen)
        { fprintf(stdout,"%s: %s.%lld: QV line does not have the same length as sequence line\n",
                         Prog_Name,file_name,lineno);
          exit (1);
        }
    }

  totbp += bp;
  nread += nr;
  if (GROUP && readgroup != NULL)
    { if (bp > gmaxbp)
        gmaxbp = bp;
      if (nr > gmaxnr)
        gmaxnr = nr; 
      if (ngroup >= maxgroup)
        { maxgroup = 1.2*ngroup + 10;
          gcount   = (int64 *) Realloc(gcount,sizeof(int64)*maxgroup,
                                       "Allocating group read counts");
        }
      gcount[ngroup++] = nr;

      len = strlen(readgroup+1);
      gtotc += len;
      if (len > gmaxc)
        gmaxc = len;
    }

  info->nread  = nread;
  info->maxlen = maxlen;
  info->totbp  = totbp;
  info->ngroup = ngroup;
  info->gmaxnr = gmaxnr;
  info->gmaxbp = gmaxbp;
  info->gsize  = gcount;
  info->gtotc  = gtotc;
  info->gmaxc  = gmaxc;

  free(readgroup);
}

  //  In a second pass, output .irp format for each pair

static void Output_Pairs(gzFile forward, gzFile reverse, int64 *gcount)
{ char *line;
  int   len;
  int   DO_QV;
  char *p, *readgroup;
  int   i, gc;

  DO_QV = ! SEQ_ONLY;

  gc = 0;
  readgroup = NULL;
  while ((line = read_line(forward,NULL)) != NULL)
    { if (GROUP)
        { p = line+1;
          for (i = 0; i < 4; i++)
            p = index(p+1,':');
          *p = '\0';

          if (readgroup == NULL || strcmp(readgroup,line) != 0)
            { free(readgroup);
              readgroup = strdup(line);
              printf("g %lld %ld %s\n",gcount[gc++],strlen(readgroup+1)+1,readgroup+1);
            } 
        }

      printf("P\n");
      line = read_line(forward,&len);
      printf("S %d %s\n",len,line);
      line = read_line(forward,NULL);
      line = read_line(forward,NULL);
      if (DO_QV)
        printf("Q %d %s\n",len,line);

      line = read_line(reverse,NULL);
      line = read_line(reverse,&len);
      printf("S %d %s\n",len,line);
      line = read_line(reverse,NULL);
      line = read_line(reverse,NULL);
      if (DO_QV)
        printf("Q %d %s\n",len,line);
    }
}

int main(int argc, char *argv[])
{ gzFile      input1, input2;
  char       *fname1, *fname2;
  Header_Info stats1, stats2;

  //  Output file type and provenance

  { int    i, len;
    char   date[26];
    time_t seconds;

    printf("1 3 seq 1 0\n");
    printf("2 3 irp\n");
    printf("# ! 1\n");
    len = -1;
    for (i = 1; i < argc; i++)
      len += strlen(argv[i]);
    printf("+ ! %d\n",len+34);
    if (len > 24)
      printf("@ ! %d\n",len);
    else
      printf("@ ! 24\n");
    printf("! 7 VGPpair 3 1.0 %d",len);
    for (i = 1; i < argc; i++)
      printf(" %s",argv[i]);
    seconds = time(NULL);
    ctime_r(&seconds,date);
    date[24] = '\0';
    printf(" 24 %s\n",date);
  }

  //  Parse command line options

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("VGPpair")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vsg")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE  = flags['v'];
    SEQ_ONLY = flags['s'];
    GROUP    = flags['g'];

    if (argc != 3)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, show progress as proceed.\n");
        fprintf(stderr,"      -s: Output sequences only, skip QVs\n");
        exit (1);
      }
  }

  //  Open input files

  { char *pwd1, *pwd2;

    fname1 = Root(argv[1],".fastq.gz");
    pwd1   = PathTo(argv[1]);
    input1 = gzopen(Catenate(pwd1,"/",fname1,".fastq"),"r");
    if (input1 == NULL)
      { fname1 = Root(argv[1],".fastq");
        input1 = gzopen(Catenate(pwd1,"/",fname1,".fastq"),"r");
        if (input1 == NULL)
          { fprintf(stdout,"%s: Cannot open %s either as a .fastq\n",Prog_Name,argv[1]);
            exit (1);
          }
      }

    fname2 = Root(argv[2],".fastq.gz");
    pwd2   = PathTo(argv[2]);
    input2 = gzopen(Catenate(pwd2,"/",fname2,".fastq"),"r");
    if (input2 == NULL)
      { fname2 = Root(argv[2],".fastq");
        input2 = gzopen(Catenate(pwd2,"/",fname2,".fastq"),"r");
        if (input2 == NULL)
          { fprintf(stdout,"%s: Cannot open %s either as a .fastq\n",Prog_Name,argv[2]);
            exit (1);
          }
      }

    if (strcmp(fname1,fname2) == 0)
      { char *aname1, *aname2;

        if (strcmp(pwd1,pwd2) == 0)
          { fprintf(stderr,"%s: Forward and reverse files are the same file?\n",Prog_Name);
            exit (1);
          }
        aname1 = Strdup(Catenate(pwd1,"/",fname1,""),NULL);
        aname2 = Strdup(Catenate(pwd2,"/",fname2,""),NULL);
        free(fname2);
        free(fname1);
        fname1 = aname1;
        fname2 = aname2;
      }
    free(pwd2);
    free(pwd1);
  }

  //  Scan 1: check syntax and count

  if (VERBOSE)
    { fprintf(stderr,"  Checking syntax of forward file %s\n",argv[1]);
      fflush(stderr);
    }

  Check_FastQ(input1,fname1,&stats1);

  if (VERBOSE)
    { fprintf(stderr,"  Checking syntax of reverse file %s\n",argv[2]);
      fflush(stderr);
    }

  Check_FastQ(input2,fname2,&stats2);

  if (stats1.nread != stats2.nread)
    { fprintf(stdout,"%s: Number of reads in two files not equal: %lld vs %lld\n",
                     Prog_Name,stats1.nread,stats2.nread);
      exit (1);
    }

  //  Output .irp sizes

  printf("# P %lld\n",stats1.nread);
  printf("# S %lld\n",2*stats1.nread);
  printf("+ S %lld\n",stats1.totbp + stats2.totbp);
  if (stats1.maxlen > stats2.maxlen)
    printf("@ S %d\n",stats1.maxlen);
  else
    printf("@ S %d\n",stats2.maxlen);

  if ( ! SEQ_ONLY)
    { printf("# Q %lld\n",2*stats1.nread);
      printf("+ Q %lld\n",stats1.totbp + stats2.totbp);
      if (stats1.maxlen > stats2.maxlen)
        printf("@ Q %d\n",stats1.maxlen);
      else
        printf("@ Q %d\n",stats2.maxlen);
    }

  if (GROUP)
    { printf("# g %d\n",stats1.ngroup);
      printf("+ g %lld\n",stats1.gtotc);
      printf("@ g %lld\n",stats1.gmaxc);
      printf("%% g # P %lld\n",stats1.gmaxnr);
      printf("%% g # S %lld\n",2*stats1.gmaxnr);
      printf("%% g + S %lld\n",stats1.gmaxbp + stats2.gmaxbp);
      if ( ! SEQ_ONLY)
        { printf("%% g # Q %lld\n",2*stats1.gmaxnr);
          printf("%% g + Q %lld\n",stats1.gmaxbp + stats2.gmaxbp);
        }
    }

  //  Scan 2: output the data in .irp format

  if (VERBOSE)
    { fprintf(stderr,"  Outputting read pairs in .irp format to stdout\n");
      fflush(stderr);
    }

  gzrewind(input1);
  gzrewind(input2);
  Output_Pairs(input1,input2,stats1.gsize);

  //  Tidy up just for good form

  gzclose(input2);
  gzclose(input1);
  free(fname2);
  free(fname1);
  exit (0);
}
