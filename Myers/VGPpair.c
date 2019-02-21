#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <time.h>

#include "gene_core.h"

int    VERBOSE;  //  Verbose mode?
int    SEQ_ONLY; //  Do not include QV strings
int    GROUP;    //  Group lanes of reads

static char *Usage = "[-vsg] <forward:fastq> <reverse:fastq>";

//  Read next line into a buffer and return a pointer to the buffer and set *plen
//    the length of the line.  NB: replaces '\n' with '\0'.

static char *read_line(FILE *input, int *plen)
{ static char *buffer;
  static int   bmax = 0;
  int len;

  if (bmax == 0)
    { bmax = 500;
      buffer = (char *) Malloc(bmax,"Allocating read buffer");
      if (buffer == NULL)
        exit (1);
    }

  if (fgets(buffer,bmax,input) == NULL)
    return (NULL);

  len = strlen(buffer);
  while (buffer[len-1] != '\n')
    { bmax = ((int) (1.4*bmax)) + 100;
      buffer = (char *) Realloc(buffer,bmax,"Reallocating read buffer");
      if (buffer == NULL)
        exit (1);
      if (fgets(buffer+len,bmax-len,input) == NULL)
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
  { int   nread;
    int   maxlen;
    int64 totbp;
    int   ngroup;
    int   gmaxnr;
    int   gmaxlen;
    int64 gmaxbp;
    int  *gsize;
  } Header_Info;

static void Check_FastQ(FILE *input, char *file_name, Header_Info *info)
{ char  *line;
  int    len, qlen, lineno;
  int    nr, nread, maxlen, ngroup, gmaxnr;
  int   *gcount, maxgroup;
  int64  bp, totbp, gmaxbp;
  char  *readgroup, *p;
  int    i;

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
  while ((line = read_line(input,NULL)) != NULL)
    { lineno += 1;
      if (line[0] != '@')
        { fprintf(stdout,"%s: %s.%d: Entry header does not start with an @-sign\n",
                         Prog_Name,file_name,lineno);
          exit (1);
        }

      if (GROUP)
        { p = line+1;
          for (i = 0; i < 4; i++)
            { p = index(p+1,':');
              if (p == NULL)
                { fprintf(stdout,"%s: %s.%d: Entry header does not appear to be an Illumina header\n",
                                 Prog_Name,file_name,lineno);
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
                      gcount   = (int *) Realloc(gcount,sizeof(int)*maxgroup,
                                                 "Allocating group read counts");
                    }
                  gcount[ngroup++] = nr;
                }
    
              readgroup = strdup(line);
              bp = 0;
              nr = 0;
            } 
        }

      line = read_line(input,&len);
      lineno += 1;
      if (line == NULL)
        { fprintf(stdout,"%s: %s.%d: Incomplete fastq entry\n",Prog_Name,file_name,lineno);
          exit (1);
        }

      nr  += 1;
      bp  += len;
      if (len > maxlen)
        maxlen = len;

      line = read_line(input,NULL);
      lineno += 1;
      if (line == NULL)
        { fprintf(stdout,"%s: %s.%d Incomplete fastq entry\n",Prog_Name,file_name,lineno);
          exit (1);
        }
      if (line[0] != '+')
        { fprintf(stdout,"%s: %s.%d Divider line does not start with a +-sign\n",
                         Prog_Name,file_name,lineno);
          exit (1);
        }
      line = read_line(input,&qlen);
      lineno += 1;
      if (line == NULL)
        { fprintf(stdout,"%s: %s.%d Incomplete fastq entry\n",Prog_Name,file_name,lineno);
          exit (1);
        }
      if (len != qlen)
        { fprintf(stdout,"%s: %s.%d: QV line does not have the same length as sequence line\n",
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
          gcount   = (int *) Realloc(gcount,sizeof(int)*maxgroup,
                                     "Allocating group read counts");
        }
      gcount[ngroup++] = nr;
    }

  info->nread  = nread;
  info->maxlen = maxlen;
  info->totbp  = totbp;
  info->ngroup = ngroup;
  info->gmaxnr = gmaxnr;
  info->gmaxbp = gmaxbp;
  info->gsize  = gcount;

  free(readgroup);
}

  //  In a second pass, output .irp format for each pair

static void Output_Pairs(FILE *forward, FILE *reverse, int *gcount)
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
              printf("g %d %ld %s\n",gcount[gc++],strlen(readgroup+1)+1,readgroup+1);
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
{ FILE       *input1, *input2;
  char       *fname1, *fname2;
  Header_Info stats1, stats2;

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

    fname1 = Root(argv[1],".fastq");
    pwd1   = PathTo(argv[1]);
    input1 = fopen(Catenate(pwd1,"/",fname1,".fastq"),"r");
    if (input1 == NULL)
      { fprintf(stdout,"%s: Cannot open %s either as a .fastq or .fasta\n",Prog_Name,argv[1]);
        exit (1);
      }

    fname2 = Root(argv[2],".fastq");
    pwd2   = PathTo(argv[2]);
    input2 = fopen(Catenate(pwd2,"/",fname2,".fastq"),"r");
    if (input2 == NULL)
      { fprintf(stdout,"%s: Cannot open %s either as a .fastq or .fasta\n",Prog_Name,argv[2]);
        exit (1);
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
    { fprintf(stdout,"%s: Cannot open %s either as a .fastq or .fasta\n",Prog_Name,argv[2]);
      exit (1);
    }

  //  Output file type and provenance

  { int    len;
    char   date[26];
    time_t seconds;

    printf("1 3 irp 1 0\n");
    printf("# ! 1\n");
    printf("! 7 VGPpair 1 0");
    if (VERBOSE+SEQ_ONLY+GROUP)
      { len = strlen(argv[1]) + strlen(argv[2]) + VERBOSE + SEQ_ONLY + GROUP + 3;
        printf(" %d -",len);
        if (VERBOSE)
          printf("v");
        if (SEQ_ONLY)
          printf("s");
        if (GROUP)
          printf("g");
      }
    else
      { len = strlen(argv[1]) + strlen(argv[2]) + 1;
        printf(" %d",len);
      }
    printf(" %s %s",argv[1],argv[2]);
    seconds = time(NULL);
    ctime_r(&seconds,date);
    date[24] = '\0';
    printf(" 24 %s\n",date);
  }

  //  Output .irp header

  printf("# P %d\n",stats1.nread);
  printf("# S %d\n",2*stats1.nread);
  printf("+ S %lld\n",stats1.totbp + stats2.totbp);
  if (stats1.maxlen > stats2.maxlen)
    printf("@ S %d\n",stats1.maxlen);
  else
    printf("@ S %d\n",stats2.maxlen);

  if ( ! SEQ_ONLY)
    { printf("# Q %d\n",2*stats1.nread);
      printf("+ Q %lld\n",stats1.totbp + stats2.totbp);
      if (stats1.maxlen > stats2.maxlen)
        printf("@ Q %d\n",stats1.maxlen);
      else
        printf("@ Q %d\n",stats2.maxlen);
    }

  if (GROUP)
    { printf("# g %d\n",stats1.ngroup);
      printf("%% g # P %d\n",stats1.gmaxnr);
      printf("%% g # S %d\n",2*stats1.gmaxnr);
      printf("%% g + S %lld\n",stats1.gmaxbp + stats2.gmaxbp);
      if ( ! SEQ_ONLY)
        { printf("%% g # Q %d\n",2*stats1.gmaxnr);
          printf("%% g + Q %lld\n",stats1.gmaxbp + stats2.gmaxbp);
        }
    }

  //  Scan 2: output the data in .irp format

  if (VERBOSE)
    { fprintf(stderr,"  Outputting read pairs in .irp format to stdout\n");
      fflush(stderr);
    }

  rewind(input1);
  rewind(input2);
  Output_Pairs(input1,input2,stats1.gsize);

  //  Tidy up just for good form

  fclose(input2);
  fclose(input1);
  free(fname2);
  free(fname1);
  exit (0);
}
