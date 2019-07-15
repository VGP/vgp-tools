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

static char *Usage = "[-vsg] <forward:fast[aq][.gz]> [<reverse:fast[aq][.gz]>";

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

static void Check_FastQ(FILE *input, int fastq, char *file_name, Header_Info *info)
{ char   *line, hsymbol;
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

  if (fastq)
    hsymbol = '@';
  else
    hsymbol = '>';

  nr = 0;
  bp = 0;
  line = read_line(input,NULL);
  while (line != NULL)
    { lineno += 1;
      if (line[0] != hsymbol)
        { fprintf(stdout,"%s: %s.%lld: Entry header does not start with an %c-sign\n",
                         Prog_Name,file_name,lineno,hsymbol);
          exit (1);
        }

      if (GROUP)
        { p = line+1;
          for (i = 0; i < 3; i++)
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

      if (fastq)
        { line = read_line(input,&len);
          lineno += 1;
          if (line == NULL)
            { fprintf(stdout,"%s: %s.%lld: Incomplete fastq entry\n",Prog_Name,file_name,lineno);
              exit (1);
            }
    
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
            { fprintf(stdout,"%s: %s.%lld: QV line doesn't have the same length as sequence line\n",
                             Prog_Name,file_name,lineno);
              exit (1);
            }

          line = read_line(input,NULL);
        }
      else
        { len = 0;
          while (1)
            { line = read_line(input,&qlen);
              if (line == NULL || line[0] == '>')
                break;
              len += qlen;
            }
          if (len == 0)
            { if (line == NULL)
                fprintf(stdout,"%s: Sequence missing after last fasta header\n",Prog_Name);
              else
                fprintf(stdout,"%s: Sequence missing between two fasta headers\n",Prog_Name);
              exit (1);
            }
        }

      nr  += 1;
      bp  += len;
      if (len > maxlen)
        maxlen = len;
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

static void Output_Sequence(int ispair, FILE *forward, FILE *reverse, Header_Info *info)
{ char  *line;
  int    len;
  int    DO_QV;
  char  *p, *readgroup;
  int    i, gc;
  int64 *gcount;
  char  *seq;
  int    rmax;

  DO_QV  = ! SEQ_ONLY;
  gcount = info->gsize;
  rmax   = info->maxlen + 2;

  seq = (char *) Malloc(rmax+1,"Allocating sequence buffer\n");
  if (seq == NULL)
    exit (1);

  gc = 0;
  readgroup = NULL;
  while ((line = read_line(forward,NULL)) != NULL)
    { if (GROUP)
        { p = line+1;
          for (i = 0; i < 3; i++)
            p = index(p+1,':');
          *p = '\0';

          if (readgroup == NULL || strcmp(readgroup,line) != 0)
            { free(readgroup);
              readgroup = strdup(line);
              printf("g %lld %ld %s\n",gcount[gc++],strlen(readgroup+1),readgroup+1);
            } 
        }

      if (ispair)
        printf("P\n");

      fgets(seq,rmax,forward);
      len = strlen(seq);
      seq[--len] = '\0';
      printf("S %d %s\n",len,seq);

      line = read_line(forward,NULL);
      line = read_line(forward,NULL);
      if (DO_QV)
        printf("Q %d %s\n",len,line);

      if (ispair)
        { line = read_line(reverse,NULL);

          fgets(seq,rmax,reverse);
          len = strlen(seq);
          seq[--len] = '\0';
          printf("S %d %s\n",len,seq);

          line = read_line(reverse,NULL);
          line = read_line(reverse,NULL);
          if (DO_QV)
            printf("Q %d %s\n",len,line);
        }
    }
}

int main(int argc, char *argv[])
{ FILE       *input1, *input2;
  char       *fname1, *fname2;
  Header_Info stats1, stats2;
  int         ispair, isfastq;

  //  Parse command line options

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("VGPseq")

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

    if (argc != 2 && argc != 3)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, show progress as proceed.\n");
        fprintf(stderr,"      -s: Output sequences only, skip QVs\n");
        fprintf(stderr,"      -g: Output group lanes based on Illumina lanes\n");
        exit (1);
      }

    ispair = (argc == 3);
  }

  //  Open input files

  { char *suffix[5] = { ".fastq.gz", ".fasta.gz", ".fastq", ".fasta", ".gz" };
    char *pwd1, *pwd2;
    int   i, this;

    pwd1 = PathTo(argv[1]);
    OPEN(argv[1],pwd1,fname1,input1,suffix,5)
    if (input1 == NULL)
      { fprintf(stderr,"%s: Cannot open %s as an .fast[aq] file\n",Prog_Name,argv[1]);
        exit (1);
      }
    if (i == 4)
      isfastq = (strcmp(".fastq",fname1+(strlen(fname1)-6)) == 0);
    else
      isfastq = (i%2 == 0);

    if (ispair)
      { pwd2 = PathTo(argv[2]);
        OPEN(argv[2],pwd2,fname2,input2,suffix,5)
        if (input2 == NULL)
          { fprintf(stderr,"%s: Cannot open %s as an .fast[aq] file\n",Prog_Name,argv[2]);
            exit (1);
          }
        if (i == 4)
          this = (strcmp(".fastq",fname2+(strlen(fname2)-6)) == 0);
        else
          this = (i%2 == 0);
        if (this != isfastq)
          { fprintf(stderr,"%s: Pair is not both .fastq or both .fasta\n",Prog_Name);
            exit (1);
          }

        if (strcmp(fname1,fname2) == 0)
          { char *aname1, *aname2;
    
/*
            if (strcmp(pwd1,pwd2) == 0)
              { fprintf(stderr,"%s: Forward and reverse files are the same file?\n",Prog_Name);
                exit (1);
              }
*/
            aname1 = Strdup(Catenate(pwd1,"/",fname1,""),NULL);
            aname2 = Strdup(Catenate(pwd2,"/",fname2,""),NULL);
            free(fname2);
            free(fname1);
            fname1 = aname1;
            fname2 = aname2;
          }
        free(pwd2);
      }
    free(pwd1);
  }

  //  Scan 1: check syntax and count

  if (VERBOSE)
    { if (argc == 2)
        fprintf(stderr,"  Checking syntax of file %s\n",argv[1]);
      else
        fprintf(stderr,"  Checking syntax of forward file %s\n",argv[1]);
      fflush(stderr);
    }

  if ( ! isfastq)
    SEQ_ONLY = 1;

  Check_FastQ(input1,isfastq,fname1,&stats1);

  if (ispair)
    { if (VERBOSE)
        { fprintf(stderr,"  Checking syntax of reverse file %s\n",argv[2]);
          fflush(stderr);
        }

      Check_FastQ(input2,isfastq,fname2,&stats2);

      if (stats1.nread != stats2.nread)
        { fprintf(stdout,"%s: Number of reads in two files not equal: %lld vs %lld\n",
                         Prog_Name,stats1.nread,stats2.nread);
          exit (1);
        }

      stats1.totbp += stats2.totbp;
      if (stats1.maxlen < stats2.maxlen)
        stats1.maxlen = stats2.maxlen;
      if (stats1.gmaxnr < stats2.gmaxnr)
        stats1.gmaxnr = stats2.gmaxnr;
      stats1.gmaxbp += stats2.gmaxbp;
    }

  //  Output .irp header

  { int    i, clen, optl;
    char   date[26];
    time_t seconds;

    optl = VERBOSE + SEQ_ONLY + GROUP;
    if (optl == 0)
      clen = -1;
    else
      clen = optl+1;
    for (i = 1; i < argc; i++)
      clen += strlen(argv[i])+1;
    seconds = time(NULL);
    ctime_r(&seconds,date);
    date[24] = '\0';

    printf("1 3 seq 1 0\n");
    if (ispair)
      printf("2 3 irp\n");
    if (GROUP)
      printf("# g %d\n",stats1.ngroup);
    if (ispair)
      printf("# P %lld\n",stats1.nread);
    printf("# S %lld\n",(1+ispair)*stats1.nread);
    if ( ! SEQ_ONLY)
      printf("# Q %lld\n",(1+ispair)*stats1.nread);

    if (GROUP)
      printf("+ g %lld\n",stats1.gtotc);
    printf("+ S %lld\n",stats1.totbp);
    if ( ! SEQ_ONLY)
      printf("+ Q %lld\n",stats1.totbp);

    if (GROUP)
      printf("@ g %lld\n",stats1.gmaxc);
    printf("@ S %d\n",stats1.maxlen);
    if ( ! SEQ_ONLY)
      printf("@ Q %d\n",stats1.maxlen);

    if (GROUP)
      { if (ispair)
          printf("%% g # P %lld\n",stats1.gmaxnr);
        printf("%% g # S %lld\n",(1+ispair)*stats1.gmaxnr);
        if ( ! SEQ_ONLY)
          printf("%% g # Q %lld\n",(1+ispair)*stats1.gmaxnr);
        printf("%% g + S %lld\n",stats1.gmaxbp);
        if ( ! SEQ_ONLY)
          printf("%% g + Q %lld\n",stats1.gmaxbp);
      }

    printf("! 6 VGPseq 3 1.0 %d",clen);
    if (optl)
      { printf(" -");
        if (VERBOSE)
          printf("v");
        if (SEQ_ONLY)
          printf("s");
        if (GROUP)
          printf("g");
      }
    for (i = 1; i < argc; i++)
      printf(" %s",argv[i]);
    printf(" 24 %s\n",date);
  }
  
  //  Scan 2: output the data in .irp format

  if (VERBOSE)
    { fprintf(stderr,"  Outputting read pairs in .irp format to stdout\n");
      fflush(stderr);
    }

  rewind(input1);
  if (ispair)
    rewind(input2);
  Output_Sequence(ispair,input1,input2,&stats1);

  //  Tidy up just for good form

  fclose(input2);
  fclose(input1);
  free(fname2);
  free(fname1);
  exit (0);
}
