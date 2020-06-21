/*******************************************************************************************
 *
 *  Utility for displaying the information in the overlaps of a .las file in a very
 *    simple to parse format.
 *
 *  Author:    Gene Myers
 *  Creation:  Feb 2019
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <zlib.h>
#include <time.h>
#include <limits.h>
#include <stdbool.h>

#include "gene_core.h"
#include "../Core/ONElib.h"

#include "VGPschema.h"

#undef  DEBUG_OUT

static int     VERBOSE;
static int     NTHREADS;
static int     HAS_QVS;

static char *Usage =
    "[-v] [-T<int(4)>] <forward:seq> <reverse:seq>";


/*******************************************************************************************
 *
 *  Parallel:  Each thread processes is a contiguous stripe across the .las input files
 *               sending the compressed binary data lines to their assigned OneFile.
 *
 ********************************************************************************************/


typedef struct
  { OneFile     *vf;       //  OneFile for output
    OneFile     *v1;       //  OneFiles for input
    OneFile     *v2;
    int          beg;      //  Range of reads to process
    int          end;
  } Thread_Arg;

  //  Write alignment records in relevant partition

static void *output_thread(void *arg)
{ Thread_Arg *parm = (Thread_Arg *) arg;
  int64       beg  = parm->beg;
  int64       end  = parm->end;
  OneFile    *vf   = parm->vf;
  OneFile    *v1   = parm->v1;
  OneFile    *v2   = parm->v2;
  int         i, j, n, t1, t2;

#define TRANSFER(vi,ti,vo)				\
{ for (j = 0; j < vi->info[ti]->nField; j++)		\
    vo->field[j] = vi->field[j];			\
  oneWriteLine(vo,ti,oneLen(vi),oneString(vi));		\
}

  oneGotoObject(v1,beg);
  oneGotoObject(v2,beg);
  t1 = oneReadLine(v1);
  t2 = oneReadLine(v2);
  if (t1 != 'S' || t2 != 'S')
    { fprintf(stderr,"%s: Fatal, goto line is not 'S' (%c,%c,%lld)\n",Prog_Name,t1,t2,beg);
      exit (1);
    }
  for (i = beg; i < end; i++)
    { oneWriteLine(vf,'P',0,NULL);

      TRANSFER(v1,t1,vf)
      n  = oneLen(v1);
      t1 = oneReadLine(v1);
      while (t1 != 'S')
        { if (t1 == 'Q')
            { if (n != oneLen(v1))
                { fprintf(stderr,"%s: Q string not same length in forward file, line %lld\n",
                                 Prog_Name,v1->line);
                  exit (1);
                }
              n = 0;
            }
          else if (t1 == 0 || t1 == 'g')
            break; 
          TRANSFER(v1,t1,vf)
          t1 = oneReadLine(v1);
        }
      if (HAS_QVS && n > 0)
        { fprintf(stderr,"%s: Missing Q-line\n",Prog_Name);
          exit (1);
        }

      TRANSFER(v2,t2,vf)
      n  = oneLen(v2);
      t2 = oneReadLine(v2);
      while (t2 != 'S')
        { if (t2 == 'Q')
            { if (n != oneLen(v2))
                { fprintf(stderr,"%s: Q string not same length in reverse file, line %lld\n",
                                 Prog_Name,v2->line);
                  exit (1);
                }
              n = 0;
            }
          else if (t2 == 0)
            break; 
          if (t2 != 'g')
            TRANSFER(v2,t2,vf)
          t2 = oneReadLine(v2);
        }
      if (HAS_QVS && n > 0)
        { fprintf(stderr,"%s: Missing Q-line\n",Prog_Name);
          exit (1);
        }

      if (t1 == 'g')
        { TRANSFER(v1,t1,vf)
          t1 = oneReadLine(v1);
          if (t1 != 'S' && t1 != 0)
            { fprintf(stderr,"%s: group line does not precede sequence line",Prog_Name);
              fprintf(stderr," in forward file, line %lld\n",v1->line);
              exit (1);
            }
        }
    }

  return (NULL);
}


/****************************************************************************************
 *                                                                                     
 *  The top-level program                                                              
 *                                                                                     
 ****************************************************************************************/

int main(int argc, char *argv[])
{ char      *fname1, *fname2;
  char      *command;
  OneSchema *schema;

  //  Capture command line for provenance & load ONE schema

  { int   n, i;
    char *c;

    n = 0;
    for (i = 1; i < argc; i++)
      n += strlen(argv[i])+1;

    command = Malloc(n+1,"Allocating command string");
    if (command == NULL)
      exit (1);

    c = command;
    if (argc >= 1)
      { c += sprintf(c,"%s",argv[1]);
        for (i = 2; i < argc; i++)
          c += sprintf(c," %s",argv[i]);
      }
    *c = '\0';

    schema = oneSchemaCreateFromText(vgpSchemaText);
  }

  //  Process options

  { int    i, j, k;
    int    flags[128];
    char  *eptr;

    ARG_INIT("VGPpair")

    NTHREADS = 4;

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
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];

    if (argc != 3)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: verbose mode, output progress as proceed\n");
        fprintf(stderr,"      -T: Number of threads to use\n");
        exit (1);
      }
  }

  //  Get proper names of both inputs

  { char *pwd, *root;
    int   i;
    FILE *input;
    char *suffix[1] = { ".seq" };

    pwd    = PathTo(argv[1]);
    OPEN(argv[1],pwd,root,input,suffix,1)
    if (input == NULL)
      { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[1]);
        exit (1);
      }
    fname1 = Strdup(Catenate(pwd,"/",root,suffix[i]),"Allocating first sequence file");
    free(root);
    free(pwd);
    fclose(input);

    pwd    = PathTo(argv[2]);
    OPEN(argv[1],pwd,root,input,suffix,1)
    if (input == NULL)
      { fprintf(stderr,"%s: Cannot open %s\n",Prog_Name,argv[2]);
        exit (1);
      }
    fname2 = Strdup(Catenate(pwd,"/",root,suffix[i]),"Allocating second sequence file");
    free(root);
    free(pwd);
    fclose(input);
  }

  { Thread_Arg  parm[NTHREADS];

#ifndef DEBUG_OUT
    pthread_t   threads[NTHREADS];
#endif

    //  Produce output in parallel threads based on partition

    { OneFile *vf, *v1, *v2;
      int      i, j;
      int64    nreads;

      if (VERBOSE)
        { fprintf(stderr,"  Splitting .seq files into %d segments for pair merging\n",NTHREADS);
          fflush(stderr);
        }

      v1 = oneFileOpenRead(fname1,schema,"seq",NTHREADS);
      v2 = oneFileOpenRead(fname2,schema,"seq",NTHREADS);

      nreads = v1->info['S']->given.count;
      if (nreads != v2->info['S']->given.count)
        { fprintf(stderr,"%s: The files do nat have the same number of sequences!\n",Prog_Name);
          exit (1);
        }

      HAS_QVS = (v1->info['Q']->given.count > 0);

      vf = oneFileOpenWriteNew("-",schema,"irp",true,NTHREADS);

      oneInheritProvenance(vf,v1);
      oneInheritProvenance(vf,v2);
      oneAddProvenance(vf,Prog_Name,"1.0",command,NULL);

      oneWriteHeader(vf);

      if (oneReadLine(v1) == 'g')
        TRANSFER(v1,'g',vf);

#ifdef DEBUG_OUT
      fprintf(stderr,"Opened\n");
      fflush(stdout);
#endif

      if (VERBOSE)
        { fprintf(stderr,"  Producing .irp segements in parallel\n");
          fflush(stderr);
        }

      for (i = 0; i < NTHREADS; i++)
        { parm[i].beg = (nreads * i) / NTHREADS;
          parm[i].end = (nreads * (i+1)) / NTHREADS;
          parm[i].v1  = v1+i;
          parm[i].v2  = v2+i;
          parm[i].vf  = vf+i;
        }
      if (parm[0].end == 0 && nreads > 0)
        parm[0].end = parm[1].beg = 1;

      //  Generate the data lines in parallel threads

#ifdef DEBUG_OUT
      for (i = 0; i < NTHREADS; i++)
        output_thread(parm+i);
#else
      for (i = 0; i < NTHREADS; i++)
        pthread_create(threads+i,NULL,output_thread,parm+i);

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);
#endif

      if (VERBOSE)
        { fprintf(stderr,"  Cat'ing .irp segments\n");
          fflush(stderr);
        }

      oneFileClose(vf);
      oneFileClose(v1);
      oneFileClose(v2);
    }

    //  Free everything as a matter of good form

    free(fname1);
    free(fname2);
    free(command);
  }

  oneSchemaDestroy(schema);

  if (VERBOSE)
    { fprintf(stderr,"  Done\n");
      fflush(stderr);
    }

  exit (0);
}
