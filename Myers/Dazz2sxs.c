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
#include <sys/types.h>
#include <sys/stat.h>
#include <zlib.h>
#include <time.h>
#include <limits.h>

#include "gene_core.h"

static char *Usage =
    " [-vidtg] <src1:.pbr[.gz]> [<src2:.pbr[.gz]] <align:las> ...";


/****************************************************************************************\
*                                                                                        *
*  Daligner codes needed to make this module stand alone                                 *
*                                                                                        *
\****************************************************************************************/

#define TRACE_XOVR  125
#define BLOCK_SYMBOL '@'

#define COMP_FLAG  0x1
#define COMP(x)   ((x) & COMP_FLAG)

typedef struct
  { void     *trace;
    int       tlen;
    int       diffs;
    int       abpos, bbpos;
    int       aepos, bepos;
  } Path;

typedef struct {
  Path    path;         /* Path: begin- and end-point of alignment + diffs    */
  uint32  flags;        /* Pipeline status and complementation flags          */
  int     aread;        /* Id # of A sequence                                 */
  int     bread;        /* Id # of B sequence                                 */
} Overlap;

static int64 PtrSize   = sizeof(void *);
static int64 OvlIOSize = sizeof(Overlap) - sizeof(void *);

static int Read_Overlap(FILE *input, Overlap *ovl)
{ if (fread( ((char *) ovl) + PtrSize, OvlIOSize, 1, input) != 1)
    return (1);
  return (0);
}

static int Read_Trace(FILE *input, Overlap *ovl, int tbytes)
{ if (tbytes > 0 && ovl->path.tlen > 0)
    { if (fread(ovl->path.trace, tbytes*ovl->path.tlen, 1, input) != 1)
        return (1);
    }
  return (0);
}

void Decompress_TraceTo16(Overlap *ovl)
{ uint16 *t16 = (uint16 *) ovl->path.trace;
  uint8  *t8  = (uint8  *) ovl->path.trace;
  int     j;

  for (j = ovl->path.tlen-1; j >= 0; j--)
    t16[j] = t8[j];
}


typedef struct
  { int first, last, next;
    char *root, *pwd, *ppnt;
    char *slice;
  } Block_Looper;

static Block_Looper *Parse_Block_LAS_Arg(char *arg)
{ Block_Looper *parse;
  char *pwd, *root;
  char *ppnt, *cpnt;
  int   first, last;

  parse = (Block_Looper *) Malloc(sizeof(Block_Looper),"Allocating parse node");
  pwd   = PathTo(arg);
  root  = Root(arg,".las");
  if (parse == NULL || pwd == NULL || root == NULL)
    exit (1);

  ppnt = index(root,BLOCK_SYMBOL);
  if (ppnt == NULL)
    first = last = -1;
  else
    { if (index(ppnt+1,BLOCK_SYMBOL) != NULL)
        { fprintf(stderr,"%s: Two or more occurences of %c-sign in source name '%s'\n",
                         Prog_Name,BLOCK_SYMBOL,root);
          exit (1);
        }
      *ppnt++ = '\0';
      first = strtol(ppnt,&cpnt,10);
      if (cpnt == ppnt)
        { first = 1;
          last  = INT_MAX;
        }
      else
        { if (first < 1)
            { fprintf(stderr,
                      "%s: Integer following %c-sigan is less than 1 in source name '%s'\n",
                      Prog_Name,BLOCK_SYMBOL,root);
              exit (1);
            }
          if (*cpnt == '-')
            { ppnt = cpnt+1;
              last = strtol(ppnt,&cpnt,10);
              if (cpnt == ppnt)
                { fprintf(stderr,"%s: Second integer must follow - in source name '%s'\n",
                                 Prog_Name,root);
                  exit (1);
                }
              if (last < first)
                { fprintf(stderr,
                          "%s: 2nd integer is less than 1st integer in source name '%s'\n",
                          Prog_Name,root);
                  exit (1);
                }
              ppnt = cpnt;
            }
          else
            { last = INT_MAX;
              ppnt = cpnt;
            }
        }
    }

  parse->pwd   = pwd;
  parse->root  = root;
  parse->ppnt  = ppnt;
  parse->first = first;
  parse->last  = last;
  parse->next  = first-1;
  parse->slice = NULL;

  return (parse);
}

FILE *Next_Block_Arg(Block_Looper *parse)
{ char *disp;
  FILE *input;

  parse->next += 1;
  if (parse->next > parse->last)
    return (NULL);

  if (parse->next < 0)
    disp  = parse->root;
  else
    disp = Numbered_Suffix(parse->root,parse->next,parse->ppnt);

  if ((input = fopen(Catenate(parse->pwd,"/",disp,".las"),"r")) == NULL)
    { if (parse->last != INT_MAX)
        { fprintf(stderr,"%s: %s.las is not present\n",Prog_Name,disp);
          exit (1);
        }
      return (NULL);
    }
  return (input);
}

static char *Block_Arg_Root(Block_Looper *parse)
{ char *name;

  if (parse->next < 0)
    name = parse->root;
  else
    name = Numbered_Suffix(parse->root,parse->next,parse->ppnt);
  return (Strdup(name,"Allocating block root"));
}

static void Free_Block_Arg(Block_Looper *parse)
{ free(parse->root);
  free(parse->pwd);
  free(parse->slice);
  free(parse);
}


/****************************************************************************************\
*                                                                                        *
*  Main code for the program proper                                                      *
*                                                                                        *
\****************************************************************************************/


FILE *fzopen(const char *path, const char *mode)
{ gzFile zfp;

  zfp = gzopen(path,mode);
  if (zfp == NULL)
    return (fopen(path,mode));

  return (funopen(zfp,
                  (int(*)(void*,char*,int))gzread,
                  (int(*)(void*,const char*,int))gzwrite,
                  (fpos_t(*)(void*,fpos_t,int))gzseek,
                  (int(*)(void*))gzclose) );
}

static int Header[128] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  };

int *Fetch_Length_Vector(FILE *input, int *nread, char *fname)
{ int *rlen, nr;
  char code, which;

  fscanf(input,"1 %d",&nr);
  { char type[nr+1];

    fscanf(input," %s",type);
    if (strcmp(type,"seq") != 0)
      { fprintf(stderr,"%s: File %s is not a .seq 1-code file\n",Prog_Name,fname);
        exit (1);
      }
  }
  fscanf(input,"%*[^\n]\n");

  while (fscanf(input,"%c",&code) == 1)
    if (Header[(int) code])
      { if (code == '#')
          { fscanf(input," %c",&which);
            if (which == 'S')
              { fscanf(input," %d",nread);
                break;
              }
          }
        fscanf(input,"%*[^\n]\n");
      }
    else
      ungetc(code,input);

  rlen = (int *) Malloc(sizeof(int)*(*nread),"Allocating read length vector");
  if (rlen == NULL)
    exit (1);

  nr = 0;
  while (fscanf(input,"%c",&code) == 1)
    { if (code == 'S')
        fscanf(input," %d",rlen+nr++);
      fscanf(input,"%*[^\n]\n");
    }

  *nread = nr;
  return (rlen);
}

int main(int argc, char *argv[])
{ int       nread1, nread2;
  int      *rlen1,  *rlen2;
  char     *fname1, *fname2;
  Overlap   _ovl, *ovl = &_ovl;

  FILE   *input;
  int     tspace, tbytes, small;
  int     trmax;
  int     ngroup, gmax;
  int    *gcount;

  int     ISTWO, VERBOSE;
  int     DOGROUP, DOCOORD, DODIFF, DOTRACE; 

  //  Process options

  { int    i, j, k;
    int    flags[128];

    ARG_INIT("LA2sxs")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vidtg")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    DOGROUP = flags['g'];
    DOCOORD = flags['i'];
    DODIFF  = flags['d'];
    DOTRACE = flags['t'];

    if (argc <= 2)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"      A #a #b  - (#a,#b) have an LA between them\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"   -i: I #ab #ae #alen #bb #be #blen - #a[#ab,#ae] aligns with");
        fprintf(stderr," #b[#bb,#be]\n");
        fprintf(stderr,"   -d: D #                           - there are # differences");
        fprintf(stderr," in the LA\n");
        fprintf(stderr,"   -t: T #n #y^#n                    - there are #n trace point");
        fprintf(stderr," intervals for the LA\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"   -g  Output la's in read pile groups\n");

        exit (1);
      }
  }

  //  Get read lengths and # of reads from sequence files

  { char *pwd, *root;
    FILE *input1, *input2;

    pwd    = PathTo(argv[1]);
    root   = Root(argv[1],".pbr.gz");
    fname1 = Catenate(pwd,"/",root,".pbr.gz");
    if ((input1 = fzopen(fname1,"r")) == NULL)
      { free(root);
        root   = Root(argv[1],".pbr");
        fname1 = Catenate(pwd,"/",root,".pbr");
        if ((input1 = fzopen(fname1,"r")) == NULL)
          { fprintf(stderr,"%s: Cannot open 1-code file %s\n",Prog_Name,argv[1]);
            exit (1);
          }
      }
    free(pwd);
    free(root);

    if (VERBOSE)
      { fprintf(stderr,"  Scanning first sequence file\n");
        fflush(stderr);
      }

    fname1 = Strdup(fname1,"Allocating first sequence file");
    rlen1  = Fetch_Length_Vector(input1,&nread1,argv[1]);
    fclose(input1);

    pwd    = PathTo(argv[2]);
    root   = Root(argv[2],".pbr.gz");
    fname2 = Catenate(pwd,"/",root,".pbr.gz");
    if ((input2 = fzopen(fname2,"r")) == NULL)
      { free(root);
        root   = Root(argv[2],".pbr");
        fname2 = Catenate(pwd,"/",root,".pbr");
        if ((input2 = fzopen(fname2,"r")) == NULL)
          ISTWO = 0;
        else
          ISTWO = 1;
      }
    else
      ISTWO = 1;
    free(pwd);
    free(root);

    if (ISTWO)
      { if (VERBOSE)
          { fprintf(stderr,"  Scanning second sequence file\n");
            fflush(stderr);
          }

        fname2 = Strdup(fname2,"Allocating second sequence file");
        rlen2  = Fetch_Length_Vector(input1,&nread2,argv[2]);
        fclose(input2);
      }
    else
      { rlen2  = rlen1;
        nread2 = nread1;
	fname2 = fname1;
      }
  }

  { int   c, j, al, ar, tlen;
    int64 novls, odeg, omax, sdeg, smax, tmax, ttot;

    //  For each record do

    novls = omax = smax = ttot = tmax = 0;
    sdeg  = odeg = 0;

    ngroup = 0;
    gmax   = 0;
    gcount = NULL;

    tspace = -1;
    for (c = 2+ISTWO; c < argc; c++)
      { Block_Looper *parse;

        parse = Parse_Block_LAS_Arg(argv[c]);

        while ((input = Next_Block_Arg(parse)) != NULL)
          { int64 no;
            int   mspace;

            //  Initiate file reading
  
            if (VERBOSE)
              { fprintf(stderr,"  Scanning %s.las for header info\n",Block_Arg_Root(parse));
                fflush(stderr);
              }

            if (fread(&no,sizeof(int64),1,input) != 1)
              { fprintf(stderr,"%s: System error, read failed!\n",Prog_Name);
                exit (1);
              }
            if (fread(&mspace,sizeof(int),1,input) != 1)
              { fprintf(stderr,"%s: System error, read failed!\n",Prog_Name);
                exit (1);
              }

            if (tspace < 0)
              { tspace = mspace; 
                if (tspace <= TRACE_XOVR && tspace != 0)
                  { small  = 1;
                    tbytes = sizeof(uint8);
                  }
                else
                  { small  = 0;
                    tbytes = sizeof(uint16);
                  }
              }
            else if (mspace != tspace)
              { fprintf(stderr,"%s: Input .las files have different trace spacing!\n",Prog_Name);
                exit (1);
              }

            al = -1;
            for (j = 0; j < no; j++)
              { Read_Overlap(input,ovl);
                tlen = ovl->path.tlen;
                fseeko(input,tlen*tbytes,SEEK_CUR);

                ar = ovl->aread;
                if (ar != al)
                  { if (sdeg > smax)
                      smax = sdeg;
                    if (odeg > omax)
                      omax = odeg;
                    novls += odeg;
                    ttot  += sdeg;
                    if (al >= 0)
                      { if (ngroup >= gmax)
                          { gmax = 1.2*ngroup + 1000;
                            gcount = Realloc(gcount,sizeof(int)*gmax,
                                            "Reallocating group count vector");
                          }
                        gcount[ngroup++] = odeg;
                      }
                    sdeg = odeg = 0;
                    al = ar;
                  }

                odeg  += 1;
                sdeg  += tlen;
                if (tlen > tmax)
                  tmax = tlen;
              }

            if (sdeg > smax)
              smax = sdeg;
            if (odeg > omax)
              omax = odeg;
            novls += odeg;
            ttot  += sdeg;
            if (ngroup > gmax)
              { gmax = 1.2*ngroup + 1000;
                gcount = Realloc(gcount,sizeof(int)*gmax,"Reallocating group count vector");
              }
            gcount[ngroup++] = odeg;
            ngroup += 1;

            fclose(input);
          }
        Free_Block_Arg(parse);
      }

    trmax = tmax;

    { int    i, clen, optl;
      char   date[26];
      time_t seconds;

      optl = DOGROUP+DOCOORD+DODIFF+DOTRACE+VERBOSE;
      if (optl > 0)
        clen = optl+1;
      else
        clen = -1;
      for (i = 1; i < argc; i++)
        clen += strlen(argv[i])+1;

      printf("1 3 sxs 1 0\n");
      printf("# ! 1\n");
      if (DOGROUP)
        printf("# g %d\n",ngroup);
      printf("# A %lld\n",novls);
      if (DOCOORD)
        printf("# I %lld\n",novls);
      if (DODIFF)
        printf("# D %lld\n",novls);
      if (DOTRACE)
        printf("# T %lld\n# Z %lld\n",novls,novls);

      printf("+ ! %d\n",clen+33);
      if (DOGROUP)
        printf("+ g %d\n",ngroup);
      if (DOTRACE)
        printf("+ T %lld\n+ Z %lld\n",ttot,ttot);

      if (clen > 24)
        printf("@ ! %d\n",clen);
      else
        printf("@ ! 24\n");
      if (DOGROUP)
        printf("@ g 1\n");
      if (DOTRACE)
        printf("@ T %lld\n@ Z %lld\n",tmax,tmax);

      if (DOGROUP)
        { printf("%% g # A %lld\n",omax);
          if (DOCOORD)
            printf("%% g # I %lld\n",omax);
          if (DODIFF)
            printf("%% g # D %lld\n",omax);
          if (DOTRACE)
            { printf("%% g # T %lld\n%% g # Z %lld\n",omax,omax);
              printf("%% g + T %lld\n%% g + Z %lld\n",smax,smax);
            }
        }

      printf("\n! 6 LA2sxs 3 1.0 %d",clen);
      if (optl > 0)
        { printf(" -");
          if (VERBOSE)
            printf("v");
          if (DOCOORD)
            printf("i");
          if (DODIFF)
            printf("d");
          if (DOTRACE)
            printf("t");
          if (DOGROUP)
            printf("g");
        }
      for (i = 1; i < argc; i++)
        printf(" %s",argv[i]);
      seconds = time(NULL);
      ctime_r(&seconds,date);
      date[24] = '\0';
      printf(" 24 %s\n",date);

      printf("\n> %ld %s S %d\n",strlen(fname1),fname1,nread1);
      printf("> %ld %s S %d\n\n",strlen(fname2),fname2,nread2);
    }
  }

  //  Read the files and display records
  
  { int     c, j, k;
    int     ar, al, ng;
    int     blen;
    uint16 *trace;

    trace = (uint16 *) Malloc(sizeof(uint16)*trmax,"Allocating trace vector");
    if (trace == NULL)
      exit (1);

    al = -1;
    ng = 0;
    for (c = 2+ISTWO; c < argc; c++)
      { Block_Looper *parse;

        parse = Parse_Block_LAS_Arg(argv[c]);

        while ((input = Next_Block_Arg(parse)) != NULL)
          { int64 no;
            int   mspace;

            //  Initiate file reading
  
            if (VERBOSE)
              { fprintf(stderr,"  Scanning %s.las for conversion\n",Block_Arg_Root(parse));
                fflush(stderr);
              }

            fread(&no,sizeof(int64),1,input);
            fread(&mspace,sizeof(int),1,input);

            //  For each record do

            for (j = 0; j < no; j++)

               //  Read it in

              { Read_Overlap(input,ovl);
                ovl->path.trace = (void *) trace;
                Read_Trace(input,ovl,tbytes);

                //  Determine if it should be displayed

                ar = ovl->aread;
                if (ar != al)
                  { if (DOGROUP)
                      printf("g %d 1 *\n",gcount[ng++]);
                    al = ar;
                  }

                //  Display it
            
                printf("A %d %d\n",ar+1,ovl->bread+1);

                blen = rlen2[ovl->bread];
                if (DOCOORD)
                  { if (COMP(ovl->flags))
                      printf("I %d %d %d %d %d %d\n",
                             ovl->path.abpos,ovl->path.aepos,rlen1[ar],
                             blen - ovl->path.bbpos, blen - ovl->path.bepos,blen);
                    else
                      printf("I %d %d %d %d %d %d\n",
                             ovl->path.abpos,ovl->path.aepos,rlen1[ar],
                             ovl->path.bbpos,ovl->path.bepos,rlen2[ovl->bread]);
                  }
        
                if (DODIFF)
                  printf("D %d\n",ovl->path.diffs);
        
                if (DOTRACE)
                  { int tlen  = ovl->path.tlen;
                    int at, bt;
        
                    if (small)
                      Decompress_TraceTo16(ovl);
        
                    at = (ovl->path.abpos / tspace) * tspace;
                    printf("U %d %d",(tlen>>1)+1,ovl->path.abpos);
                    for (k = 0; k < tlen-2; k += 2)
                      { at += tspace;
                        printf(" %d",at);
                      }
                    printf(" %d\n",ovl->path.aepos);
                 
                    bt = ovl->path.bbpos;
                    printf("V %d %d",(tlen>>1)+1,ovl->path.bbpos);
                    for (k = 0; k < tlen; k += 2)
                      { bt += trace[k+1];
                        printf(" %d",bt);
                      }
                    printf("\n");
        
                    printf("Z %d",tlen>>1);
                    for (k = 0; k < tlen; k += 2)
                      printf(" %3d",trace[k]);
                    printf("\n");
                  }
              }
            fclose(input);
          }
        Free_Block_Arg(parse);
      }
    free(trace);
  }

  exit (0);
}
