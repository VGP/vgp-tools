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

#undef  DEBUG_FIND
#undef  DEBUG_OUT

#define IO_BLOCK  100000

static int     VERBOSE;
static int     NTHREADS;
static int     DOGROUP;
static int     DOCOORD;
static int     DODIFF;
static int     DOTRACE; 

static int     NREAD1;
static int     NREAD2;
static int     RMAX1;
static int     RMAX2;
static int    *RLEN1;
static int    *RLEN2;

static int     TSPACE;
static int     TBYTES;

static char *Usage =
    " [-vidtg] [-T<int(4)>] <src1:.pbr> [<src2:.pbr>] <align:las> ...";

typedef struct
  { char      *fname;   //  Full path name of file
    int64      fsize;   //  size of file
  } File_Object;


/****************************************************************************************
 *                                                                                      
 *  Daligner codes needed to make this module stand alone                               
 *                                                                                      
 ****************************************************************************************/

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
static int64 OvlSize   = sizeof(Overlap);
static int64 OvlIOSize = sizeof(Overlap) - sizeof(void *);

static int64 Read_Overlap(FILE *input, Overlap *ovl)
{ if (fread( ((char *) ovl) + PtrSize, OvlIOSize, 1, input) != 1)
    return (-1);
  return (OvlIOSize);
}

static int64 Read_Trace(FILE *input, Overlap *ovl, int tbytes)
{ if (tbytes > 0 && ovl->path.tlen > 0)
    { if (fread(ovl->path.trace, tbytes*ovl->path.tlen, 1, input) != 1)
        return (-1);
    }
  return (tbytes*ovl->path.tlen);
}

void Decompress_Trace8(Overlap *ovl, int64 *bdels, int64 *diffs)
{ uint8  *t8  = (uint8 *) ovl->path.trace;
  int     i, j;

  for (i = j = 0; j < ovl->path.tlen; j += 2, i += 1)
    { diffs[i] = t8[j];
      bdels[i] = t8[j+1];
    }
}

void Decompress_Trace16(Overlap *ovl, int64 *bdels, int64 *diffs)
{ uint16 *t16 = (uint16 *) ovl->path.trace;
  int     i, j;

  for (i = j = 0; j < ovl->path.tlen; j += 2, i += 1)
    { diffs[i] = t16[j];
      bdels[i] = t16[j+1];
    }
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

static char *Block_Arg_Name(Block_Looper *parse)
{ char *name, *disp;

  if (parse->next < 0)
    disp = parse->root;
  else
    disp = Numbered_Suffix(parse->root,parse->next,parse->ppnt);
  name = Catenate(parse->pwd,"/",disp,".las");
  return (Strdup(name,"Allocating block root"));
}

static void Free_Block_Arg(Block_Looper *parse)
{ free(parse->root);
  free(parse->pwd);
  free(parse->slice);
  free(parse);
}


/****************************************************************************************
 *                                                                                        
 *  Code for creating a read length vector                                                
 *                                                                                        
 *****************************************************************************************/

typedef struct
  { OneFile     *vf;       //  OneFile for input
    int          beg;      //  Range of reads to process
    int          end;
    int         *len;
  } Vector_Arg;

  //  Display each read (and/or QV streams) in the active DB according to the
  //    range pairs in pts[0..reps) and according to the display options.

static void *fetch_thread(void *arg)
{ Vector_Arg *parm = (Vector_Arg *) arg;
  int         end  = parm->end;
  OneFile    *vf   = parm->vf;
  int        *len  = parm->len;
  int         i;

  for (i = parm->beg; i < end; i++) 
    { oneGotoObject(vf,i);
      oneReadLine(vf);
      len[i] = oneLen(vf);
    }

  return (NULL);
}
       
int *Fetch_Length_Vector(OneFile *vf, int *nread, int *rmax)
{ int       *rlen, i;
  int64      nreads;
  Vector_Arg parm[NTHREADS];
  pthread_t  threads[NTHREADS];

  nreads = vf->info['S']->given.count;
  rlen   = (int *) Malloc(sizeof(int)*nreads,"Allocating read length vector");
  if (rlen == NULL)
    exit (1);

  for (i = 0; i < NTHREADS; i++)
    { parm[i].beg = (nreads * i) / NTHREADS;
      parm[i].end = (nreads * (i+1)) / NTHREADS;
      parm[i].len = rlen;
      parm[i].vf  = vf+i;
      pthread_create(threads+i,NULL,fetch_thread,parm+i);
    }

  for (i = 0; i < NTHREADS; i++)
    pthread_join(threads[i],NULL);

  *rmax  = vf->info['S']->given.max;
  *nread = (int) nreads;

  oneFileClose(vf);

  return (rlen);
}



/****************************************************************************************
 *                                                                                        
 *  Find the nearest overlap record after point beg, return the file offset or -1 if not found.
 *                                                                                        
 ****************************************************************************************/

static int64 find_nearest(int fid, uint8 *buf, int64 beg, int *alast)
{ Overlap *ovl;
  Path    *path;
  uint8   *ouf, *trace8;
  uint16  *trace16;
  int64    pos, rb, x;
  int      bact, dact, tln;
  int      i, blen, diff;

#ifdef DEBUG_FIND
  fprintf(stderr,"\nFind starting at %lld\n",beg);
#endif

  lseek(fid,beg,SEEK_SET);
  rb  = read(fid,buf,IO_BLOCK);
  pos = OvlIOSize;
  ouf = buf-OvlSize;   
  while (pos < rb)
    { ovl  = (Overlap *) (ouf+pos);
      path = (Path *) ovl;
      tln  = path->tlen;

      if (ovl->aread >= NREAD1 || ovl->bread >= NREAD2 || tln <= 0)
        tln = 1;
      else if (path->aepos > RMAX1 || path->bepos > RMAX2)
        tln = 1;
      else if (path->abpos >= path->aepos || path->bbpos >= path->bepos)
        tln = 1;
      else if (tln != 2*((path->aepos + (TSPACE-1))/TSPACE - path->abpos/TSPACE))
        tln = 1;

      if (pos + tln*TBYTES > rb)
        { x = pos - OvlIOSize;
          beg += x;
          rb  -= x;
          memcpy(buf,buf+x,rb);
          rb += read(fid,buf+rb,IO_BLOCK-rb);
          pos = OvlIOSize;

          if (tln == 1)
            { if (pos > rb)
                return (-1);
            }
          else 
            { if (pos + tln*TBYTES > rb)
                { pos += 1;
                  if (pos > rb)
                    return (-1);
                }
            }

          continue;
        }

      if (tln == 1)
        { pos += 1;
          continue;
        }

#ifdef DEBUG_FIND
      fprintf(stderr,"\n%lld: %d[%d..%d] vs %d[%d..%d] tln = %d\n",pos,
                     ovl->aread,path->abpos,path->aepos,ovl->bread,path->bbpos,path->bepos,tln);
#endif

      bact = path->bepos - path->bbpos;
      dact = path->diffs;

      if (TBYTES == 1)
        { trace8 = (uint8 *) buf+pos;
          diff  = 0;
          blen  = 0;
          for (i = 0; i < tln; i += 2)
            { if (blen > bact)
                break;
              if (diff > dact)
                break;
              diff += trace8[i];
              blen += trace8[i+1];
            }
        }
      else
        { trace16 = (uint16 *) buf+pos;
          diff  = 0;
          blen  = 0;
          for (i = 0; i < tln; i += 2)
            { if (blen > bact)
                break;
              if (diff > dact)
                break;
              diff += trace16[i];
              blen += trace16[i+1];
            }
        }

#ifdef DEBUG_FIND
      fprintf(stderr,"   trace %d of %d: %d vs %d and %d vs %d\n",i,tln,blen,bact,diff,dact);
#endif

      if (i >= tln && blen == bact && diff == dact)
        if (path->aepos <= RLEN1[ovl->aread] && path->bepos <= RLEN2[ovl->bread])
          { *alast = ovl->aread;
            return ((beg+pos)+tln*TBYTES);
          }

      pos += 1;
    }

  return (-1);
}


/*******************************************************************************************
 *
 *  Parallel:  Each thread processes is a contiguous stripe across the .las input files
 *               sending the compressed binary data lines to their assigned OneFile.
 *
 ********************************************************************************************/

typedef struct
  { File_Object *fobj;     //  Reference to array of file objects
    uint8       *buf;      //  IO buffer for this thread
    OneFile     *vf;       //  OneFile for output
    int          bidx;     //  Scan range is [bidx:beg,eidx:end)
    int64        beg;
    int          eidx;
    int64        end;
    int          alast;
  } Thread_Arg;

  //  Write alignment records in relevant partition

static void *output_thread(void *arg)
{ Thread_Arg  *parm  = (Thread_Arg *) arg;
  File_Object *fobj  = parm->fobj;
  OneFile     *vf    = parm->vf;

  void        (*decon)(Overlap *, int64 *, int64 *);
  int          trmax;
  int64        epos, pos;
  FILE        *fid;
  int          f, al;
  Overlap     _ovl, *ovl = &_ovl;
  int64       *trace, *bdels, *diffs;

  if (TBYTES == sizeof(uint8))
    decon  = Decompress_Trace8;
  else
    decon  = Decompress_Trace16;

  trmax = 10000;
  trace = Malloc(3*trmax*sizeof(int64),"Allocating trace vector");
  if (trace == NULL)
    exit (1);
  bdels = trace + trmax;
  diffs = bdels + trmax;

  //  Do relevant section of each file assigned to this thread in sequence

  al = parm->alast;
  for (f = parm->bidx; f <= parm->eidx; f++)
    { fid = fopen(fobj[f].fname,"r");
      if (f < parm->eidx)
        epos = fobj[f].fsize;
      else
        epos = parm->end;
      if (f > parm->bidx || parm->beg == 0)
        parm->beg = sizeof(int) + sizeof(int64);

#ifdef DEBUG_OUT
      fprintf(stderr,"Block: %12lld to %12lld --> %8lld\n",parm->beg,epos,epos - parm->beg);
      fflush(stdout);
#endif

      pos = parm->beg;
      fseek(fid,pos,SEEK_SET); 
      while (pos < epos)
        { int      ar;
          int      blen;
          char     groupno[25];

          //  Read it in

          pos += Read_Overlap(fid,ovl);
          if (ovl->path.tlen > trmax)
            { trmax = 1.2*ovl->path.tlen + 1000;
              trace = Realloc(trace,3*trmax*sizeof(int64),"Reallocating trace vector");
              if (trace == NULL)
                exit (1);
              bdels = trace + trmax;
              diffs = bdels + trmax;
            }
          ovl->path.trace = (void *) trace;
          pos += Read_Trace(fid,ovl,TBYTES);

          //  Write requested lines to VgpFile

          ar = ovl->aread;
          if (ar != al)
            { if (DOGROUP)
                { int glen = sprintf(groupno,"%d",ar+1);

                  oneInt(vf,0) = 0;
                  oneInt(vf,1) = glen;
                  oneWriteLine(vf,'g',glen,groupno);
                }
              al = ar;
            }

          oneInt(vf,0) = ar+1;
          oneInt(vf,1) = ovl->bread+1;
          oneWriteLine(vf,'A',0,NULL);

          blen = RLEN2[ovl->bread];
          if (DOCOORD)
            { oneInt(vf,0) = ovl->path.abpos;
              oneInt(vf,1) = ovl->path.aepos;
              oneInt(vf,2) = RLEN1[ar];
              if (COMP(ovl->flags))
                { oneInt(vf,3) = ovl->path.bepos;
                  oneInt(vf,4) = ovl->path.bbpos;
                }
              else
                { oneInt(vf,3) = ovl->path.bbpos;
                  oneInt(vf,4) = ovl->path.bepos;
                }
              oneInt(vf,5) = blen;
              oneWriteLine(vf,'I',0,NULL);
            }
  
          if (DODIFF)
            { oneInt(vf,0) = ovl->path.diffs;
              oneWriteLine(vf,'D',0,NULL);
            }
  
          if (DOTRACE)
            { int tlen = (ovl->path.tlen >> 1);
  
              decon(ovl,bdels,diffs);
           
              oneInt(vf,0) = tlen;
              oneWriteLine(vf,'W',tlen,bdels);
              oneWriteLine(vf,'X',tlen,diffs);
            }
        }
    }

  free(trace);
  return (NULL);
}


/****************************************************************************************
 *                                                                                     
 *  The top-level program                                                              
 *                                                                                     
 ****************************************************************************************/

int main(int argc, char *argv[])
{ OneSchema *schema;
  char      *fname1, *fname2;
  char      *command;
  int        nfiles;
  int        ISTWO;

  //  Capture command line for provenance

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

    ARG_INIT("Dazz2sxs")

    NTHREADS = 4;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vidtg")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
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
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
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
        fprintf(stderr,"   -v: verbose mode, output progress as proceed\n");

        exit (1);
      }
  }

  //  Get read lengths and # of reads from sequence files in fnameA, rlenA, nreadA where A = 1,2

  { char    *pwd, *root;
    int      i;
    OneFile *vf;
    FILE    *input;
    char    *suffix[1] = { ".pbr" };

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

    if (VERBOSE)
      { fprintf(stderr,"  Scanning .pbr file %s\n",fname1);
        fflush(stderr);
      }

    vf = oneFileOpenRead(fname1,schema,"seq",NTHREADS);
    RLEN1 = Fetch_Length_Vector(vf,&NREAD1,&RMAX1);

    pwd = PathTo(argv[2]);
    OPEN(argv[2],pwd,root,input,suffix,1)
    if (input == NULL)
      { ISTWO  = 0;
        fname2 = fname1;
      }
    else
      { ISTWO  = 1;
        fname2 = Strdup(Catenate(pwd,"/",root,suffix[i]),"Allocating second sequence file");
        free(root);
        fclose(input);
      }
    free(pwd);

    if (ISTWO)
      { if (VERBOSE)
          { fprintf(stderr,"  Scanning .pbr file %s\n",fname2);
            fflush(stderr);
          }
        vf = oneFileOpenRead(fname2,schema,"seq",NTHREADS);
        RLEN2 = Fetch_Length_Vector(vf,&NREAD2,&RMAX2);
      }
    else
      { RLEN2  = RLEN1;
        NREAD2 = NREAD1;
        RMAX2  = RMAX1;
	fname2 = fname1;
     }

  }

  //  Determine number of .las files, nfiles

  { int           c;
    Block_Looper *parse;
    FILE         *input;

    nfiles = 0;
    for (c = 2+ISTWO; c < argc; c++)
      { parse = Parse_Block_LAS_Arg(argv[c]);
        while ((input = Next_Block_Arg(parse)) != NULL)
          { nfiles += 1;
            fclose(input);
          }
        Free_Block_Arg(parse);
      }

    if (VERBOSE)
      { fprintf(stderr,"  Partitioning %d .las files into %d parts\n",nfiles,NTHREADS);
        fflush(stderr);
      }
  }

  //  Find partition points dividing data in all files into NTHREADS roughly equal parts
  //    and then in parallel threads produce the output for each part.

  { File_Object fobj[nfiles];
    Thread_Arg  parm[NTHREADS];

#ifndef DEBUG_OUT
    pthread_t   threads[NTHREADS];
#endif

    { uint8        *bf;
      int           f, i, c;
      int64         b, work, wper;
      int           mspace;
      Block_Looper *parse;
      FILE         *input;

      //  Allocate IO buffer space for threads

      bf = Malloc(NTHREADS*IO_BLOCK,"Allocating IO_Buffer\n");
      if (bf == NULL)
        exit (1);

      //  Get name and size of each file in 'fobj[]', determin trace spacing and
      //    byte size of elements.-        

      TSPACE = -1;
      work   = 0;
      f      = 0;
      for (c = 2+ISTWO; c < argc; c++)
        { parse = Parse_Block_LAS_Arg(argv[c]);
          while ((input = Next_Block_Arg(parse)) != NULL)
            { struct stat info;
              int64       no;

              fobj[f].fname = Block_Arg_Name(parse);

              if (stat(fobj[f].fname, &info) == -1)
                { fprintf(stderr,"%s: Cannot get stats for %s\n",Prog_Name,fobj[f].fname);
                  exit (1);
                }

              fobj[f].fsize = info.st_size;
              work += info.st_size - (sizeof(int) + sizeof(int64));

              if (fread(&no,sizeof(int64),1,input) != 1)
                { fprintf(stderr,"%s: System error, read failed!\n",Prog_Name);
                  exit (1);
                }
              if (fread(&mspace,sizeof(int),1,input) != 1)
                { fprintf(stderr,"%s: System error, read failed!\n",Prog_Name);
                  exit (1);
                }

              if (TSPACE < 0)
                TSPACE = mspace;
              else if (mspace != TSPACE)
                { fprintf(stderr,"%s: Input .las files have different trace spacing!\n",Prog_Name);
                  exit (1);
                }

              f += 1;
              fclose(input);
            }
          Free_Block_Arg(parse);
        }

      if (TSPACE <= TRACE_XOVR && TSPACE != 0)
        TBYTES = sizeof(uint8);
      else
        TBYTES = sizeof(uint16);

      //  Allocate work evenly amongst threads, setting up search start
      //    point for each thread.  Also find the beginning of data in
      //    each file that a thread will start in (place in end.fpos)

      wper = work / NTHREADS;

      f = 0;
      b = 0;
      work = fobj[f].fsize;
      for (i = 0; i < NTHREADS; i++)
        { parm[i].fobj   = fobj;
          parm[i].buf    = bf + i*IO_BLOCK;

          if (b != 0)
            parm[i].end = sizeof(int) + sizeof(int64);
          else
            parm[i].end = 0;

#ifdef DEBUG_FIND
          fprintf(stderr," %2d: %1d %10lld (%10lld)\n",i,f,b,parm[i].end);
          fflush(stdout);
#endif

          parm[i].beg  = b;
          parm[i].bidx = f;

          work -= wper;
          while (work < .01*IO_BLOCK)
            { if (f == nfiles-1)
                { NTHREADS = i+1;
                  break;
                }
              work += fobj[++f].fsize;
            }
          b = fobj[f].fsize - work;
          if (b < 0)
            { work += b;
              b = 0;
            }
        }
    }

    { int i, f, fid;

      //  For each non-zero start point find synchronization point in
      //    bam/sam file.  If can't find it then start at beginning of
      //    next file, and if at or before first data line then signal
      //    start at beginning by zero'ing the synch point.

      for (i = 0; i < NTHREADS; i++)
        { if (parm[i].beg != 0)
            { fid = open(fobj[parm[i].bidx].fname,O_RDONLY);
              parm[i].beg = find_nearest(fid,parm[i].buf,parm[i].beg,&(parm[i].alast));
              if (parm[i].beg < 0 || parm[i].beg >= fobj[parm[i].bidx].fsize)
                { parm[i].beg   = 0;
                  parm[i].bidx += 1;
                }
              else if (parm[i].beg <= parm[i].end)
                parm[i].beg = 0;
              close(fid);
            }
          if (parm[i].beg == 0)
            parm[i].alast = -1;
        }

      //  Paranoid: if one thread's synch point overtakes the next one (will almost
      //    certainly never happen unless files very small and threads very large),
      //    remove the redundant threads.

      f = 0;
      for (i = 1; i < NTHREADS; i++)
        if (parm[i].bidx > parm[f].bidx || parm[i].beg > parm[f].beg)
          parm[++f] = parm[i];
      NTHREADS = f+1;

      //  Develop end points of each threads work using the start point of the next thread

      for (i = 1; i < NTHREADS; i++)
        if (parm[i].beg == 0)
          { parm[i-1].end  = fobj[parm[i].bidx-1].fsize;
            parm[i-1].eidx = parm[i].bidx-1;
          }
        else
          { parm[i-1].end  = parm[i].beg;
            parm[i-1].eidx = parm[i].bidx;
          }
      parm[NTHREADS-1].end  = fobj[nfiles-1].fsize;
      parm[NTHREADS-1].eidx = nfiles-1;

#if defined(DEBUG_FIND) || defined(DEBUG_OUT)
      fprintf(stderr,"\nPartition:\n");
      for (i = 0; i < NTHREADS; i++)
        { fprintf(stderr," %2d: %2d / %12lld",i,parm[i].bidx,parm[i].beg);
          fprintf(stderr,"  -  %2d / %12lld\n",parm[i].eidx,parm[i].end);
        }
      fflush(stdout);
#endif
    }

    //  Produce output in parallel threads based on partition

    { OneFile *vf;
      int      i;

      vf = oneFileOpenWriteNew("-",schema,"sxs",true,NTHREADS);

      oneAddProvenance(vf,Prog_Name,"1.0",command,NULL);
      oneAddReference(vf,fname1,NREAD1);
      oneAddReference(vf,fname2,NREAD2);

      oneWriteHeader(vf);

      oneInt(vf,0) = TSPACE;
      oneWriteLine(vf,'T',0,NULL);

#ifdef DEBUG_OUT
      fprintf(stderr,"Opened\n");
      fflush(stdout);
#endif

      if (VERBOSE)
        { fprintf(stderr,"  Producing .sxs segements in parallel\n");
          fflush(stderr);
        }

      //  Generate the data lines in parallel threads

#ifdef DEBUG_OUT
      for (i = 0; i < NTHREADS; i++)
        { parm[i].vf = vf+i;
          output_thread(parm+i);
        }
#else
      for (i = 0; i < NTHREADS; i++)
        { parm[i].vf = vf+i;
          pthread_create(threads+i,NULL,output_thread,parm+i);
        }

      for (i = 0; i < NTHREADS; i++)
        pthread_join(threads[i],NULL);
#endif

      oneFileClose(vf);
    }

    //  Free everything as a matter of good form

    { int f;

      free(fname1);
      free(RLEN1);
      if (ISTWO)
        { free(fname2);
          free(RLEN2);
        }
      for (f = 0; f < nfiles; f++)
        free(fobj[f].fname);
      free(parm[0].buf);
      free(command);
    }
  }

  oneSchemaDestroy(schema);

  if (VERBOSE)
    { fprintf(stderr,"  Done\n");
      fflush(stderr);
    }

  exit (0);
}
