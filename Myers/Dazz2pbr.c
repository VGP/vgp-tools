/*******************************************************************************************
 *
 *  Display a portion of the data-base and selected information in 1-code format.
 *
 *  Author:  Gene Myers
 *  Date  :  November 2015
 *
 ********************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdbool.h>

#include "gene_core.h"
#include "../Core/ONElib.h"

#include "VGPschema.h"

#define PATHSEP "/."

#define IO_BLOCK 10000000

static char *Usage = "[-vagu] [-T<int(4)>] <path:db>";


/*******************************************************************************************
 *
 *  Dazzler system codes ripped out of context to make Dazz2pbr stand-alone
 *
 ********************************************************************************************/

#define MAX_NAME 10000      //  Longest file name or fasta header line

#define DB_NFILE  "files = %9d\n"   //  number of files
#define DB_FDATA  "  %9d %s %s\n"   //  last read index + 1, fasta prolog, file name
#define DB_NBLOCK "blocks = %9d\n"  //  number of blocks
#define DB_PARAMS "size = %11lld cutoff = %9d all = %1d\n"  //  block size, len cutoff, all in well
#define DB_BDATA  " %9d %9d\n"      //  First read index (untrimmed), first read index (trimmed)

#define DB_QV   0x03ff   //  Mask for 3-digit quality value
#define DB_CSS  0x0400   //  This is the second or later of a group of subreads from a given insert
#define DB_BEST 0x0800   //  This is the "best" subread of a given insert (may be the only 1)

#define DB_ARROW 0x2     //  DB is an arrow DB
#define DB_ALL   0x1     //  all wells are in the trimmed DB

//  Fields have different interpretations if a .db versus a .dam

typedef struct
  { int     origin; //  Well # (DB), Contig # (DAM)
    int     rlen;   //  Length of the sequence (Last pulse = fpulse + rlen)
    int     fpulse; //  First pulse (DB), left index of contig in scaffold (DAM)
    int64   boff;   //  Offset (in bytes) of compressed read in 'bases' file, or offset of
                    //    uncompressed bases in memory block
    int64   coff;   //  Offset (in bytes) of compressed quiva streams in '.qvs' file (DB),
                    //  Offset (in bytes) of scaffold header string in '.hdr' file (DAM)
                    //  4 compressed shorts containing snr info if an arrow DB.
    int     flags;  //  QV of read + flags above (DB only)
  } DAZZ_READ;

typedef struct _track
  { struct _track *next;   //  Link to next track
    char          *name;   //  Symbolic name of track
    int            size;   //  Size in bytes of anno records
    int            nreads; //  Number of reads in track
    void          *anno;   //  over [0,nreads]: read i annotation: int, int64, or 'size' records
    int           *alen;   //  length of track data for read i (if data != NULL)
    void          *data;   //  data[anno[i] .. anno[i]+alen[i[) is data for read i (if data != NULL)
    int            loaded; //  Is track data loaded in memory?
    int64          dmax;   //  Largest read data segment in bytes
  } DAZZ_TRACK;

typedef struct
  { struct _track *next;
    char          *name;
    int64         *aoff;    //  offset in file or memory of arrow vector for read i
    void          *arrow;   //  FILE * to the .arw file if not loaded, memory block otherwise
    int            loaded;  //  Are arrow vectors loaded in memory?
  } DAZZ_ARROW;

typedef struct
  { int         ureads;     //  Total number of reads in untrimmed DB
    int         treads;     //  Total number of reads in trimmed DB
    int         cutoff;     //  Minimum read length in block (-1 if not yet set)
    int         allarr;     //  DB_ALL | DB_ARROW
    float       freq[4];    //  frequency of A, C, G, T, respectively

    //  Set with respect to "active" part of DB (all vs block, untrimmed vs trimmed)

    int         maxlen;     //  length of maximum read (initially over all DB)
    int64       totlen;     //  total # of bases (initially over all DB)

    int         nreads;     //  # of reads in actively loaded portion of DB
    int         trimmed;    //  DB has been trimmed by cutoff/all
    int         part;       //  DB block (if > 0), total DB (if == 0)
    int         ufirst;     //  Index of first read in block (without trimming)
    int         tfirst;     //  Index of first read in block (with trimming)

       //  In order to avoid forcing users to have to rebuild all thier DBs to accommodate
       //    the addition of fields for the size of the actively loaded trimmed and untrimmed
       //    blocks, an additional read record is allocated in "reads" when a DB is loaded into
       //    memory (reads[-1]) and the two desired fields are crammed into the first two
       //    integer spaces of the record.

    char       *path;       //  Root name of DB for .bps, .qvs, and tracks
    int         loaded;     //  Are reads loaded in memory?
    void       *bases;      //  file pointer for bases file (to fetch reads from),
                            //    or memory pointer to uncompressed block of all sequences.
    DAZZ_READ  *reads;      //  Array [-1..nreads] of DAZZ_READ
    DAZZ_TRACK *tracks;     //  Linked list of loaded tracks
  } DAZZ_DB;

static char       *atrack_name = ".@arw";
static DAZZ_DB    *Arrow_DB = NULL;         //  Last db/arw used by "Load_Arrow"
static DAZZ_ARROW *Arrow_Ptr;               //    Becomes invalid after closing

static int Open_DB(char* path, DAZZ_DB *db)
{ DAZZ_DB dbcopy;
  char   *root, *pwd, *bptr, *fptr;
  int     nreads;
  FILE   *index, *dbvis, *bases;
  int     status, plen, isdam;
  int     part, cutoff, all;
  int     ufirst, tfirst, ulast, tlast;
  char   *name;

  status = -1;
  dbcopy = *db;

  plen = strlen(path);
  if (strcmp(path+(plen-4),".dam") == 0)
    { root = Root(path,".dam");
      isdam = 1;
    }
  else
    { if (strcmp(path+(plen-3),".db") == 0)
        isdam = -1;
      else
        isdam = 0;
      root = Root(path,".db");
    }
  pwd = PathTo(path);

  bptr = rindex(root,'.');
  if (bptr != NULL && bptr[1] != '\0' && bptr[1] != '-')
    { part = strtol(bptr+1,&fptr,10);
      if (*fptr != '\0' || part == 0)
        part = 0;
      else
        *bptr = '\0';
    }
  else
    part = 0;

  name = Malloc(strlen(pwd)+strlen(root)+10,"Allocating name string");
  if (name == NULL)
    exit (1);

  if (isdam > 0)
    sprintf(name,"%s/%s.dam",pwd,root);
  else
    sprintf(name,"%s/%s.db",pwd,root);
  if (name == NULL)
    return (-1);
  if ((dbvis = fopen(name,"r")) == NULL)
    { if (isdam < 0)
        { fprintf(stderr,"%s: Could not open DB %s\n",Prog_Name,path);
          exit (1);
        }
      if (isdam > 0)
        { fprintf(stderr,"%s: Could not open DAM %s\n",Prog_Name,path);
          exit (1);
        }
      sprintf(name,"%s/%s.dam",pwd,root);
      if (name == NULL)
        return (-1);
      if ((dbvis = fopen(name,"r")) == NULL)
        { fprintf(stderr,"%s: Could not open %s as a DB or a DAM\n",Prog_Name,path);
          exit (1);
        }
      isdam = 1;
    }
  if (isdam < 0)
    isdam = 0;

  sprintf(name,"%s/.%s.idx",pwd,root);
  if ((index = fopen(name,"r")) == NULL)
    { fprintf(stderr,"%s: Cannot open .%s.idx for reading\n",Prog_Name,root);
      exit (1);
    }
  if (fread(db,sizeof(DAZZ_DB),1,index) != 1)
    { fprintf(stderr,"%s: Index file (.idx) of %s is junk\n",Prog_Name,root);
      exit (1);
    }

  { int   p, nblocks, nfiles;
    int64 size;
    char  fname[MAX_NAME], prolog[MAX_NAME];

    nblocks = 0;
    if (fscanf(dbvis,DB_NFILE,&nfiles) != 1)
      { fprintf(stderr,"%s: Stub file (.db) of %s is junk\n",Prog_Name,root);
        exit (1);
      }
    for (p = 0; p < nfiles; p++)
      if (fscanf(dbvis,DB_FDATA,&tlast,fname,prolog) != 3)
        { fprintf(stderr,"%s: Stub file (.db) of %s is junk\n",Prog_Name,root);
          exit (1);
        }
    if (fscanf(dbvis,DB_NBLOCK,&nblocks) != 1)
      if (part == 0)
        { cutoff = 0;
          all    = DB_ALL;
        }
      else
        { fprintf(stderr,"%s: DB %s has not yet been partitioned, cannot request a block !\n",
                         Prog_Name,root);
          exit (1);
        }
    else
      { if (fscanf(dbvis,DB_PARAMS,&size,&cutoff,&all) != 3)
          { fprintf(stderr,"%s: Stub file (.db) of %s is junk\n",Prog_Name,root);
            exit (1);
          }
        if (part > nblocks)
          { fprintf(stderr,"%s: DB %s has only %d blocks\n",Prog_Name,root,nblocks);
            exit (1);
          }
      }

    if (part > 0)
      { for (p = 1; p <= part; p++)
          if (fscanf(dbvis,DB_BDATA,&ufirst,&tfirst) != 2)
            { fprintf(stderr,"%s: Stub file (.db) of %s is junk\n",Prog_Name,root);
              exit (1);
            }
        if (fscanf(dbvis,DB_BDATA,&ulast,&tlast) != 2)
          { fprintf(stderr,"%s: Stub file (.db) of %s is junk\n",Prog_Name,root);
            exit (1);
          }
      }
    else
      { ufirst = tfirst = 0;
        ulast  = db->ureads;
        tlast  = db->treads;
      }
  }

  db->trimmed = 0;
  db->tracks  = NULL;
  db->part    = part;
  db->cutoff  = cutoff;
  db->allarr |= all;
  db->ufirst  = ufirst;
  db->tfirst  = tfirst;

  nreads = ulast-ufirst;
  if (part <= 0)
    { db->reads = (DAZZ_READ *) Malloc(sizeof(DAZZ_READ)*(nreads+2),"Allocating Open_DB index");
      if (db->reads == NULL)
        exit (1);

      db->reads += 1;
      if (fread(db->reads,sizeof(DAZZ_READ),nreads,index) != (size_t) nreads)
        { fprintf(stderr,"%s: Index file (.idx) of %s is junk\n",Prog_Name,root);
          exit (1);
        }
    }
  else
    { DAZZ_READ *reads;
      int        i, r, maxlen;
      int64      totlen;

      reads = (DAZZ_READ *) Malloc(sizeof(DAZZ_READ)*(nreads+2),"Allocating Open_DB index");
      if (reads == NULL)
        exit (1);
      reads += 1;

      fseeko(index,sizeof(DAZZ_READ)*ufirst,SEEK_CUR);
      if (fread(reads,sizeof(DAZZ_READ),nreads,index) != (size_t) nreads)
        { fprintf(stderr,"%s: Index file (.idx) of %s is junk\n",Prog_Name,root);
          exit (1);
        }

      totlen = 0;
      maxlen = 0;
      for (i = 0; i < nreads; i++)
        { r = reads[i].rlen;
          totlen += r;
          if (r > maxlen)
            maxlen = r;
        }

      db->maxlen = maxlen;
      db->totlen = totlen;
      db->reads  = reads;
    }

  ((int *) (db->reads))[-1] = ulast - ufirst;   //  Kludge, need these for DB part
  ((int *) (db->reads))[-2] = tlast - tfirst;

  db->nreads = nreads;
  sprintf(name,"%s/.%s",pwd,root);
  db->path   = Strdup(name,"Allocating Open_DB path");
  if (db->path == NULL)
    exit (1);
  sprintf(name,"%s/.%s.bps",pwd,root);
  bases = fopen(name,"r");
  if (bases == NULL)
    { fprintf(stderr,"%s: Cannot open %s/.%s.bps for reading\n",Prog_Name,pwd,root);
      exit (1);
    }
  db->bases = (void *) bases;
  db->loaded = 0;

  status = isdam;

  if (bptr != NULL)
    *bptr = '.';

  free(name);
  free(pwd);
  free(root);

  if (status < 0)
    *db = dbcopy;

  return (status);
}

static void Trim_DB(DAZZ_DB *db)
{ int         i, j, r, f;
  int         allflag, cutoff, css;
  int64       totlen;
  int         maxlen, nreads;
  DAZZ_TRACK *record;
  DAZZ_READ  *reads;

  if (db->trimmed) return;

  if (db->cutoff <= 0 && (db->allarr & DB_ALL) != 0) return;

  { int load_error;

    load_error = db->loaded;
    for (record = db->tracks; record != NULL; record = record->next)
      if (record->name == atrack_name)
        { if (((DAZZ_ARROW *) record)->loaded)
            load_error = 1;
        }
    if (load_error)
      { fprintf(stderr,"%s: Cannot load anything before trim (Trim_DB)\n",Prog_Name);
        exit (1);
      }
  }

  cutoff = db->cutoff;
  if ((db->allarr & DB_ALL) != 0)
    allflag = 0;
  else
    allflag = DB_BEST;

  reads  = db->reads;
  nreads = db->nreads;

  for (record = db->tracks; record != NULL; record = record->next)
    if (record->name == atrack_name)
      { DAZZ_ARROW *atrack = (DAZZ_ARROW *) record;
        int64      *aoff   = atrack->aoff;

        for (j = i = 0; i < nreads; i++)
          if ((reads[i].flags & DB_BEST) >= allflag && reads[i].rlen >= cutoff)
            aoff[j++] = aoff[i];
        atrack->aoff = Realloc(aoff,sizeof(int64)*j,NULL);
      }

  css    = 0;
  totlen = maxlen = 0;
  for (j = i = 0; i < nreads; i++)
    { f = reads[i].flags;
      if ((f & DB_CSS) == 0)
        css = 0;
      r = reads[i].rlen;
      if ((f & DB_BEST) >= allflag && r >= cutoff)
        { totlen += r;
          if (r > maxlen)
            maxlen = r;
          reads[j] = reads[i];
          if (css)
            reads[j++].flags |= DB_CSS;
          else
            reads[j++].flags &= ~DB_CSS;
          css = 1;
        }
    }
 
  db->totlen  = totlen;
  db->maxlen  = maxlen;
  db->nreads  = j;
  db->trimmed = 1;

  if (j < nreads)
    { db->reads = Realloc(reads-1,sizeof(DAZZ_READ)*(j+2),NULL);
      db->reads += 1;
    }
}

static void Close_Arrow(DAZZ_DB *db)
{ DAZZ_ARROW *atrack;

  Arrow_DB = NULL;
  if (db->tracks != NULL && db->tracks->name == atrack_name)
    { atrack = (DAZZ_ARROW *) db->tracks;
      if (atrack->loaded)
        free(atrack->arrow);
      else
        fclose((FILE *) atrack->arrow);
      free(atrack->aoff);
      db->tracks = db->tracks->next;
      free(atrack);
    }
}

static void Close_DB(DAZZ_DB *db)
{ if (db->loaded)
    free(((char *) (db->bases)) - 1);
  else if (db->bases != NULL)
    fclose((FILE *) db->bases);
  if (db->reads != NULL)
    free(db->reads-1);
  free(db->path);

  Close_Arrow(db);
}

static char *New_Read_Buffer(DAZZ_DB *db)
{ char *read;

  read = (char *) Malloc(db->maxlen+4,"Allocating New Read Buffer");
  if (read == NULL)
    exit (1);
  return (read+1);
}

static void Load_Read(DAZZ_DB *db, int i, char *read)
{ FILE      *bases  = (FILE *) db->bases;
  int64      off;
  int        len, clen;
  DAZZ_READ *r = db->reads;

  if (i < 0 || i >= db->nreads)
    { fprintf(stderr,"%s: Index out of bounds (Load_Read)\n",Prog_Name);
      exit (1);
    }

  if (db->loaded)
    { len = r[i].rlen;
      strncpy(read,(char *) bases + r[i].boff,len);
      read[len] = 4;
      Lower_Read(read);
      read[-1] = '\0';
      return;
    }

  off = r[i].boff;
  len = r[i].rlen;

  if (ftello(bases) != off)
    fseeko(bases,off,SEEK_SET);
  clen = COMPRESSED_LEN(len);
  if (clen > 0)
    { if (fread(read,clen,1,bases) != 1)
        { fprintf(stderr,"%s: Failed read of .bps file (Load_Read)\n",Prog_Name);
          exit (1);
        } 
    }   
  Uncompress_Read(len,read);
  Lower_Read(read);
  read[-1] = '\0';
} 

static void Open_Arrow(DAZZ_DB *db)
{ int64      *avector;
  DAZZ_ARROW *atrack;
  FILE       *afile;
  DAZZ_READ  *reads;
  int         i, nreads;
  char       *name;

  if (db->tracks != NULL && db->tracks->name == atrack_name)
    return;

  name = Malloc(strlen(db->path)+10,"Allocating name string");
  if (name == NULL)
    exit (1);

  if ((db->allarr & DB_ARROW) == 0)
    { fprintf(stderr,"%s: The DB is not an Arrow database (Open_Arrow)\n",Prog_Name);
      exit (1);
    }
  if (db->loaded)
    { fprintf(stderr,"%s: Cannot open Arrow vectors after loading all reads (Open_Arrow)\n",
                        Prog_Name);
      exit (1);
    }

  sprintf(name,"%s.arw",db->path);
  afile = fopen(name,"r");
  if (afile == NULL)
    { fprintf(stderr,"%s: Cannot open %s.arw for reading\n",Prog_Name,db->path);
      exit (1);
    }

  nreads  = db->nreads;
  avector = (int64 *) Malloc(sizeof(int64)*nreads,"Allocating Arrow index");
  atrack  = (DAZZ_ARROW *) Malloc(sizeof(DAZZ_ARROW),"Allocating Arrow track");
  if (avector == NULL || atrack == NULL)
    { fclose(afile);
      if (avector != NULL)
        free(avector);
      exit (1);
    }
  db->tracks     = (DAZZ_TRACK *) atrack;
  atrack->next   = NULL;
  atrack->name   = atrack_name;
  atrack->aoff   = avector;
  atrack->arrow  = (void *) afile;
  atrack->loaded = 0;


  reads = db->reads;
  for (i = 0; i < nreads; i++)
    avector[i] = reads[i].boff;
}

static int Load_Arrow(DAZZ_DB *db, int i, char *arrow)
{ FILE      *afile;
  int64      off;
  int        len, clen;

  if (db != Arrow_DB)
    { if (db->tracks == NULL || db->tracks->name != atrack_name)
        { fprintf(stderr,"%s: Arrow data is not available (Load_Arrow)\n",Prog_Name);
          exit (1);
        }
      Arrow_Ptr = (DAZZ_ARROW *) db->tracks;
      Arrow_DB  = db;
    }

  if (i < 0 || i >= db->nreads)
    { fprintf(stderr,"%s: Index out of bounds (Load_Arrow)\n",Prog_Name);
      exit (1);
    }

  afile = (FILE *) Arrow_Ptr->arrow;
  off   = Arrow_Ptr->aoff[i];
  len   = db->reads[i].rlen;

  if (ftello(afile) != off)
    fseeko(afile,off,SEEK_SET);
  clen = COMPRESSED_LEN(len);
  if (clen > 0)
    { if (fread(arrow,clen,1,afile) != 1)
        { fprintf(stderr,"%s: Failed read of .arw file (Load_Arrow)\n",Prog_Name);
          exit (1);
        }
    }
  Uncompress_Read(len,arrow);
  Letter_Arrow(arrow);
  arrow[-1] = '\0';
  return (0);
}


/*******************************************************************************************
 *
 *  Data output thread
 *
 ********************************************************************************************/

static int         DOARW, DOGRP, TRIM;
static char      **fhead = NULL;
static int        *findx = NULL;

typedef struct
  { OneFile     *vf;       //  OneFile for output
    int          beg;      //  Range of reads to output
    int          end;
    char        *arg;
  } Thread_Arg;

  //  Output reads [beg,end) as per the options into the OneFile for
  //     the thread.

static void *output_thread(void *arg)
{ Thread_Arg *parm = (Thread_Arg *) arg;
  OneFile    *vf   = parm->vf;
  int         beg  = parm->beg;
  int         end  = parm->end;

  DAZZ_DB    _db, *db = &_db;
  DAZZ_READ  *reads;
  char       *read, *arrow;
  int         i;
  int         map;

  Open_DB(parm->arg,db);

  if (DOARW)
    Open_Arrow(db);

  if (TRIM)
    Trim_DB(db);

  read = New_Read_Buffer(db);
  if (DOARW)
    arrow = New_Read_Buffer(db);
  else
    arrow = NULL;

  map = -1;
  while (findx[map+1] == 0)
    map += 1;
  while (beg > findx[map])
    map += 1;
  reads  = db->reads;

  for (i = beg; i < end; i++)
    { int         len;
      DAZZ_READ  *r;

      r   = reads + i;
      len = r->rlen;

      if (DOGRP && i == findx[map])
        { map += 1;
          len  = strlen(fhead[map]);
          oneInt(vf,0) = 0;
          oneInt(vf,1) = len;
          oneWriteLine(vf,'g',len,fhead[map]);
        }

      Load_Read(db,i,read);
      if (DOARW)
        Load_Arrow(db,i,arrow);

      oneInt(vf,0) = len;
      oneWriteLine(vf,'S',len,read);
      
      oneInt(vf,0) = r->origin;
      oneInt(vf,1) = r->fpulse;
      oneInt(vf,2) = r->fpulse+len;
      oneReal(vf,3) = (r->flags & DB_QV) / 1000.;
      oneWriteLine(vf,'W',0,NULL);

      if (DOARW)
        { int   j, snr[4];
          int64 big;

          big = *((uint64 *) &(r->coff));
          for (j = 0; j < 4; j++)
            { snr[3-j] = (big & 0xffff);
              big >>= 16;
            }
          oneReal(vf,0) = snr[0] / 100.;
          oneReal(vf,1) = snr[1] / 100.;
          oneReal(vf,2) = snr[2] / 100.;
          oneReal(vf,3) = snr[3] / 100.;
          oneWriteLine(vf,'N',0,NULL);

          oneInt(vf,0) = len;
          oneWriteLine(vf,'A',len,arrow);
        }
    }

  Close_DB(db);

  return (NULL);
}


/*******************************************************************************************
 *
 *  Main
 *
 ********************************************************************************************/

int main(int argc, char *argv[])
{ OneSchema  *schema;
  int         hasArrow;
  char       *command;
  int         nfiles;
  int64       nreads;
  DAZZ_DB    _db, *db = &_db;

  int         NTHREADS;
  int         VERBOSE;

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

  //  Process arguments

  { int  i, j, k;
    int  flags[128];
    char *eptr;

    ARG_INIT("Dazz2pbr")

    NTHREADS = 4;

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vagu")
            break;
          case 'T':
            ARG_POSITIVE(NTHREADS,"Number of threads")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE = flags['v'];
    TRIM    = 1-flags['u'];
    DOARW   = flags['a'];
    DOGRP   = flags['g'];

    if (argc <= 1)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"          W # # # #    - well, pulse start, end, and qv\n");
        fprintf(stderr,"          S #          - quality of read (#/1000)\n");
        fprintf(stderr,"      -a: N # # # #    - SNR of ACGT channels (#/100)\n");
        fprintf(stderr,"          A # string   - arrow pulse-width string\n");
        fprintf(stderr,"      -g: g # # string - cell size and name\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: verbose mode, output progress as proceed\n");
        fprintf(stderr,"      -u: Export untrimmed DB (default is trimmed DB).\n");
        fprintf(stderr,"      -T: Number of threads to use\n");
        exit (1);
      }
  }


  //  Open DB

  { int   status;

    status = Open_DB(argv[1],db);
    if (status < 0)
      exit (1);
    if (status == 1)
      { fprintf(stderr,"%s: Cannot export a .dam to .pbr format\n",Prog_Name);
        exit (1);
      }

    hasArrow = ((db->allarr & DB_ARROW) != 0);
    if (DOARW)
      { if (!hasArrow)
          { fprintf(stderr,"%s: -a option set but no Arrow data in DB\n",Prog_Name);
            exit (1);
          }
      }
  }

  if (VERBOSE)
    { fprintf(stderr,"  Analyzing contents of DB %s\n",argv[1]);
      fflush(stderr);
    }

  //  Load QVs if requested

  if (DOARW)
    Open_Arrow(db);

  //  Get SMRT cell names from the .db prolog

  { char *pwd, *root;
    FILE *dstub;
    char *dstub_name;
    int   i;

    root = Root(argv[1],".db");
    pwd = PathTo(argv[1]);
    if (db->part > 0)
      *rindex(root,'.') = '\0';
    dstub_name = Strdup(Catenate(pwd,"/",root,".db"),"Allocating db file name");
    if (dstub_name == NULL)
      exit (1);
    dstub = fopen(dstub_name,"r");
    if (dstub == NULL)
      { fprintf(stderr,"%s: Cannot open %s for reading\n",Prog_Name,dstub_name);
        exit (1);
      }
    free(pwd);
    free(root);

    fscanf(dstub,DB_NFILE,&nfiles);

    fhead = (char **) Malloc(sizeof(char *)*nfiles,"Allocating file list");
    findx = (int *) Malloc(sizeof(int *)*(nfiles+1),"Allocating file index");
    if (fhead == NULL || findx == NULL)
      exit (1);

    findx += 1;
    findx[-1] = 0;

    for (i = 0; i < nfiles; i++)
      { char prolog[MAX_NAME], fname[MAX_NAME];

        fscanf(dstub,DB_FDATA,findx+i,fname,prolog);
        if ((fhead[i] = Strdup(prolog,"Adding to file list")) == NULL)
          exit (1);
      }

    free(dstub_name);
    fclose(dstub);

    //  If TRIM (the default) then "trim" prolog ranges and the DB

    if (TRIM)
      { int        nid, oid, lid;
        int        cutoff, allflag;
        DAZZ_READ *reads;

        reads  = db->reads - db->ufirst;
        cutoff = db->cutoff;
        if ((db->allarr & DB_ALL) != 0)
          allflag = 0;
        else
          allflag = DB_BEST;
        
        nid = 0;
        oid = db->ufirst;
        lid = oid + db->nreads;
        for (i = 0; i < nfiles; i++)
          { while (oid < findx[i] && oid < lid)
              { if ((reads[oid].flags & DB_BEST) >= allflag && reads[oid].rlen >= cutoff)
                  nid++;
                oid += 1;
              }
            findx[i] = nid;
          }
        nreads = nid;
      }

    else if (db->part > 0)
      { for (i = 0; i < nfiles; i++)
          findx[i] -= db->ufirst;
        nreads = db->nreads;
      }

    else
      nreads = db->nreads;
  }

  Close_DB(db);

  //  Setup a OneFile for each thread, put the header in the first one

  { OneFile   *vf;
    Thread_Arg parm[NTHREADS];
    pthread_t  threads[NTHREADS];
    int        i;

    if (VERBOSE)
      { fprintf(stderr,"  Partitioning Dazzler DB into %d parts\n",NTHREADS);
        fflush(stderr);
      }

    vf = oneFileOpenWriteNew("-",schema,"pbr",true,NTHREADS);
    oneAddProvenance(vf,Prog_Name,"1.0",command,NULL);
    oneWriteHeader(vf);

    for (i = 0; i < NTHREADS; i++)
      { parm[i].vf = vf+i;
        parm[i].beg = (nreads * i) / NTHREADS;
        parm[i].end = (nreads * (i+1)) / NTHREADS;
        parm[i].arg = argv[1];
      }

    if (VERBOSE)
      { fprintf(stderr,"  Producing .pbr segements in parallel\n");
        fflush(stderr);
      }

    //  Generate the data lines in parallel threads

    for (i = 0; i < NTHREADS; i++)
      pthread_create(threads+i,NULL,output_thread,parm+i);

    for (i = 0; i < NTHREADS; i++)
      pthread_join(threads[i],NULL);

    oneFileClose(vf);

    for (i = 0; i < nfiles; i++)
      free(fhead[i]);
    free(fhead);
    free(findx-1);
  }

  oneSchemaDestroy(schema);

  if (VERBOSE)
    { fprintf(stderr,"  Done\n");
      fflush(stderr);
    }

  exit (0);
}
