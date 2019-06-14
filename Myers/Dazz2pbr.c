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

#include "gene_core.h"

#define HIDE_FILES

#ifdef HIDE_FILES
#define PATHSEP "/."
#else
#define PATHSEP "/"
#endif

static char *Usage = "[-aguU] <path:db|dam>";


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
  char   *root, *pwd, *bptr, *fptr, *cat;
  int     nreads;
  FILE   *index, *dbvis, *bases;
  int     status, plen, isdam;
  int     part, cutoff, all;
  int     ufirst, tfirst, ulast, tlast;

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

  if (isdam > 0)
    cat = Catenate(pwd,"/",root,".dam");
  else
    cat = Catenate(pwd,"/",root,".db");
  if (cat == NULL)
    return (-1);
  if ((dbvis = fopen(cat,"r")) == NULL)
    { if (isdam < 0)
        { fprintf(stderr,"%s: Could not open DB %s\n",Prog_Name,path);
          exit (1);
        }
      if (isdam > 0)
        { fprintf(stderr,"%s: Could not open DAM %s\n",Prog_Name,path);
          exit (1);
        }
      cat = Catenate(pwd,"/",root,".dam");
      if (cat == NULL)
        return (-1);
      if ((dbvis = fopen(cat,"r")) == NULL)
        { fprintf(stderr,"%s: Could not open %s as a DB or a DAM\n",Prog_Name,path);
          exit (1);
        }
      isdam = 1;
    }
  if (isdam < 0)
    isdam = 0;

  if ((index = fopen(Catenate(pwd,PATHSEP,root,".idx"),"r")) == NULL)
    { fprintf(stderr,"%s: Cannot open %s for reading\n",Prog_Name,Catenate(".",root,".idx",""));
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
  db->path   = Strdup(Catenate(pwd,PATHSEP,root,""),"Allocating Open_DB path");
  if (db->path == NULL)
    exit (1);
  bases = fopen(Catenate(db->path,"","",".bps"),"r");
  if (bases == NULL)
    { fprintf(stderr,"%s: Cannot open %s for reading\n",Prog_Name,Catenate(db->path,".bps","",""));
      exit (1);
    }
  db->bases = (void *) bases;
  db->loaded = 0;

  status = isdam;

  if (bptr != NULL)
    *bptr = '.';

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

static void Load_Read(DAZZ_DB *db, int i, char *read, int ascii)
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
      if (ascii == 0)
        { if (*read < 4)
            read[-1] = read[len] = 4;
          else
            { read[len] = '\0';
              Number_Read(read);
              read[-1] = 4;
            }
        }
      else
        { if (*read < 4)
            { read[len] = 4;
              if (ascii == 1)
                Lower_Read(read);
              else
                Upper_Read(read);
              read[-1] = '\0';
            }
          else
            { read[len] = '\0';
              if ((ascii == 1) != islower(*read))
                Change_Read(read);
            }
          read[-1] = '\0';
        }
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
  if (ascii == 1)
    { Lower_Read(read);
      read[-1] = '\0';
    } 
  else if (ascii == 2)
    { Upper_Read(read);
      read[-1] = '\0';
    } 
  else
    read[-1] = 4;
} 

static void Open_Arrow(DAZZ_DB *db)
{ int64      *avector;
  DAZZ_ARROW *atrack;
  FILE       *afile;
  DAZZ_READ  *reads;
  int         i, nreads;

  if (db->tracks != NULL && db->tracks->name == atrack_name)
    return;

  if ((db->allarr & DB_ARROW) == 0)
    { fprintf(stderr,"%s: The DB is not an Arrow database (Open_Arrow)\n",Prog_Name);
      exit (1);
    }
  if (db->loaded)
    { fprintf(stderr,"%s: Cannot open Arrow vectors after loading all reads (Open_Arrow)\n",
                        Prog_Name);
      exit (1);
    }

  afile = fopen(Catenate(db->path,"","",".arw"),"r");
  if (afile == NULL)
    { fprintf(stderr,"%s: Cannot open %s for reading\n",Prog_Name,Catenate(db->path,".arw","",""));
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

static int Load_Arrow(DAZZ_DB *db, int i, char *arrow, int ascii)
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
  if (ascii == 1)
    { Letter_Arrow(arrow);
      arrow[-1] = '\0';
    }
  else
    arrow[-1] = 4;
  return (0);
}


/*******************************************************************************************
 *
 *  The main code
 *
 ********************************************************************************************/

int main(int argc, char *argv[])
{ DAZZ_DB    _db, *db = &_db;
  int         Arrow_DB;
  int         FirstRead;
  int         maxgrp = 0;
  int64      *gcount = NULL;

  int         nfiles;
  char      **fhead = NULL;
  int        *findx = NULL;

  int         TRIM, UPPER;
  int         DOARW, DOGRP;

  //  Process arguments

  { int  i, j, k;
    int  flags[128];

    ARG_INIT("Dazz2pbr")

    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("aguU")
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    TRIM  = 1-flags['u'];
    UPPER = flags['U'];
    DOARW = flags['a'];
    DOGRP = flags['g'];

    if (argc <= 1)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage);
        fprintf(stderr,"\n");
        fprintf(stderr,"          W # # # #    - well, pulse start, end, and qv\n");
        fprintf(stderr,"          S #          - quality of read (#/1000)\n");
        fprintf(stderr,"      -a: N # # # #    - SNR of ACGT channels (#/100)\n");
        fprintf(stderr,"          A # string   - arrow pulse-width string\n");
        fprintf(stderr,"      -g: g # # string - cell size and name\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -u: Export untrimmed DB (default is trimmed DB).\n");
        fprintf(stderr,"      -U: Output base pairs in upper case letters\n");
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

    Arrow_DB = ((db->allarr & DB_ARROW) != 0);
    if (DOARW)
      { if (!Arrow_DB)
          { fprintf(stderr,"%s: -a option set but no Arrow data in DB\n",Prog_Name);
            exit (1);
          }
      }
  }

  //  Load QVs if requested

  if (DOARW)
    Open_Arrow(db);

  //  If get prolog and file names and index ranges from the .db or .dam file 

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
      }

    else if (db->part > 0)
      { for (i = 0; i < nfiles; i++)
          findx[i] -= db->ufirst;
      }
  }

  if (TRIM)
    { Trim_DB(db);
      FirstRead = db->tfirst;
    }
  else
    FirstRead = db->ufirst;

  //  Scan to count the size of things

  { DAZZ_READ  *reads;
    int         i, map;
    int64       noreads, ngroup;
    int64       seqmax, seqtot;
    int64       gmxread, gmxtot;
    int64       gtotc, gmaxc;
    int64       gread, gtotbp;

    map    = 0;
    reads  = db->reads;

    noreads = 0;
    seqmax  = 0;
    seqtot  = 0;
    ngroup  = 0;
    gmxread = 0;
    gmxtot  = 0;
    gtotc   = 0;
    gmaxc   = 0;

    gread   = 0;
    gtotbp  = 0;
    map     = -1;

    for (i = 0; i < db->nreads; i++)
      { int         len;
        DAZZ_READ  *r;

        r   = reads + i;
        len = r->rlen;

        noreads += 1;
        if (len > seqmax)
          seqmax = len;
        seqtot += len;

        if (DOGRP && i+FirstRead >= findx[map])
          { if (gread > 0)
              { if (gread > gmxread)
                  gmxread = gread;
                if (gtotbp > gmxtot)
                  gmxtot = gtotbp;
                if (ngroup >= maxgrp)
                  { maxgrp = 1.2*ngroup + 100;
                    gcount = Realloc(gcount,sizeof(int64)*maxgrp,"Allocating count vector");
                  }
                gcount[ngroup++] = gread;
              }

            gread  = 1;
            gtotbp = len;
            map   += 1;

            len = strlen(fhead[map]);
            gtotc += len;
            if (len > gmaxc)
              gmaxc = len; 
          }
        else
          { gread  += 1;
            gtotbp += len;
          }
      }
    if (gread > 0)
      { if (gread > gmxread)
          gmxread = gread;
        if (gtotbp > gmxtot)
          gmxtot = gtotbp;
        if (ngroup >= maxgrp)
          { maxgrp = 1.2*ngroup + 100;
            gcount = Realloc(gcount,sizeof(int64)*maxgrp,"Allocating count vector");
          }
        gcount[ngroup++] = gread;
      }

    //  Output file type and provenance

    { int    clen, optl;
      char   date[26];
      time_t seconds;

      printf("1 3 seq 1 0\n");
      printf("2 3 pbr\n");
      // printf("# ! 1\n");
      if (DOGRP)
        printf("# g %lld\n",ngroup);
      printf("# S %lld\n",noreads);
      printf("# W %lld\n",noreads);
      if (DOARW)
        { printf("# A %lld\n",noreads);
          printf("# N %lld\n",noreads);
        }

      optl = (1-TRIM) + UPPER + DOGRP + DOARW;
      if (optl == 0)
        clen = -1;
      else
        clen = optl+1;
      for (i = 1; i < argc; i++)
        clen += strlen(argv[i])+1;

      // printf("+ ! %d\n",clen+35);
      if (DOGRP)
        printf("+ g %lld\n",gtotc);
      printf("+ S %lld\n",seqtot);
      if (DOARW)
        printf("+ A %lld\n",seqtot);

      // if (clen > 24)
        // printf("@ ! %d\n",clen);
      // else
        // printf("@ ! 24\n");
      if (DOGRP)
        printf("@ g %lld\n",gmaxc);
      printf("@ S %lld\n",seqmax);
      if (DOARW)
        printf("@ A %lld\n",seqmax);

      if (DOGRP)
        { printf("%% g # S %lld\n",gmxread);
          printf("%% g # W %lld\n",gmxread);
          if (DOARW)
            printf("%% g # A %lld\n",gmxread);
          printf("%% g + S %lld\n",gmxtot);
          if (DOARW)
            printf("%% g + A %lld\n",gmxtot);
        }

      printf("! 8 Dazz2pbr 3 1.0 %d",clen);
      if (optl > 0)
        { printf(" -");
          if (DOARW)
            printf("a");
          if (DOGRP)
            printf("g");
          if (!TRIM)
            printf("u");
          if (UPPER)
            printf("U");
        }
      for (i = 1; i < argc; i++)
        printf(" %s",argv[i]);
      seconds = time(NULL);
      ctime_r(&seconds,date);
      date[24] = '\0';
      printf(" 24 %s\n",date);
    }
  }

  //  Display each read (and/or QV streams) in the active DB according to the
  //    range pairs in pts[0..reps) and according to the display options.

  { DAZZ_READ  *reads;
    char       *read, *arrow;
    int         i;
    int         map, ngroup;

    read = New_Read_Buffer(db);
    if (DOARW)
      arrow = New_Read_Buffer(db);
    else
      arrow = NULL;
    UPPER += 1;

    map    = -1;
    ngroup = 0;
    reads  = db->reads;

    for (i = 0; i < db->nreads; i++)
      { int         len;
        int         flags, qv;
        DAZZ_READ  *r;

        r   = reads + i;
        len = r->rlen;
        flags = r->flags;
        qv    = (flags & DB_QV);

        if (DOGRP && i+FirstRead >= findx[map])
          { printf("g %lld %ld %s\n",gcount[ngroup++],strlen(fhead[map+1]),fhead[map+1]);
            map += 1;
          }

        Load_Read(db,i,read,UPPER);
        if (DOARW)
          Load_Arrow(db,i,arrow,1);

        printf("S %d ",len);
        printf("%.*s\n",len,read);
        printf("W %d %d %d 0.%03d\n",r->origin,r->fpulse,r->fpulse+len,qv);

        if (DOARW)
          { int   j, snr[4];
            int64 big;

            big = *((uint64 *) &(r->coff));
            for (j = 0; j < 4; j++)
              { snr[3-j] = (big & 0xffff);
                big    >>= 16;
              }
            printf("N %d %d %d %d\n",snr[0],snr[1],snr[2],snr[3]);
            printf("A %d ",len);
            printf("%.*s\n",len,arrow);
          }
      }
  }

  fclose(stdout);

  { int i;

    for (i = 0; i < nfiles; i++)
      free(fhead[i]);
    free(fhead);
    free(findx-1);
  }

  Close_DB(db);

  exit (0);
}
