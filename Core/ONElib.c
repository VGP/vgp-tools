/*****************************************************************************************
 *
 *  File: ONElib.c
 *    implementation for ONElib.h
 *
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University and Eugene Myers 2019-
 *
 * HISTORY:
 * Last edited: Apr 23 02:46 2020 (rd109)
 * * Apr 23 00:31 2020 (rd109): global rename of VGP to ONE, Vgp to One, vgp to one
 * * Apr 20 11:27 2020 (rd109): added VgpSchema to make schema dynamic
 * * Dec 27 09:46 2019 (gene): style edits + compactify code
 * * Jul  8 04:28 2019 (rd109): refactored to use info[]
 * * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *
 ****************************************************************************************/

#include <assert.h>
#include <sys/errno.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <stdarg.h>
#include <time.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/uio.h>
#include <math.h>

#include "ONElib.h"

// set major and minor code versions

#define MAJOR 1
#define MINOR 1

//  utilities with implementation at the end of the file

static void  die(char *format, ...);                  //  print message to stderr and exit -1
static void *myalloc(size_t size);                    //  allocate block, die if malloc fails
static void *mycalloc(size_t number, size_t size);    //  allocate & zero # objects of size
#define new(n,type)  (type *) myalloc((n)*sizeof(type))    // actually use these not myalloc
#define new0(n,type)  (type *) mycalloc((n),sizeof(type))

// global required for parallelisation

static pthread_mutex_t mutexInit = PTHREAD_MUTEX_INITIALIZER;

// forward declarations of serialisation functions lower in the file

OneCodec *vcCreate();
void      vcAddToTable(OneCodec *vc, int len, char *bytes);
void      vcAddHistogram(OneCodec *vc, OneCodec *vh);
void      vcCreateCodec(OneCodec *vc, int partial);
void      vcDestroy(OneCodec *vc);
int       vcMaxSerialSize();
int       vcSerialize(OneCodec *vc, void *out);
OneCodec *vcDeserialize(void *in);
int       vcEncode(OneCodec *vc, int ilen, char *ibytes, char *obytes);
int       vcDecode(OneCodec *vc, int ilen, char *ibytes, char *obytes);

/***********************************************************************************
 *
 *    ONE_FILE CREATION & DESTRUCTION
 *
 **********************************************************************************/

/******************* OneInfo ********************/

static OneInfo *infoCreate (int nField)
{ OneInfo *vi = new0 (1, OneInfo) ;
  vi->nField = nField ;
  if (nField) vi->fieldType = new (nField, OneType) ;
  return vi;
}

static OneInfo *infoDeepCopy (OneInfo *vi0)
{ OneInfo *vi = new (1, OneInfo) ;
  *vi = *vi0 ;
  if (vi->nField)
    { vi->fieldType = new (vi->nField, OneType) ;
      memcpy (vi->fieldType, vi0->fieldType, vi->nField*sizeof(OneType)) ;
    }
  if (vi->fieldCodec) vi->fieldCodec = vcCreate() ;
  if (vi->listCodec && vi->listCodec != DNAcodec) vi->listCodec = vcCreate() ;
  return vi ;
}

static void infoDestroy (OneInfo *vi)
{ if (vi->buffer && ! vi->isUserBuf) free (vi->buffer) ;
  if (vi->fieldCodec) vcDestroy (vi->fieldCodec) ;
  if (vi->listCodec) vcDestroy (vi->listCodec) ;
  if (vi->fieldType) free (vi->fieldType) ;
  free(vi);
}

/******************* OneSchema ********************/

// a bit of sugar to set the OneInfo list information
static int listEltSize[8] = { 0, 0, 0, 0, 1, sizeof(I64), sizeof(double), 1 };
#define VL(vi,i,t,vs)    if (vi->listEltSize) \
    die ("multiple list types for linetype %c in schema for filetype %s", t, vs->primary) ; \
    vi->listEltSize = listEltSize[vi->fieldType[i]] ; \
    vi->listField = i

static OneSchema *schemaLoadRecord (OneSchema *vs, OneFile *vf)
{
  // parse a schema specfication line from vf and add into vs
  // return value is vs unless a new primary type is declared, in which case vs->nxt
  char t ;

  switch (vf->lineType)
    {
    case 'X':  // ignore - blank or comment line in schema file
      break ;
    case 'P':
      if (*vs->primary && !vs->objectType)
	die ("schema: file type %s has no object type", vs->primary) ;
      if (oneLen(vf) != 3) die ("schema: primary name %s is not 3 letters", oneString(vf)) ;
      OneSchema *vsNxt = new0 (1, OneSchema) ;
      vsNxt->major = vs->major ; vsNxt->minor = vs->minor ;
      vsNxt->nBinaryHeader = vs->nBinaryHeader ; // transfer the header count
      vsNxt->nBinary = vsNxt->nBinaryHeader ; // start count for this schema at header count
      vs->nxt = vsNxt ;
      vs = vsNxt ;
      strcpy (vs->primary, oneString(vf)) ;
      break ;
    case 'S':
      if (oneLen(vf) != 3) die ("schema: secondary name %s is not 3 letters", oneString(vf)) ;
      if (vs->nSecondary)
	{ char **temp = vs->secondary ;
	  vs->secondary = new (vs->nSecondary+1, char*) ;
	  memcpy (vs->secondary, temp, vs->nSecondary*sizeof(char*)) ;
	  free (temp) ;
	}
      else
	vs->secondary = new (1, char*) ;
      vs->secondary[vs->nSecondary] = new0 (4, char) ;
      strcpy (vs->secondary[vs->nSecondary++], oneString(vf)) ;
      break ;
    case 'V':
      vs->major = oneInt(vf,0) ; vs->minor = oneInt(vf,1) ;
      break ;
    case 'A': case 'L': case 'C': case 'F': case 'E':
      t = oneChar(vf,0) ;
      if (vs->info[t])
	die ("duplicate schema specification for linetype %c in filetype %s", t, vs->primary) ;
      if (t >= 'a' && t <= 'z') // the group type
	{ if (vs->groupType) die ("second group type in schema for filetype %s", vs->primary) ;
	  vs->groupType = t ;
	}
      else if (!vs->objectType && t >= 'A' && t <= 'Z')
	vs->objectType = t ;
      else if ((t < 'A' || t > 'Z') && *vs->primary) // allow non-alphabetic lines in header
	die ("non-alphabetic linetype %c (ascii %d) in schema for filetype %s",t,t,vs->primary) ;
      OneInfo *vi = vs->info[t] = infoCreate (oneInt(vf,1)) ;
      int i = 0 ; char *s = oneString(vf) ;
      for (i = 0 ; i < vi->nField ; ++i, s = oneNextString(vf,s))
	if (!strcmp (s, "INT")) vi->fieldType[i] = vINT ;
	else if (!strcmp (s, "REAL")) vi->fieldType[i] = vREAL ;
	else if (!strcmp (s, "CHAR")) vi->fieldType[i] = vCHAR ;
	else if (!strcmp (s, "STRING")) vi->fieldType[i] = vSTRING ;
	else if (!strcmp (s, "INT_LIST")) { vi->fieldType[i] = vINT_LIST ; VL(vi,i,t,vs) ; }
	else if (!strcmp (s, "REAL_LIST")) { vi->fieldType[i] = vREAL_LIST ; VL(vi,i,t,vs) ; }
	else if (!strcmp (s, "STRING_LIST")) { vi->fieldType[i] = vSTRING_LIST; VL(vi,i,t,vs) ; }
	else if (!strcmp (s, "DNA"))
	  { vi->fieldType[i] = vSTRING ; VL(vi,i,t,vs) ;
	    vi->listCodec = DNAcodec ; vi->isUseListCodec = TRUE ;
	    if (vf->lineType != 'L') die ("linetype in schema for DNA lines must be L") ;
	  }
	else
	  die ("ONE schema error: bad field %d of %d type %s in line %d",
	       i, vi->nField, s, vf->line) ;
      vi->isIntListDiff = TRUE ; // harmless if not INT_LIST (and will go)
      if (vf->lineType == 'C' || vf->lineType == 'B') vi->listCodec = vcCreate () ;
      if (vf->lineType == 'F' || vf->lineType == 'B') vi->fieldCodec = vcCreate () ;
      if (vf->lineType != 'A') // create binary encoding
	vi->binaryTypePack = (vs->nBinary++ << 2) | 0x80;
      if (vs->nBinary >= 32)
	die ("ONE schema error file type %s: too many line specs >= 32", vs->primary) ;
      if (vi->nField > vs->nFieldMax)
	vs->nFieldMax = vi->nField ;
      break ;
    default:
      die ("unrecognized schema line %d starting with %c", vf->line, vf->lineType) ;
    }

  return vs ;
}

static void oneFileDestroy (OneFile *vf) ; // need a forward declaration here

OneSchema *oneSchemaCreateFromFile (char *filename)
{
  FILE *fs = fopen (filename, "r") ;
  if (!fs) return 0 ;
  OneSchema *vs = new0 (1, OneSchema) ;

  OneFile *vf = new0 (1, OneFile) ;      // shell object to support bootstrap
  // bootstrap specification of linetypes to read schemas
  vf->info['A'] = infoCreate (2) ;       // ASCII only linetype - used for many header line types
  vf->info['A']->fieldType[0] = vCHAR ;
  vf->info['A']->fieldType[1] = vSTRING_LIST ; vf->info['A']->listField = 1 ;
  vf->info['L'] = infoCreate (2) ;       // general line type, which can be binary
  vf->info['L']->fieldType[0] = vCHAR ;
  vf->info['L']->fieldType[1] = vSTRING_LIST ; vf->info['L']->listField = 1 ;
  vf->info['P'] = infoCreate (1) ;       // to define the schema for parsing a .def file
  vf->info['P']->fieldType[0] = vSTRING ; vf->info['P']->listField = 0 ;
  vf->info['/'] = infoCreate (0) ;       // to store comments
  vf->field = new (2, OneField) ;

  // first load the universal header and footer (non-alphabetic) line types 
  // do this by writing their schema into a temporary file and parsing it into the base schema
  //  vf->f = tmpfile () ;                   // unique temp file destroyed on program exit
  vf->f = fopen ("TEST","w+") ;              // unique temp file destroyed on program exit
  fprintf (vf->f, "A 1 3 6 STRING 3 INT 3 INT         first line: 3-letter type, major, minor version\n") ;
  fprintf (vf->f, "A 2 1 6 STRING                     subtype: 3-letter subtype\n") ;
  fprintf (vf->f, "A # 2 4 CHAR 3 INT                 linetype, count\n") ;
  fprintf (vf->f, "A @ 2 4 CHAR 3 INT                 linetype, list max\n") ;
  fprintf (vf->f, "A + 2 4 CHAR 3 INT                 linetype, list total\n") ;
  fprintf (vf->f, "A %% 4 4 CHAR 4 CHAR 4 CHAR 3 INT  group, #/+, linetype, value\n") ;
  fprintf (vf->f, "A ! 1 11 STRING_LIST               provenance: program, version, command, date\n") ;
  fprintf (vf->f, "A < 2 6 STRING 3 INT               reference: filename, object count\n") ;
  fprintf (vf->f, "A > 1 6 STRING                     deferred: filename\n") ;
  fprintf (vf->f, "A $ 1 3 INT                        binary file - goto footer: isBigEndian\n") ;
  fprintf (vf->f, "A ^ 0                              binary file: end of footer designation\n") ;
  fprintf (vf->f, "A - 1 3 INT                        binary file: offset of start of footer\n") ;
  fprintf (vf->f, "L & 1 8 INT_LIST                   binary file: object index\n") ;
  fprintf (vf->f, "L * 1 8 INT_LIST                   binary file: group index\n") ;
  fprintf (vf->f, "L : 1 6 STRING                     binary file: field codec\n") ;
  fprintf (vf->f, "L ; 1 6 STRING                     binary file: list codec\n") ;
  fprintf (vf->f, "A / 1 6 STRING                     binary file: comment\n") ;
  if (fseek (vf->f, 0, SEEK_SET)) die ("ONE schema failure: cannot rewind tmp file") ;
  while (oneReadLine (vf))
    schemaLoadRecord (vs, vf) ;
  vs->major = MAJOR ; vs->minor = MINOR ;

  // next reuse the temp file to load the schema for reading schemas
  if (fseek (vf->f, 0, SEEK_SET)) die ("ONE schema failure: cannot rewind tmp file") ;
  fprintf (vf->f, "P 3 def                      this is the primary file type for schemas\n") ;
  fprintf (vf->f, "L P 1 6 STRING               primary type name\n") ;
  fprintf (vf->f, "L S 1 6 STRING               secondary type names\n") ;
  fprintf (vf->f, "L V 2 3 INT 3 INT            versions: major, minor - group field\n") ;
  fprintf (vf->f, "L L 2 4 CHAR 11 STRING_LIST  basic linetype no compression\n") ;
  fprintf (vf->f, "L C 2 4 CHAR 11 STRING_LIST  linetype with list compression\n") ;
  fprintf (vf->f, "L F 2 4 CHAR 11 STRING_LIST  linetype with field compression\n") ;
  fprintf (vf->f, "L B 2 4 CHAR 11 STRING_LIST  linetype with both list and field compression\n") ;
  fprintf (vf->f, "L X 0                        blank lines, for spacers or comments\n") ;
  fprintf (vf->f, "\n") ; // terminator
  if (fseek (vf->f, 0, SEEK_SET)) die ("ONE schema failure: cannot rewind tmp file") ;
  OneSchema *vs0 = vs ;  // need this because loadInfo() updates vs on reading P lines
  vf->line = 0 ;
  while (oneReadLine (vf))
    vs = schemaLoadRecord (vs, vf) ;
  OneSchema *vsDef = vs ; // will need this to destroy it once the true schema is read
  oneFileDestroy (vf) ;   // this will also effectively remove the temp file on closing

  // finally read the schema itself
  if (!(vf = oneFileOpenRead (filename, vs0, "def", 1)))
    return 0 ;
  vs = vs0 ; // set back to vs0, so next filetype spec will replace vsScm
  while (oneReadLine (vf))
    vs = schemaLoadRecord (vs, vf) ;
  oneFileDestroy (vf) ;
  oneSchemaDestroy (vsDef) ; // no longer need this, and can destroy because unlinked from vs0
  
  return vs0 ;
}

void oneSchemaDestroy (OneSchema *vs)
{ int i ;
  while (vs)
    { for (i = 0 ; i < 128 ; ++i) if (vs->info[i]) infoDestroy (vs->info[i]) ;
      if (vs->nSecondary)
	{ for (i = 0 ; i < vs->nSecondary ; ++i) free (vs->secondary[i]) ;
	  free (vs->secondary) ;
	}
      OneSchema *t = vs->nxt ;
      free (vs) ;
      vs = t ;
    }
}

/*************************************/

static OneFile *oneFileCreate (OneSchema *vs, char *type)
{ // searches through the linked list of vs to find type, either as primary or a secondary
  // if found fills and returns vf, else returns 0
  
  int         i, j ;
  OneFile    *vf = new0 (1, OneFile) ;
  char       *secondary = 0 ;

  // transfer header info
  for (i = 0 ; i < 128 ; ++i)
    if (vs->info[i]) vf->info[i] = infoDeepCopy (vs->info[i]) ;

  // find type in schema 
  while ((vs = vs->nxt))
    if (!strcmp (type, vs->primary))
      break ;
    else if (vs->nSecondary)
      { for (j = 0 ; j < vs->nSecondary ; ++j)
	  if (!strcmp (type, vs->secondary[j])) break ;
	if (j < vs->nSecondary) { secondary = vs->secondary[j] ; break ; }
      }
  if (!vs)
    return 0 ; // failed to find a match
  
  // transfer info from matched schema
  for (i = 0 ; i < 128 ; ++i)
    if (vs->info[i]) vf->info[i] = infoDeepCopy (vs->info[i]) ;

  // set other information
  vf->major = vs->major ; vf->minor = vs->minor ;
  vf->objectType = vs->objectType ;
  vf->groupType = vs->groupType ;
  strcpy (vf->fileType, vs->primary) ;
  if (secondary) strcpy (vf->subType, secondary) ;
  vf->codecTrainingSize = 100000;
  vf->nFieldMax = vs->nFieldMax ;
  vf->field = new (vf->nFieldMax, OneField) ;

  // determine endian of machine
  { int   t = 1;
    char *b = (char *) (&t);
    vf->isBig = (*b == 0);
  }

  return vf ;
}

static void provRefDefCleanup (OneFile *vf)
{ int n ;

  if (vf->provenance)
    { OneProvenance *p = vf->provenance ;
      for (n = vf->info['!']->accum.count ; n-- ; p++)
	{ free (p->program) ;
	  free (p->version) ;
	  free (p->command) ;
	  free (p->date) ; }
      free (vf->provenance) ;
    }
  if (vf->reference)
    { OneReference *r = vf->reference ;
      for (n = vf->info['<']->accum.count ; n-- ; r++) 
	free (r->filename) ;
      free (vf->reference) ;
    }
  if (vf->deferred)
    { OneReference *r = vf->deferred ;
      for (n = vf->info['>']->accum.count ; n-- ; r++) 
	free (r->filename) ;
      free (vf->deferred) ;
    }
}

static void oneFileDestroy (OneFile *vf)
{ int       i, j;
  OneInfo *li, *lx;

  if (vf->share)
    { for (i = 0; i < 128 ; i++)
        { lx = vf->info[i];
          if ((i == '&' || i == '*') && ! vf->isWrite)
            continue;
          if (lx != NULL)
            { for (j = 1; j < vf->share; j++)
                { li = vf[j].info[i];
                  if (li->fieldCodec == lx->fieldCodec)
                    li->fieldCodec = NULL;
                  if (li->listCodec == lx->listCodec)
                    li->listCodec  = NULL;
                  infoDestroy(li);
                }
            }
        }

      for (j = 1; j < vf->share; j++)
        { provRefDefCleanup (&vf[j]) ;
          if (vf[j].codecBuf   != NULL) free (vf[j].codecBuf);
          if (vf[j].f          != NULL) fclose (vf[j].f);
        }
    }

  provRefDefCleanup (vf) ;
  if (vf->codecBuf != NULL) free (vf->codecBuf);
  if (vf->f != NULL && vf->f != stdout) fclose (vf->f);

  for (i = 0; i < 128 ; i++)
    if (vf->info[i] != NULL)
      infoDestroy (vf->info[i]);

  free (vf->field) ;

  free (vf) ;
}

/***********************************************************************************
 *
 *    ASCII PARSING UTILITIES: error reporting, lexical level
 *
 **********************************************************************************/

void parseError (OneFile *vf, char *format, ...)
{ va_list args;

  fprintf (stderr, "ONE PARSE ERROR ");

  va_start (args, format);
  vfprintf (stderr, format, args);
  va_end (args);

  vf->lineBuf[vf->linePos] = '\0';
  fprintf (stderr, ", line %lld: %s\n", vf->line, vf->lineBuf);

  exit (1);
}

static char inline vfGetc(OneFile *vf)
{ char c = getc(vf->f);
  if (vf->linePos < 127)
    vf->lineBuf[vf->linePos++] = c;
  return c;
}

static inline void eatWhite (OneFile *vf)
{ char x = vfGetc(vf);
  if (x == ' ') // 200414: removed option to have tab instead of space
    return;
  parseError (vf, "failed to find expected space separation character");
}

static inline char readChar(OneFile *vf)
{ eatWhite(vf);
  return vfGetc(vf);
}

static inline char *readBuf(OneFile *vf)
{ char x, *cp, *endBuf;

  eatWhite (vf);
  endBuf = vf->numberBuf + 32;
  for (cp = vf->numberBuf; cp < endBuf ; cp++)
    { x = vfGetc(vf);
      if (isspace(x) || x == '\0' || x == EOF)
        break;
      *cp = x;
    }
  if (cp >= endBuf)
    { cp[-1] = 0;
      parseError (vf, "overlong item %s", vf->numberBuf);
    }
  else
    { ungetc (x, vf->f);
      vf->linePos -= 1;
      *cp = 0;
    }
  return vf->numberBuf;
}

static inline I64 readInt(OneFile *vf)
{ char *e, *b;
  I64   x;

  b = readBuf(vf);
  x = strtoll(b, &e, 10);
  if (e == b)
    parseError (vf, "empty int field");
  if (*e != '\0')
    parseError (vf, "bad int");
  return x;
}

static inline double readReal(OneFile *vf)
{ char  *e, *b;
  double x;

  b = readBuf(vf);
  x = strtod (b, &e);
  if (e == b)
    parseError (vf, "empty real field");
  if (*e != '\0')
    parseError (vf, "bad real");
  return (x);
}

static inline void readString(OneFile *vf, char *buf, I64 n)
{ eatWhite (vf);
  if (vf->isCheckString)
    { char *cp = buf;
      --cp;
      while (n-- && (*++cp = vfGetc (vf)))
        if (*cp == '\n' || *cp == EOF)
      break;
      if (++n)
        parseError (vf, "line too short %d", buf);
      *++cp = 0;
    }
  else if ((I64) fread (buf, 1, n, vf->f) != n)
    die ("ONE parse error: failed to read %d byte string", n);
}

static inline void readFlush (OneFile *vf) // reads to the end of the line and stores as comment
{ char       x;
  int        n = 0;
  OneInfo *li = vf->info['/'] ;

  // check the first character - if it is newline then done
  x = getc (vf->f) ; 
  if (x == '\n')
    return ;
  else if (x != ' ')
    parseError (vf, "comment not separated by a space") ;

  // else the remainder of the line is a comment
  if (!li->bufSize)
    { li->bufSize = 1024 ;
      li->buffer = new (1024, char) ;
    }
  while ((x = getc (vf->f)) && x != '\n')
    if (x == EOF)
      parseError (vf, "premature end of file");
    else
      { if ((n+1) >= li->bufSize)
	  { char *s = new (2*li->bufSize, char) ;
	    memcpy (s, li->buffer, li->bufSize) ;
	    free (li->buffer) ;
	    li->buffer = s ;
	    li->bufSize *= 2 ;
	  }
	((char*)li->buffer)[n] = x ;
	++n ;
      }
  ((char*)li->buffer)[n] = 0 ; // string terminator
}


/***********************************************************************************
 *
 *    LIST BUFFER & COUNT MANAGEMENT: error reporting, lexical level
 *
 **********************************************************************************/

  //  Ensure line type t buffer can handles size+nStrings, and accumulate counts

static inline void updateCountsAndBuffer (OneFile *vf, char t, I64 size, I64 nStrings)
{ OneInfo *li;

  li = vf->info[(int) t];
  li->accum.total += size;
  if (size > li->accum.max)
    li->accum.max = size;
  size += nStrings;             // need to allocate space for terminal 0s
  if ( ! li->isUserBuf && size > li->bufSize)   // expand buffer
    { if (li->buffer != NULL) free (li->buffer);
      li->bufSize = size;
      li->buffer  = new (size*li->listEltSize, void);
    }
}

  //  Called when a new group starts or eof, accumulate group counts since last group start

static inline void updateGroupCount(OneFile *vf, BOOL isGroupLine)
{ int        i;
  OneInfo   *li;
  OneCounts  *ci;

  for (i = 'A'; i <= 'Z' ; i++)
    { li = vf->info[i];
      if (li != NULL)
        { ci = &(li->accum);
          if (vf->inGroup)
            { if (ci->groupCount < ci->count - li->gCount)
                ci->groupCount = ci->count - li->gCount;
              if (ci->groupTotal < ci->total - li->gTotal)
                ci->groupTotal = ci->total - li->gTotal;
            }
          else
            { li->oCount = ci->count;
              li->oTotal = ci->total;
            }
          li->gCount = ci->count;
          li->gTotal = ci->total;
        }
    }
  if (isGroupLine)
    { vf->group  += 1;
      vf->inGroup = TRUE;
    }
}


/***********************************************************************************
 *
 *   BINARY INT LIST COMPACTION & UNCOMPACTION
 *
 **********************************************************************************/

  //  EWM: it did not look like these have been fully debugged.  Logic does not
  //    accommodate negative numbers, other bits of code looked incomplete.  I
  //    have heavily rewritten and these need debug.

static char *compactIntList (OneFile *vf, OneInfo *li, I64 len, char *buf)
{ char *y;
  int   d, k;
  I64   z, i, mask, *ibuf;

  ibuf = (I64 *) buf;

  if (li->isIntListDiff)
    for (i = len-1; i > 0; i--)
      ibuf[i] -= ibuf[i-1];

  mask = 0;                    // find how many top bytes can be skipped
  for (i = 0; i < len; i++)
    if (ibuf[i] >= 0) 
      mask |= ibuf[i];
    else
      mask |= -(ibuf[i]+1);

  k = li->listEltSize;
  mask >>= 7;
  for (d = 1; d < k; d++)
    { if (mask == 0)
        break;
      mask >>= 8;
    }
  z = k - d;   // number of 0 bytes

  if (z == 0)
    return (buf);
  
  if (buf != li->buffer && ! li->isUserBuf && (I64) (li->bufSize*sizeof(I64)) < d*len)
    { if (li->buffer != NULL)
        free (li->buffer);
      li->bufSize = ((d*len) / sizeof(I64)) + 1;
      li->buffer = new (li->bufSize * sizeof(I64), void);
    }

  y = li->buffer;
  if (vf->isBig)     // copy d bytes per I64, ignoring z before or after depending on isBig
    while (len--)
      { buf += z;
        for (k = 0; k < d; k++)
          *y++ = *buf++;
      }
  else
    while (len--)
      { for (k = 0; k < d; k++)
          *y++ = *buf++;
        buf += z;
      }

  // finally record the number of zero bytes in the top bits of the len field
  oneInt(vf,li->listField) |= (z << 56);
  
  return (li->buffer);
}

static void decompactIntList (OneFile *vf, OneInfo *li, I64 len, char *buf)
{ I64   i, *x;
  int   d, z, k;
  char *s, *t;

  z   = (oneInt (vf, li->listField) >> 56);

  if (z > 0)                      // decompacts in place
    { d = li->listEltSize - z;
      s = buf + d*len;
      t = s + z*len; 
      if (vf->isBig)
        do
          { for (k = 0; k < d; k++)
              *--t = *--s;
            if (*s & 0x80)
              for (k = 0; k < z; k++)
                *--t = 0xff;
            else
              for (k = 0; k < z; k++)
                *--t = 0x0;
          }
        while (s > buf);
      else
        do
          { if (s[-1] & 0x80)
              for (k = 0; k < z; k++)
                *--t = 0xff;
            else
              for (k = 0; k < z; k++)
                *--t = 0;
            for (k = 0; k < d; k++)
              *--t = *--s;
          }
        while (s > buf);
    }
  
  if (li->isIntListDiff)
    { x = (I64 *) buf;
      for (i = 1; i < len; i++)
        x[i] += x[i-1];
    }
}


/***********************************************************************************
 *
 *  ONE_READ_LINE:
 *      Reads the next line and returns FALSE at end of file or on error. The line is
 *      parsed according to its linetype and contents accessed by macros that follow.
 *      The top bit of the first character determines whether the line is binary or ascii
 *
 **********************************************************************************/

  //  Read a string list, first into new allocs, then into sized line buffer.
  //    Annoyingly inefficient, but we don't use it very much.

static void readStringList(OneFile *vf, char t, I64 len)
{ int    j;
  I64    totLen, sLen;
  char **string, *buf;

  totLen = 0;
  string = new (len, char *);
  for (j = 0; j < len ; ++j)
    { sLen = readInt (vf);
      totLen += sLen;
      string[j] = new (sLen+1, char);
      readString (vf, string[j], sLen);
    }

  updateCountsAndBuffer (vf, t, totLen, len);

  buf = (char *) vf->info[(int) t]->buffer;
  for (j = 0; j < len ; ++j)
    { strcpy (buf, string[j]);
      buf += strlen(buf) + 1;
      free (string[j]);
    }
  free (string);
}

BOOL oneReadLine (OneFile *vf)
{ BOOL      isAscii;
  U8        x;
  char      t;
  OneInfo *li;

  if (vf->isWrite)
    die ("ONE read error: trying to read a line from a file open for writing");
  if (vf->isFinal)
    die ("ONE read error: cannot read more data after counts are finalized");

  vf->linePos = 0;                 // must come before first vfGetc()
  x = vfGetc (vf);                 // read first char
  if (feof (vf->f) || x == '\n')   // blank line (x=='\n') is end of records marker before footer
    { vf->lineType = 0 ;           // additional marker of end of file
      return (FALSE);
    }

  vf->line += 1;      // otherwise assume this is a good line, and die if not
  if (x & 0x80)
    { isAscii = FALSE;
      t = vf->binaryTypeUnpack[x];
    }
  else
    { isAscii = TRUE;
      t = x;
    }
  vf->lineType = t;

  li = vf->info[(int) t];
  if (li == NULL)
    parseError (vf, "unknown line type %c(%d was %d) line %d", t, t, x, (int)vf->line);
  li->accum.count += 1;
  if (t == vf->objectType)
    vf->object += 1;
  if (t == vf->groupType)
    updateGroupCount (vf, TRUE);

  // fprintf (stderr, "reading line %lld type %c - expecting %d fields\n", vf->line, t, li->nField) ;

  if (isAscii)           // read field by field according to ascii spec
    { int     i, j;
      I64    *ilst, len;
      double *rlst;

      for (i = 0; i < li->nField; i++)
        switch (li->fieldType[i])
        { case vINT:
            vf->field[i].i = readInt (vf);
	    //	    printf ("  field %d int %d\n", i, (int)oneInt(vf,i)) ; 
	    break;
          case vREAL:
            vf->field[i].r = readReal (vf);
            break;
          case vCHAR:
            vf->field[i].c = readChar (vf);
	    //	    printf ("  field %d char %c\n", i, (int)oneChar(vf,i)) ;
            break;
          case vSTRING:
            len = readInt (vf);
            vf->field[i].len = len;
            updateCountsAndBuffer (vf, t, len, 1);
            readString (vf, (char*) li->buffer, len);
            break;
          case vINT_LIST:
            len = readInt (vf);
            vf->field[i].len = len;
            updateCountsAndBuffer (vf, t, len, 0);
            ilst = (I64 *) li->buffer;
            for (j = 0; j < len; ++j)
              ilst[j] = readInt(vf);
            break;
          case vREAL_LIST:
            len = readInt (vf);
            vf->field[i].len = len;
            updateCountsAndBuffer (vf, t, len, 0);
            rlst = (double *) li->buffer;
            for (j = 0; j < len; ++j)
              rlst[j] = readReal (vf);
            break;
          case vSTRING_LIST: // vSTRING_LIST - inefficient for now - also used for binary
            len = readInt (vf);
            vf->field[i].len = len;
	    //	    printf ("  field %d string list len %d\n", i, (int)oneLen(vf)) ;
            readStringList (vf, t, len);
            break;
        }
      readFlush (vf);
    }

  else        // binary - block read fields and list, potentially compressed
    { int nField;
      I64 nBits, listLen, usedBytes, listSize;

      nField = li->nField;
      if (nField > 0)
        { if (x & 0x1)                       // fields are compressed
            { nBits = (U8) getc (vf->f);     // NB only fields under 255 bits are compressed
              if (fread (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1)
                die ("ONE read error: fail to read compressed fields");

              vcDecode (li->fieldCodec, nBits, vf->codecBuf, (char *) vf->field);
            }
          else
            { if (fread (vf->field, sizeof(OneField), nField, vf->f) != (unsigned long) nField)
                die ("ONE read error: fail to read binary fields");
            }
        }

      if (t == vf->groupType)
        { I64 *groupIndex = (I64 *) vf->info['*']->buffer;
          oneInt(vf,0)    = groupIndex[vf->group] - groupIndex[vf->group-1];
        }

      if (li->listEltSize > 0) // there is a list
        { listLen = oneLen(vf);
          li->accum.total += listLen;
          if (listLen > li->accum.max)
            li->accum.max = listLen;

          if (listLen > 0)
            { if (li->fieldType[li->listField] == vSTRING_LIST) // handle as ASCII
                readStringList (vf, t, listLen);

              else if (x & 0x2)     // list is compressed
                { if (fread (&nBits, sizeof(I64), 1, vf->f) != 1)
                    die ("ONE read error: fail to read list nBits");
                  if (fread (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1)
                    die ("ONE read error: fail to read compressed list");
                  vcDecode (li->listCodec, nBits, vf->codecBuf, li->buffer);
                }

              else
                { usedBytes = li->listEltSize - (oneInt(vf,li->listField) >> 56); // data bytes
                  listSize  = listLen * usedBytes;

                  if ((I64) fread (li->buffer, 1, listSize, vf->f) != listSize)
                    die ("ONE read error: list read %lld not %lld", x, listSize);
                }

              if (li->fieldType[li->listField] == vINT_LIST)
                decompactIntList (vf, li, listLen, li->buffer);
            }

          if (li->fieldType[li->listField] == vSTRING)
            ((char *) li->buffer)[listLen] = '\0'; // 0 terminate
        }

      { U8 peek = getc(vf->f) ; // check if next line is a comment - if so then read it
	ungetc(peek, vf->f) ;
	if (peek & 0x80)
	  peek = vf->binaryTypeUnpack[peek];
	if (peek == '/') // a comment
	  { OneField keepField0 = vf->field[0] ;
	    oneReadLine (vf) ; // read comment line into vf->info['/']->buffer
	    vf->lineType = t ;
	    vf->field[0] = keepField0 ;
	  }
      }
    }

  return t;
}

char *oneReadComment (OneFile *vf)
{ char *comment = (char*)(vf->info['/']->buffer) ;

  if (comment && *comment != 0)
    return comment ;
  else
    return 0 ;
}

/***********************************************************************************
 *
 *   ONE_FILE_OPEN_READ:
 *     Opens file for reading and reads all lines of header if it exists.
 *     If there is no header then fileType must be given (i.e. non-zero),
 *     otherwise if fileType is non-zero then it must match the type in the header.
 *     If the file begins with a $-line then it is binary and the routine reads
 *     the footer that includes decompressors and object index.
 *
 **********************************************************************************/

OneFile *oneFileOpenRead (const char *path, OneSchema *vs, char *fileType, int nthreads)
{
  OneFile *vf;
  FILE    *f;
  int      curLine = 0;
  int      i, j;
  off_t    startOff, footOff;
  U8       c ;

  if (strcmp (path, "-") == 0)
    f = stdin;
  else
    { f = fopen (path, "r");
      if (f == NULL)
        return NULL;
    }
    
#define OPEN_ERROR1(x) { fprintf (stderr,"ONE file error: %s\n",x); fclose(f) ; return NULL; }
#define OPEN_ERROR3(x,y,z) { fprintf (stderr,"ONE file error: ") ; \
      fprintf (stderr,x,y,z) ; fprintf (stderr, "\n") ; fclose(f) ; return NULL ; }
    
  c = getc(f);
  if (feof(f))
    OPEN_ERROR1("file is empty") ;

  if (c == '1')
    { int  major, minor, slen;
      char name[4];
      
      if (fscanf (f, " %d", &slen) != 1)
	OPEN_ERROR1("line 1: failed to read type name length") ;
      if (slen != 3)
	OPEN_ERROR1("line 1: type name is not three letters") ;
      if (fscanf (f, " %s %d %d", name, &major, &minor) != 3)
	OPEN_ERROR1("line 1: failed to read remainder of line") ;
      while (getc(f) != '\n')
	if (feof(f))
	  OPEN_ERROR1("end of file before end of line 1") ;
      ++curLine ;

      vf = oneFileCreate (vs, name) ;
      if (!vf)
	OPEN_ERROR3("unknown primary file type %s in header line %d", name, 1);

      if (major != vf->major)
	{ oneFileDestroy (vf) ;
	  OPEN_ERROR3("major version file %d != code %d", major, (int) vf->major) ;
	}
      if (minor > vf->minor)
	{ oneFileDestroy (vf) ;
	  OPEN_ERROR3("minor version file %d > code %d", minor, (int) vf->minor) ;
	}
      if (fileType && strcmp (fileType, vf->fileType))
	{ oneFileDestroy (vf) ;
	  OPEN_ERROR3("primary fileType mismatch file %s != %s", fileType, vf->fileType) ;
	}
    }
  else
    { if (fileType == NULL)
	OPEN_ERROR1("fileType not defined in file or code") ;
      ungetc(c,f);
      vf = oneFileCreate (vs, fileType) ;
    }

  vf->f = f;
  vf->line = curLine;

  if (nthreads > 1)
    { OneFile *v, *vf0 = vf ;

      if (strcmp (path, "-") == 0)
        die ("ONE error: parallel input incompatible with stdin as input");

      vf->share = nthreads ;
      vf = new (nthreads, OneFile);
      vf[0] = *vf0 ;
      free (vf0) ; // NB free() not oneFileDestroy because don't want deep destroy

      for (i = 1; i < nthreads; i++)
        { v = oneFileCreate(vs, vf->fileType);
	  v->share = -i;
          v->f = fopen (path, "r");

          vf[i] = *v;
          free (v);
        }
    }

  // read header and (optionally) footer
  // recognise end of header by peeking at the first char to check if alphabetic 
 
  f = vf->f;
  vf->isCheckString = TRUE;   // always check strings while reading header
  while (TRUE)
    { U8 peek = getc(f);

      if (feof(f))       // loop exit at end of file
        break;
      ungetc(peek, f);

      if (peek & 0x80)
        peek = vf->binaryTypeUnpack[peek];

      if (isalpha(peek))
        break;    // loop exit at standard data line

      else if (peek == '!')  // hack to insert a count of 4 for vSTRING_LIST
        { getc(f);
          ungetc('4',f);
          ungetc(' ',f);
          ungetc('!',f);
        }
      
      vf->f = f;
      oneReadLine(vf);  // can't fail because we checked file eof already
      f = vf->f;

      switch (vf->lineType)
	{
	case '1':
          parseError(vf, "1 should be first line in header");
          break;

        case '2':
          { char   *sub = oneString(vf) ;
	    OneSchema *os = vs ;
	    while ((os = os->nxt))
	      { for (j = 0 ; j < os->nSecondary ; ++j)
		  if (!strcmp (sub, os->secondary[j])) break ;
		if (j < os->nSecondary) break ;
	      }
            if (!os)
              parseError (vf, "unknown secondary subType %s", sub);
            if (strcmp (os->primary, vf->fileType))
              parseError (vf, "subtype %s not compatible with primary type %d",
			  sub, vf->fileType);
            strcpy (vf->subType, sub) ;
            break;
          }

        case '#':
        case '@':
        case '+':
        case '%':
          { char      c = oneChar(vf,0);
            OneInfo *li = vf->info[(int) c];

            if (li == NULL)
              parseError (vf, "unknown line type %c", c);
            switch (vf->lineType)
            { case '#':
                li->given.count = oneInt(vf,1);
                if (c == vf->objectType && vf->isBinary) // allocate space for object index
                  { vf->info['&']->bufSize = li->given.count;
                    vf->info['&']->buffer  = new (li->given.count, I64);
                  }
                if (c == vf->groupType && vf->isBinary) // allocate space for group index
                  { vf->info['*']->bufSize = li->given.count+1; // +1 for end value
                    vf->info['*']->buffer  = new (vf->info['*']->bufSize, I64);
                  }
                break;
              case '@':
                li->given.max = oneInt(vf,1);
                li->bufSize = li->given.max + 1; // allow for string terminators
                li->buffer = new (li->bufSize*li->listEltSize, void);
                break;
              case '+':
                li->given.total = oneInt(vf,1);
                break;
              case '%':
                c  = oneChar(vf,2);
                li = vf->info[(int) c];
                if (li == NULL)
                  parseError (vf, "unknown line type %c", c);
                c = oneChar(vf,1);
                if (c == '#')
                  li->given.groupCount = oneInt(vf,3);
                else if (c == '+')
                  li->given.groupTotal = oneInt(vf,3);
                else
                  parseError (vf, "unrecognised symbol %c", c);
                break;
            }
            break;
          }

        case '!':     // NB need to copy the strings
          { char *prog     = oneString(vf);
            char *version  = prog + strlen(prog) + 1;
            char *command  = version + strlen(version) + 1;
            char *date = command + strlen(command) + 1;

            vf->info['!']->accum.count -= 1; // to avoid double counting
            oneAddProvenance (vf, prog, version, command, date);
            break;
          }

        case '<':
          vf->info['<']->accum.count -= 1; // to avoid double counting
          oneAddReference (vf, oneString(vf), oneInt(vf,1));
          break;

        case '>':
          vf->info['>']->accum.count -= 1; // to avoid double counting
          oneAddDeferred (vf, oneString(vf));
          break;

        // Below here are binary file header types - requires given.count/given.max first

        case '$':  // read footer - goto end, find offset to start of footer and go there
          if (oneInt(vf,0) != vf->isBig)
            die ("ONE file error: endian mismatch - convert file to ascii");
          vf->isBinary = TRUE;

          startOff = ftello (f);
          if (fseek (f, -sizeof(off_t), SEEK_END) != 0)
            die ("ONE file error: can't seek to final line");

          if (fread (&footOff, sizeof(off_t), 1, f) != 1)
            die ("ONE file error: can't read footer offset");

          if (fseeko (f, footOff, SEEK_SET) != 0)
            die ("ONE file error: can't seek to start of footer");
          break;

        case '^':    // end of footer - return to where we jumped from header
          if (fseeko (f, startOff, SEEK_SET) != 0)
            die ("ONE file error: can't seek back");
          break;

        case '&':
          vf->isIndexIn = TRUE;
          break;

        case '*':
          break;

        case 1:
          vf->info[(int) oneChar(vf,0)]->fieldCodec = vcDeserialize (oneString(vf));
          break;

        case 2:
          vf->info[(int) oneChar(vf,0)]->listCodec = vcDeserialize (oneString(vf));
          break;

        default:
          parseError (vf, "unknown header line type %c", vf->lineType);
          break;
      }
    }
  vf->isCheckString = FALSE;   // user can set this back to TRUE if they wish

  // allocate codec buffer - always allocate enough to handle fields of all line types

  { I64      size;
    OneInfo *li, *l0;
    OneFile *v;

    size = vf->nFieldMax * sizeof(OneField) ;
    for (i = 0; i < 128; ++i)
      if (vf->info[i])
	{ li = vf->info[i];
	  if (li->listCodec && size < li->given.max * li->listEltSize)
	    size = li->given.max * li->listEltSize;
	}
    vf->codecBuf     = new (size+1, void);  // add one for worst case codec usage
    vf->codecBufSize = size+1;

    startOff = ftello (f);
    for (i = 1; i < nthreads; i++)
      { v = vf+i;

        infoDestroy (v->info['&']);
        infoDestroy (v->info['*']);
        v->info['&'] = NULL;
        v->info['*'] = NULL;

        for (j = 0; j < 128; j++)
          { li = v->info[j];
            if (li != NULL)
              { l0 = vf->info[j];
                li->fieldCodec = l0->fieldCodec;
                li->listCodec  = l0->listCodec;
                if (li->listEltSize > 0)
                  { li->bufSize = l0->bufSize;
                    li->buffer  = new (l0->bufSize*l0->listEltSize, void);
                  }
                li->given = l0->given;
              }
          }

        v->codecBuf = new (size+1, void);
        v->codecBufSize = size+1;
        if (fseeko (v->f, startOff, SEEK_SET) != 0)
          die ("ONE file error: can't seek to start of data");

        v->info['&'] = vf->info['&'];
        v->info['*'] = vf->info['*'];

        v->isIndexIn = vf->isIndexIn;
        strcpy (v->subType, vf->subType) ;
      }
  }
  
  vf->f = f;
  return (vf);
}

/***********************************************************************************
 *
 *   ONE_USER_BUFFER / CLOSE / GOTO
 *
 **********************************************************************************/

  // This lets the user reassign the buffer that lists in a particular line type are read into.
  //   If this is not set, a default buffer is provided.  If buffer == NULL then the package
  //   reverts to the default buffer.  This routine can be called repeatedly.
  // NB the package doesn't check the size of a user supplied buffer - the user must allocate
  //   enough memory for all forthcoming list data.

void oneUserBuffer (OneFile *vf, char lineType, void *buffer)
{ OneInfo *li;

  li = vf->info[(int) lineType];
  if (buffer != NULL)
    { if ( ! li->isUserBuf && li->buffer != NULL)
        { free (li->buffer);
          li->bufSize = 0;
        }
      li->buffer    = buffer;
      li->isUserBuf = TRUE;
    }
  else
    { if (li->isUserBuf)
        { li->bufSize = li->given.max + 1;
          li->buffer  = new (li->given.max*li->listEltSize, void);
        }
      li->isUserBuf = FALSE;
    }
}

BOOL oneGotoObject (OneFile *vf, I64 i)
{ if (vf != NULL && vf->isIndexIn && vf->objectType)
    if (0 <= i && i < vf->info[(int) vf->objectType]->given.count)
      if (fseek (vf->f, ((I64 *) vf->info['&']->buffer)[i], SEEK_SET) == 0)
        { vf->object = i;
          return (TRUE);
        }
  return (FALSE);
}

I64 oneGotoGroup (OneFile *vf, I64 i)
{ if (vf != NULL && vf->isIndexIn && vf->groupType)
    if (0 <= i && i < vf->info[(int) vf->groupType]->given.count)
      { I64 *groupIndex = (I64 *) vf->info['*']->buffer;
        if (!oneGotoObject(vf,groupIndex[i]))
	  return (0);
        return (groupIndex[i+1] - groupIndex[i]);
      }
  return (0);
}

/***********************************************************************************
 *
 *   ONE_OPEN_WRITE_(NEW | FROM)
 *
 **********************************************************************************/

OneFile *oneFileOpenWriteNew (const char *path, OneSchema *vs, char *fileType,
                              BOOL isBinary, int nthreads)
{ OneFile *vf;
  FILE    *f;
  
  if (strcmp (path, "-") == 0)
    f = stdout;
  else
    { f = fopen (path, "w");
      if (f == NULL)
        return NULL ;
    }

  vf = oneFileCreate (vs, fileType);
  if (!vf)
    return NULL ;
  
  vf->f = f;
  vf->isWrite  = TRUE;
  vf->isBinary = isBinary;
  vf->isLastLineBinary = TRUE; // we don't want to add a newline before the first true line
  
  vf->codecBufSize = vf->nFieldMax*sizeof(OneField) + 1;
  vf->codecBuf     = new (vf->codecBufSize, void); 

  if (nthreads > 1)
    { OneFile *v, *vf0 = vf ;
      int      i ;
      char     name[100] ;
      int      pid = getpid() ;

      vf->share = nthreads ;
      vf->fieldLock = mutexInit;
      vf->listLock  = mutexInit;
      vf = new (nthreads, OneFile);
      vf[0] = *vf0 ;
      free (vf0) ; // NB free() not oneFileDestroy because don't want deep destroy
      
      for (i = 1; i < nthreads; i++)
	{ v = oneFileCreate (vs, fileType);

	  v->isWrite  = TRUE;
	  v->isBinary = isBinary;
          v->isLastLineBinary = isBinary;
	  
	  v->codecBufSize = vf->nFieldMax*sizeof(OneField) + 1;
	  v->codecBuf     = new (v->codecBufSize, void);
	  v->codecTrainingSize /= 3*nthreads;
          v->share = -i;

          sprintf(name,".part.%d.%d",pid,i);
          f = fopen (name, "w");
          if (f == NULL)
            die ("ONE file error: cannot create temporary file %d for parallel write", i);
	  v->f = f;

	  vf[i] = *v;
	  free (v);
	}
    }

  return vf;
}

OneFile *oneFileOpenWriteFrom (const char *path, OneSchema *vs, OneFile *vfIn,
			       BOOL useAccum, BOOL isBinary, int nthreads)
{ OneFile  *vf;
  OneInfo *li;
  int       i;
  I64       size, sz;

  vf = oneFileOpenWriteNew (path, vs, *vfIn->subType ? vfIn->subType : vfIn->fileType,
			    isBinary, nthreads);

  oneInheritProvenance (vf, vfIn);
  oneInheritReference  (vf, vfIn);
  oneInheritDeferred   (vf, vfIn);

  if (useAccum)
    { for (i = 0; i < 128 ; ++i)
        if (vfIn->info[i])
          vf->info[i]->given = vfIn->info[i]->accum;
    }
  else
    { for (i = 0; i < 128 ; ++i)
        if (vfIn->info[i])
          vf->info[i]->given = vfIn->info[i]->given;
    }

  // allocate codec buffer - always allocate enough to handle fields of all line types

  size = vf->codecBufSize;
  for (i = 0; i < 128 ; ++i)
    { li = vf->info[i];
      if (li != NULL && li->listCodec != NULL)
        { sz = li->given.max * li->listEltSize;
          if (sz >= size)
            size = sz+1;
        }
    }
  if (size > vf->codecBufSize)
    { free (vf->codecBuf);
      vf->codecBufSize = size;
      vf->codecBuf     = new (size, void);
    }
  
  return (vf);
}


/***********************************************************************************
 *
 *    SETTING UP PROVENANCE, REFERENCES, & DEFERRALS
 *
 **********************************************************************************/

static BOOL addProvenance(OneFile *vf, OneProvenance *from, int n)
{ I64 i ;
  OneInfo   *l = vf->info['!'];
  I64         o = l->accum.count;
  OneProvenance *p;

  if (n == 0)
    return (FALSE);
  if (vf->isHeaderOut)
    die("ONE error: can't addProvenance after writing header");

  l->accum.count += n;

  p = new(o+n, OneProvenance);
  if (o > 0)
    memcpy (p, vf->provenance, o*sizeof(OneProvenance));
  memcpy (p+o, from, n*sizeof(OneProvenance));
  free (vf->provenance);
  vf->provenance = p;

  // finally create self-owned copy of all fields

  p = p+o ;
  for (i = 0 ; i < n ; ++i, ++p)
    { p->program = strdup(p->program) ;
      p->version = strdup(p->version) ;
      p->command = strdup(p->command) ;
      p->date = strdup(p->date) ;
    }

  return (TRUE);
}

BOOL oneInheritProvenance(OneFile *vf, OneFile *source)
{ return (addProvenance(vf, source->provenance, source->info['!']->accum.count)); }

BOOL oneAddProvenance(OneFile *vf, char *prog, char *version, char *command, char *date)
{ OneProvenance p;

  p.program = prog;
  p.version = version;
  p.command = command;
  if (date != NULL)
    p.date = date;
  else
    { time_t t = time(NULL);
      p.date = new (20, char);
      strftime(p.date, 20, "%F_%T", localtime(&t));
    }
  addProvenance (vf, &p, 1);
  if (date == NULL)
    free (p.date) ;
  return TRUE ; // always added something
}

static BOOL addReference(OneFile *vf, OneReference *from, int n, BOOL isDeferred)
{ I64        o;
  OneInfo  *l;
  OneReference *r, **t;
  I64 i ;

  if (n == 0)
    return FALSE;
  if (vf->isHeaderOut)
    die ("ONE error: can't addReference after writing header");

  if (isDeferred)
    { l = vf->info['>'];
      t = &(vf->deferred);
    }
  else
    { l = vf->info['<'];
      t = &(vf->reference);
    }
  o = l->accum.count;
  l->accum.count += n;

  r = new (o+n, OneReference);
  if (o > 0)
    memcpy (r, *t, o*sizeof(OneReference));
  memcpy (r+o, from, n*sizeof(OneReference));
  free (*t);
  *t = r;

  r += o ; // make self-owned copy of filename strings
  for (i = 0 ; i < n ; ++i, ++r)
    r->filename = strdup (r->filename) ;

  return TRUE;
}

BOOL oneInheritReference(OneFile *vf, OneFile *source)
{ return (addReference(vf, source->reference, source->info['<']->accum.count, FALSE)); }

BOOL oneAddReference(OneFile *vf, char *filename, I64 count)
{ OneReference ref;
  ref.filename = filename;
  ref.count    = count;
  return (addReference(vf, &ref, 1, FALSE));
}

BOOL oneInheritDeferred (OneFile *vf, OneFile *source)
{ return (addReference (vf, source->deferred, source->info['>']->accum.count, TRUE)); }

BOOL oneAddDeferred (OneFile *vf, char *filename)
{ OneReference ref;
  ref.filename = filename;
  return (addReference (vf, &ref, 1, TRUE));
}

/***********************************************************************************
 *
 *   ONE_WRITE_HEADER / FOOTER
 *
 **********************************************************************************/

void oneWriteHeader (OneFile *vf)
{ int         i,n;
  OneReference  *r;
  OneProvenance *p;
  OneInfo   *li;

  if ( ! vf->isWrite)
    die ("ONE error: trying to write header to a file open for reading");
  if (vf->line > 0)
    die ("ONE error: cannot write header after writing one or more data lines");
  if (vf->info[(int) vf->objectType]->given.count == 0 && ! vf->isBinary)
    die ("ONE error: information for ASCII header is not present, use oneFileOpenWriteFrom");

  fprintf (vf->f, "1 %lu %s %lld %lld", strlen(vf->fileType), vf->fileType, vf->major, vf->minor);
  vf->line += 1;

  if (*vf->subType)
    { fprintf (vf->f, "\n2 %lu %s", strlen(vf->subType), vf->subType);
      vf->line += 1;
    }

  r = vf->reference;
  n = vf->info['<']->accum.count;
  for (i = 0; i < n; i++, r++)
    { fprintf (vf->f, "\n< %lu %s %lld", strlen(r->filename), r->filename, r->count);
      vf->line += 1;
    }

  r = vf->deferred;
  n = vf->info['>']->accum.count;
  for (i = 0; i < n; i++, r++)
    { fprintf (vf->f, "\n> %lu %s", strlen(r->filename), r->filename);
      vf->line += 1;
    }
  
  p = vf->provenance; 
  n = vf->info['!']->accum.count;
  for (i = 0; i < n; i++, p++)
    { fprintf (vf->f, "\n! %lu %s %lu %s %lu %s %lu %s",
                     strlen (p->program), p->program, strlen (p->version), p->version,
                     strlen (p->command), p->command, strlen (p->date), p->date);
      vf->line += 1;
    }

  if (vf->isBinary)         // defer writing rest of header
    { fprintf (vf->f, "\n$ %d", vf->isBig);
      vf->line += 1;
    }
  else             // write counts based on those supplied in input header
    { for (i = 'A'; i <= 'Z'+1 ; i++)
	{ if (i == 'Z'+1)
	    { if (vf->groupType) // NB group types are all lower case so > 'Z'+1
		i = vf->groupType ;
	      else
		break ;
	    }
	  li = vf->info[i];
          if (li != NULL && li->given.count > 0)
            { fprintf (vf->f, "\n# %c %lld", i, li->given.count);
              vf->line += 1;
              if (li->given.max > 0)
                { fprintf (vf->f, "\n@ %c %lld", i, li->given.max);
                  vf->line += 1;
                }
              if (li->given.total > 0)
                { fprintf (vf->f, "\n+ %c %lld", i, li->given.total);
                  vf->line += 1;
                }
              if (li->given.groupCount > 0)
                { fprintf (vf->f, "\n%% %c # %c %lld", vf->groupType, i, li->given.groupCount);
                  vf->line += 1;
                }
              if (li->given.groupTotal > 0)
                { fprintf (vf->f, "\n%% %c + %c %lld", vf->groupType, i, li->given.groupTotal);
                  vf->line += 1;
                }
            }
        }
    }
  fflush (vf->f);

  vf->isLastLineBinary = FALSE;
  vf->isHeaderOut = TRUE;
}


/***********************************************************************************
 *
 *   ONE_WRITE_LINE
 *
 **********************************************************************************/

static int writeStringList (OneFile *vf, char t, int len)
{ OneInfo *li;
  int       j, nByteWritten = 0;
  I64       sLen, totLen;
  char     *buf;

  buf = (char *) vf->info[(int) t]->buffer;

  totLen = 0;
  for (j = 0; j < len; j++)
    { sLen = strlen (buf);
      totLen += sLen;
      nByteWritten += fprintf (vf->f, " %lld %s", sLen, buf);
      buf += sLen + 1;
    }

  li = vf->info[(int) t];
  li->accum.total += totLen;
  if (li->accum.max < totLen)
    li->accum.max = totLen;

  return nByteWritten ;
}

// process is to fill fields by assigning to macros, then call - list contents are in buf
// NB adds '\n' before writing line not after, so user fprintf() can add extra material
// first call will write initial header, allowing space for count sizes to expand on close

void oneWriteLine (OneFile *vf, char t, I64 listLen, void *listBuf)
{ I64       i, j;
  OneInfo *li;

  //  fprintf (stderr, "write line type %c\n", t) ;
  
  if ( ! vf->isWrite)
    die ("ONE write error: trying to write a line to a file open for reading");
  if (vf->isFinal && isalpha(t))
    die ("ONE write error: annot write more data after counts are finalized %c", t);

  li = vf->info[(int) t];
  if (li == NULL)
    die ("ONE write error: line type %c not present in file spec %s ", t, vf->fileType);

  if (listBuf == NULL)
    listBuf = li->buffer;

  if ( ! vf->isLastLineBinary)      // terminate previous ascii line
    fputc ('\n', vf->f);

  vf->line  += 1;
  li->accum.count += 1;
  if (t == vf->groupType)
    updateGroupCount(vf, TRUE);

  if (li->listEltSize > 0)  // need to write the list
    { if (listLen >= 0)
	vf->field[li->listField].len = listLen ;
      else
	die ("ONE write error: listLen %lld must be non-negative", listLen) ;
    }

  // BINARY - block write and optionally compress
  
  if (vf->isBinary)
    { U8  x, cBits;
      int nField;
      I64 fieldSize, nBits, listBytes, listSize;

      if (!vf->isLastLineBinary)
	vf->byte = ftello (vf->f) ;

      if (t == vf->objectType) // update index and increment object count
        { OneInfo *lx = vf->info['&'];

          if (vf->object >= lx->bufSize) // first ensure enough space
            { I64  ns = (lx->bufSize << 1) + 0x20000;
              I64 *nb = new (ns, I64);
	      
              memcpy(nb, lx->buffer, lx->bufSize*sizeof(I64));
              free (lx->buffer);
              lx->buffer  = nb;
              lx->bufSize = ns;
            }
          ((I64 *) lx->buffer)[vf->object] = vf->byte;
#define CHECK_INDEX
#ifdef CHECK_INDEX
          { if (ftello (vf->f) != ((I64 *) lx->buffer)[vf->object])
	      die ("ONE write error: byte offset index error") ;
	  }
#endif
          vf->object += 1;
        }
      if (t == vf->groupType)
        { OneInfo *lx = vf->info['*'];
	  
          if (vf->group >= lx->bufSize) // still room for final value because one ahead here
            { I64  ns, *nb;

              ns = (lx->bufSize << 1) + 0x20000;
              nb = new (ns, I64);
              memcpy(nb, lx->buffer, lx->bufSize*sizeof(I64));
              free (lx->buffer);
              lx->buffer  = nb;
              lx->bufSize = ns;
            }
       
          ((I64 *) lx->buffer)[vf->group-1] = vf->object; // group # already advanced
        }

      nField    = li->nField;
      fieldSize = nField*sizeof(OneField);

      if (li->listEltSize > 0 && li->fieldType[li->listField] == vINT_LIST)
	listBuf = compactIntList (vf, li, listLen, listBuf);

      x = li->binaryTypePack;   //  Binary line code + compression flags
      if (li->isUseListCodec)
        x |= 0x02;

      if (li->isUseFieldCodec)
        { nBits = vcEncode (li->fieldCodec, fieldSize, (char *) vf->field, vf->codecBuf);
          if (nBits < 256)
            { x |= 0x01;
              cBits = nBits;

              fputc (x, vf->f);
              fputc (cBits, vf->f) ;
              if (fwrite (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1)
                die("ONE write error: fail to write compressed fields");
	      vf->byte += 2 + ((nBits+7) >> 3) ;
            }
          else
            { fputc (x, vf->f);
              if (nField > 0)
                { if (fwrite (vf->field, fieldSize, 1, vf->f) != 1)
                    die ("ONE write error: write fields: t %c, nField %d, fieldSize %lld",
			 t, nField, fieldSize);
                }
	      vf->byte += 1 + fieldSize ;
            }
        }
      else
        { fputc (x, vf->f);
          if (nField > 0)
            if (fwrite (vf->field, fieldSize, 1, vf->f) != 1)
              die ("ONE write error: write fields: t %c, nField %d, fieldSize %lld",
		   t, nField, fieldSize);
	  vf->byte += 1 + fieldSize ;

          if (li->fieldCodec != NULL)
            { vcAddToTable (li->fieldCodec, fieldSize, (char *) vf->field);
              li->fieldTack += fieldSize;

              if (li->fieldTack > vf->codecTrainingSize)
                { if (vf->share == 0)
                    { vcCreateCodec (li->fieldCodec, 1);
                      li->isUseFieldCodec = TRUE;
                    }
                  else
                    { OneFile  *ms;
                      OneInfo *lx;

                      if (vf->share < 0)
                        { ms = vf + vf->share;
                          lx = ms->info[(int) t]; 
                        }
                      else
                        { ms = vf;
                          lx = li;
                        }

                      pthread_mutex_lock(&ms->fieldLock);

                      if ( ! li->isUseFieldCodec)

                        { if (vf->share < 0)
                            { lx->fieldTack += li->fieldTack;
                              li->fieldTack = 0;
                            }
                          if (lx->fieldTack > ms->codecTrainingSize)
                            { for (i = 1; i < ms->share; i++)
                                vcAddHistogram (lx->fieldCodec,
                                                ms[i].info[(int) t]->fieldCodec);
			      vcCreateCodec (lx->fieldCodec, 1);
                              for (i = 1; i < ms->share; i++)
                                { OneCodec *m = ms[i].info[(int) t]->fieldCodec;
                                  ms[i].info[(int) t]->fieldCodec = lx->fieldCodec;
                                  vcDestroy (m);
                                }
                              lx->isUseFieldCodec = TRUE;
                              for (i = 1; i < ms->share; i++)
                                ms[i].info[(int) t]->isUseFieldCodec = TRUE;
                            }
                        }

                      pthread_mutex_unlock(&ms->fieldLock);
                    }
                }
            }
        }

      // Write the list if there is one

      if (li->listEltSize)
        { li->accum.total += listLen;
          if (listLen > li->accum.max)
            li->accum.max = listLen;

          if (listLen > 0)
            { listBytes = li->listEltSize - (oneInt(vf, li->listField) >> 56); // data bytes
              listSize  = listLen * listBytes;

              if (li->fieldType[li->listField] == vSTRING_LIST) // handle as ASCII
                vf->byte += writeStringList (vf, t, listLen);

              else if (x & 0x2)
                { if (listSize >= vf->codecBufSize)
                    { free (vf->codecBuf);
                      vf->codecBufSize = listSize+1;
                      vf->codecBuf     = new (vf->codecBufSize, void);
                    }
                  nBits = vcEncode (li->listCodec, listSize, listBuf, vf->codecBuf);
                  if (fwrite (&nBits, sizeof(I64), 1, vf->f) != 1)
                    die ("ONE write error: failed to write list nBits");
                  if (fwrite (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1)
                    die ("ONE write error: failed to write compressed list");
		  vf->byte += sizeof(I64) + ((nBits+7) >> 3);
                }

              else
                { if (fwrite (listBuf, listSize, 1, vf->f) != 1)
                    die ("ONE write error: failed to write list field %d listLen %lld listSize %lld listBuf %lx",
			 li->listField, listLen, listSize, listBuf);
		  vf->byte += listSize;
                  if (li->listCodec != NULL)
                    { vcAddToTable (li->listCodec, listSize, listBuf);
                      li->listTack += listSize;

                      if (li->listTack > vf->codecTrainingSize)
                        { if (vf->share == 0)
                            { vcCreateCodec (li->listCodec, 1);
                              li->isUseListCodec = TRUE;
                            }
                          else
                            { OneFile  *ms;
                              OneInfo *lx;

                              if (vf->share < 0)
                                { ms = vf + vf->share;
                                  lx = ms->info[(int) t]; 
                                }
                              else
                                { ms = vf;
                                  lx = li;
                                }

                              pthread_mutex_lock(&ms->listLock);

                              if ( ! li->isUseListCodec)

                                { if (vf->share < 0)
                                    { lx->listTack += li->listTack;
                                      li->listTack = 0;
                                    }
                                  if (lx->listTack > ms->codecTrainingSize)
                                    { for (i = 1; i < ms->share; i++)
                                        vcAddHistogram (lx->listCodec,
                                                        ms[i].info[(int) t]->listCodec);
                                      vcCreateCodec (lx->listCodec, 1);
                                      for (i = 1; i < ms->share; i++)
                                        { OneCodec *m = ms[i].info[(int) t]->listCodec;
                                          ms[i].info[(int) t]->listCodec = lx->listCodec;
                                          vcDestroy (m);
                                        }
                                      lx->isUseListCodec = TRUE;
                                      for (i = 1; i < ms->share; i++)
                                        ms[i].info[(int) t]->isUseListCodec = TRUE;
                                    }
                                }
 
                              pthread_mutex_unlock(&ms->listLock);
                            }
                        }
                    }
                }
            }
        }
      vf->isLastLineBinary = TRUE;
    }

  // ASCII - write field by field

  else
    { fputc (t, vf->f);
      for (i = 0; i < li->nField; i++)
        switch (li->fieldType[i])
        { case vINT:
            fprintf (vf->f, " %lld", vf->field[i].i);
            break;
          case vREAL:
            fprintf (vf->f, " %f", vf->field[i].r);
            break;
          case vCHAR:
            fprintf (vf->f, " %c", vf->field[i].c);
            break;
          case vSTRING:
          case vINT_LIST:
          case vREAL_LIST:
          case vSTRING_LIST:
            li->accum.total += listLen;
            if (listLen > li->accum.max)
              li->accum.max = listLen;

            fprintf (vf->f, " %lld", listLen);
            if (li->fieldType[i] == vSTRING)
              { if (listLen > INT_MAX)
                  die ("ONE write error: string length %lld > current max %d", listLen, INT_MAX);
                fprintf (vf->f, " %.*s", (int) listLen, (char *) listBuf);
              }
            else if (li->fieldType[i] == vINT_LIST)
              { I64 *b = (I64 *) listBuf;
                for (j = 0; j < listLen ; ++j)
                  fprintf (vf->f, " %lld", b[j]);
              }
            else if (li->fieldType[i] == vREAL_LIST)
              { double *b = (double *) listBuf;
                for (j = 0; j < listLen ; ++j)
                  fprintf (vf->f, " %f", b[j]);
              }
            else // vSTRING_LIST
              writeStringList (vf, t, listLen);
            break;
        }
      vf->isLastLineBinary = FALSE;
    }
}

void oneWriteComment (OneFile *vf, char *comment)
{
  if (vf->isLastLineBinary)
    oneWriteLine (vf, '/', strlen(comment), comment) ;
  else
    fprintf (vf->f, " %s", comment) ;
}

/***********************************************************************************
 *
 *    MERGING, FOOTER HANDLING, AND CLOSE
 *
 **********************************************************************************/

static void oneWriteFooter (OneFile *vf)
{ int      i,n;
  off_t    footOff;
  OneInfo *li;
  char    *codecBuf ;
  
  footOff = ftello (vf->f);
  if (footOff < 0)
    die ("ONE write error: failed footer ftell");

  //  first the per-linetype information
  codecBuf = new (vcMaxSerialSize()+1, char) ; // +1 for added up unused 0-terminator
  for (i = 'A'; i <= 'Z'+1 ; i++)
    { if (i == 'Z'+1)
	{ if (vf->groupType) // NB group types are all lower case so > 'Z'+1
	    i = vf->groupType ;
	  else
	    break ;
	}
      li = vf->info[i];
      if (li != NULL && li->accum.count > 0)
        { fprintf (vf->f, "# %c %lld\n", i, li->accum.count);
	  if (li->listEltSize)
            { fprintf (vf->f, "@ %c %lld\n", i, li->accum.max);
	      fprintf (vf->f, "+ %c %lld\n", i, li->accum.total);
	    }
	  if (vf->groupType && i != vf->groupType && vf->group > 0)
	    { fprintf (vf->f, "%% %c # %c %lld\n", vf->groupType, i, li->accum.groupCount);
	      if (li->listEltSize)
		fprintf (vf->f, "%% %c + %c %lld\n", vf->groupType, i, li->accum.groupTotal);
	    }
          if (li->isUseFieldCodec)
            { oneChar(vf,0) = i;
              n = vcSerialize (li->fieldCodec, codecBuf);
              oneWriteLine (vf, 1, n, codecBuf);
            }
          if (li->isUseListCodec && li->listCodec != DNAcodec)
            { oneChar(vf,0) = i;
              n = vcSerialize (li->listCodec, codecBuf);
              oneWriteLine (vf, 2, n, codecBuf);
            }
        }
    }
  free (codecBuf) ;

  oneWriteLine (vf, '&', vf->object, NULL); // number of objects in file = length of index
  // NB NULL here and below for '*' defaults writing info->buffer, which contains the index

  if (vf->groupType > 0 && vf->group > 0)
    { ((I64 *) vf->info['*']->buffer)[vf->group] = vf->object;
      oneWriteLine (vf, '*', vf->group+1, NULL); // number of groups in file + 1 = length of index
    }

  fprintf (vf->f, "^\n");

  if (fwrite (&footOff, sizeof(off_t), 1, vf->f) != 1)
    die ("ONE write error: failed writing footer offset");
}

void oneFinalizeCounts(OneFile *vf)
{ int       i, j, n, k, len;
  OneInfo *li, *ln;

  if (vf->share < 0)
    die ("ONE write error: cannot call oneFileClose on a slave OneFile");

  vf->isFinal = TRUE;

  if (vf->share == 0)
    { updateGroupCount(vf,FALSE);
      return;
    }

  len = vf->share;
  
  //  Close current groups at the end of each part (if any)

  if (vf->groupType > 0)
    for (i = 'A'; i <= 'Z'; i++)
      if (vf->info[i] != NULL)
        for (j = 0; j < len; j++)
          if (vf[j].inGroup)
            { I64 oc, ot;

              ot = oc = 0;
              for (k = j+1; k < len; k++)
                if (vf[k].inGroup)
                  { oc += vf[k].info[i]->oCount;
                    ot += vf[k].info[i]->oTotal; 
                    break;
                  }
                else
                  { oc += vf[k].info[i]->accum.count;
                    ot += vf[k].info[i]->accum.total;
                  }

              li = vf[j].info[i];
              if ((li->accum.count - li->gCount) + oc > li->accum.groupCount)
                li->accum.groupCount = (li->accum.count - li->gCount) + oc;
              if ((li->accum.total - li->gTotal) + ot > li->accum.groupTotal)
                li->accum.groupTotal = (li->accum.total - li->gTotal) + ot;
            }

  //  first the per-linetype information

  n = vf->groupType;
  if (n == 0)
    n = 'Z';

  for (i = 'A'; i <= n; i++)
    { ln = vf->info[i];
      for (j = 1; j < len; j++)
        { li = (vf+j)->info[i];
          if (li != NULL && li->accum.count > 0)
            { ln->accum.count += li->accum.count;
              if (li->accum.max > ln->accum.max)
                ln->accum.max = li->accum.max;
              ln->accum.total += li->accum.total;
              if (li->accum.groupCount > ln->accum.groupCount)
                ln->accum.groupCount = li->accum.groupCount;
              if (li->accum.groupTotal > ln->accum.groupTotal)
                ln->accum.groupTotal = li->accum.groupTotal;
            }
        }
    }

  if ( ! vf->isBinary)
    return;

  //  Stitch the group index together

  if (vf->groupType > 0)
    { I64 *gb, *gi, off;
      int  ns;

      ns = 0; 
      for (j = 0; j < len; j++)
        ns += vf[j].group;
      gb = new (ns+1, I64);

      ns = 0;
      off = 0;
      for (j = 0; j < len; j++)
        { li = vf[j].info['*'];
          gi = (I64 *) (li->buffer);
          for (i = 0; i < vf[j].group; i++)
            gb[ns++] = gi[i] + off;
          off += vf[j].object;
        }
      gb[ns] = off;
      li = vf->info['*'];
      free(li->buffer);
      li->buffer  = gb;
      li->bufSize = ns+1;
      vf->group = ns;
    }

  //  Stitch the object index together

  { int  ns;
    I64 *gb, *gi, off;

    ns = 0;
    for (j = 0; j < len; j++)
      ns += vf[j].object;
    gb = new (ns, I64);

    ns = 0;
    off = 0;
    for (j = 0; j < len; j++)
      { li = vf[j].info['&'];
        gi = (I64 *) (li->buffer);
        for (i = 0; i < vf[j].object; i++)
          gb[ns++] = gi[i] + off;
        off += ftello(vf[j].f);
      }

    li = vf->info['&'];
    free(li->buffer);
    li->buffer  = gb;
    li->bufSize = ns;
    vf->object  = ns;
  }
}

// automatically rewrites header if allowed when writing

void oneFileClose (OneFile *vf)
{
  if (vf->share < 0)
    die ("ONE file error: cannot call oneFileClose on a slave OneFile");

  if (vf->isWrite)
    {
      if (!vf->isFinal) // RD moved this here from above - surely only needed if isWrite
	oneFinalizeCounts (vf);
      
      if (vf->share > 0)
        { int  i, pid, fid, nread;
          char name[100], *buf;

          buf = new (10000000, char);
          pid = getpid();
          for (i = 1; i < vf->share; i++)
            { fclose (vf[i].f);
              vf[i].f = NULL;
              sprintf(name,".part.%d.%d",pid,i);
              fid = open(name,O_RDONLY);
              while ((nread = read(fid,buf,10000000)) > 0)
                if ((int) fwrite(buf,1,nread,vf->f) != nread)
                  die ("ONE write error: while cat'ing thread bits (oneFileClose)");
              if (unlink(name) < 0)
                die ("ONE write error: could not delete thread file %s", name);
            }
          free(buf);
        }
      fputc ('\n', vf->f);  // end of file if ascii, end of data marker if binary
      if (vf->isBinary) // write the footer
        oneWriteFooter (vf);
    }
  
  oneFileDestroy (vf);
}

/***********************************************************************************
 *
 *  Length limited Huffman Compressor/decompressor with special 2-bit compressor for DNA
 *  Author:          Gene Myers
 *  Creation date:   June 27, 2019
 *
 *  inline both compression.h and compression.c here
 *
 **********************************************************************************/

#undef  DEBUG
#undef  TEST

  //  To create a compressor, get an initially empty object with vcCreate, then
  //    add a significant corpus of the byte data to be compressed with vcAddToTable,
  //    and finally create a Huffman codec based on this corpus by calling
  //    vcCreateCodec.  The parameter "partial" should be set if not all the data
  //    to be compressed has been scanned.  At this point you have a compressor ready
  //    to operate.  You can destroy/free it with vcDestroy.

OneCodec *vcCreate();
void      vcAddToTable(OneCodec *vc, int len, char *bytes);
void      vcCreateCodec(OneCodec *vc, int partial);
void      vcDestroy(OneCodec *vc);

  //  In the instance of accumulating data over multiple threads, vcAddHistogram, will
  //    add the counts in the table for vh, to the table for vc.

void      vcAddHistogram(OneCodec *vc, OneCodec *vh);

  //  A diagnostic routine: shows you the compression scheme and if the distribution
  //    of the scanned corpus is available, it shows you that too.  Output to file 'to'.

void      vcPrint(OneCodec *vc, FILE *to);

  //  You can encode and decode where ibytes/ilen are the input and the output
  //    is placed at obytes and the length of the compressed/decompressed result
  //    is returned as the value of the function.  For vcEncode, ilen is the size
  //    of the uncompressed input in bytes, and the return value is the size of
  //    the compressed output in **bits**.  The converse is true for vcDecode, i.e
  //    ilen is the number of bits in the compressed input, and the return value
  //    is the number of bytes in the uncompressed output.  The routines are endian safe.

int       vcEncode(OneCodec *vc, int ilen, char *ibytes, char *obytes);
int       vcDecode(OneCodec *vc, int ilen, char *ibytes, char *obytes);

  //  Rather than directly reading or writing an encoding of a compressor, the routines
  //    below serialize or deserialize the compressor into/outof a user-supplied buffer.
  //    vcMaxSerialSize gives the maximum size of a serialized compressor so the user
  //    can arrange a buffer of the appropriate size.  vcSerialize serializes vc into
  //    buffer 'out' and returns the # of bytes in the encoding.  vcDeserialize will reverse
  //    the process given a serialization.  The routines are endian-safe.

int       vcMaxSerialSize();
int       vcSerialize(OneCodec *vc, void *out);
OneCodec *vcDeserialize(void *in);

typedef unsigned long long uint64;
typedef unsigned int       uint32;
typedef unsigned short     uint16;
typedef unsigned char      uint8;

#define HUFF_CUTOFF  12     //  This cannot be larger than 16 !

  //  Endian flipping macros

#define FLIP64(p)	\
{ uint8 x = p[0];	\
  p[0] = p[7];		\
  p[7] = x;		\
  x     = p[1];		\
  p[1] = p[6];		\
  p[6] = x;		\
  x     = p[2];		\
  p[2] = p[5];		\
  p[5] = x;		\
  x     = p[3];		\
  p[3] = p[4];		\
  p[4] = x;		\
}

#define FLIP32(p)	\
{ uint8 x = p[0];	\
  p[0] = p[3];		\
  p[3] = x;		\
  x     = p[1];		\
  p[1] = p[2];		\
  p[2] = x;		\
}

#define FLIP16(p)	\
{ uint8 x = p[0];	\
  p[0] = p[1];		\
  p[1] = x;		\
}

/*******************************************************************************************
 *
 *  Routines for computing a length-limited Huffman Encoding Scheme
 *
 ********************************************************************************************/

#define EMPTY        0      //  Compressor just created, histogram zero'd
#define FILLED       1      //  Compressor histogram being filled, no codec
#define CODED_WITH   2      //  Compressor has a codec (can no longer accumulate histogram)
#define CODED_READ   3      //  Compressor has codec but no histogram as was created by read

typedef struct
  { int    state;            //  1 of the 4 states immediately above
    int    isbig;            //  endian of the current machine
    uint16 codebits[256];    //  Code esc_code is the special code for
    uint8  codelens[256];    //    non-Huffman exceptions
    char   lookup[0x10000];  //  Lookup table (just for decoding)
    int    esc_code;         //  The special escape code (-1 if not partial)
    int    esc_len;          //  The length in bits of the special code (if present)
    uint64 hist[256];        //  Byte distribution for codec
  } _OneCodec;

  //  The special "predefined" DNA compressor

static _OneCodec _DNAcodec = { .state = CODED_READ };
OneCodec  *DNAcodec = (OneCodec *) &_DNAcodec;

  //  Create an EMPTY compressor object with zero'd histogram and determine machine endian

OneCodec *vcCreate()
{ _OneCodec *v;
  int i;

  v = (_OneCodec *) malloc(sizeof(_OneCodec));
  if (v == NULL)
    { fprintf(stderr,"vcCreate: Could not allocate compressor\n");
      exit (1);
    }

  v->state = EMPTY;
  for (i = 0; i < 256; i++)
    v->hist[i] = 0;

  { uint32 t;
    uint8 *b;

    t = 1;
    b = (uint8 *) (&t);
    v->isbig = (b[0] == 0);
  }

  return ((OneCodec *) v);
}

  //  Free a compressor object

void vcDestroy(OneCodec *vc)
{ _OneCodec *v = (_OneCodec *) vc;
  if (vc != DNAcodec)
    free(v);
}

  //  Add the frequencies of bytes in bytes[0..len) to vc's histogram
  //    State becomes FILLED

void vcAddToTable(OneCodec *vc, int len, char *bytes)
{ _OneCodec *v = (_OneCodec *) vc;
  uint8 *data = (uint8 *) bytes;
  int i;

  for (i = 0; i < len; i++)
    v->hist[(int) data[i]] += 1;
  if (v->state < FILLED)
    v->state = FILLED;
}

  //  Add the frequencies of bytes in bytes[0..len) to vc's histogram
  //    State becomes FILLED

void vcAddHistogram(OneCodec *vc, OneCodec *vh)
{ _OneCodec *v = (_OneCodec *) vc;
  _OneCodec *h = (_OneCodec *) vh;
  int i;

  if (v->state >= CODED_WITH)
    { fprintf(stderr,"vcAddHistogram: Compressor already has a codec\n");
      exit (1);
    }
  if (h->state == CODED_READ)
    { fprintf(stderr,"vcAddHistogram: Source compressor doesn't have a histogram\n");
      exit (1);
    }

  for (i = 0; i < 256; i++)
    v->hist[i] += h->hist[i];
  v->state = FILLED;
}

  //  Check vc has a non-empty distribution histogram and if so then build
  //    length-limited Huffman tables for the bytes that occur in the histogram,
  //    plus a special escape code if partial is set and there is at least one byte
  //    with a zero count in the histogram.  The algorithm is by Larmore & Hirschberg,
  //    JACM 73, 3 (1990).

uint64 *HIST;

int HSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);
  return (HIST[x] - HIST[y]);
}

void vcCreateCodec(OneCodec *vc, int partial)
{ _OneCodec *v = (_OneCodec *) vc;

  uint64  *hist;
  char    *look;
  uint8   *lens;
  uint16  *bitv;

  int      code[256];
  int      leng[256];
  uint16   bits[256];
  int      ncode, dcode, ecode;

  int      i;

  if (v->state >= CODED_WITH)
    { fprintf(stderr,"vcCreateCoder: Compressor already has a codec\n");
      exit (1);
    }
  if (v->state == EMPTY)
    { fprintf(stderr,"vcCreateCoder: Compressor has no byte distribution data\n");
      exit (1);
    }

  hist  = v->hist;
  look  = v->lookup;
  lens  = v->codelens;
  bitv  = v->codebits;

  ecode = -partial;
  ncode = 0;
  for (i = 0; i < 256; i++)
    if (hist[i] > 0)
      code[ncode++] = i;
    else if (ecode < 0)
      { ecode = i;
        code[ncode++] = i;
      }
  dcode = 2*ncode;

  if (ecode < 0)
    partial = 0;

  HIST = hist;
  qsort(code,ncode,sizeof(int),HSORT);

#ifdef DEBUG
  fprintf(stderr,"\nSorted Codes %d:\n",ncode);
  for (i = 0; i < ncode; i++)
    fprintf(stderr," %3d: %3d %10llu\n",i,code[i],hist[code[i]]);
#endif

  { uint8   matrix[HUFF_CUTOFF][dcode];
    uint64  count1[dcode], count2[dcode], countb[ncode];
    uint64 *lcnt, *ccnt, *swp;
    int     llen, span;
    int     j, k, n, L;

    for (n = 0; n < ncode; n++)
      { count1[n] = countb[n] = hist[code[n]];
        leng[n] = 0;
      }

#ifdef DEBUG
    fprintf(stderr,"\nCoin Filter:\n");
    fprintf(stderr,"  Row %2d:",HUFF_CUTOFF);
    for (n = 0; n < ncode; n++)
      fprintf(stderr," %lld*",countb[n]);
    fprintf(stderr,"\n");
#endif

    lcnt = count1;
    ccnt = count2;
    llen = ncode-1;
    for (L = HUFF_CUTOFF-1; L > 0; L--)
      { j = 0;
        k = 0;
        for (n = 0; j < ncode || k < llen; n++)
          { if (k >= llen || (j < ncode && countb[j] <= lcnt[k] + lcnt[k+1]))
              { ccnt[n] = countb[j];
                matrix[L][n] = 1;
                j += 1;
              }
            else
              { ccnt[n] = lcnt[k] + lcnt[k+1];
                matrix[L][n] = 0;
                k += 2;
              }
          }
        llen = n-1;
        swp  = lcnt;
        lcnt = ccnt;
        ccnt = swp;

#ifdef DEBUG
        fprintf(stderr,"  Row %2d:",L);
        for (n = 0; n <= llen; n++)
          fprintf(stderr," %lld%c",lcnt[n],matrix[L][n]?'*':'+');
        fprintf(stderr,"\n");
#endif
      }

    span = 2*(ncode-1);
    for (L = 1; L < HUFF_CUTOFF; L++)
      { j = 0;
        for (n = 0; n < span; n++)
          { if (matrix[L][n])
              leng[j++] += 1;
          }
        span = 2*(span-j);
      }
    for (n = 0; n < span; n++)
      leng[n] += 1;

#ifdef DEBUG
    fprintf(stderr,"\nBack Trace:\n");
    span = 2*(ncode-1);
    for (L = 1; L < HUFF_CUTOFF; L++)
      { j = 0;
        fprintf(stderr,"  Row %2d:",L);
        for (n = 0; n < span; n++)
          { if (matrix[L][n])
              j += 1;
            fprintf(stderr," %c",matrix[L][n]?'*':'+');
          }
        fprintf(stderr,"\n");
        span = 2*(span-j);
      }
    fprintf(stderr,"  Length:");
    for (n = 0; n < ncode; n++)
      fprintf(stderr," %d",leng[n]);
    fprintf(stderr,"\n");
#endif
  }

  { int    n, llen;
    uint16 lbits;

    llen  = leng[0];
    lbits = bits[0] = (1 << llen) - 1;
    for (n = 1; n < ncode; n++)
      { while ((lbits & 0x1) == 0)
          { lbits >>= 1;
            llen -= 1;
          }
        lbits -= 1;
        while (llen < leng[n])
          { lbits = (lbits << 1) | 0x1;
            llen += 1;
          }
        bits[n] = lbits;
      }

#ifdef DEBUG
    { int j;

      fprintf(stderr,"\nCodes:\n");
      for (n = 0; n < ncode; n++)
        { fprintf(stderr,"   %3d: %2d ",code[n],leng[n]);
          for (j = leng[n]-1; j >= 0; j--)
            fprintf(stderr,"%x",(bits[n]>>j)&0x1);
          fprintf(stderr,"\n");
        }
    }
#endif
  }

  for (i = 0; i < 256; i++)
    { lens[i] = 0;
      bitv[i] = 0;
    }

  for (i = 0; i < ncode; i++)
    { lens[code[i]] = leng[i];
      bitv[code[i]] = bits[i];
    }

  { int    j, powr;    //  Fill in a decoder table giving the next Huffman code
    uint16 base;       //    that is a prefix of the next 16 bits

    for (i = 0; i < 256; i++)
      { if (lens[i] > 0)
          { base = (bitv[i] << (16-lens[i]));
            powr = (1 << (16-lens[i]));
            for (j = 0; j < powr; j++)
              look[base+j] = i;
          }
      }
  }

  if (partial)
    { v->esc_code = ecode;
      v->esc_len  = lens[ecode];
      lens[ecode] = 0;
    }
  else
    v->esc_code = -1;
  v->state = CODED_WITH;
}

  //  For debug, give a nice print out of the distribution histogram (if present)
  //     and the Huffman codec

void vcPrint(OneCodec *vc, FILE *to)
{ _OneCodec *v = (_OneCodec *) vc;

  uint64  total_bits, ucomp_bits, count;
  uint16  mask, code, *bits;
  uint64 *hist;
  uint8  *lens;
  int     clen;
  int     hashist;
  int     i, k;

  if (vc == DNAcodec)
    { fprintf(to,"    DNAcompressor\n");
      return;
    }

  if (v->state < CODED_WITH)
    { fprintf(stderr,"vcPrint: Compressor has no codec\n");
      exit (1);
    }
  hashist = (v->state == CODED_WITH);

  bits = v->codebits;
  lens = v->codelens;

  if (hashist)
    { hist = v->hist;
      total_bits = 0;
      ucomp_bits = 0;

      count = 0;
      for (i = 0; i < 256; i++)
        count += hist[i];

      fprintf(to,"\nHistogram:\n");
      for (i = 0; i < 256; i++)
        if (hist[i] > 0)
          { if (isprint(i))
              fprintf(to,"      %c: %12llu %5.1f%%\n",i,hist[i],(hist[i]*100.)/count);
            else
              fprintf(to,"    %3d: %12llu %5.1f%%\n",i,hist[i],(hist[i]*100.)/count);
          }
    }

  fprintf(to,"\nCode Table:\n");
  for (i = 0; i < 256; i++)
    { clen = lens[i];
      if (i == v->esc_code)
        clen = v->esc_len;
      if (clen > 0)
        { mask = (1 << clen);
          code = bits[i];
          if (isprint(i))
            fprintf(to,"   %c: %2d ",i,clen);
          else
            fprintf(to," %3d: %2d ",i,clen);
          for (k = 0; k < clen; k++)
            { mask >>= 1;
              if (code & mask)
                fprintf(to,"1");
              else
                fprintf(to,"0");
            }
          if (i == v->esc_code)
            fprintf(to," ***\n");
          else
            { fprintf(to,"\n");
              if (hashist)
                { total_bits += clen*hist[i];
                  ucomp_bits += (hist[i]<<3);
                }
            }
        }
    }
  if (hashist)
    fprintf(to,"\nTotal Bytes = %lld (%.2f%%)\n",(total_bits-1)/8+1,(100.*total_bits)/ucomp_bits);
}


/*******************************************************************************************
 *
 *  Read and Write Huffman Schemes (actually just (de)serialize)
 *
 ********************************************************************************************/

  //  Maximum # of bytes in a serialized compressor code

int vcMaxSerialSize()
{ return (257 + 2*sizeof(int) + 256*sizeof(uint16)); }

  //  Code the compressor into blob 'out' and return number of bytes in the code

int vcSerialize(OneCodec *vc, void *out)
{ _OneCodec *v = (_OneCodec *) vc;
  
  int     i;
  uint16 *bits;
  uint8  *lens, *o;

  if (vc == DNAcodec)
    return (0);

  if (v->state < CODED_WITH)
    { fprintf(stderr,"vcWrite: Compressor does not have a codec\n");
      exit (1);
    }

  lens = v->codelens;
  bits = v->codebits;
  o    = (uint8 *) out;

  //  Only need to record endian, escape code, code lengths, and codes for those
  //    with non-zero length

  *o++ = v->isbig;
  memcpy(o,&(v->esc_code),sizeof(int));
  o += sizeof(int);
  memcpy(o,&(v->esc_len),sizeof(int));
  o += sizeof(int);
  for (i = 0; i < 256; i++)
    { *o++ = lens[i];
      if (lens[i] > 0 || i == v->esc_code)
        { memcpy(o,bits+i,sizeof(uint16));
          o += sizeof(uint16);
        }
    }
  return (o - (uint8 *) out);
}

  //  Create a compressor object from the serialized code in blob 'in'.
  //    The compressor does not have the original histogram from which
  //    its codec was created.  If the endian of the current machine and
  //    the one that serialized the compressor don't match, then all relevant
  //    items are byte-flipped.

OneCodec *vcDeserialize(void *in)
{ _OneCodec *v;

  char    *look;
  uint8   *lens, *ip;
  uint16  *bits, base;
  int      i, j, powr;

  v = (_OneCodec *) malloc(sizeof(_OneCodec));
  if (v == NULL)
    { fprintf(stderr,"vcRead: Could not allocate compressor\n");
      exit (1);
    }

  v->state = CODED_READ;
  lens = v->codelens;
  bits = v->codebits;
  look = v->lookup;
  ip   = (uint8 *) in;

  { uint32 t;
    uint8 *b;

    t = 1;
    b = (uint8 *) (&t);
    v->isbig = (b[0] == 0);
  }

  if (v->isbig != *ip++)  // If endians out and in don't match then flip item bytes as needed
    { FLIP32(ip)
      memcpy(&(v->esc_code),ip,sizeof(int));
      ip += sizeof(int);
      FLIP32(ip)
      memcpy(&(v->esc_len),ip,sizeof(int));
      ip += sizeof(int);
      for (i = 0; i < 256; i++)
        { lens[i] = *ip++;
          if (lens[i] > 0 || i == v->esc_code)
            { FLIP16(ip)
              memcpy(bits+i,ip,sizeof(uint16));
              ip += sizeof(uint16);
            }
          else
            bits[i] = 0;
        }
    }
  else
    { memcpy(&(v->esc_code),ip,sizeof(int));
      ip += sizeof(int);
      memcpy(&(v->esc_len),ip,sizeof(int));
      ip += sizeof(int);
      for (i = 0; i < 256; i++)
        { lens[i] = *ip++;
          if (lens[i] > 0 || i == v->esc_code)
            { memcpy(bits+i,ip,sizeof(uint16));
              ip += sizeof(uint16);
            }
          else
            bits[i] = 0;
        }
    }

  if (v->esc_code >= 0)
    lens[v->esc_code] = v->esc_len;
  for (i = 0; i < 256; i++)
    { if (lens[i] > 0)
        { base = (bits[i] << (16-lens[i]));
          powr = (1 << (16-lens[i]));
          for (j = 0; j < powr; j++)
            look[base+j] = i;
        }
    }
  if (v->esc_code >= 0)
    lens[v->esc_code] = 0;

  return ((OneCodec *) v);
}


/*******************************************************************************************
 *
 *  Encoders and Decoders
 *
 ********************************************************************************************/

static uint8 Number[128] =
    { 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 3, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
    };

  //  Compress DNA into 2-bits per base

int Compress_DNA(int len, char *s, char *t)
{ int    i, j;
  uint8 *s0, *s1, *s2, *s3;

  s0 = (uint8 *) s;
  s1 = s0+1;
  s2 = s1+1;
  s3 = s2+1;

  len -= 3;
  for (i = j = 0; i < len; i += 4)
    t[j++] = (Number[s0[i]] << 6) | (Number[s1[i]] << 4) | (Number[s2[i]] << 2) | Number[s3[i]];
  switch (i-len)
  { case 0:
      t[j++] = (Number[s0[i]] << 6) | (Number[s1[i]] << 4) | (Number[s2[i]] << 2);
      break;
    case 1:
      t[j++] = (Number[s0[i]] << 6) | (Number[s1[i]] << 4);
      break;
    case 2:
      t[j++] = (Number[s0[i]] << 6);
      break;
    default:
      break;
  }

  return ((len+3)<<1);
}

  //  Encode ibytes[0..ilen) according to compressor vc and place in obytes
  //  Return the # of bits used.

int vcEncode(OneCodec *vc, int ilen, char *ibytes, char *obytes)
{ _OneCodec *v = (_OneCodec *) vc;

  uint64  c, ocode, *ob;
  int     n, k, rem, tbits, ibits, esc, elen;
  uint8  *clens, x, *bcode, *bb;
  uint16 *cbits;

  if (vc == DNAcodec)
    return (Compress_DNA(ilen,ibytes,obytes));

  if (v->state < CODED_WITH)
    { fprintf(stderr,"vcEncode: Compressor does not have a codec\n");
      exit (1);
    }

  esc   = v->esc_code;
  elen  = v->esc_len;
  clens = v->codelens;
  cbits = v->codebits;
  ibits = (ilen << 3);
  bcode = (uint8 *) &ocode;

#define OCODE(L,C)				\
{ rem -= L;					\
  if (rem <= 0)					\
    { ocode |= (C >> (-rem));			\
      *ob++ = ocode;				\
      if (rem < 0)				\
        { rem   += 64;				\
          ocode = (C << rem);			\
        }					\
      else					\
        { rem   = 64;				\
          ocode = 0;				\
        }					\
    } 						\
  else						\
    ocode |= (C << rem);			\
}

  ob    = (uint64 *) obytes;
  tbits = 2;
  rem   = 62;
  if (v->isbig)
    ocode = 0x4000000000000000llu;
  else
    ocode = 0;
  for (k = 0; k < ilen; k++)
    { x = ibytes[k];
      n = clens[x];
      if (n == 0)
        { if (esc < 0)
            { fprintf(stderr,"Compression lib: No code for %c(%x) and no escape code\n",x,x);
              exit (1);
            }
          c = cbits[esc];
          tbits += 8+elen;
          if (tbits > ibits)
            break;
          OCODE(elen,c);
          c = x;
          OCODE(8,c);
        }
      else
        { tbits += n;
          if (tbits > ibits)
            break;
          c = cbits[x];
          OCODE(n,c);
        }
    }
  
  if (k < ilen)
    { *obytes = 0xff;
      memcpy(obytes+1,ibytes,ilen);
      return (ibits+8);
    }

  bb = (uint8 *) ob;
  if (v->isbig)
    { rem = ((71-rem)>>3);
      for (k = 0; k < rem; k++)
        *bb++ = bcode[k];
    }
  else
    { rem = 7 - ((63-rem)>>3);
      for (k = 7; k >= rem; k--)
        *bb++ = bcode[k];
    }

  if (tbits >= 64 && !v->isbig)
    { x = obytes[7];
      obytes[7] = obytes[0];
      obytes[0] = x;
    }

  return (tbits);
}

  //  Uncompress read from 2-bits per base into [0-3] per byte representation

static char Base[4] = { 'a', 'c', 'g', 't' };

int Uncompress_DNA(char *s, int len, char *t)
{ int   i, tlen, byte;
  char *t0, *t1, *t2, *t3;

  t0 = t;
  t1 = t0+1;
  t2 = t1+1;
  t3 = t2+1;

  tlen = len-3;
  for (i = 0; i < tlen; i += 4)
    { byte = *s++;
      t0[i] = Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 4) & 0x3];
      t2[i] = Base[(byte >> 2) & 0x3];
      t3[i] = Base[byte & 0x3];
    }

  switch (i-tlen)
  { case 0:
      byte = *s++;
      t0[i] = Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 4) & 0x3];
      t2[i] = Base[(byte >> 2) & 0x3];
      break;
    case 1:
      byte = *s++;
      t0[i] = Base[(byte >> 6) & 0x3];
      t1[i] = Base[(byte >> 4) & 0x3];
      break;
    case 2:
      byte = *s++;
      t0[i] = Base[(byte >> 6) & 0x3];
      break;
    default:
      break;
  }

  return (len);
}

  //  Decode ilen bits in ibytes, into obytes according to vc's codec
  //  Return the number of bytes decoded.

int vcDecode(OneCodec *vc, int ilen, char *ibytes, char *obytes)
{ _OneCodec *v = (_OneCodec *) vc;

  char   *look;
  uint8  *lens, *q;
  uint64  icode, ncode, *p;
  int     rem, nem;
  uint8   c, *o;
  int     n, k, elen, inbig, esc;

  if (vc == DNAcodec)
    return (Uncompress_DNA(ibytes,ilen>>1,obytes));

  if (v->state < CODED_WITH)
    { fprintf(stderr,"vcDecode: Compressor does not have a codec\n");
      exit (1);
    }

  if (*((uint8 *) ibytes) == 0xff)
    { int olen = (ilen>>3)-1;
      memcpy(obytes,ibytes+1,olen);
      return (olen);
    }

  p = (uint64 *) ibytes;

  inbig = (*ibytes & 0x40);
  if (!inbig && ilen >= 64)
    { uint8 x = ibytes[7];
      ibytes[7] = ibytes[0];
      ibytes[0] = x;
    }

  if (inbig != v->isbig)
    { q = (uint8 *) ibytes;
      for (k = 64; k <= ilen; k += 64)
        { FLIP64(q)
          q += 8;
        }
    }

  lens = v->codelens;
  look = v->lookup;
  esc  = v->esc_code;
  elen = v->esc_len;

#define GET(n)						\
  ilen  -= n;						\
  icode <<= n;						\
  rem   -= n;						\
  while (rem < 16)					\
    { int z = 64-rem;					\
      icode |= (ncode >> rem);				\
      if (nem > z)					\
        { nem -= z;					\
          ncode <<= z;					\
          rem = 64;					\
          break;					\
        }						\
      else						\
        { rem += nem; 					\
          if (rem >= ilen)				\
            break;					\
          else if (ilen-rem < 64)			\
            { nem = ilen-rem;				\
              q = (uint8 *) p;				\
              ncode = 0;				\
              for (k = 0; k < nem; k += 8)		\
                ncode |= (((uint64) (*q++)) << (56-k));	\
            }						\
          else						\
            { ncode = *p++;				\
              nem   = 64;				\
            }						\
	}						\
    }
 
  if (ilen < 64)
    { q = (uint8 *) ibytes;
      icode = 0;
      for (k = 0; k < ilen; k += 8)
        icode |= (((uint64) (*q++)) << (56-k));
    }
  else
    icode = *p++;
  o = (uint8 *) obytes;
  icode <<= 2;
  ilen -= 2;
  rem   = 62;
  if (rem > ilen)
    rem = ilen;
  ncode = 0;
  nem   = 0;
  while (ilen > 0)
    { c = look[icode >> 48];
      if (c == esc)
        { GET(elen)
          c = (icode >> 56);
          GET(8);
        }
      else
        { n = lens[(int) c];
          GET(n)
        }
      *o++ = c;
    }

  return (o - (uint8 *) obytes);
}
  
/***********************************************************************************
 *
 *    UTILITIES: memory allocation, file opening, timer
 *    adapted from Richard's utilities
 *
 **********************************************************************************/

static void die(char *format, ...)
{ va_list args;

  va_start (args, format);
  fprintf (stderr, "FATAL ERROR: ");
  vfprintf (stderr, format, args);
  fprintf (stderr, "\n");
  va_end (args);
  exit (-1);
}

static I64 nAlloc = 0;
static I64 totalAlloc = 0;

static void *myalloc(size_t size)
{ void *p;

  p = malloc(size);
  if (p == NULL) die("myalloc failure requesting %d bytes", size);
  nAlloc     += 1;
  totalAlloc += size;
  return (p);
}

static void *mycalloc(size_t number, size_t size)
{ void *p;

  p = calloc(number,size);
  if (p == NULL) die("mycalloc failure requesting %d objects of size %d", number, size);
  nAlloc     += 1;
  totalAlloc += size*number;
  return p;
}

/********************* end of file ***********************/
