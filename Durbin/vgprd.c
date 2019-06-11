/*  File: vgprd.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description: implementation for vgprd.h
 * Exported functions:
 * HISTORY:
 * Last edited: May  1 08:37 2019 (rd109)
 * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "vgprd.h"

#include "vgpformat_0_1.h"

#include <assert.h>

#define WITH_ZLIB		/* this allows fzopen() to open gzip files - undef if problems */

#include <stdlib.h>		/* strcmp etc. */
#include <string.h>		/* strcmp etc. */

static void writeHeader (VgpFile *vf) ;

void vgpFileClose (VgpFile *vf) /* automatically rewrites header if allowed when writing */
{
  if (vf->isWrite && vf->isHeader && !fseek (vf->f, 0, SEEK_SET)) /* 0 is success for fseek() */
    writeHeader (vf) ;		/* rewrite the header */
  fclose (vf->f) ;
  if (vf->provenance) free (vf->provenance) ; /* can't know who owns subfields to free them */
  if (vf->reference) free (vf->reference) ;
  if (vf->deferred) free (vf->deferred) ;
  int i ;
  for (i = 0 ; i < 128 ; ++i) if (vf->buffer[i] && !vf->isUserBuf[i]) free (vf->buffer[i]) ;
  free (vf) ;
}

VgpFile *vgpFileOpenRead (const char *path, FileType type)
/* opens file and reads header if it exists
   if there is no header then type must be given
   if there is a header and type is non-zero then it must match 
*/
{
  int i ;

  static FileSpecification *formatSpec = 0 ;
  if (!formatSpec) formatSpec = vgpDefineFormat () ; /* initialise from included file */
  
  VgpFile *vf = new0 (1, VgpFile) ;
  if (!(vf->f = fzopen (path, "r"))) { free(vf) ; return 0 ; }

  /* read header - do this by peeking at the first char to check if not alphabetic */
  BOOL isFirst = TRUE ;
  while (!feof (vf->f))
    { char peek = getc(vf->f) ; ungetc(peek, vf->f) ;
      if ((peek >= 'A' && peek <= 'Z') || (peek >= 'a' && peek <= 'z'))
	{ if (isFirst)
	    { if (!type)	/* have to define type in function call or file */
		{ fprintf (stderr, "VGP file open error: type not defined in file or code\n") ;
		  fclose (vf->f) ; free (vf) ; return 0 ;
		}
	      vf->type = type ; vf->spec = &formatSpec[type] ;
	      isFirst = FALSE ;
	    }
	  break ;
	}
      if (isFirst) /* parse the file type and version line */
	{ vf->spec = &formatSpec[0] ; /* this contains the first line specification */
	  if (!vgpReadLine (vf)) die ("badly formatted VGP first header line %d", vf->line) ;
	  { char *s = vgpString (vf) ; FileType t = 0 ;
	    for (i = 1 ; i < MAX_FILE ; ++i) if (!strcmp (s, fileTypeName[i])) t = i ;
	    if (!t) die ("unknown primary file type %s in header line 1", s) ;
	    if (type && t != type)
	      die ("primary file type mismatch file %s != requested %s", s, fileTypeName[type]) ;
	    vf->type = t ; vf->spec = &formatSpec[t] ;
	  }
	  if ((vf->major = vgpInt(vf,1)) != vf->spec->major)
	    die ("major version file %d != parser %d", vf->major, vf->spec->major) ;
	  if ((vf->minor = vgpInt(vf,2)) > vf->spec->minor)
	    die ("minor version file %d > parser %d", vf->minor, vf->spec->minor) ;
	  isFirst = FALSE ;
	}
      else
	{ if (!vgpReadLine (vf)) die ("badly formatted VGP file line %d", vf->line) ;
	  switch (vf->lineType)
	    {
	    case '2':
	      { char *s = vgpString (vf) ; SubType t = 0 ;
		for (i = 1 ; i < MAX_SUB ; ++i) if (!strcmp (s, subTypeName[i])) t = i ;
		if (!t) die ("unknown secondary subtype %s line %d", s, vf->line) ;
		if (subPrimary[t] != vf->type)
		  die ("subtype %s not compatible with primary type %d line %d",
		       s, fileTypeName[vf->type], vf->line) ;
		vf->sub = t ;
		break ;
	      }
	    case '#': case '@': case '+': case '%':
	      { char lt = vgpChar(vf,0) ;
		if (!vf->spec->line[lt]) die ("unknown line type %c line %d", lt, vf->line) ;
		if (vf->lineType == '#') vf->expectCount[lt] = vgpInt(vf,1) ;
		else if (vf->lineType == '@')
		  { vf->expectMax[lt] = vgpInt(vf,1) ;
		    vf->buffer[lt] = new (vf->expectMax[lt]+1, char) ;
		  }
		else if (vf->lineType == '+') vf->expectTotal[lt] = vgpInt(vf,1) ;
		else /* must be % */
		  { // if (vgpChar(vf,1) != vf->spec->objectType)
		    //   die ("object type %c does not match expected %c line %d",
		    //     vgpChar(vf,1), vf->spec->objectType, vf->line) ;
		    lt = vgpChar(vf,2) ;
		    if (vgpChar(vf,1) == '#') vf->expectGroupCount[lt] = vgpInt(vf,3) ;
		    else if (vgpChar(vf,1) == '+') vf->expectGroupTotal[lt] = vgpInt(vf,3) ;
		    else die ("unrecognised symbol %c line %d", vgpChar(vf,1), vf->line) ;
		  }
		break ;
	      }
	    case '!': 
	      { char *prog = vgpString(vf) ;
		char *version = prog + strlen(prog) + 1 ;
		char *command = version + strlen(version) + 1 ;
		char *dateTime = command + strlen(command) + 1 ;
		vgpAddProvenance (vf, prog, version, command, dateTime) ;
	      }
	      break ;
	    case '<': vgpAddReference (vf, vgpString(vf), vgpInt(vf,1)) ; break ;
	    case '>': vgpAddReference (vf, vgpString(vf), vgpInt(vf,1)) ; break ;
	    default: die ("unknown header line type %c line %d", vf->lineType, vf->line) ;
	    }
	}
    }
  
  return vf ;
}

static inline char readChar (VgpFile *vf) ;
static inline I64 readInt (VgpFile *vf) ;
static inline double readReal (VgpFile *vf) ;
static inline void readString (VgpFile *vf, char *buf, I64 n) ;
static inline void readFlush (VgpFile *vf) ; /* reads to the end of the line */

static inline void confirmBufferSize (VgpFile *vf, char t, I64 size)
{
  vf->total[t] += size ;
  int byteSize = vf->spec->line[t]->listByteSize ;
  BOOL isNewMax = FALSE ; if (size > vf->max[t]) { vf->max[t] = size ; isNewMax = TRUE ; }
  if (!vf->isUserBuf[t]) /* check there is space to write */
    { if (!vf->buffer[t])  /* make it */
	{ size = (vf->expectMax[t] > vf->max[t]) ? vf->expectMax[t] : vf->max[t] ;
	  vf->buffer[t] = new (size*byteSize, char) ;
	}
      else if (isNewMax && vf->max[t] > vf->expectMax[t])
	{ free (vf->buffer[t]) ;
	  vf->buffer[t] = new (size*byteSize, char) ;
	}
    }
}

BOOL vgpReadLine (VgpFile *vf)
/* this reads the next line and returns FALSE at end of file or on error
   the line is parsed according to its linetype and contents accessed by macros that follow
*/
{ char t = getc(vf->f) ; /* read first char */
  if (feof (vf->f)) return FALSE ; /* after getting char to deal with final char */
  /* otherwise  presume we are at a good line, and die if not */
  vf->lineType = t ; ++vf->line ;
  LineSpecification *ls = vf->spec->line[t] ;
  if (!ls) die ("unknown line type %c line %d", t, vf->line) ;
  int i, j ;
  for (i = 0 ; ls->field[i] && i < MAX_FIELD ; ++i)
    switch (ls->field[i])
      {
      case INT: vf->field[i].i = readInt (vf) ;  break ;
      case REAL: vf->field[i].r = readReal (vf) ; break ;
      case CHAR: vf->field[i].c = readChar (vf) ; break ;
      case STRING: case INT_LIST: case REAL_LIST: case STRING_LIST:
	{ I64 len = readInt (vf) ;
	  vf->field[i].len = len ;
	  if (ls->field[i] == STRING)
	    { confirmBufferSize (vf, t, len+1) ;
	      char *buf = (char*) vf->buffer[t] ;
	      readString (vf, buf, len) ;
	    }
	  else if (ls->field[i] == INT_LIST)
	    { confirmBufferSize (vf, t, len) ;
	      I64 *buf = (I64*) vf->buffer[t] ;
	      for (j = 0 ; j < len ; ++j) *buf++ = readInt (vf) ;
	    }
	  else if (ls->field[i] == REAL_LIST)
	    { confirmBufferSize (vf, t, len) ;
	      double *buf = (double*) vf->buffer[t] ;
	      for (j = 0 ; j < len ; ++j) *buf++ = readReal (vf) ;
	    }
	  else	/* STRING_LIST - inefficient for now */
	    { I64 totLen = 0 ;
	      char **string = new (len, char*) ;
	      for (j = 0 ; j < len ; ++j)
		{ I64 sLen = readInt (vf) ;
		  totLen += sLen + 1 ;
		  string[j] = new (sLen, char) ;
		  readString (vf, string[j], sLen) ;
		}
	      confirmBufferSize (vf, t, totLen) ;
	      char *buf = (char*) vf->buffer[t] ;
	      for (j = 0 ; j < len ; ++j)
		{ strcpy (buf, string[j]) ;
		  buf += strlen(buf) + 1 ;
		  free (string[j]) ;
		}
	      free (string) ;
	    }
	}
      }
  readFlush (vf) ;
  return TRUE ;
}

void vgpUserBuffer (VgpFile *vf, char lineType, void* buffer)
/* this lets the user reassign the buffer that lists are read into
   if this is not set a default is provided
   in either case, the buffer can be 
   this can be called repeatedly, so the location can be changed, e.g. for each line, or group
   NB the package doesn't check the size - the user must allocate enough memory
   if buffer == 0 then revert to the package default
*/
{ assert (vf && vf->spec && vf->spec->line[lineType]) ;
  if (buffer)
    { if (!vf->isUserBuf[lineType] && vf->buffer[lineType]) free (vf->buffer[lineType]) ;
      vf->buffer[lineType] = buffer ;
      vf->isUserBuf[lineType] = TRUE ;
    }
  else
    vf->isUserBuf[lineType] = FALSE ;
}

VgpFile *vgpFileOpenWrite (const char *path, FileType type, SubType sub, BOOL isGz)
{
  int i, j, k ;

  if (sub && subPrimary[sub] != type)
    die ("subtype %s is not secondary for type %s", subTypeName[sub], fileTypeName[type]) ;
    
  FileSpecification *formatSpec = 0 ;
  if (!formatSpec) formatSpec = vgpDefineFormat () ; /* initialisation */
  
  VgpFile *vf = new0 (1, VgpFile) ;
  if (isGz) vf->f = fzopen (path, "w") ; else vf->f = fopen (path, "w") ;
  if (!vf->f) { free(vf) ; return 0 ; }
  vf->isWrite = TRUE ;
  vf->type = type ;
  vf->sub = sub ;
  vf->spec = &formatSpec[type] ;
  vf->major = vf->spec->major ;
  vf->minor = vf->spec->minor ;
			
  return vf ;
}

BOOL vgpWriteLine (VgpFile *vf, char lineType, void *buf)
/* process is to fill fields by assigning to macros, then call - list contents are in buf */
/* NB adds '\n' before writing line not after, so user fprintf() can add extra material */
/* first call will write initial header, allowing space for count sizes to expand on close */
{
  I64 i, j, len ;
  LineSpecification *ls = vf->spec->line[lineType] ;
  if (!ls) die ("vgpWriteLine error: line type %c not present in file spec %s",
		lineType, fileTypeName[vf->type]) ;
  fprintf (vf->f, "\n%c", lineType) ;
  ++vf->count[lineType] ;
  if (lineType >= 'a' && lineType <= 'z') /* it is a group */
    for (i = 'A' ; i <= 'Z' ; ++i)	  /* update the group info for the previous group */
      { if (vf->count[i] - vf->gCount[i] > vf->groupCount[i])
	  vf->groupCount[i] = vf->count[i] - vf->gCount[i] ;
	if (vf->total[i] - vf->gTotal[i] > vf->groupTotal[i])
	  vf->groupTotal[i] = vf->total[i] - vf->gTotal[i] ;
	vf->gCount[i] = vf->count[i] ; vf->gTotal[i] = vf->total[i] ; /* and reset counters */
      }
  for (i = 0 ; i < MAX_FIELD ; ++i)
   { if (!ls->field[i]) break ;
     switch (ls->field[i])
       {
       case INT: fprintf (vf->f, " %lld", vf->field[i].i) ; break ;
       case REAL: fprintf (vf->f, " %f", vf->field[i].r) ; break ;
       case CHAR: fprintf (vf->f, " %c", vf->field[i].c) ; break ;
       case STRING: case INT_LIST: case REAL_LIST: case STRING_LIST:
	 if (ls->field[i] == STRING) len = strlen(buf) ; else len = vf->field[i].len ;
	 vf->total[lineType] += len ; if (len > vf->max[lineType]) vf->max[lineType] = len ;
	 fprintf (vf->f, " %lld", len) ;
	 if (ls->field[i] == STRING)
	   fprintf (vf->f, " %s", (char*)buf) ;
	 else if (ls->field[i] == INT_LIST)
	   for (j = 0 ; j < len ; ++j)
	     { fprintf (vf->f, " %lld", *((I64*)buf)) ; buf += sizeof(I64) ; }
	 else if (ls->field[i] == REAL_LIST)
	   for (j = 0 ; j < len ; ++j)
	     { fprintf (vf->f, " %f", *((double*)buf)) ; buf += sizeof(double) ; }
	 else			/* STRING_LIST */
	   for (j = 0 ; j < len ; ++j)
	     { fprintf (vf->f, " %lu %s", strlen((char*)buf), (char*)buf) ;
	       buf += strlen((char*)buf) + 1 ;
	     }
       }
   }
  return TRUE ;
}

static void writeHeader (VgpFile *vf)
{
  int i ;
  int N = 0 ;			/* number of chars printed */
  char groupChar = 0 ;
 
  N = fprintf (vf->f, "1 %lu %s %lld %lld",
	       strlen(fileTypeName[vf->type]), fileTypeName[vf->type], vf->major, vf->minor) ;
  if (vf->sub)
    N += fprintf (vf->f, "\n2 %lu %s", strlen(subTypeName[vf->sub]), subTypeName[vf->sub]) ;

  for (i = 0 ; i < 128 ; ++i)
    { if (vf->count[i])
	{ N += fprintf (vf->f, "\n# %c %lld", i, vf->count[i]) ;
	  if (i >= 'a' && i <= 'z')
	    { if (groupChar) die ("code does not support two group types %c %c", groupChar, i) ;
	      groupChar = i ;
	    }
	}
      if (vf->max[i])
	{ N += fprintf (vf->f, "\n@ %c %lld", i, vf->max[i]) ;
	  N += fprintf (vf->f, "\n+ %c %lld", i, vf->total[i]) ;
	}
    }
  if (groupChar)
    for (i = 'A' ; i <= 'Z' ; ++i)
      { if (vf->groupCount[i])
	  N += fprintf (vf->f, "\n%% %c # %c %lld", groupChar, i, vf->groupCount[i]) ;
	if (vf->groupTotal[i])
	  N += fprintf (vf->f, "\n%% %c + %c %lld", groupChar, i, vf->groupTotal[i]) ;
      }
  
  Provenance *p = vf->provenance ; 
  for (i = vf->count['!'] ; i-- ; ++p)
    N += fprintf (vf->f, "\n! %s %s %s %s", p->program, p->version, p->command, p->date) ;
  Reference *r = vf->reference ;
  for (i = vf->count['<'] ; i-- ; ++r)
    N += fprintf (vf->f, "\n< %s %lld", r->filename, r->count) ;
  r = vf->deferred ;
  for (i = vf->count['>'] ; i-- ; ++r)
    N += fprintf (vf->f, "\n> %s %lld", r->filename, r->count) ;

  if (!vf->isHeader)
    { vf->isHeader = TRUE ;
      int nLineType = 0, nListType = 0 ;
      for (i = 0 ; i < 128 ; ++i)
	if (vf->spec->line[i])
	  { ++nLineType ; if (vf->spec->line[i]->listByteSize) ++nListType ; }
      vf->headerSize = N + (25+27)*nLineType + 25*nListType + (25+27)*nListType ;
	/* N covers lines 1, 2, provenance and references, then rest for #, @, + and % lines */
      char *space = new0 (N - vf->headerSize, char) ;
      fwrite (space, 1, N - vf->headerSize, vf->f) ; /* write buffer space */
    }
  else if (N > vf->headerSize)
    die ("new header size %d is bigger than space allowed for %d", N, vf->headerSize) ;
}

/************ the next set of functions manage provenance and reference lines ***************/

static BOOL addProvenance (VgpFile *vf, Provenance *from, int n)
{
  if (!n) return FALSE ;
  if (vf->isHeader) die ("can't addProvenance after writing header") ;
  Provenance *p = new (vf->count['!'] + n, Provenance) ;
  I64 *count = &vf->count['!'], *max = &vf->max['!'], *total = &vf->total['!'] ;
  if (vf->provenance) memcpy (p, vf->provenance, *count*sizeof(Provenance)) ;
  memcpy (p+*count*sizeof(Provenance), from, n*sizeof(Provenance)) ;
  if (vf->provenance) free (vf->provenance) ;
  vf->provenance = p ;
  *count += n ;
  while (n--)
    { int len = strlen (from->program) ; *total += len ; if (len > *max) *max = len ;
      len = strlen (from->version) ; *total += len ; if (len > *max) *max = len ;
      len = strlen (from->command) ; *total += len ; if (len > *max) *max = len ;
      len = strlen (from->date) ; *total += len ; if (len > *max) *max = len ;
      ++from ;
    }
  return TRUE ;
}

BOOL vgpInheritProvenance (VgpFile *vf, VgpFile *source)
{ return addProvenance (vf, source->provenance, source->count['!']) ; }

#include <time.h>		/* for time functions */

BOOL vgpAddProvenance (VgpFile *vf, char *prog, char *version, char *command, char *date)
{ Provenance p ; p.program = prog ; p.version = version ; p.command = command ;
  if (date) p.date = date ;
  else
    { p.date = new (32, char) ;
      time_t tim = time(0) ;
      strftime (p.date, 32, "%F_%T", localtime (&tim)) ;
    }
  return addProvenance (vf, &p, 1) ;
}

static BOOL addReference (VgpFile *vf, Reference *from, int n, BOOL isDeferred)
{
  if (!n) return FALSE ;
  if (vf->isHeader) die ("can't addReference after writing header") ;
  Reference **target = isDeferred ? &vf->deferred : &vf->reference ;
  I64 *count, *max, *total ;
  if (isDeferred) { count = &vf->count['>'] ; max = &vf->max['>'] ; total = &vf->total['>'] ; }
  else { count = &vf->count['<'] ; max = &vf->max['<'] ; total = &vf->total['<'] ; }
  Reference *r = new (*count + n, Reference) ;
  if (*target) memcpy (r, *target, *count*sizeof(Reference)) ;
  memcpy (r+*count*sizeof(Reference), from, n*sizeof(Reference)) ;
  if (*target) free (*target) ;
  *target = r ;
  *count += n ;
  while (n--)
    { int len = strlen (from->filename) ; *total += len ; if (len > *max) *max = len ;
      ++from ;
    }
  return TRUE ;
}

BOOL vgpInheritReference (VgpFile *vf, VgpFile *source) /* as for provenance */
{ return addReference (vf, source->reference, source->count['<'], FALSE) ; }

BOOL vgpAddReference (VgpFile *vf, char *filename, I64 count)
{ Reference ref ; ref.filename = filename ; ref.count = count ;
  return addReference (vf, &ref, 1, FALSE) ;
}

BOOL vgpInheritDeferred (VgpFile *vf, VgpFile *source)
{ return addReference (vf, source->deferred, source->count['>'], TRUE) ; }

BOOL vgpAddDeferred (VgpFile *vf, char *filename, I64 count)
{ Reference ref ; ref.filename = filename ; ref.count = count ;
  return addReference (vf, &ref, 1, TRUE) ;
}

static inline void eatWhite (VgpFile *vf)
{ char x = getc (vf->f) ;
  if (x != ' ' && x != '\t') die ("failed to find expected whitespace line %d", vf->line) ;
}

static inline char readChar (VgpFile *vf) { eatWhite (vf) ; return getc (vf->f) ; }

static inline char* readBuf (VgpFile *vf)
{ eatWhite (vf) ;
  static char buf[32] ; char *endBuf = buf+32 ;
  char *cp = buf ; --cp ;
  while ((++cp < endBuf) && (*cp = getc(vf->f)))
    if (*cp == ' ' || *cp == '\t' || *cp == '\n' || *cp == EOF) break ;
  if (cp == endBuf) { *--cp = 0 ; die ("overlong item %s line %d", buf, vf->line) ; }
  ungetc (*cp, vf->f) ; *cp = 0 ;
  return buf ;
}

static inline I64 readInt (VgpFile *vf)
{ char *ep, *buf = readBuf (vf) ;
  I64 x = strtoll (buf, &ep, 10) ; if (*ep) die ("bad int line %d", vf->line) ;
  return x ;
}

static inline double readReal (VgpFile *vf)
{ char *ep, *buf = readBuf (vf) ;
  double x = strtod (buf, &ep) ; if (*ep) die ("bad real line %d", vf->line) ;
  return x ;
}

static inline void readString (VgpFile *vf, char *buf, I64 n)
{ eatWhite (vf) ;
  char *cp = buf ; --cp ;
  while (n-- && (*++cp = getc(vf->f)))
    if (*cp == '\n' || *cp == EOF) break ;
  if (++n) { die ("line too short %d", buf, vf->line) ; }
  *++cp = 0 ;
}

static inline void readFlush (VgpFile *vf) /* reads to the end of the line */
{ char x ;
  while ((x = getc (vf->f)) && x != '\n')
    if (x == EOF) die ("premature end of file line %d", vf->line) ;
}

/*********** utilities inlined here so we don't need to import RD package *************/

#include <stdarg.h>

void die (char *format, ...)
{
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "FATAL ERROR: ") ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, "\n") ;
  va_end (args) ;

  exit (-1) ;
}

long totalAllocated = 0 ;

void *myalloc (size_t size)
{
  void *p = (void*) malloc (size) ;
  if (!p) die ("myalloc failure requesting %d bytes", size) ;
  totalAllocated += size ;
  return p ;
}

void *mycalloc (size_t number, size_t size)
{
  void *p = (void*) calloc (number, size) ;
  if (!p) die ("mycalloc failure requesting %d objects of size %d", number, size) ;
  totalAllocated += size*number ;
  return p ;
}

#ifdef WITH_ZLIB
#include <zlib.h>
#endif

FILE *fzopen(const char *path, const char *mode)
{  /* very cool from https://stackoverflow.com/users/3306211/fernando-mut */
#ifdef WITH_ZLIB
  gzFile zfp;			/* fernando said *zfp - makes me worry.... */

  /* try gzopen */
  zfp = gzopen(path,mode);
  if (zfp == NULL)
    return fopen(path,mode);

  /* open file pointer */
  return funopen(zfp,
                 (int(*)(void*,char*,int))gzread,
                 (int(*)(void*,const char*,int))gzwrite,
                 (fpos_t(*)(void*,fpos_t,int))gzseek,
                 (int(*)(void*))gzclose);
#else
  return fopen(path,mode);
#endif
}

FILE *fopenTag (char* root, char* tag, char* mode)
{
  if (strlen (tag) > 30) die ("tag %s in fopenTag too long - should be < 30 chars", tag) ;
  char *fileName = new (strlen (root) + 32, char) ;
  strcpy (fileName, root) ;
  strcat (fileName, ".") ;
  strcat (fileName, tag) ;
  FILE *f = fzopen (fileName, mode) ;
  free (fileName) ;
  return f ;
}

/***************** rusage for timing information ******************/

#include <sys/resource.h>
#ifndef RUSAGE_SELF     /* to prevent "RUSAGE_SELF redefined" gcc warning, fixme if this is more intricate */
#define RUSAGE_SELF 0
#endif

#ifdef RUSAGE_STRUCTURE_DEFINITIONS
struct rusage {
  struct timeval ru_utime; /* user time used */
  struct timeval ru_stime; /* system time used */
  long ru_maxrss;          /* integral max resident set size */
  long ru_ixrss;           /* integral shared text memory size */
  long ru_idrss;           /* integral unshared data size */
  long ru_isrss;           /* integral unshared stack size */
  long ru_minflt;          /* page reclaims */
  long ru_majflt;          /* page faults */
  long ru_nswap;           /* swaps */
  long ru_inblock;         /* block input operations */
  long ru_oublock;         /* block output operations */
  long ru_msgsnd;          /* messages sent */
  long ru_msgrcv;          /* messages received */
  long ru_nsignals;        /* signals received */
  long ru_nvcsw;           /* voluntary context switches */
  long ru_nivcsw;          /* involuntary context switches */
};

struct timeval {
  time_t       tv_sec;   /* seconds since Jan. 1, 1970 */
  suseconds_t  tv_usec;  /* and microseconds */
} ;
#endif /* RUSAGE STRUCTURE_DEFINITIONS */

static struct rusage rOld, rFirst ;

void timeUpdate (FILE *f)
{
  static BOOL isFirst = 1 ;
  struct rusage rNew ;
  int secs, usecs ;

  getrusage (RUSAGE_SELF, &rNew) ;
  if (!isFirst)
    { secs = rNew.ru_utime.tv_sec - rOld.ru_utime.tv_sec ;
      usecs =  rNew.ru_utime.tv_usec - rOld.ru_utime.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "user\t%d.%06d", secs, usecs) ;
      secs = rNew.ru_stime.tv_sec - rOld.ru_stime.tv_sec ;
      usecs =  rNew.ru_stime.tv_usec - rOld.ru_stime.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "\tsystem\t%d.%06d", secs, usecs) ;
      fprintf (f, "\tmax_RSS\t%ld", rNew.ru_maxrss - rOld.ru_maxrss) ;
      fprintf (f, "\tmemory\t%li", totalAllocated) ;   
      fputc ('\n', f) ;
    }
  else
    { rFirst = rNew ;
      isFirst = FALSE ;
    }

  rOld = rNew ;
}

void timeTotal (FILE *f) { rOld = rFirst ; timeUpdate (f) ; }

/********************* end of file ***********************/
