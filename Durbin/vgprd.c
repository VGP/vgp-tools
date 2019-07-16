/*  File: vgprd.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description: implementation for vgprd.h
 * Exported functions:
 * HISTORY:
 * Last edited: Jul 16 16:17 2019 (rd109)
 * * Jul  8 04:28 2019 (rd109): refactored to use lineInfo[]
 * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "vgprd.h"

#include "vgpformat_1_0.h"	/* includes the format specification */

#include <assert.h>
#include <sys/errno.h>		/* for errno */
#include <sys/types.h>		/* for off_t */

//#define WITH_ZLIB		/* this allows fzopen() to open gzip files - undef if problems */

#include <stdlib.h>		/* strcmp etc. */
#include <string.h>		/* strcmp etc. */

/****************** creators and destructors *****************/

static VgpFile *vgpFileCreate (FileType fileType)
{
  VgpFile *vf = new0 (1, VgpFile) ;

  /* define the header lines */
  vf->lineInfo['1'] = vgpDefineLine (STRING, INT, INT, 0, 0, 0) ; /* type, major, minor */
  vf->lineInfo['2'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0) ;     /* subtype */
  vf->lineInfo['#'] = vgpDefineLine (CHAR, INT, 0, 0, 0, 0) ;     /* linetype count */
  vf->lineInfo['@'] = vgpDefineLine (CHAR, INT, 0, 0, 0, 0) ;     /* linetype max */
  vf->lineInfo['+'] = vgpDefineLine (CHAR, INT, 0, 0, 0, 0) ;     /* linetype total */
  vf->lineInfo['%'] = vgpDefineLine (CHAR, CHAR, CHAR, INT, 0, 0) ; /* group #/+ linetype value */
  vf->lineInfo['!'] = vgpDefineLine (STRING_LIST, 0, 0, 0, 0, 0) ; /* name version command date */
  vf->lineInfo['<'] = vgpDefineLine (STRING, INT, 0, 0, 0, 0) ;   /* filename objectcount */
  vf->lineInfo['>'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0) ;     /* filename */
  /* below here is magic for binary files */
  vf->lineInfo['$'] = vgpDefineLine (0, 0, 0, 0, 0, 0) ;  /* goto 32 bytes from end of file */
  vf->lineInfo[1] = vgpDefineLine (CHAR, STRING, 0, 0, 0, 0) ; /* field codec (binary only) */
  vf->lineInfo[2] = vgpDefineLine (CHAR, STRING, 0, 0, 0, 0) ; /* list codec (binary only) */
  vf->lineInfo['&'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0) ;   /* object index */
  vf->lineInfo['^'] = vgpDefineLine (0, 0, 0, 0, 0, 0) ;       /* end of footer - return to top */
  vf->lineInfo['-'] = vgpDefineLine (INT, 0, 0, 0, 0, 0) ; /* seek offset to footer */

  defineFormat (vf, fileType) ;		/* from vgpformat_1_0.h */

      /* fill in the nField and listByteSize information, and perform format consistency checks */
  int j, k, n = 0 ;
  int listSize[8] = { 0, 0, 0, 0, 1, sizeof(I64), sizeof(double), 1 } ;

  for (j = 0 ; j < 128 ; ++j)
    if (vf->lineInfo[j])
      { for (k = 0 ; k < MAX_FIELD ; k++) /*  */
	  if (listSize[vf->lineInfo[j]->fieldSpec[k]])
	    { if (vf->lineInfo[j]->listByteSize)
		die ("VGP spec error file type %s: two list types in record %d",
		     fileTypeName[vf->fileType], j) ;
	      vf->lineInfo[j]->listByteSize = listSize[vf->lineInfo[j]->fieldSpec[k]] ;
	      vf->lineInfo[j]->listField = k+1 ;
	    }
	for (k = 0 ; k < MAX_FIELD && vf->lineInfo[j]->fieldSpec[k] ; ++k) ;
	vf->lineInfo[j]->nField = k ;
	if (j >= 'a' && j <= 'z' && j != vf->groupType)
	  die ("VGP spec error file type %s: extra group type %c not permitted",
	       fileTypeName[vf->fileType], j) ;
      }
  if (vf->objectType && !vf->lineInfo[vf->objectType])
    die ("VGP spec error file type %s: missing line spec for object type %c",
	 fileTypeName[vf->fileType], vf->objectType) ;
  if (vf->groupType && !vf->lineInfo[vf->groupType])
    die ("VGP spec error file type %s: missing line spec for group type %c",
	 fileTypeName[vf->fileType], vf->groupType) ;

  /* make the binaryType pack and unpack lookups */
  for (j = 0 ; j < 128 ; ++j)
    if (vf->lineInfo[j])
      { unsigned char x = (n++ << 2) | 0x80 ;
	vf->lineInfo[j]->binaryTypePack = x ;
	vf->binaryTypeUnpack[x] = j ;
	vf->binaryTypeUnpack[x | 0x01] = j ;
	vf->binaryTypeUnpack[x | 0x02] = j ;
	vf->binaryTypeUnpack[x | 0x03] = j ;
      }
  if (n >= 32) die ("VGP spec error file type %s: too many line specs %d >= 32",
		    fileTypeName[vf->fileType], n) ;

  vf->codecTrainingSize = 100000 ;
  vf->lineInfo[1]->bufSize = vcMaxSerialSize() + 1 ; /* +1 for added but unused 0-terminator */
  vf->lineInfo[1]->buffer = new(vf->lineInfo[1]->bufSize,char) ;
  vf->lineInfo[2]->bufSize = vcMaxSerialSize() + 1 ;
  vf->lineInfo[2]->buffer = new(vf->lineInfo[2]->bufSize,char) ;

  return vf ;
}

static void lineInfoDestroy (LineInfo *li)
{
  if (li->buffer && !li->isUserBuf) free (li->buffer) ;
  if (li->fieldCodec) vcDestroy (li->fieldCodec) ;
  if (li->listCodec) vcDestroy (li->listCodec) ;
}

void vgpFileDestroy (VgpFile *vf)
{
  if (vf->provenance) free (vf->provenance) ; /* can't know who owns subfields to free them */
  if (vf->reference) free (vf->reference) ;
  if (vf->deferred) free (vf->deferred) ;
  if (vf->codecBuf) free (vf->codecBuf) ;
  if (vf->f) fclose (vf->f) ;
  int i ; for (i = 0 ; i < 128 ; ++i) if (vf->lineInfo[i]) lineInfoDestroy (vf->lineInfo[i]) ;
  free (vf) ;
}

/************ top level read functions ***********/

static void parseError (VgpFile *vf, char *format, ...) ;
static char inline vfGetc (VgpFile *vf) { char c = getc (vf->f) ; if (vf->linePos < 127) vf->lineBuf[vf->linePos++] = c ; return c ; }

VgpFile *vgpFileOpenRead (const char *path, FileType fileType)
/* opens file and reads header if it exists
   if there is no header then fileType must be given
   if there is a header and fileType is non-zero then it must match 
   if there is $ line (file is binary) then read footer including decompressors and index
*/
{
  int i ;
  
  VgpFile *vf = vgpFileCreate (fileType) ;
  if (!(vf->f = fopen (path, "r"))) { free(vf) ; return 0 ; }

  /* read header and (optionally) footer
     recognise end of header by peeking at the first char to check if alphabetic 
  */
  while (TRUE)
    { unsigned char peek = getc(vf->f) ; if (feof (vf->f)) break ; /* loop exit at end of file */
      ungetc(peek, vf->f) ;
      if (peek & 0x80) peek = vf->binaryTypeUnpack[peek] ;
      if ((peek >= 'A' && peek <= 'Z') || (peek >= 'a' && peek <= 'z'))
	{ if (!vf->fileType)    /* have to define fileType in function call or file */
	    { fprintf (stderr, "VGP file open error: fileType not defined in file or code\n") ;
	      vgpFileDestroy (vf) ;
	      return 0 ;
	    }
	  break ;		/* loop exit at standard data line */
	}
      else			/* this is a header line */
	{ if (peek == '!')	/* hack to insert a count of 4 for STRING_LIST */
	    { getc(vf->f) ; ungetc('4',vf->f) ; ungetc(' ',vf->f) ; ungetc('!',vf->f) ; }
	  vgpReadLine (vf) ;	/* can't fail because we checked file end already */
	  if (vf->line == 1 && vf->lineType != '1') die ("first header line must start with 1") ;
	  switch (vf->lineType)
	    {
	    case '1':
	      if (vf->line != 1) parseError (vf, "1 line not first") ;
	      if (vgpInt(vf,1) != vf->major)
		die ("major version file %d != code %d", vgpInt(vf,1), vf->major) ;
	      if (vgpInt(vf,2) != vf->minor)
		die ("minor version file %d > code %d", vgpInt(vf,2), vf->minor) ;
	      { char *s = vgpString (vf) ; FileType t = 0 ;
		for (i = 1 ; i < MAX_FILE ; ++i) if (!strcmp (s, fileTypeName[i])) t = i ;
		if (!t) die ("unknown primary fileType %s in header line 1", s) ;
		if (!fileType)
		  { FILE *f = vf->f ; vf->f = 0 ; vgpFileDestroy (vf) ;
		    vf = vgpFileCreate (t) ; vf->f = f ; vf->line = 1 ;
		  }
		else if (fileType && t != fileType)
		  die ("primary fileType mismatch file %s != %s", s, fileTypeName[fileType]) ;
	      }
	      break ;
	    case '2':
	      { char *s = vgpString (vf) ; SubType t = 0 ;
		for (i = 1 ; i < MAX_SUB ; ++i) if (!strcmp (s, subTypeName[i])) t = i ;
		if (!t) parseError (vf, "unknown secondary subType %s", s) ;
		if (subPrimary[t] != vf->fileType)
		  parseError (vf, "subtype %s not compatible with primary type %d",
			      s, fileTypeName[vf->fileType]) ;
		vf->subType = t ;
		break ;
	      }
	    case '#': case '@': case '+': case '%':
	      { char c = vgpChar(vf,0) ;
		LineInfo *li = vf->lineInfo[c] ;
		if (!li) parseError (vf, "unknown line type %c", c) ;
		if (vf->lineType == '#')
		  { li->expectCount = vgpInt(vf,1) ;
		    if (c == vf->objectType)
		      { vf->lineInfo['&']->bufSize = li->expectCount ;
			vf->lineInfo['&']->buffer = new(li->expectCount*sizeof(I64), char) ;
		      }
		  }
		else if (vf->lineType == '@')
		  { li->expectMax = vgpInt(vf,1) ;
		    li->bufSize = li->expectMax + 1 ; /* allow for string terminators */
		    li->buffer = new(li->bufSize * li->listByteSize, char) ;
		  }
		else if (vf->lineType == '+')
		  li->expectTotal = vgpInt(vf,1) ;
		else /* must be % */
		  { c = vgpChar(vf,2) ;
		    li = vf->lineInfo[c] ;
		    if (!li) parseError (vf, "unknown line type %c", c) ;
		    if (vgpChar(vf,1) == '#') li->expectGroupCount = vgpInt(vf,3) ;
		    else if (vgpChar(vf,1) == '+') li->expectGroupTotal = vgpInt(vf,3) ;
		    else parseError (vf, "unrecognised symbol %c", vgpChar(vf,1)) ;
		  }
		break ;
	      }
	    case '!': 		/* NB need to copy the strings - small memory leak */
	      { char *prog = vgpString(vf) ;
		char *version = prog + strlen(prog) + 1 ;
		char *command = version + strlen(version) + 1 ;
		char *dateTime = strdup(command + strlen(command) + 1) ;
		command = strdup(command) ; version = strdup(version) ; prog = strdup(prog) ;
		--vf->lineInfo['!']->count ; /* to avoid double counting */
		vgpAddProvenance (vf, prog, version, command, dateTime) ;
	      }
	      break ;
	    case '<':
	      --vf->lineInfo['<']->count ; /* to avoid double counting */
	      vgpAddReference (vf, vgpString(vf), vgpInt(vf,1)) ;
	      break ;
	    case '>':
	      --vf->lineInfo['>']->count ; /* to avoid double counting */
	      vgpAddDeferred (vf, vgpString(vf)) ;
	      break ;
	    /* below here are binary file header types - require expectCount/expectMax first */
	    case '$':  /* read footer - goto end, find offset to start of footer and go there */
	      errno = 0 ; off_t startOff = ftello (vf->f) ; fseek (vf->f, -32, SEEK_END) ;
	      if (errno) die ("can't seek to final line") ;
	      vgpReadLine (vf) ;
	      fseeko (vf->f, (off_t) vgpInt(vf,0), SEEK_SET) ;
	      if (errno) die ("can't seek to start of footer") ;
	      break ;
	    case '^': /* end of footer - return to where we left header */
	      errno = 0 ; fseeko (vf->f, startOff, SEEK_SET) ; if (errno) die ("can't seek back") ;
	      break ;
	    case '&':
	      vf->isIndex = TRUE ;
	      break ;
	    case 1:
	      vf->lineInfo[vgpChar(vf,0)]->fieldCodec = vcDeserialize (vgpString(vf)) ;
	      break ;
	    case 2:
	      vf->lineInfo[vgpChar(vf,0)]->listCodec = vcDeserialize (vgpString(vf)) ;
	      break ;
	    default: parseError (vf, "unknown header line type %c", vf->lineType) ;
	    }
	}
    }

  /* allocate codec buffer - always allocate enough to handle fields of all line types */
  I64 size = MAX_FIELD*sizeof(Field) ;
  for (i = 0 ; i < 128 ; ++i)
    if (vf->lineInfo[i] && vf->lineInfo[i]->listCodec &&
	size < vf->lineInfo[i]->expectMax * vf->lineInfo[i]->listByteSize)
      size = vf->lineInfo[i]->expectMax * vf->lineInfo[i]->listByteSize ;
  vf->codecBuf = new (++size, char) ; /* add one for worst case codec usage */
  
  return vf ;
}

static inline char readChar (VgpFile *vf) ;
static inline I64 readInt (VgpFile *vf) ;
static inline double readReal (VgpFile *vf) ;
static inline void readString (VgpFile *vf, char *buf, I64 n) ;
static inline void readFlush (VgpFile *vf) ; /* reads to the end of the line */

static inline void confirmBufferSize (VgpFile *vf, char t, I64 size, I64 nStrings)
{ LineInfo *li = vf->lineInfo[t] ;
  li->total += size ; if (size > li->max) li->max = size ;
  size += nStrings ;		/* need to allocate space for terminal 0s */
  if (!li->isUserBuf && size > li->bufSize)   /* expand buffer */
    { if (li->buffer) free (li->buffer) ;
      li->bufSize = size ; li->buffer = new (size*li->listByteSize, char) ;
    }
}

static inline void updateGroupCount (VgpFile *vf, BOOL isGroupLine)
{ int i ;
  for (i = 'A' ; i <= 'Z' ; ++i)
    if (vf->lineInfo[i] && vf->lineInfo[i]->count)
      { LineInfo *li = vf->lineInfo[i] ;
	if (vf->inGroup && li->groupCount < li->count - li->gCount)
	  li->groupCount = li->count - li->gCount ;
	if (vf->inGroup && li->groupTotal < li->total - li->gTotal)
	  li->groupTotal = li->total - li->gTotal ;
	li->gCount = li->count ; li->gTotal = li->total ;
      }
  if (isGroupLine)
    { ++vf->group ;
      vf->inGroup = TRUE ;
    }
}

static void readStringList (VgpFile *vf, char t, I64 len) /* annoying */
{ int j ;
  I64 totLen = 0 ;
  char **string = new (len, char*) ;
  for (j = 0 ; j < len ; ++j)
    { I64 sLen = readInt (vf) ;
      totLen += sLen ;
      string[j] = new (sLen+1, char) ;
      readString (vf, string[j], sLen) ;
    }
  confirmBufferSize (vf, t, totLen, len) ;
  char *buf = (char*) vf->lineInfo[t]->buffer ;
  for (j = 0 ; j < len ; ++j)
    { strcpy (buf, string[j]) ;
      buf += strlen(buf) + 1 ;
      free (string[j]) ;
    }
  free (string) ;
}

BOOL vgpReadLine (VgpFile *vf)
/* this reads the next line and returns FALSE at end of file or on error
   the line is parsed according to its linetype and contents accessed by macros that follow
   the top bit of the first character determines whether the line is binary or ascii
*/
{
  BOOL isAscii ;
  char t ;
  
  vf->linePos = 0 ;		/* must come before first vfGetc() */
  unsigned char x = vfGetc (vf) ;	/* read first char */
  if (feof (vf->f) || x == '\n') /* blank line (x=='\n') is end of records marker before footer */
    { updateGroupCount (vf, FALSE) ;
      return FALSE ;
    }

  ++vf->line ;			/* otherwise assume this is a good line, and die if not */
  if (x & 0x80) { isAscii = FALSE ; t = vf->binaryTypeUnpack[x] ; }
  else { isAscii = TRUE ; t = x ; }
  vf->lineType = t ;
  LineInfo *li = vf->lineInfo[t] ;
  if (!li) parseError (vf, "unknown line type %c", t) ;
  ++li->count ;
  if (t == vf->objectType) ++vf->object ;
  if (t == vf->groupType) updateGroupCount (vf, TRUE) ;

  // if (vf->line < 50) printf ("reading line %lld type %d = %c\n", vf->line, t, t) ;

  if (isAscii)			/* read field by field according to spec */
    { int i, j ;
      I64 len ;
      for (i = 0 ; i < MAX_FIELD && li->fieldSpec[i] ; ++i)
	switch (li->fieldSpec[i])
	  {
	  case INT: vf->field[i].i = readInt (vf) ;  break ;
	  case REAL: vf->field[i].r = readReal (vf) ; break ;
	  case CHAR: vf->field[i].c = readChar (vf) ; break ;
	  case STRING:
	    len = readInt (vf) ; vf->field[i].len = len ;
	    confirmBufferSize (vf, t, len, 1) ;
	    readString (vf, (char*) li->buffer, len) ;
	    break ;
	  case INT_LIST:
	    len = readInt (vf) ; vf->field[i].len = len ;
	    confirmBufferSize (vf, t, len, 0) ;
	    { I64 *buf = (I64*) li->buffer ;
	      for (j = 0 ; j < len ; ++j) *buf++ = readInt(vf) ;
	    }
	    break ;
	  case REAL_LIST:
	    len = readInt (vf) ; vf->field[i].len = len ;
	    confirmBufferSize (vf, t, len, 0) ;
	    { double *buf = (double*) li->buffer ;
	      for (j = 0 ; j < len ; ++j) *buf++ = readReal (vf) ;
	    }
	    break ;
	  case STRING_LIST: /* STRING_LIST - inefficient for now - also used for binary */
	    len = readInt (vf) ; vf->field[i].len = len ;
	    readStringList (vf, t, len) ;
	    break ;
	  }
      readFlush (vf) ;
    }
  else				/* binary - block read fields and list, potentially compressed */
   { int nField = li->nField ;
     if (nField)
       { if (x & 0x1) 		/* fields are compressed */
	   { int nBits = (unsigned char) getc (vf->f) ; /* NB only compress fields if under 255 bits */
	     if (fread (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1) die ("compressed fields") ;
	     vcDecode (li->fieldCodec, nBits, vf->codecBuf, (char*) vf->field) ;
	   }
	 else
	   if (fread (vf->field, sizeof(Field), nField, vf->f) != nField) die ("fields") ;
       }
     int listField = li->listField ;
     if (listField-- && vgpLen(vf, listField)) 	/* NB subtract one to set actual field number */
       { I64 listLen = vgpLen(vf, listField) ;
	 li->total += listLen ; if (listLen > li->max) li->max = listLen ;
	 if (x & 0x2)		/* list is compressed */
	   { I64 nBits ;
	     if (fread (&nBits, sizeof(I64), 1, vf->f) != 1) die ("list nBits") ;
	     if (fread (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1) die ("compressed list") ;
	     vcDecode (li->listCodec, nBits, vf->codecBuf, li->buffer) ;
	   }
	 else if (li->fieldSpec[listField] == STRING_LIST) /* handle as ASCII */
	   readStringList (vf, t, vgpLen(vf,listField)) ;
	 else
	   { I64 size = vgpLen(vf,listField) * li->listByteSize ;
	     if (fread (li->buffer, size, 1, vf->f) != 1) die ("list") ;
	     if (li->fieldSpec[listField] == STRING)
	       ((char*)li->buffer)[size] = 0 ; /* 0 terminate */
	   }
       }
   }

  return TRUE ;
}

void vgpUserBuffer (VgpFile *vf, char lineType, void* buffer)
/* this lets the user reassign the buffer that lists are read into
   if this is not set a default is provided
   this can be called repeatedly, so the location can be changed, e.g. for each line, or group
   NB the package doesn't check the size - the user must allocate enough memory
   if buffer == 0 then revert to the package default
*/
{ LineInfo *li = vf->lineInfo[lineType] ;
  if (buffer)
    { if (!li->isUserBuf && li->buffer) { free (li->buffer) ; li->bufSize = 0 ; }
      li->buffer = buffer ;
      li->isUserBuf = TRUE ;
    }
  else
    { if (li->isUserBuf)
	{ li->bufSize = li->expectMax ;
	  li->buffer = new (li->expectMax * li->listByteSize, char) ;
	}
      li->isUserBuf = FALSE ;
    }
}

void vgpClose (VgpFile *vf) /* automatically rewrites header if allowed when writing */
{
  updateGroupCount (vf, FALSE) ; /* count last group */

  if (vf->isWrite)
    fputc ('\n', vf->f) ;	/* end of file if ascii, end of data marker if binary */

  if (vf->isWrite && vf->isBinary) /* write the footer */
    { int i, j ;
      errno = 0 ; off_t footOff = ftello (vf->f) ; if (errno) die ("failed footer ftell") ;
      
      /*  first the content information, which has to be ascii because we don't have codecs yet */
      char g = vf->groupType ; if (!g) g = 'Z' ; /* if no group then just loop over upper case */
      char *buf = new (vcMaxSerialSize (), char) ;
      for (i = 'A' ; i <= g ; ++i)
	if (vf->lineInfo[i] && vf->lineInfo[i]->count)
	  { LineInfo *li = vf->lineInfo[i] ;
	    fprintf (vf->f, "# %c %lld\n", i, li->count) ;
	    if (li->max) fprintf (vf->f, "@ %c %lld\n", i, li->max) ;
	    if (li->total) fprintf (vf->f, "+ %c %lld\n", i, li->total) ;
	    if (li->groupCount) fprintf (vf->f, "%% %c # %c %lld\n", g, i, li->groupCount) ;
	    if (li->groupTotal) fprintf (vf->f, "%% %c + %c %lld\n", g, i, li->groupTotal) ;
	    if (li->isUseFieldCodec)
	      { vgpChar(vf,0) = i ;
		vgpLen(vf,1) = vcSerialize (li->fieldCodec, buf) ;
		vgpWriteLine (vf, 1, buf) ;
	      } /* list codec must come after count */
	    if (li->isUseListCodec && li->listCodec != DNAcodec)
	      { vgpChar(vf,0) = i ;
		vgpLen(vf,1) = vcSerialize (li->listCodec, buf) ;
		vgpWriteLine (vf, 2, buf) ;
	      }
	  }
      free (buf) ;

      /* then the index, which might be compressed by a codec */
      vgpLen(vf,0) = vf->object ; /* number of objects in file = length of index */
      vgpWriteLine (vf, '&', vf->lineInfo['&']->buffer) ;

      /* finally we must write the termination line for the footer, then the last pointer line */
      fprintf (vf->f, "^\n") ;
      { errno = 0 ; off_t lineOff = ftello (vf->f) ; 
	vf->isBinary = FALSE ;
	vgpInt(vf,0) = footOff ; vgpWriteLine (vf, '-', 0) ;
	vf->isBinary = TRUE ;
	lineOff = ftello(vf->f) - lineOff ; if (errno) die ("failed final line ftells") ;
	while (lineOff < 31) { fputc (' ', vf->f) ; ++lineOff ; } fputc ('\n', vf->f) ;
      }
    }

  vgpFileDestroy (vf) ;
}

BOOL vgpGotoObject (VgpFile *vf, I64 i)
{
  errno = 0 ;
  if (vf && vf->isIndex && i >= 0 && i < vf->lineInfo[vf->objectType]->expectCount
      && !fseek (vf->f, ((I64*)vf->lineInfo['&']->buffer)[i], SEEK_SET))
    { vf->object = i ;
      return TRUE ;
    }
  else
    return FALSE ;
}

/***************** top level write functions *******************/

VgpFile *vgpFileOpenWriteNew (const char *path, FileType fileType, SubType subType, BOOL isBinary)
{
  int i, j, k ;

  if (subType && subPrimary[subType] != fileType)
    die ("subtype %s is not secondary for filetype %s",
	 subTypeName[subType], fileTypeName[fileType]) ;
  
  VgpFile *vf = vgpFileCreate (fileType) ;
  if (!strcmp (path, "-")) vf->f = stdout ;
  else if (!(vf->f = fopen (path, "w"))) { free(vf) ; return 0 ; }
  vf->isWrite = TRUE ;
  vf->subType = subType ;
  vf->isBinary = isBinary ;

  vf->codecBufSize = MAX_FIELD*sizeof(Field) + 1 ;
  vf->codecBuf = new (vf->codecBufSize, char) ; 

  return vf ;
}

VgpFile *vgpFileOpenWriteFrom (const char *path, VgpFile *vfIn, BOOL isBinary)
{
  VgpFile *vf = vgpFileOpenWriteNew (path, vfIn->fileType, vfIn->subType, isBinary) ;
  vgpInheritProvenance (vf, vfIn) ;
  vgpInheritReference (vf, vfIn) ;
  vgpInheritDeferred (vf, vfIn) ;
  int i ;
  for (i = 0 ; i < 128 ; ++i)
    if (vfIn->lineInfo[i])
      { LineInfo *li = vf->lineInfo[i], *liIn = vfIn->lineInfo[i] ;
	li->expectCount = liIn->expectCount ;
	li->expectMax = liIn->expectMax ;
	li->expectTotal = liIn->expectTotal ;
	li->expectGroupCount = liIn->expectGroupCount ;
	li->expectGroupTotal = liIn->expectGroupTotal ;
      }

  /* allocate codec buffer - always allocate enough to handle fields of all line types */
  I64 size = vf->codecBufSize ;
  for (i = 0 ; i < 128 ; ++i)
    if (vf->lineInfo[i] && vf->lineInfo[i]->listCodec &&
	size <= vf->lineInfo[i]->expectMax * vf->lineInfo[i]->listByteSize)
      size = vf->lineInfo[i]->expectMax * vf->lineInfo[i]->listByteSize + 1 ;
  if (size > vf->codecBufSize)
    { free (vf->codecBuf) ; vf->codecBufSize = size ; vf->codecBuf = new (size, char) ; }
  
  return vf ;
}

static void writeStringList (VgpFile *vf, char t, char *buf, int len)
{
  int j ;
  I64 totLen = 0 ;
  for (j = 0 ; j < len ; ++j)
    { I64 sLen = strlen ((char*)buf) ;
      totLen += sLen ;
      fprintf (vf->f, " %lld %s", sLen, (char*)buf) ;
      buf += sLen + 1 ;
    }
  LineInfo *li = vf->lineInfo[t] ;
  li->total += totLen - len ;
  if (li->max < totLen) li->max = totLen ;
}


static void writeBinaryIntList (VgpFile *vf, LineInfo *li, char *buf, int len)
{
  I64 nBits ;
  I64 i, j, mask = 0, *x = (I64*)buf ;

  for (i = 0 ; i < len ; ++i) mask |= *x++ ; /* find how many top bytes can be skipped */
  for (i = 0 ; i < 64 ; i += 8) if (!(mask >> (i*8))) break ;
  i >>= 3 ; j = 8 - i ; 	/* number of data bytes, number of 0 bytes */

  I64 listSize = i * len ;
  if (j)
    { 
    }
  //  else
  //    buffer = buf ;
  
  //  nBits = vcEncode (li->listCodec, listSize, buffer, vf->codecBuf) ;
  nBits |= (i << 60) ;		/* use top bits to store how many data bytes */
  if (fwrite (&nBits, sizeof(I64), 1, vf->f) != 1) die ("list nBits") ;
  if (fwrite (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1) die ("compressed list") ;
}

void vgpWriteLine (VgpFile *vf, char t, void *buf)
/* process is to fill fields by assigning to macros, then call - list contents are in buf */
/* NB adds '\n' before writing line not after, so user fprintf() can add extra material */
/* first call will write initial header, allowing space for count sizes to expand on close */
{
  I64 i, j, len ;
  LineInfo *li = vf->lineInfo[t] ;
  if (!li) die ("vgpWriteLine error: line type %c not present in file spec %s ",
		t, fileTypeName[vf->fileType]) ;

  if (!vf->isLastLineBinary) fputc ('\n', vf->f) ; /* terminate previous ascii line */

  ++vf->line ;
  ++li->count ;
  if (t == vf->groupType) updateGroupCount (vf, TRUE) ;
  if (t == vf->objectType)
    { LineInfo *lindex = vf->lineInfo['&'] ;
      if (vf->object >= lindex->bufSize)
	{ I64 newSize = (lindex->bufSize + (1 << 16)) << 1 ;
	  I64* a = lindex->buffer ;
	  lindex->buffer = new(newSize, I64) ;
	  memcpy (lindex->buffer, a, lindex->bufSize*sizeof(I64)) ;
	  free (a) ;
	  lindex->bufSize = newSize ;
	}
      ((I64*)(lindex->buffer))[vf->object] = ftello (vf->f) ; /* must come after writing newline */
      ++vf->object ;
    }

  //  printf ("line %lld type %d = %c\n", vf->line, t, t) ;
  
  if (!vf->isBinary)		/* ASCII - write field by field */
    { fputc (t, vf->f) ;
      for (i = 0 ; i < MAX_FIELD && li->fieldSpec[i] ; ++i)
	switch (li->fieldSpec[i])
	  {
	  case INT: fprintf (vf->f, " %lld", vf->field[i].i) ; break ;
	  case REAL: fprintf (vf->f, " %f", vf->field[i].r) ; break ;
	  case CHAR: fprintf (vf->f, " %c", vf->field[i].c) ; break ;
	  case STRING: case INT_LIST: case REAL_LIST: case STRING_LIST:
	    len = vf->field[i].len ; fprintf (vf->f, " %lld", len) ;
	    li->total += len ; if (len > li->max) li->max = len ;
	    if (li->fieldSpec[i] == STRING)
	      { if (len > INT_MAX)
		  die ("write problem: string length %lld > current max %d", len, INT_MAX) ;
		fprintf (vf->f, " %.*s", (int)len, (char*)buf) ;
	      }
	    else if (li->fieldSpec[i] == INT_LIST)
	      for (j = 0 ; j < len ; ++j)
		{ fprintf (vf->f, " %lld", *((I64*)buf)) ; buf += sizeof(I64) ; }
	    else if (li->fieldSpec[i] == REAL_LIST)
	      for (j = 0 ; j < len ; ++j)
		{ fprintf (vf->f, " %f", *((double*)buf)) ; buf += sizeof(double) ; }
	    else		/* STRING_LIST */
	      writeStringList (vf, t, buf, len) ;
	  }
      vf->isLastLineBinary = FALSE ;
    }
  else				/* binary - block write and optionally compress */
    { unsigned char x = li->binaryTypePack ;
      I64 fieldsSize = li->nField*sizeof(Field), nBits ;
      if (li->isUseFieldCodec)
	{ nBits = vcEncode (li->fieldCodec, fieldsSize, (char*)vf->field, vf->codecBuf) ;
	  if (nBits < 256) x |= 0x01 ;
	}
      if (li->isUseListCodec) x |= 0x02 ; /* need to do this before writing out x */
      fputc (x, vf->f) ;

      int nField = li->nField ;
      if (x & 0x01)		/* fields compress */
	{ unsigned char cBits = nBits ; fputc (cBits, vf->f) ;
	  if (fwrite (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1) die ("compressed fields") ;
	}
      else
	{ if (fwrite (vf->field, fieldsSize, 1, vf->f) != 1)
	    die ("write fields: t %c, nField %d, fieldsSize %lld", t, nField, fieldsSize) ;
	  if (li->fieldCodec)
	    { vcAddToTable (li->fieldCodec, fieldsSize, (char*)vf->field) ;
	      if (li->count*fieldsSize > vf->codecTrainingSize)
		{ vcCreateCodec (li->fieldCodec, 1) ; /* 1 makes escape code for unseen chars */
		  li->isUseFieldCodec = TRUE ;
		}
	    }
	}

      int listField = li->listField ;
      if (listField--)		/* NB subtract one to set actual field number */
	{ I64 listLen = vgpLen(vf, listField) ;
	  if (listLen)
	    { li->total += listLen ; if (listLen > li->max) li->max = listLen ;
	      I64 listSize = listLen * li->listByteSize ;
	      if (x & 0x2)		/* list is compressed */
		{ if (listSize >= vf->codecBufSize)
		    { free (vf->codecBuf) ; vf->codecBufSize = listSize + 1 ;
		      vf->codecBuf = new (vf->codecBufSize, char) ;
		    }
		  //		  if (li->fieldSpec[listField] == INT_LIST)
		  //		    writeBinaryIntList (vf, li, buf, listLen) ;
		  else
		    { nBits = vcEncode (li->listCodec, listSize, buf, vf->codecBuf) ;
		      if (fwrite (&nBits, sizeof(I64), 1, vf->f) != 1) die ("list nBits") ;
		      if (fwrite (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1)
			die ("compressed list") ;
		    }
		}
	      else if (li->fieldSpec[listField] == STRING_LIST) /* handle as ASCII */
		writeStringList (vf, t, buf, vgpLen(vf,listField)) ;
	      else
		{ if (fwrite (buf, listSize, 1, vf->f) != 1) die ("write list") ;
		  if (li->listCodec)
		    { vcAddToTable (li->listCodec, listSize, buf) ;
		      if (li->total > vf->codecTrainingSize)
			{ vcCreateCodec (li->listCodec, 1) ; /* 1 makes escape code for unseen chars */
			  li->isUseListCodec = TRUE ;
			}
		    }
		}
	    }
	}
      vf->isLastLineBinary = TRUE ;
    }
}

void vgpWriteHeader (VgpFile *vf)
{
  int i, j ;
  char groupChar = 0 ;
 
  fprintf (vf->f, "1 %lu %s %lld %lld",
	   strlen(fileTypeName[vf->fileType]), fileTypeName[vf->fileType], vf->major, vf->minor) ;
  ++vf->line ;
  if (vf->subType)
    { fprintf (vf->f, "\n2 %lu %s", strlen(subTypeName[vf->subType]), subTypeName[vf->subType]) ;
      ++vf->line ;
    }

  Reference *r = vf->reference ;
  for (i = vf->lineInfo['<']->count ; i-- ; ++r)
    { fprintf (vf->f, "\n< %lu %s %lld", strlen(r->filename), r->filename, r->count) ;
      ++vf->line ;
    }
  r = vf->deferred ;
  for (i = vf->lineInfo['>']->count ; i-- ; ++r)
    { fprintf (vf->f, "\n> %lu %s", strlen(r->filename), r->filename) ; ++vf->line ; }
  
  Provenance *p = vf->provenance ; 
  for (i = vf->lineInfo['!']->count ; i-- ; ++p)
    { fprintf (vf->f, "\n! %lu %s %lu %s %lu %s %lu %s",
	       strlen (p->program), p->program, strlen (p->version), p->version,
	       strlen (p->command), p->command, strlen (p->date), p->date) ;
      ++vf->line ;
    }

  if (vf->isBinary)		     /* defer writing rest of header */
    { fprintf (vf->f, "\n$") ; ++vf->line ; }
  else				     /* write counts */
    { char g = vf->groupType ; if (!g) g = 'Z' ; /* if no group then upper case only */
      for (i = 'A' ; i <= g ; ++i)
	if (vf->lineInfo[i] && vf->lineInfo[i]->expectCount)
	  { LineInfo *li = vf->lineInfo[i] ;
	    fprintf (vf->f, "\n# %c %lld", i, li->expectCount) ; ++vf->line ;
	    if (li->expectMax)
	      { fprintf (vf->f, "\n@ %c %lld", i, li->expectMax) ; ++vf->line ; }
	    if (li->expectTotal)
	      { fprintf (vf->f, "\n+ %c %lld", i, li->expectTotal) ; ++vf->line ; }
	    if (li->expectGroupCount)
	      { fprintf (vf->f, "\n%% %c # %c %lld", g, i, li->expectGroupCount) ; ++vf->line ; }
	    if (li->expectGroupTotal)
	      { fprintf (vf->f, "\n%% %c + %c %lld", g, i, li->expectGroupTotal) ; ++vf->line ; }
	  }
    }
  fflush (vf->f) ;

  vf->isHeader = TRUE ;
}

/************ the next set of functions manage provenance and reference lines ***************/

static BOOL addProvenance (VgpFile *vf, Provenance *from, int n)
{
  if (!n) return FALSE ;
  if (vf->isHeader) die ("can't addProvenance after writing header") ;
  I64 *count = &vf->lineInfo['!']->count ;
  Provenance *p = new (*count + n, Provenance) ;
  if (vf->provenance) memcpy (p, vf->provenance, *count*sizeof(Provenance)) ;
  memcpy (p + *count, from, n*sizeof(Provenance)) ;
  if (vf->provenance) free (vf->provenance) ;
  vf->provenance = p ;
  *count += n ;
  return TRUE ;
}

BOOL vgpInheritProvenance (VgpFile *vf, VgpFile *source)
{ return addProvenance (vf, source->provenance, source->lineInfo['!']->count) ; }

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
  I64 *count = isDeferred ? &vf->lineInfo['>']->count : &vf->lineInfo['<']->count ;
  Reference *r = new (*count + n, Reference) ;
  if (*target) memcpy (r, *target, *count*sizeof(Reference)) ;
  memcpy (r + *count, from, n*sizeof(Reference)) ;
  if (*target) free (*target) ;
  *target = r ;
  *count += n ;
  return TRUE ;
}

BOOL vgpInheritReference (VgpFile *vf, VgpFile *source) /* as for provenance */
{ return addReference (vf, source->reference, source->lineInfo['<']->count, FALSE) ; }

BOOL vgpAddReference (VgpFile *vf, char *filename, I64 count)
{ Reference ref ; ref.filename = filename ; ref.count = count ;
  return addReference (vf, &ref, 1, FALSE) ;
}

BOOL vgpInheritDeferred (VgpFile *vf, VgpFile *source)
{ return addReference (vf, source->deferred, source->lineInfo['>']->count, TRUE) ; }

BOOL vgpAddDeferred (VgpFile *vf, char *filename)
{ Reference ref ; ref.filename = filename ;
  return addReference (vf, &ref, 1, TRUE) ;
}

static inline void eatWhite (VgpFile *vf)
{ char x = vfGetc (vf) ;
  if (x != ' ' && x != '\t') parseError (vf, "failed to find expected whitespace") ;
}

static inline char readChar (VgpFile *vf) { eatWhite (vf) ; return vfGetc (vf) ; }

static inline char* readBuf (VgpFile *vf)
{ eatWhite (vf) ;
  char *endBuf = vf->numberBuf + 32 ;
  char *cp = vf->numberBuf ; --cp ;
  while ((++cp < endBuf) && (*cp = vfGetc (vf)))
    if (*cp == ' ' || *cp == '\t' || *cp == '\n' || *cp == EOF) break ;
  if (cp == endBuf) { *--cp = 0 ; parseError (vf, "overlong item %s", vf->numberBuf) ; }
  ungetc (*cp, vf->f) ; --vf->linePos ; *cp = 0 ;
  return vf->numberBuf ;
}

static inline I64 readInt (VgpFile *vf)
{ char *ep, *buf = readBuf (vf) ;
  BOOL isMinus = (*buf == '-') ; if (isMinus) ++buf ;
  if (!*buf) parseError (vf, "empty int field") ;
  I64 x = strtoll (buf, &ep, 10) ; if (*ep) parseError (vf, "bad int") ;
  if (isMinus) x = -x ;
  return x ;
}

static inline double readReal (VgpFile *vf)
{ char *ep, *buf = readBuf (vf) ;
  if (!*buf) parseError (vf, "empty real field") ;
  double x = strtod (buf, &ep) ; if (*ep) parseError (vf, "bad real line") ;
  return x ;
}

static inline void readString (VgpFile *vf, char *buf, I64 n)
{ eatWhite (vf) ;
  char *cp = buf ; --cp ;
  while (n-- && (*++cp = vfGetc (vf)))
    if (*cp == '\n' || *cp == EOF) break ;
  if (++n) { parseError (vf, "line too short %d", buf) ; }
  *++cp = 0 ;
}

static inline void readFlush (VgpFile *vf) /* reads to the end of the line */
{ char x ;
  while ((x = getc (vf->f)) && x != '\n')
    if (x == EOF) parseError (vf, "premature end of file") ;
}

#include <stdarg.h>

void parseError (VgpFile *vf, char *format, ...)
{
  va_list args ;

  fprintf (stderr, "PARSE ERROR ") ;

  va_start (args, format) ;
  vfprintf (stderr, format, args) ;
  va_end (args) ;

  vf->lineBuf[vf->linePos] = 0 ; /* terminate */
  fprintf (stderr, ", line %lld: %s\n", vf->line, vf->lineBuf) ;

  exit (-1) ;
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

long nAlloc = 0 ;
long totalAlloc = 0 ;

void *myalloc (size_t size)
{
  void *p = (void*) malloc (size) ;
  if (!p) die ("myalloc failure requesting %d bytes", size) ;
  ++nAlloc ;
  totalAlloc += size ;
  return p ;
}

void *mycalloc (size_t number, size_t size)
{
  void *p = (void*) calloc (number, size) ;
  if (!p) die ("mycalloc failure requesting %d objects of size %d", number, size) ;
  ++nAlloc ;
  totalAlloc += size*number ;
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

/* NB this is not threadsafe - do not call inside threads */

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
      fprintf (f, "\tnnew\t%li", nAlloc) ;   
      fprintf (f, "\ttotnew\t%li", totalAlloc) ;   
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
