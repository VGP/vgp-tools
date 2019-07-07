/*  File: vgprd.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description: implementation for vgprd.h
 * Exported functions:
 * HISTORY:
 * Last edited: Jul  7 22:37 2019 (rd109)
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

static void parseError (VgpFile *vf, char *format, ...) ;
static char inline vfGetc (VgpFile *vf) { char c = getc (vf->f) ; if (vf->linePos < 127) vf->lineBuf[vf->linePos++] = c ; return c ; }
static inline void updateGroupCount (VgpFile *vf, BOOL isGroupLine) ;
  
void vgpClose (VgpFile *vf) /* automatically rewrites header if allowed when writing */
{
  updateGroupCount (vf, FALSE) ; /* count last group */

  if (vf->isWrite)
    fputc ('\n', vf->f) ;	/* end of file if ascii, end of data marker if binary */

  if (vf->isWrite && vf->isBinary) /* write the footer */
    { int i, j ;
      errno = 0 ; off_t footOff = ftello (vf->f) ; if (errno) die ("failed footer ftell") ;
      
      /*  first the content information, which has to be ascii because we don't have codecs yet */
      for (i = 'A' ; i <= 'Z' ; ++i)	/* don't write metadata for header symbols */
	if (vf->count[i])
	  { fprintf (vf->f, "# %c %lld\n", i, vf->count[i]) ;
	    if (vf->max[i]) fprintf (vf->f, "@ %c %lld\n", i, vf->max[i]) ;
	    if (vf->total[i]) fprintf (vf->f, "+ %c %lld\n", i, vf->total[i]) ;
	  }
      char g = vf->spec->groupType ;
      if (g && vf->count[g])
	{ fprintf (vf->f, "# %c %lld\n", g, vf->count[g]) ;
	  if (vf->max[g]) fprintf (vf->f, "@ %c %lld\n", g, vf->max[g]) ;
	  if (vf->total[g]) fprintf (vf->f, "+ %c %lld\n", g, vf->total[g]) ;
	  for (j = 'A' ; j <= 'Z' ; ++j)
	    { if (vf->groupCount[j]) fprintf (vf->f, "%% %c # %c %lld\n",g,j,vf->groupCount[j]) ;
	      if (vf->groupTotal[j]) fprintf (vf->f, "%% %c + %c %lld\n",g,j,vf->groupTotal[j]) ;
	    }
	}

      /* next the codecs, which need to come after the counts so we know buffer sizes */
      char *buf = new (vcMaxSerialSize (), char) ;
      for (i = 0 ; i < 128 ; ++i)
	{ if (vf->isUseFieldCodec[i])
	    { vgpChar(vf,0) = vcSerialize (vf->fieldCodec[i], buf) ;
	      vgpWriteLine (vf, 1, buf) ;
	    }
	  if (vf->isUseListCodec[i] && vf->listCodec[i] != DNAcodec)
	    { vgpChar(vf,0) = vcSerialize (vf->listCodec[i], buf) ;
	      vgpWriteLine (vf, 2, buf) ;
	    }
	}
      free (buf) ;

      /* then the index, which might be compressed by a codec */
      vgpLen(vf,0) = vf->object ; /* number of objects in file = length of index */
      vgpWriteLine (vf, '&', vf->buffer['&']) ;

      /* finally we must write the end of footer, then the last pointer line */
      fprintf (vf->f, "^\n") ;
      { errno = 0 ; off_t lineOff = ftello (vf->f) ; 
	FileSpecification *keepSpec = vf->spec ; vf->spec = &lastLineSpec ; vf->isBinary = FALSE ;
	vgpInt(vf,0) = footOff ; vgpWriteLine (vf, '-', 0) ;
	vf->spec = keepSpec ; vf->isBinary = TRUE ;
	lineOff = ftello(vf->f) - lineOff ; if (errno) die ("failed final line ftells") ;
	while (lineOff < 31) { fputc (' ', vf->f) ; ++lineOff ; } fputc ('\n', vf->f) ;
      }
    }

  fclose (vf->f) ;

  if (vf->provenance) free (vf->provenance) ; /* can't know who owns subfields to free them */
  if (vf->reference) free (vf->reference) ;
  if (vf->deferred) free (vf->deferred) ;
  if (vf->codecBuf) free (vf->codecBuf) ;
  int i ;
  for (i = 0 ; i < 128 ; ++i)
    { if (vf->buffer[i] && !vf->isUserBuf[i]) free (vf->buffer[i]) ;
      if (vf->fieldCodec[i]) vcDestroy (vf->fieldCodec[i]) ;
      if (vf->listCodec[i]) vcDestroy (vf->listCodec[i]) ;
    }
  free (vf) ;
}

VgpFile *vgpFileOpenRead (const char *path, FileType fileType)
/* opens file and reads header if it exists
   if there is no header then fileType must be given
   if there is a header and fileType is non-zero then it must match 
   if there is $ line (file is binary) then read footer including decompressors and index
*/
{
  int i ;

  FileSpecification *formatSpec = vgpDefineFormat () ; /* initialise from included file */
  
  VgpFile *vf = new0 (1, VgpFile) ;
  if (!(vf->f = fzopen (path, "r"))) { free(vf) ; return 0 ; }
  vf->spec = &firstLineSpec ;	/* this contains the first line specification, only */

  /* read header and (optionally) footer
     recognise end of header by peeking at the first char to check if alphabetic 
  */
  while (TRUE)
    { unsigned char peek = getc(vf->f) ; if (feof (vf->f)) break ; /* loop exit at end of file */
      ungetc(peek, vf->f) ;
      if (peek & 0x80) peek = vf->spec->binaryTypeUnpack[peek] ;
      if ((peek >= 'A' && peek <= 'Z') || (peek >= 'a' && peek <= 'z'))
	{ if (vf->spec == &firstLineSpec) /* have not read any lines - missing header */
	    { if (!fileType)	/* have to define fileType in function call or file */
		{ fprintf (stderr, "VGP file open error: fileType not defined in file or code\n") ;
		  fclose (vf->f) ; free (vf) ; return 0 ;
		}
	      vf->fileType = fileType ; vf->spec = &formatSpec[fileType] ;
	      vf->major = vf->spec->major ; vf->minor = vf->spec->minor ;
	    }
	  break ;		/* loop exit at standard data line */
	}
      else			/* this is a header line */
	{ if (peek == '!')	/* hack to insert a count of 4 for STRING_LIST */
	    { getc(vf->f) ; ungetc('4',vf->f) ; ungetc(' ',vf->f) ; ungetc('!',vf->f) ; }
	  vgpReadLine (vf) ;	/* can't fail because we checked file end already */
	  switch (vf->lineType)
	    {
	    case '1':
	      { char *s = vgpString (vf) ; FileType t = 0 ;
		for (i = 1 ; i < MAX_FILE ; ++i) if (!strcmp (s, fileTypeName[i])) t = i ;
		if (!t) die ("unknown primary fileType %s in header line 1", s) ;
		if (fileType && t != fileType)
		  die ("primary fileType mismatch file %s != %s", s, fileTypeName[fileType]) ;
		vf->fileType = t ; vf->spec = &formatSpec[t] ;
	      }
	      if ((vf->major = vgpInt(vf,1)) != vf->spec->major)
		die ("major version file %d != code %d", vf->major, vf->spec->major) ;
	      if ((vf->minor = vgpInt(vf,2)) > vf->spec->minor)
		die ("minor version file %d > code %d", vf->minor, vf->spec->minor) ;
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
	      { char lt = vgpChar(vf,0) ;
		if (!vf->spec->line[lt]) parseError (vf, "unknown line type %c", lt) ;
		if (vf->lineType == '#')
		  { vf->expectCount[lt] = vgpInt(vf,1) ;
		    if (lt == vf->spec->objectType)
		      { vf->bufSize['&'] = vf->expectCount[lt] ; /* allow for string terminators */
			vf->buffer['&'] = new(vf->bufSize['&']*sizeof(I64), char) ;
		      }
		  }
		else if (vf->lineType == '@')
		  { vf->expectMax[lt] = vgpInt(vf,1) ;
		    vf->bufSize[lt] = vf->expectMax[lt] + 1 ; /* allow for string terminators */
		    vf->buffer[lt] = new(vf->bufSize[lt] * vf->spec->line[lt]->listByteSize, char) ;
		  }
		else if (vf->lineType == '+') vf->expectTotal[lt] = vgpInt(vf,1) ;
		else /* must be % */
		  { lt = vgpChar(vf,2) ;
		    if (!vf->spec->line[lt]) parseError (vf, "unknown line type %c", lt) ;
		    if (vgpChar(vf,1) == '#') vf->expectGroupCount[lt] = vgpInt(vf,3) ;
		    else if (vgpChar(vf,1) == '+') vf->expectGroupTotal[lt] = vgpInt(vf,3) ;
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
		--vf->count['!'] ; /* to avoid double counting */
		vgpAddProvenance (vf, prog, version, command, dateTime) ;
	      }
	      break ;
	    case '<':
	      --vf->count['<'] ; /* to avoid double counting */
	      vgpAddReference (vf, vgpString(vf), vgpInt(vf,1)) ;
	      break ;
	    case '>':
	      --vf->count['>'] ; /* to avoid double counting */
	      vgpAddDeferred (vf, vgpString(vf)) ;
	      break ;
	    /* below here are binary file header types - require expectCount/expectMax first */
	    case '$':  /* read footer - goto end, find offset to start of footer and go there */
	      errno = 0 ; off_t startOff = ftello (vf->f) ; fseek (vf->f, -32, SEEK_END) ;
	      if (errno) die ("can't seek to final line") ;
	      FileSpecification *fsKeep = vf->spec ; vf->spec = &lastLineSpec ;
	      vgpReadLine (vf) ; vf->spec = fsKeep ;
	      fseeko (vf->f, (off_t) vgpInt(vf,0), SEEK_SET) ;
	      if (errno) die ("can't seek to start of footer") ;
	      break ;
	    case '^': /* end of footer - return to where we left header */
	      errno = 0 ; fseeko (vf->f, startOff, SEEK_SET) ; if (errno) die ("can't seek back") ;
	      break ;
	    case '&': /* index - nothing to do - the index is buffer['&']!*/
	      break ;
	    case 1:
	      vf->fieldCodec[vgpChar(vf,0)] = vcDeserialize (vgpString(vf)) ;
	      break ;
	    case 2:
	      vf->listCodec[vgpChar(vf,0)] = vcDeserialize (vgpString(vf)) ;
	      break ;
	    default: parseError (vf, "unknown header line type %c", vf->lineType) ;
	    }
	}
    }

  /* allocate codec buffer - always allocate enough to handle fields */
  I64 size = MAX_FIELD*sizeof(Field) ;
  for (i = 0 ; i < 128 ; ++i)
    if (vf->listCodec[i] &&
	size < vf->expectMax[vf->lineType]*vf->spec->line[vf->lineType]->listByteSize)
      size = vf->expectMax[vf->lineType]*vf->spec->line[vf->lineType]->listByteSize ;
  vf->codecBuf = new (++size, char) ; /* add one for worst case codec usage */
  
  return vf ;
}

static inline char readChar (VgpFile *vf) ;
static inline I64 readInt (VgpFile *vf) ;
static inline double readReal (VgpFile *vf) ;
static inline void readString (VgpFile *vf, char *buf, I64 n) ;
static inline void readFlush (VgpFile *vf) ; /* reads to the end of the line */

static inline void confirmBufferSize (VgpFile *vf, char t, I64 size, I64 nStrings)
{
  vf->total[t] += size ;
  if (size > vf->max[t]) vf->max[t] = size ;
  size += nStrings ;		/* need to allocate space for terminal 0s */
  if (!vf->isUserBuf[t] && size > vf->bufSize[t])   /* expand buffer */
    { if (vf->buffer[t]) free (vf->buffer[t]) ;
      vf->buffer[t] = new (size*vf->spec->line[t]->listByteSize, char) ;
      vf->bufSize[t] = size ;
    }
}

static inline void updateGroupCount (VgpFile *vf, BOOL isGroupLine)
{ int i ;
  for (i = 'A' ; i <= 'Z' ; ++i)
    if (vf->count[i])
      { if (vf->inGroup && vf->count[i] - vf->gCount[i] > vf->groupCount[i])
	  vf->groupCount[i] = vf->count[i] - vf->gCount[i] ;
	if (vf->inGroup && vf->total[i] - vf->gTotal[i] > vf->groupTotal[i])
	  vf->groupTotal[i] = vf->total[i] - vf->gTotal[i] ;
	vf->gCount[i] = vf->count[i] ; vf->gTotal[i] = vf->total[i] ;
      }
  if (!isGroupLine)
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
  char *buf = (char*) vf->buffer[t] ;
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
  if (x & 0x80) { isAscii = FALSE ; t = vf->spec->binaryTypeUnpack[x] ; }
  else { isAscii = TRUE ; t = x ; }
  vf->lineType = t ;
  LineSpecification *ls = vf->spec->line[t] ;
  if (!ls) parseError (vf, "unknown line type %c", t) ;
  ++vf->count[t] ;
  if (t == vf->spec->objectType) ++vf->object ;
  if (t >= 'a' && t <= 'z') updateGroupCount (vf, TRUE) ;

  if (isAscii)			/* read field by field according to spec */
    { int i, j ;
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
		{ confirmBufferSize (vf, t, len, 1) ;
		  char *buf = (char*) vf->buffer[t] ;
		  readString (vf, buf, len) ;
		}
	      else if (ls->field[i] == INT_LIST)
		{ confirmBufferSize (vf, t, len, 0) ;
		  I64 *buf = (I64*) vf->buffer[t] ;
		  for (j = 0 ; j < len ; ++j) *buf++ = readInt (vf) ;
		}
	      else if (ls->field[i] == REAL_LIST)
		{ confirmBufferSize (vf, t, len, 0) ;
		  double *buf = (double*) vf->buffer[t] ;
		  for (j = 0 ; j < len ; ++j) *buf++ = readReal (vf) ;
		}
	      else              /* STRING_LIST - inefficient for now - also used for binary */
		readStringList (vf, t, len) ;
	    }
	  }
      readFlush (vf) ;
    }
  else				/* binary - block read fields and list, potentially compressed */
   { int nField = ls->nField ;
     if (nField)
       { if (x & 0x1) 		/* fields are compressed */
	   { int nBits = (uint8_t) getc (vf->f) ; /* NB only compress fields if under 255 bits */
	     if (fread (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1) die ("compressed fields") ;
	     vcDecode (vf->fieldCodec[t], nBits, vf->codecBuf, (char*) vf->field) ;
	   }
	 else
	   if (fread (vf->field, sizeof(Field), nField, vf->f) != nField) die ("fields") ;
       }
     
     int listField = ls->listField ;
     if (listField-- && vgpLen(vf, listField)) 	/* NB subtract one to set actual field number */
       { I64 listLen = vgpLen(vf, listField) ;
	 vf->total[t] += listLen ; if (listLen > vf->max[t]) vf->max[t] = listLen ;
	 if (x & 0x2)		/* list is compressed */
	   { I64 nBits ;
	     if (fread (&nBits, sizeof(I64), 1, vf->f) != 1) die ("list nBits") ;
	     if (fread (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1) die ("compressed list") ;
	     vcDecode (vf->listCodec[t], nBits, vf->codecBuf, vf->buffer[t]) ;
	   }
	 else if (ls->field[listField] == STRING_LIST) /* handle as ASCII */
	   readStringList (vf, t, vgpLen(vf,listField)) ;
	 else
	   { I64 size = vgpLen(vf,listField) * ls->listByteSize ;
	     if (fread (vf->buffer[t], size, 1, vf->f) != 1) die ("list") ;
	     if (ls->field[listField] == STRING)
	       ((char*)vf->buffer[t])[size] = 0 ; /* 0 terminate */
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
{ assert (vf && vf->spec && vf->spec->line[lineType]) ;
  if (buffer)
    { if (!vf->isUserBuf[lineType] && vf->buffer[lineType]) free (vf->buffer[lineType]) ;
      vf->buffer[lineType] = buffer ;
      vf->isUserBuf[lineType] = TRUE ;
    }
  else
    { if (vf->isUserBuf[lineType]) { vf->buffer[lineType] = 0 ; vf->bufSize[lineType] = 0 ; }
      vf->isUserBuf[lineType] = FALSE ;
    }
}

VgpFile *vgpFileOpenWriteNew (const char *path, FileType fileType, SubType subType, BOOL isBinary)
{
  int i, j, k ;

  if (subType && subPrimary[subType] != fileType)
    die ("subtype %s is not secondary for filetype %s", subTypeName[subType], fileTypeName[fileType]) ;
    
  FileSpecification *formatSpec = vgpDefineFormat () ; /* initialisation */
  
  VgpFile *vf = new0 (1, VgpFile) ;
  if (!strcmp (path, "-")) vf->f = stdout ;
  else if (!(vf->f = fopen (path, "w"))) { free(vf) ; return 0 ; }
  vf->isWrite = TRUE ;
  vf->fileType = fileType ;
  vf->subType = subType ;
  vf->spec = &formatSpec[fileType] ;
  vf->major = vf->spec->major ;
  vf->minor = vf->spec->minor ;
  vf->isBinary = isBinary ;
			
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
    { vf->expectCount[i] = vfIn->expectCount[i] ;
      vf->expectMax[i] = vfIn->expectMax[i] ;
      vf->expectTotal[i] = vfIn->expectTotal[i] ;
      vf->expectGroupCount[i] = vfIn->expectGroupCount[i] ;
      vf->expectGroupTotal[i] = vfIn->expectGroupTotal[i] ;
    }

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
  vf->total[t] += totLen - len ;
  if (totLen > vf->max[t]) vf->max[t] = totLen ;
}

void vgpWriteLine (VgpFile *vf, char t, void *buf)
/* process is to fill fields by assigning to macros, then call - list contents are in buf */
/* NB adds '\n' before writing line not after, so user fprintf() can add extra material */
/* first call will write initial header, allowing space for count sizes to expand on close */
{
  I64 i, j, len ;
  LineSpecification *ls = vf->spec->line[t] ;
  if (!ls) die ("vgpWriteLine error: line type %c not present in file spec %s ",
		t, fileTypeName[vf->fileType]) ;

  if (!vf->isLastLineBinary) fputc ('\n', vf->f) ; /* terminate previous ascii line */

  ++vf->line ;
  ++vf->count[t] ;
  if (t >= 'a' && t <= 'z') updateGroupCount (vf, TRUE) ;
  if (t == vf->spec->objectType)
    { if (vf->object >= vf->bufSize['&'])
	{ I64 newSize = (vf->bufSize['&'] + (1 << 16)) << 1 ;
	  I64* a = vf->buffer['&'] ;
	  vf->buffer['&'] = new(newSize, I64) ;
	  memcpy (vf->buffer['&'], a, vf->bufSize['&']) ;
	  vf->bufSize['&'] = newSize ;
	}
      ((I64*)vf->buffer['&'])[vf->object] = ftello (vf->f) ; /* must come after writing newline */
      ++vf->object ;
    }
  
  if (!vf->isBinary)		/* ASCII - write field by field */
    { fputc (t, vf->f) ;
      for (i = 0 ; i < MAX_FIELD ; ++i)
	{ if (!ls->field[i]) break ;
	  switch (ls->field[i])
	    {
	    case INT: fprintf (vf->f, " %lld", vf->field[i].i) ; break ;
	    case REAL: fprintf (vf->f, " %f", vf->field[i].r) ; break ;
	    case CHAR: fprintf (vf->f, " %c", vf->field[i].c) ; break ;
	    case STRING: case INT_LIST: case REAL_LIST: case STRING_LIST:
	      len = vf->field[i].len ; fprintf (vf->f, " %lld", len) ;
	      vf->total[t] += len ; if (len > vf->max[t]) vf->max[t] = len ;
	      if (ls->field[i] == STRING)
		{ if (len > INT_MAX)
		    die ("write problem: string length %lld > current max %d", len, INT_MAX) ;
		  fprintf (vf->f, " %.*s", (int)len, (char*)buf) ;
		}
	      else if (ls->field[i] == INT_LIST)
		for (j = 0 ; j < len ; ++j)
		  { fprintf (vf->f, " %lld", *((I64*)buf)) ; buf += sizeof(I64) ; }
	      else if (ls->field[i] == REAL_LIST)
		for (j = 0 ; j < len ; ++j)
		  { fprintf (vf->f, " %f", *((double*)buf)) ; buf += sizeof(double) ; }
	      else		/* STRING_LIST */
		writeStringList (vf, t, buf, len) ;
	    }
	}
      vf->isLastLineBinary = FALSE ;
    }
  else				/* binary - block write and optionally compress */
    { unsigned char x = vf->spec->binaryTypePack[t] ;
      I64 fieldsSize = ls->nField*sizeof(Field), nBits ;
      if (vf->isUseFieldCodec[t])
	{ nBits = vcEncode (vf->fieldCodec[t], fieldsSize, (char*)vf->field, vf->codecBuf) ;
	  if (nBits < 256) x |= 0x01 ;
	}
      if (t == 'S') die ("line %lld t %c isUseListCodec %d listCodec %lx\n",
			 vf->line, t, vf->isUseListCodec[t], vf->listCodec[t]) ;
      if (vf->isUseListCodec[t])
	{ x |= 0x02 ; /* need to do this before writing out x */
	  die ("line %lld using list codec for type %c\n", vf->line, t) ;
	}
      fputc (x, vf->f) ;

      int nField = ls->nField ;
      if (x & 0x01)		/* fields compress */
	{ unsigned char cBits = nBits ; fputc (cBits, vf->f) ;
	  if (fwrite (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1) die ("compressed fields") ;
	}
      else
	if (fwrite (vf->field, fieldsSize, 1, vf->f) != 1)
	  die ("write fields - t %c, nField %d, fieldsSize %lld", t, ls->nField, fieldsSize) ;

      int listField = ls->listField ;
      if (listField--)		/* NB subtract one to set actual field number */
	{ I64 listLen = vgpLen(vf, listField) ;
	  if (listLen)
	    { vf->total[t] += listLen ; if (listLen > vf->max[t]) vf->max[t] = listLen ;
	      I64 listSize = listLen * ls->listByteSize ;
	      if (x & 0x2)		/* list is compressed */
		{ nBits = vcEncode (vf->listCodec[t], listSize, buf, vf->codecBuf) ;
		  if (fwrite (&nBits, sizeof(I64), 1, vf->f) != 1) die ("list nBits") ;
		  if (fwrite (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1) die ("compressed list") ;
		}
	      else if (ls->field[listField] == STRING_LIST) /* handle as ASCII */
		writeStringList (vf, t, buf, vgpLen(vf,listField)) ;
	      else
		if (fwrite (buf, listSize, 1, vf->f) != 1) die ("write list") ;
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
  for (i = vf->count['<'] ; i-- ; ++r)
    { fprintf (vf->f, "\n< %lu %s %lld", strlen(r->filename), r->filename, r->count) ;
      ++vf->line ;
    }
  r = vf->deferred ;
  for (i = vf->count['>'] ; i-- ; ++r)
    { fprintf (vf->f, "\n> %lu %s", strlen(r->filename), r->filename) ; ++vf->line ; }
  
  Provenance *p = vf->provenance ; 
  for (i = vf->count['!'] ; i-- ; ++p)
    { fprintf (vf->f, "\n! %lu %s %lu %s %lu %s %lu %s",
	       strlen (p->program), p->program, strlen (p->version), p->version,
	       strlen (p->command), p->command, strlen (p->date), p->date) ;
      ++vf->line ;
    }

  if (vf->isBinary)		     /* defer writing rest of header */
    { fprintf (vf->f, "\n$") ; ++vf->line ; }
  else				     /* write counts */
    { for (i = 'A' ; i <= 'Z' ; ++i)	/* don't write metadata for header symbols */
	if (vf->expectCount[i])
	  { fprintf (vf->f, "\n# %c %lld", i, vf->expectCount[i]) ; ++vf->line ;
	    if (vf->expectMax[i])
	      { fprintf (vf->f, "\n@ %c %lld", i, vf->expectMax[i]) ; ++vf->line ; }
	    if (vf->expectTotal[i])
	      { fprintf (vf->f, "\n+ %c %lld", i, vf->expectTotal[i]) ; ++vf->line ; }
	  }
      char g = vf->spec->groupType ;
      if (g && vf->expectCount[g])
	{ fprintf (vf->f, "\n# %c %lld", g, vf->expectCount[g]) ; ++vf->line ;
	  if (vf->expectMax[g])
	    { fprintf (vf->f, "\n@ %c %lld", g, vf->expectMax[g]) ; ++vf->line ; }
	  if (vf->expectTotal[g])
	    { fprintf (vf->f, "\n+ %c %lld", g, vf->expectTotal[g]) ; ++vf->line ; }
	  for (j = 'A' ; j <= 'Z' ; ++j)
	    { if (vf->expectGroupCount[j])
		{ fprintf (vf->f, "\n%% %c # %c %lld", g, j, vf->expectGroupCount[j]) ; ++vf->line ; }
	      if (vf->expectGroupTotal[j])
		{ fprintf (vf->f, "\n%% %c + %c %lld", g, j, vf->expectGroupTotal[j]) ; ++vf->line ; }
	    }
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
  Provenance *p = new (vf->count['!'] + n, Provenance) ;
  I64 *count = &vf->count['!'] ;
  if (vf->provenance) memcpy (p, vf->provenance, *count*sizeof(Provenance)) ;
  memcpy (p + *count, from, n*sizeof(Provenance)) ;
  if (vf->provenance) free (vf->provenance) ;
  vf->provenance = p ;
  *count += n ;
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
  I64 *count = isDeferred ? &vf->count['>'] : &vf->count['<'] ;
  Reference *r = new (*count + n, Reference) ;
  if (*target) memcpy (r, *target, *count*sizeof(Reference)) ;
  memcpy (r + *count, from, n*sizeof(Reference)) ;
  if (*target) free (*target) ;
  *target = r ;
  *count += n ;
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
