/*****************************************************************************************
 *
 *  File: VGPlib.c
 *    implementation for VGPlib.h
 *
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University and Eugene Myers 2019-
 *
 * HISTORY:
 * Last edited: Feb  7 15:01 2020 (rd109)
 *   * Dec 27 09:46 2019 (gene): style edits + compactify code
 *   * Jul  8 04:28 2019 (rd109): refactored to use lineInfo[]
 *   * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *
 ****************************************************************************************/

#undef WITH_ZLIB         // this allows fzopen() to open gzip files - undef if problems

#include <assert.h>
#include <sys/errno.h>
#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <stdarg.h>
#include <time.h>
#include <ctype.h>
#include <pthread.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/uio.h>

#ifdef WITH_ZLIB
#include <zlib.h>
#endif

#include "VGPlib.h"

#include "compression.c" // directly include compression package - cleaner because only used here

static pthread_mutex_t mutexInit = PTHREAD_MUTEX_INITIALIZER;

/***********************************************************************************
 *
 *    VGP_FILE CREATION & DESTRUCTION
 *
 **********************************************************************************/

#include "VGPformat_1_0.h"   //  Inserts code for data line formats

static int listSize[8] = { 0, 0, 0, 0, 1, sizeof(I64), sizeof(double), 1 };

static VgpFile *vgpFileCreate(FileType fileType)
{ VgpFile   *vf;
  LineInfo **info;
  LineInfo  *line;
  int        j, k, n;
  U8         x;

  vf   = new0 (1, VgpFile);
  info = vf->lineInfo;

  // only binary format lines written in binary

  info[1]   = vgpDefineLine (CHAR, STRING, 0, 0, 0, 0); // field codec (binary only)
  info[2]   = vgpDefineLine (CHAR, STRING, 0, 0, 0, 0); // list codec (binary only)
  info['&'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0);  // object index
    info['&']->isIntListDiff = TRUE;
  info['*'] = vgpDefineLine (INT_LIST, 0, 0, 0, 0, 0);  // group index
    info['*']->isIntListDiff = TRUE;
  info['/'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0); // comment

  defineFormat (vf, fileType);      // from vgpformat_1_0.h

  // make the binaryType pack and unpack lookups
  //   EWM: previously for all line types, wasteful, only needed for those set above

  n = 0;
  for (j = 0; j < 128; ++j)
    if (info[j] != NULL)
      { x = (n++ << 2) | 0x80;
        info[j]->binaryTypePack = x;
        vf->binaryTypeUnpack[x  ] = j;
        vf->binaryTypeUnpack[x+1] = j;
        vf->binaryTypeUnpack[x+2] = j;
        vf->binaryTypeUnpack[x+3] = j;
      }
  if (n >= 32)
    die ("VGP spec error file type %s: too many line specs %d >= 32",
           fileTypeName[fileType], n);

  // define the header lines - always ASCII

  info['1'] = vgpDefineLine (STRING, INT, INT, 0, 0, 0);   // type, major, minor
  info['2'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0);       // subtype
  info['#'] = vgpDefineLine (CHAR, INT, 0, 0, 0, 0);       // linetype count
  info['@'] = vgpDefineLine (CHAR, INT, 0, 0, 0, 0);       // linetype max
  info['+'] = vgpDefineLine (CHAR, INT, 0, 0, 0, 0);       // linetype total
  info['%'] = vgpDefineLine (CHAR, CHAR, CHAR, INT, 0, 0); // group #/+ linetype value
  info['!'] = vgpDefineLine (STRING_LIST, 0, 0, 0, 0, 0);  // name version command date
  info['<'] = vgpDefineLine (STRING, INT, 0, 0, 0, 0);     // filename objectcount
  info['>'] = vgpDefineLine (STRING, 0, 0, 0, 0, 0);       // filename

  // binary format lines that are always ASCII

  info['$'] = vgpDefineLine (INT, 0, 0, 0, 0, 0);       // isBig - binary file
  info['^'] = vgpDefineLine (0, 0, 0, 0, 0, 0);         // end of footer - return to top

  // fill in the nField and listByteSize information, and perform format consistency checks

  for (j = 0; j < 128; ++j)
    { line = info[j];
      if (line != NULL)
        { line->nField = 0;
          for (k = 0; k < MAX_FIELD && line->fieldSpec[k]; k++)
            { if (listSize[line->fieldSpec[k]])
                { if (line->listByteSize > 0)
                    die ("VGP spec error file type %s: two list types in record %d",
                           fileTypeName[fileType], j);
                  line->listByteSize = listSize[line->fieldSpec[k]];
                  line->listField    = k+1;
                }
              line->nField += 1;
            }
          if (islower(j) && j != vf->groupType)
            die ("VGP spec error file type %s: extra group type %c not permitted",
                   fileTypeName[fileType], j);
        }
    }
  if (vf->objectType > 0 && info[(int) vf->objectType] == NULL)
    die ("VGP spec error file type %s: missing line spec for object type %c",
           fileTypeName[fileType], vf->objectType);
  if (vf->groupType > 0 && info[(int) vf->groupType] == NULL)
    die ("VGP spec error file type %s: missing line spec for group type %c",
           fileTypeName[fileType], vf->groupType);

  // setup for compression

  vf->codecTrainingSize = 100000;
  info[1]->bufSize = vcMaxSerialSize() + 1; // +1 for added but unused 0-terminator
  info[1]->buffer  = new (info[1]->bufSize, void);
  info[2]->bufSize = vcMaxSerialSize() + 1;
  info[2]->buffer  = new (info[2]->bufSize, void);

  // determine endian of machine

  { int   t = 1;
    char *b = (char *) (&t);
    vf->isBig = (*b == 0);
  }

  return (vf);
}

static void lineInfoDestroy(LineInfo *li)
{ if (li->buffer != NULL && ! li->isUserBuf) free (li->buffer);
  if (li->fieldCodec != NULL) vcDestroy (li->fieldCodec);
  if (li->listCodec  != NULL) vcDestroy (li->listCodec);
  free(li);
}

static void vgpFileDestroy(VgpFile *vf)
{ int       i, j;
  LineInfo *li, *lx;

  if (vf->share)
    { for (i = 0; i < 128 ; i++)
        { lx = vf->lineInfo[i];
          if ((i == '&' || i == '*') && ! vf->isWrite)
            continue;
          if (lx != NULL)
            { for (j = 1; j < vf->share; j++)
                { li = vf[j].lineInfo[i];
                  if (li->buffer != NULL && ! li->isUserBuf)
                    free (li->buffer);
                  free(li);
                }
            }
        }

      for (j = 1; j < vf->share; j++)
        { if (vf[j].provenance != NULL) free (vf[j].provenance);
          if (vf[j].reference  != NULL) free (vf[j].reference);
          if (vf[j].deferred   != NULL) free (vf[j].deferred);
          if (vf[j].codecBuf   != NULL) free (vf[j].codecBuf);
          if (vf[j].f          != NULL) fclose (vf[j].f);
        }
    }

  if (vf->provenance != NULL) free (vf->provenance);   //  EWM: Memory leak
  if (vf->reference  != NULL) free (vf->reference);    //  EWM: Memory leak
  if (vf->deferred   != NULL) free (vf->deferred);     //  EWM: Memory leak
  if (vf->codecBuf   != NULL) free (vf->codecBuf);
  if (vf->f != NULL && vf->f != stdout) fclose (vf->f);

  for (i = 0; i < 128 ; i++)
    if (vf->lineInfo[i] != NULL)
      lineInfoDestroy (vf->lineInfo[i]);

  free(vf);
}


/***********************************************************************************
 *
 *    ASCII PARSING UTILITIES: error reporting, lexical level
 *
 **********************************************************************************/

void parseError (VgpFile *vf, char *format, ...)
{ va_list args;

  fprintf (stderr, "PARSE ERROR ");

  va_start (args, format);
  vfprintf (stderr, format, args);
  va_end (args);

  vf->lineBuf[vf->linePos] = '\0';
  fprintf (stderr, ", line %lld: %s\n", vf->line, vf->lineBuf);

  exit (1);
}

static char inline vfGetc(VgpFile *vf)
{ char c = getc(vf->f);
  if (vf->linePos < 127)
    vf->lineBuf[vf->linePos++] = c;
  return (c);
}

static inline void eatWhite (VgpFile *vf)
{ char x = vfGetc(vf);
  if (x == ' ' || x == '\t')
    return;
  parseError (vf, "failed to find expected whitespace");
}

static inline char readChar(VgpFile *vf)
{ eatWhite(vf);
  return (vfGetc(vf));
}

static inline char *readBuf(VgpFile *vf)
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
  return (vf->numberBuf);
}

static inline I64 readInt(VgpFile *vf)
{ char *e, *b;
  I64   x;

  b = readBuf(vf);
  x = strtoll(b, &e, 10);
  if (e == b)
    parseError (vf, "empty int field");
  if (*e != '\0')
    parseError (vf, "bad int");
  return (x);
}

static inline double readReal(VgpFile *vf)
{ char  *e, *b;
  double x;

  b = readBuf(vf);
  x = strtod (b, &e);
  if (e == b)
    parseError (vf, "empty real field");
  if (*e != '\0')
    parseError (vf, "bad real line");
  return (x);
}

static inline void readString(VgpFile *vf, char *buf, I64 n)
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
    die ("failed to read %d byte string", n);
}

static inline void readFlush (VgpFile *vf) // reads to the end of the line and stores as comment
{ char       x;
  int        n = 0;
  LineInfo *li = vf->lineInfo['/'] ;

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

static inline void updateCountsAndBuffer (VgpFile *vf, char t, I64 size, I64 nStrings)
{ LineInfo *li;

  li = vf->lineInfo[(int) t];
  li->accum.total += size;
  if (size > li->accum.max)
    li->accum.max = size;
  size += nStrings;             // need to allocate space for terminal 0s
  if ( ! li->isUserBuf && size > li->bufSize)   // expand buffer
    { if (li->buffer != NULL) free (li->buffer);
      li->bufSize = size;
      li->buffer  = new (size*li->listByteSize, void);
    }
}

  //  Called when a new group starts or eof, accumulate group counts since last group start

static inline void updateGroupCount(VgpFile *vf, BOOL isGroupLine)
{ int       i;
  LineInfo *li;
  Counts   *ci;

  for (i = 'A'; i <= 'Z' ; i++)
    { li = vf->lineInfo[i];
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

static char *compactIntList (VgpFile *vf, LineInfo *li, I64 len, char *buf)
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

  k = li->listByteSize;
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
  vgpInt(vf,li->listField-1) |= (z << 56);
  
  return (li->buffer);
}

static void decompactIntList (VgpFile *vf, LineInfo *li, I64 len, char *buf)
{ I64   i, *x;
  int   ix, d, z, k;
  char *s, *t;

  ix = li->listField - 1 ;
  z   = (vgpInt (vf, ix) >> 56);

  if (z > 0)                      // decompacts in place
    { d = li->listByteSize - z;
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
 *  VGP_READ_LINE:
 *      Reads the next line and returns FALSE at end of file or on error. The line is
 *      parsed according to its linetype and contents accessed by macros that follow.
 *      The top bit of the first character determines whether the line is binary or ascii
 *
 **********************************************************************************/

  //  Read a string list, first into new allocs, then into sized line buffer.
  //    Annoyingly inefficient, but we don't use it very much.

static void readStringList(VgpFile *vf, char t, I64 len)
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

  buf = (char *) vf->lineInfo[(int) t]->buffer;
  for (j = 0; j < len ; ++j)
    { strcpy (buf, string[j]);
      buf += strlen(buf) + 1;
      free (string[j]);
    }
  free (string);
}

BOOL vgpReadLine (VgpFile *vf)
{ BOOL      isAscii;
  U8        x;
  char      t;
  LineInfo *li;

  if (vf->isWrite)
    die ("Trying to read a line from a file open for writing");
  if (vf->isFinal)
    die ("Cannot read more data after counts are finalized");

  vf->linePos = 0;                 // must come before first vfGetc()
  x = vfGetc (vf);                 // read first char
  if (feof (vf->f) || x == '\n')   // blank line (x=='\n') is end of records marker before footer
    return (FALSE);

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

  li = vf->lineInfo[(int) t];
  if (li == NULL)
    parseError (vf, "unknown line type %c(%d was %d) %d", t, t, x, vf->linePos);
  li->accum.count += 1;
  if (t == vf->objectType)
    vf->object += 1;
  if (t == vf->groupType)
    updateGroupCount (vf, TRUE);

  // fprintf (stderr, "reading line %lld type %c\n", vf->line, t) ;

  if (isAscii)           // read field by field according to ascii spec
    { int     i, j;
      I64    *ilst, len;
      double *rlst;

      for (i = 0; i < li->nField; i++)
        switch (li->fieldSpec[i])
        { case INT:
            vf->field[i].i = readInt (vf);
            break;
          case REAL:
            vf->field[i].r = readReal (vf);
            break;
          case CHAR:
            vf->field[i].c = readChar (vf);
            break;
          case STRING:
            len = readInt (vf);
            vf->field[i].len = len;
            updateCountsAndBuffer (vf, t, len, 1);
            readString (vf, (char*) li->buffer, len);
            break;
          case INT_LIST:
            len = readInt (vf);
            vf->field[i].len = len;
            updateCountsAndBuffer (vf, t, len, 0);
            ilst = (I64 *) li->buffer;
            for (j = 0; j < len; ++j)
              ilst[j] = readInt(vf);
            break;
          case REAL_LIST:
            len = readInt (vf);
            vf->field[i].len = len;
            updateCountsAndBuffer (vf, t, len, 0);
            rlst = (double *) li->buffer;
            for (j = 0; j < len; ++j)
              rlst[j] = readReal (vf);
            break;
          case STRING_LIST: // STRING_LIST - inefficient for now - also used for binary
            len = readInt (vf);
            vf->field[i].len = len;
            readStringList (vf, t, len);
            break;
        }
      readFlush (vf);
    }

  else        // binary - block read fields and list, potentially compressed
    { int nField, ix, nBits;
      I64 listLen, usedBytes, listSize;

      nField = li->nField;
      if (nField > 0)
        { if (x & 0x1)                       // fields are compressed
            { nBits = (U8) getc (vf->f);     // NB only fields under 255 bits are compressed
              if (fread (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1)
                die ("fail to read compressed fields");

              vcDecode (li->fieldCodec, nBits, vf->codecBuf, (char *) vf->field);
            }
          else
            if (fread (vf->field, sizeof(Field), nField, vf->f) != (unsigned long) nField)
              die ("fail to read fields");
        }

      if (t == vf->groupType)
        { I64 *groupIndex = (I64 *) vf->lineInfo['*']->buffer;
          vf->field[0].i = groupIndex[vf->group] - groupIndex[vf->group-1];
        }

      ix =  li->listField-1;
      if (ix >= 0)
        { listLen = vgpLen(vf);
          li->accum.total += listLen;
          if (listLen > li->accum.max)
            li->accum.max = listLen;

          if (listLen > 0)
            { if (li->fieldSpec[ix] == STRING_LIST) // handle as ASCII
                readStringList (vf, t, listLen);

              else if (x & 0x2)     // list is compressed
                { if (fread (&nBits, sizeof(I64), 1, vf->f) != 1)
                    die ("fail to read list nBits");
                  if (fread (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1)
                    die ("fail to read compressed list");
                  vcDecode (li->listCodec, nBits, vf->codecBuf, li->buffer);
                }

              else
                { usedBytes = li->listByteSize - (vgpInt(vf,ix) >> 56); // data bytes
                  listSize  = listLen * usedBytes;

                  if ((I64) fread (li->buffer, 1, listSize, vf->f) != listSize)
                    die ("list read %lld not %lld", x, listSize);
                }

              if (li->fieldSpec[ix] == INT_LIST)
                decompactIntList (vf, li, listLen, li->buffer);
            }

          if (li->fieldSpec[ix] == STRING)
            ((char *) li->buffer)[listLen] = '\0'; // 0 terminate
        }

      { U8 peek = getc(vf->f) ; // check if next line is a comment - if so then read it
	ungetc(peek, vf->f) ;
	if (peek & 0x80)
	  peek = vf->binaryTypeUnpack[peek];
	if (peek == '/') // a comment
	  { Field keepField0 = vf->field[0] ;
	    vgpReadLine (vf) ; // read comment line into vf->lineInfo['/']->buffer
	    vf->lineType = t ;
	    vf->field[0] = keepField0 ;
	  }
      }
    }

  return (t);
}

char *vgpReadComment (VgpFile *vf)
{
  char *comment = (char*)(vf->lineInfo['/']->buffer) ;
  if (comment && *comment != 0)
    return comment ;
  else
    return 0 ;
}

/***********************************************************************************
 *
 *   VGP_FILE_OPEN_READ:
 *     Opens file for reading and reads all lines of header if it exists.
 *     If there is no header then fileType must be given (i.e. non-zero),
 *     otherwise if fileType is non-zero then it must match the type in the header.
 *     If the file begins with a $-line then it is binary and the routine reads
 *     the footer that includes decompressors and object index.
 *
 **********************************************************************************/

VgpFile *vgpFileOpenRead (const char *path, FileType fileType, int nthreads)
{ VgpFile *vf;
  FILE    *f;
  int      curLine;
  int      i, j;
  off_t    startOff, footOff;

  if (strcmp (path, "-") == 0)
    f = stdin;
  else
    { f = fopen (path, "r");
      if (f == NULL)
        return (NULL);
    }

  // EWM: New logic, either has '1' line or doesn't.  If it does then get it now and
  //   determine or check fileType.  Avoid uncessary checks in core loop and especially
  //   creating and then destroying a VgpFile object of type 0.

  { U8       c;
    FileType t;
    int      major, minor, slen;
    char     name[10];

    c = getc(f);
    if (feof(f))
      die ("File is empty?");

    if (c != '1')
      { if (fileType == 0)
          { fprintf (stderr, "VGP file open error: fileType not defined in file or code\n");
            return (NULL);
          }
        ungetc(c,f);

        curLine = 0;
      }

    else
      { if (fscanf(f," %d %s %d %d",&slen,name,&major,&minor) != 4)  // EWM: better error report?
          { fprintf(stderr, "PARSE ERROR, line 1: 1-line syntax\n"); 
            exit (1);
          }
        while (getc(f) != '\n')
	  if (feof(f))
	    die ("end of file before end of line 1") ;

        if (major != MAJOR_NUMBER)
          die ("major version file %d != %d", major, MAJOR_NUMBER);
        if (minor != MINOR_NUMBER)
          die ("minor version file %d != %d", minor, MINOR_NUMBER);

        t = 0;
        for (i = 1; i < MAX_FILE; ++i)
          if (strcmp(name,fileTypeName[i]) == 0)
            t = i;
        if (t == 0)
          die ("unknown primary file type %s in header line 1", name);

        if (fileType > 0 && t != fileType)
          die("primary fileType mismatch file %s != %s", name, fileTypeName[fileType]);

        fileType = t;
        curLine  = 1;
     }
  }

  if (nthreads == 1)
    { vf = vgpFileCreate(fileType);
      vf->f = f;
      vf->line = curLine;
    }

  else
    { VgpFile *v;

      if (strcmp (path, "-") == 0)
        die ("Parallel input incompatible with stdin as input");

      vf = new (nthreads, VgpFile);

      for (i = 0; i < nthreads; i++)
        { v = vgpFileCreate(fileType);
          if (i > 0)
            { v->line  = 0;
              v->share = -i;
              f = fopen (path, "r");
            }
          else
            { v->line  = curLine;
              v->share = nthreads;
            }
          v->f = f;
          vf[i] = *v;
          free (v);
        }
    }

  // read header and (optionally) footer
  // recognise end of header by peeking at the first char to check if alphabetic 
 
  f = vf->f;
  vf->isCheckString = TRUE;   // always check strings while reading header
  while (TRUE)

    { U8 peek;

      peek = getc(f);
      if (feof(f))       // loop exit at end of file
        break;

      ungetc(peek, f);
      if (peek & 0x80)
        peek = vf->binaryTypeUnpack[peek];

      if (isalpha(peek))
        break;    // loop exit at standard data line

      else if (peek == '!')  // hack to insert a count of 4 for STRING_LIST
        { getc(f);
          ungetc('4',f);
          ungetc(' ',f);
          ungetc('!',f);
        }
      
      vf->f = f;
      vgpReadLine(vf);  // can't fail because we checked file eof already
      f = vf->f;

      switch (vf->lineType)

      { case '1':
          parseError(vf, "1 should be first line in header");
          break;

        case '2':
          { char   *s;
            SubType t;
	    
            t = 0;
            s = vgpString(vf);
            for (i = 1; i < MAX_SUB; ++i)
              if (!strcmp (s, subTypeName[i]))
                t = i;

            if (t == 0)
              parseError(vf, "unknown secondary subType %s", s);
            if (subPrimary[t] != vf->fileType)
              parseError(vf, "subtype %s not compatible with primary type %d",
                                s, fileTypeName[vf->fileType]);
            vf->subType = t;
            break;
          }

        case '#':
        case '@':
        case '+':
        case '%':
          { char      c = vgpChar(vf,0);
            LineInfo *li = vf->lineInfo[(int) c];

            if (li == NULL)
              parseError (vf, "unknown line type %c", c);
            switch (vf->lineType)
            { case '#':
                li->given.count = vgpInt(vf,1);
                if (c == vf->objectType && vf->isBinary) // allocate space for object index
                  { vf->lineInfo['&']->bufSize = li->given.count;
                    vf->lineInfo['&']->buffer  = new (li->given.count, I64);
                  }
                if (c == vf->groupType && vf->isBinary) // allocate space for group index
                  { vf->lineInfo['*']->bufSize = li->given.count+1; // +1 for end value
                    vf->lineInfo['*']->buffer  = new (vf->lineInfo['*']->bufSize, I64);
                  }
                break;
              case '@':
                li->given.max = vgpInt(vf,1);
                li->bufSize = li->given.max + 1; // allow for string terminators
                li->buffer = new (li->bufSize*li->listByteSize, void);
                break;
              case '+':
                li->given.total = vgpInt(vf,1);
                break;
              case '%':
                c  = vgpChar(vf,2);
                li = vf->lineInfo[(int) c];
                if (li == NULL)
                  parseError (vf, "unknown line type %c", c);
                c = vgpChar(vf,1);
                if (c == '#')
                  li->given.groupCount = vgpInt(vf,3);
                else if (c == '+')
                  li->given.groupTotal = vgpInt(vf,3);
                else
                  parseError (vf, "unrecognised symbol %c", c);
                break;
            }
            break;
          }

        case '!':     // NB need to copy the strings - small memory leak
          { char *prog     = vgpString(vf);
            char *version  = prog + strlen(prog) + 1;
            char *command  = version + strlen(version) + 1;
            char *dateTime = command + strlen(command) + 1;

            dateTime = strdup(dateTime);
            command  = strdup(command);
            version  = strdup(version);
            prog     = strdup(prog);

            vf->lineInfo['!']->accum.count -= 1; // to avoid double counting
            vgpAddProvenance (vf, prog, version, command, dateTime);
            break;
          }

        case '<':
          vf->lineInfo['<']->accum.count -= 1; // to avoid double counting
          vgpAddReference (vf, vgpString(vf), vgpInt(vf,1));
          break;

        case '>':
          vf->lineInfo['>']->accum.count -= 1; // to avoid double counting
          vgpAddDeferred (vf, vgpString(vf));
          break;

        // Below here are binary file header types - requires given.count/given.max first

        case '$':  // read footer - goto end, find offset to start of footer and go there
          if (vgpInt(vf,0) != vf->isBig)
            die ("endian mismatch - convert to ascii");
          vf->isBinary = TRUE;

          startOff = ftello (f);
          if (fseek (f, -sizeof(off_t), SEEK_END) != 0)
            die ("can't seek to final line");

          if (fread (&footOff, sizeof(off_t), 1, f) != 1)
            die ("can't read footer offset");

          if (fseeko (f, footOff, SEEK_SET) != 0)
            die ("can't seek to start of footer");
          break;

        case '^':    // end of footer - return to where we jumped from header
          if (fseeko (f, startOff, SEEK_SET) != 0)
            die ("can't seek back");
          break;

        case '&':
          vf->isIndexIn = TRUE;
          break;

        case '*':
          break;

        case 1:
          vf->lineInfo[(int) vgpChar(vf,0)]->fieldCodec = vcDeserialize (vgpString(vf));
          break;

        case 2:
          vf->lineInfo[(int) vgpChar(vf,0)]->listCodec = vcDeserialize (vgpString(vf));
          break;

        default:
          parseError (vf, "unknown header line type %c", vf->lineType);
          break;
      }
    }
  vf->isCheckString = FALSE;   // user can set this back to TRUE if they wish

  // allocate codec buffer - always allocate enough to handle fields of all line types

  { I64       size;
    LineInfo *li, *l0;
    VgpFile  *v;

    size = MAX_FIELD*sizeof(Field);
    for (i = 0; i < 128; ++i)
      { li = vf->lineInfo[i];
        if (li != NULL && li->listCodec)
          if (size < li->given.max * li->listByteSize)
            size = li->given.max * li->listByteSize;
       }
    vf->codecBuf     = new (size+1, void);  // add one for worst case codec usage
    vf->codecBufSize = size+1;

    startOff = ftello (f);
    for (i = 1; i < nthreads; i++)
      { v = vf+i;

        lineInfoDestroy (v->lineInfo['&']);
        lineInfoDestroy (v->lineInfo['*']);
        v->lineInfo['&'] = NULL;
        v->lineInfo['*'] = NULL;

        for (j = 0; j < 128; j++)
          { li = v->lineInfo[j];
            if (li != NULL)
              { l0 = vf->lineInfo[j];
                li->fieldCodec = l0->fieldCodec;
                li->listCodec  = l0->listCodec;
                if (li->listField > 0)
                  { li->bufSize = l0->bufSize;
                    li->buffer  = new (l0->bufSize*l0->listByteSize, void);
                  }
                li->given = l0->given;
              }
          }

        v->codecBuf = new (size+1, void);
        v->codecBufSize = size+1;
        if (fseeko (v->f, startOff, SEEK_SET) != 0)
          die ("can't seek to start of data");

        v->lineInfo['&'] = vf->lineInfo['&'];
        v->lineInfo['*'] = vf->lineInfo['*'];

        v->isIndexIn = vf->isIndexIn;
        v->subType   = vf->subType;
      }
  }
  
  vf->f = f;
  return (vf);
}


/***********************************************************************************
 *
 *   VGP_USER_BUFFER / CLOSE / GOTO
 *
 **********************************************************************************/

  // This lets the user reassign the buffer that lists in a particular line type are read into.
  //   If this is not set, a default buffer is provided.  If buffer == NULL then the package
  //   reverts to the default buffer.  This routine can be called repeatedly.
  // NB the package doesn't check the size of a user supplied buffer - the user must allocate
  //   enough memory for all forthcoming list data.

void vgpUserBuffer (VgpFile *vf, char lineType, void *buffer)
{ LineInfo *li;

  li = vf->lineInfo[(int) lineType];
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
          li->buffer  = new (li->given.max*li->listByteSize, void);
        }
      li->isUserBuf = FALSE;
    }
}

BOOL vgpGotoObject (VgpFile *vf, I64 i)
{ if (vf != NULL && vf->isIndexIn)
    if (0 <= i && i < vf->lineInfo[(int) vf->objectType]->given.count)
      if (fseek (vf->f, ((I64 *) vf->lineInfo['&']->buffer)[i], SEEK_SET) == 0)
        { vf->object = i;
          return (TRUE);
        }
  return (FALSE);
}

I64 vgpGotoGroup (VgpFile *vf, I64 i)
{ if (vf != NULL && vf->isIndexIn)
    if (0 <= i && i < vf->lineInfo[(int) vf->groupType]->given.count)
      { I64 *groupIndex = (I64 *) vf->lineInfo['*']->buffer;
        if (!vgpGotoObject(vf,groupIndex[i]))
	  return 0 ;
        return (groupIndex[i+1] - groupIndex[i]);
      }
  return (0);
}


/***********************************************************************************
 *
 *   VGP_OPEN_WRITE_(NEW | FROM)
 *
 **********************************************************************************/

VgpFile *vgpFileOpenWriteNew (const char *path, FileType fileType, SubType subType,
                              BOOL isBinary, int nthreads)
{ VgpFile *vf, *v;
  FILE    *f;
  int      i, pid;
  char     name[100];

  if (subType && subPrimary[subType] != fileType)
    die ("subtype %s is not secondary for filetype %s",
           subTypeName[subType], fileTypeName[fileType]);
  
  if (strcmp (path, "-") == 0)
    f = stdout;
  else
    { f = fopen (path, "w");
      if (f == NULL)
        return (NULL);
    }

  if (nthreads == 1)
    { vf = vgpFileCreate (fileType);
      vf->f = f;

      vf->isWrite  = TRUE;
      vf->subType  = subType;
      vf->isBinary = isBinary;
      vf->isLastLineBinary = TRUE; // we don't want to add a newline before the first true line

      vf->codecBufSize = MAX_FIELD*sizeof(Field) + 1;
      vf->codecBuf     = new (vf->codecBufSize, void); 

      return (vf);
    }

  pid = getpid();

  vf = new (nthreads, VgpFile);
  for (i = 0; i < nthreads; i++)
    { v = vgpFileCreate (fileType);

      v->isWrite  = TRUE;
      v->subType  = subType;
      v->isBinary = isBinary;
      
      v->codecBufSize = MAX_FIELD*sizeof(Field) + 1;
      v->codecBuf     = new (v->codecBufSize, void);

      if (i > 0)
        { v->codecTrainingSize /= 3*nthreads;
          v->share = -i;
          v->isLastLineBinary = isBinary;
          sprintf(name,".part.%d.%d",pid,i);
          f = fopen (name, "w");
          if (f == NULL)
            die ("Cannot create unique temporary file");
        }
      else
        { v->share = nthreads;
          v->isLastLineBinary = TRUE;
          v->fieldLock = mutexInit;
          v->listLock  = mutexInit;
        }
      v->f = f;
      vf[i] = *v;
      free (v);
    }

  return (vf);
}

VgpFile *vgpFileOpenWriteFrom (const char *path, VgpFile *vfIn, BOOL useAccum, BOOL isBinary, int nthreads)
{ VgpFile  *vf;
  LineInfo *li;
  int       i;
  I64       size, sz;

  vf = vgpFileOpenWriteNew (path, vfIn->fileType, vfIn->subType, isBinary, nthreads);

  vgpInheritProvenance (vf, vfIn);
  vgpInheritReference  (vf, vfIn);
  vgpInheritDeferred   (vf, vfIn);

  if (useAccum)
    { for (i = 0; i < 128 ; ++i)
        if (vfIn->lineInfo[i])
          vf->lineInfo[i]->given = vfIn->lineInfo[i]->accum;
    }
  else
    { for (i = 0; i < 128 ; ++i)
        if (vfIn->lineInfo[i])
          vf->lineInfo[i]->given = vfIn->lineInfo[i]->given;
    }

  // allocate codec buffer - always allocate enough to handle fields of all line types

  size = vf->codecBufSize;
  for (i = 0; i < 128 ; ++i)
    { li = vf->lineInfo[i];
      if (li != NULL && li->listCodec != NULL)
        { sz = li->given.max * li->listByteSize;
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

static BOOL addProvenance(VgpFile *vf, Provenance *from, int n)
{ I64         o;
  LineInfo   *l;
  Provenance *p;

  if (n == 0)
    return (FALSE);
  if (vf->isHeaderOut)
    die("can't addProvenance after writing header");

  l = vf->lineInfo['!'];
  o = l->accum.count;
  l->accum.count += n;

  p = new(o+n, Provenance);
  if (o > 0)
    memcpy (p, vf->provenance, o*sizeof(Provenance));
  memcpy (p+o, from, n*sizeof(Provenance));

  free (vf->provenance);
  vf->provenance = p;

  return (TRUE);
}

BOOL vgpInheritProvenance(VgpFile *vf, VgpFile *source)
{ return (addProvenance(vf, source->provenance, source->lineInfo['!']->accum.count)); }

BOOL vgpAddProvenance(VgpFile *vf, char *prog, char *version, char *command, char *date)
{ Provenance p;

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
  return addProvenance (vf, &p, 1);
}

static BOOL addReference(VgpFile *vf, Reference *from, int n, BOOL isDeferred)
{ I64        o;
  LineInfo  *l;
  Reference *r, **t;

  if (n == 0)
    return (FALSE);
  if (vf->isHeaderOut)
    die ("can't addReference after writing header");

  if (isDeferred)
    { l = vf->lineInfo['>'];
      t = &(vf->deferred);
    }
  else
    { l = vf->lineInfo['<'];
      t = &(vf->reference);
    }
  o = l->accum.count;
  l->accum.count += n;

  r = new (o+n, Reference);
  if (o > 0)
    memcpy (r, *t, o*sizeof(Reference));
  memcpy (r+o, from, n*sizeof(Reference));

  free (*t);
  *t = r;

  return (TRUE);
}

BOOL vgpInheritReference(VgpFile *vf, VgpFile *source)
{ return (addReference(vf, source->reference, source->lineInfo['<']->accum.count, FALSE)); }

BOOL vgpAddReference(VgpFile *vf, char *filename, I64 count)
{ Reference ref;
  ref.filename = filename;
  ref.count    = count;
  return (addReference(vf, &ref, 1, FALSE));
}

BOOL vgpInheritDeferred (VgpFile *vf, VgpFile *source)
{ return (addReference (vf, source->deferred, source->lineInfo['>']->accum.count, TRUE)); }

BOOL vgpAddDeferred (VgpFile *vf, char *filename)
{ Reference ref;
  ref.filename = filename;
  return (addReference (vf, &ref, 1, TRUE));
}


/***********************************************************************************
 *
 *   VGP_WRITE_HEADER / FOOTER
 *
 **********************************************************************************/

void vgpWriteHeader (VgpFile *vf)
{ int         i,n;
  Reference  *r;
  Provenance *p;
  LineInfo   *li;

  if ( ! vf->isWrite)
    die ("Trying to write header to a file open for reading");
  if (vf->line > 0)
    die ("Cannot write header after writing one or more data lines");
  if (vf->lineInfo[(int) vf->objectType]->given.count == 0 && ! vf->isBinary)
    die ("Information for ASCII header is not present, use vgpFileOpenWriteFrom");

  fprintf (vf->f, "1 %lu %s %lld %lld",
                 strlen(fileTypeName[vf->fileType]), fileTypeName[vf->fileType],
                 vf->major, vf->minor);
  vf->line += 1;

  if (vf->subType)
    { fprintf (vf->f, "\n2 %lu %s", strlen(subTypeName[vf->subType]), subTypeName[vf->subType]);
      vf->line += 1;
    }

  r = vf->reference;
  n = vf->lineInfo['<']->accum.count;
  for (i = 0; i < n; i++, r++)
    { fprintf (vf->f, "\n< %lu %s %lld", strlen(r->filename), r->filename, r->count);
      vf->line += 1;
    }

  r = vf->deferred;
  n = vf->lineInfo['>']->accum.count;
  for (i = 0; i < n; i++, r++)
    { fprintf (vf->f, "\n> %lu %s", strlen(r->filename), r->filename);
      vf->line += 1;
    }
  
  p = vf->provenance; 
  n = vf->lineInfo['!']->accum.count;
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
	  li = vf->lineInfo[i];
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
 *   VGP_WRITE_LINE
 *
 **********************************************************************************/

static int writeStringList (VgpFile *vf, char t, int len)
{ LineInfo *li;
  int       j, nByteWritten = 0;
  I64       sLen, totLen;
  char     *buf;

  buf = (char *) vf->lineInfo[(int) t]->buffer;

  totLen = 0;
  for (j = 0; j < len; j++)
    { sLen = strlen (buf);
      totLen += sLen;
      nByteWritten += fprintf (vf->f, " %lld %s", sLen, buf);
      buf += sLen + 1;
    }

  li = vf->lineInfo[(int) t];
  li->accum.total += totLen;
  if (li->accum.max < totLen)
    li->accum.max = totLen;

  return nByteWritten ;
}

// process is to fill fields by assigning to macros, then call - list contents are in buf
// NB adds '\n' before writing line not after, so user fprintf() can add extra material
// first call will write initial header, allowing space for count sizes to expand on close

void vgpWriteLine (VgpFile *vf, char t, I64 listLen, void *listBuf)
{ I64       i, j, ix;
  LineInfo *li;

  //  fprintf (stderr, "write line type %c\n", t) ;
  
  if ( ! vf->isWrite)
    die ("Trying to write a line to a file open for reading");
  if (vf->isFinal && isalpha(t))
    die ("Cannot write more data after counts are finalized %c" ,t);

  li = vf->lineInfo[(int) t];
  if (li == NULL)
    die ("vgpWriteLine error: line type %c not present in file spec %s ",
         t,fileTypeName[vf->fileType]);

  if (listBuf == NULL)
    listBuf = li->buffer;

  if ( ! vf->isLastLineBinary)      // terminate previous ascii line
    fputc ('\n', vf->f);

  vf->line  += 1;
  li->accum.count += 1;
  if (t == vf->groupType)
    updateGroupCount(vf, TRUE);

  ix = li->listField-1;
  if (ix >= 0)
    { if (listLen >= 0)
	vf->field[ix].len = listLen ;
      else
	die ("listLen %lld must be positive", listLen) ;
    }

  // BINARY - block write and optionally compress
  
  if (vf->isBinary)
    { U8  x, cBits;
      int nField;
      I64 fieldSize, nBits, listBytes, listSize;

      if (!vf->isLastLineBinary)
	vf->byte = ftello (vf->f) ;

      if (t == vf->objectType) // update index and increment object count
        { LineInfo *lx = vf->lineInfo['&'];

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
	      die ("byte offset index error") ;
	  }
#endif
          vf->object += 1;
        }
      if (t == vf->groupType)
        { LineInfo *lx = vf->lineInfo['*'];
	  
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
      fieldSize = nField*sizeof(Field);

      if (ix >= 0 && li->fieldSpec[ix] == INT_LIST)
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
                die("compressed fields");
	      vf->byte += 2 + ((nBits+7) >> 3) ;
            }
          else
            { fputc (x, vf->f);
              if (nField > 0)
                { if (fwrite (vf->field, fieldSize, 1, vf->f) != 1)
                    die ("write fields: t %c, nField %d, fieldSize %lld", t, nField, fieldSize);
                }
	      vf->byte += 1 + fieldSize ;
            }
        }
      else
        { fputc (x, vf->f);
          if (nField > 0)
            if (fwrite (vf->field, fieldSize, 1, vf->f) != 1)
              die ("write fields: t %c, nField %d, fieldSize %lld", t, nField, fieldSize);
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
                    { VgpFile  *ms;
                      LineInfo *lx;

                      if (vf->share < 0)
                        { ms = vf + vf->share;
                          lx = ms->lineInfo[(int) t]; 
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
                                                ms[i].lineInfo[(int) t]->fieldCodec);
			      vcCreateCodec (lx->fieldCodec, 1);
                              for (i = 1; i < ms->share; i++)
                                { vcDestroy (ms[i].lineInfo[(int) t]->fieldCodec);
                                  ms[i].lineInfo[(int) t]->fieldCodec = lx->fieldCodec;
                                }
                              lx->isUseFieldCodec = TRUE;
                              for (i = 1; i < ms->share; i++)
                                ms[i].lineInfo[(int) t]->isUseFieldCodec = TRUE;
                            }
                        }

                      pthread_mutex_unlock(&ms->fieldLock);
                    }
                }
            }
        }

      // Write the list if there is one

      if (ix >= 0)
        { li->accum.total += listLen;
          if (listLen > li->accum.max)
            li->accum.max = listLen;

          if (listLen > 0)
            { listBytes = li->listByteSize - (vgpInt(vf, ix) >> 56); // data bytes
              listSize  = listLen * listBytes;

              if (li->fieldSpec[ix] == STRING_LIST) // handle as ASCII
                vf->byte += writeStringList (vf, t, listLen);

              else if (x & 0x2)
                { if (listSize >= vf->codecBufSize)
                    { free (vf->codecBuf);
                      vf->codecBufSize = listSize+1;
                      vf->codecBuf     = new (vf->codecBufSize, void);
                    }
                  nBits = vcEncode (li->listCodec, listSize, listBuf, vf->codecBuf);
                  if (fwrite (&nBits, sizeof(I64), 1, vf->f) != 1)
                    die ("failed to write list nBits");
                  if (fwrite (vf->codecBuf, ((nBits+7) >> 3), 1, vf->f) != 1)
                    die ("failed to write compressed list");
		  vf->byte += sizeof(I64) + ((nBits+7) >> 3);
                }

              else
                { if (fwrite (listBuf, listSize, 1, vf->f) != 1)
                    die ("failed to write list ix %d listLen %lld listSize %lld listBuf %lx",
			 ix, listLen, listSize, listBuf);
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
                            { VgpFile  *ms;
                              LineInfo *lx;

                              if (vf->share < 0)
                                { ms = vf + vf->share;
                                  lx = ms->lineInfo[(int) t]; 
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
                                                        ms[i].lineInfo[(int) t]->listCodec);
                                      vcCreateCodec (lx->listCodec, 1);
                                      lx->isUseListCodec = TRUE;
                                      for (i = 1; i < ms->share; i++)
                                        { vcDestroy (ms[i].lineInfo[(int) t]->listCodec);
                                          ms[i].lineInfo[(int) t]->listCodec = lx->listCodec;
                                          ms[i].lineInfo[(int) t]->isUseListCodec = TRUE;
                                        }
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
        switch (li->fieldSpec[i])
        { case INT:
            fprintf (vf->f, " %lld", vf->field[i].i);
            break;
          case REAL:
            fprintf (vf->f, " %f", vf->field[i].r);
            break;
          case CHAR:
            fprintf (vf->f, " %c", vf->field[i].c);
            break;
          case STRING:
          case INT_LIST:
          case REAL_LIST:
          case STRING_LIST:
            li->accum.total += listLen;
            if (listLen > li->accum.max)
              li->accum.max = listLen;

            fprintf (vf->f, " %lld", listLen);
            if (li->fieldSpec[i] == STRING)
              { if (listLen > INT_MAX)
                  die ("write problem: string length %lld > current max %d", listLen, INT_MAX);
                fprintf (vf->f, " %.*s", (int) listLen, (char *) listBuf);
              }
            else if (li->fieldSpec[i] == INT_LIST)
              { I64 *b = (I64 *) listBuf;
                for (j = 0; j < listLen ; ++j)
                  fprintf (vf->f, " %lld", b[j]);
              }
            else if (li->fieldSpec[i] == REAL_LIST)
              { double *b = (double *) listBuf;
                for (j = 0; j < listLen ; ++j)
                  fprintf (vf->f, " %f", b[j]);
              }
            else // STRING_LIST
              writeStringList (vf, t, listLen);
            break;
        }
      vf->isLastLineBinary = FALSE;
    }
}

void vgpWriteComment (VgpFile *vf, char *comment)
{
  if (vf->isLastLineBinary)
    vgpWriteLine (vf, '/', strlen(comment), comment) ;
  else
    fprintf (vf->f, " %s", comment) ;
}

/***********************************************************************************
 *
 *    MERGING, FOOTER HANDLING, AND CLOSE
 *
 **********************************************************************************/

static void vgpWriteFooter (VgpFile *vf)
{ int       i,n;
  off_t     footOff;
  LineInfo *li;
  
  footOff = ftello (vf->f);
  if (footOff < 0)
    die ("failed footer ftell");

  //  first the per-linetype information

  for (i = 'A'; i <= 'Z'+1 ; i++)
    { if (i == 'Z'+1)
	{ if (vf->groupType) // NB group types are all lower case so > 'Z'+1
	    i = vf->groupType ;
	  else
	    break ;
	}
      li = vf->lineInfo[i];
      if (li != NULL && li->accum.count > 0)
        { fprintf (vf->f, "# %c %lld\n", i, li->accum.count);
	  if (li->listField)
            { fprintf (vf->f, "@ %c %lld\n", i, li->accum.max);
	      fprintf (vf->f, "+ %c %lld\n", i, li->accum.total);
	    }
	  if (vf->groupType && i != vf->groupType && vf->group > 0)
	    { fprintf (vf->f, "%% %c # %c %lld\n", vf->groupType, i, li->accum.groupCount);
	      if (li->listField)
		fprintf (vf->f, "%% %c + %c %lld\n", vf->groupType, i, li->accum.groupTotal);
	    }
          if (li->isUseFieldCodec)
            { vgpChar(vf,0) = i;
              n = vcSerialize (li->fieldCodec, vf->lineInfo[1]->buffer);
              vgpWriteLine (vf, 1, n, vf->lineInfo[i]->buffer);
            }
          if (li->isUseListCodec && li->listCodec != DNAcodec)
            { vgpChar(vf,0) = i;
              n = vcSerialize (li->listCodec, vf->lineInfo[2]->buffer);
              vgpWriteLine (vf, 2, n, vf->lineInfo[2]->buffer);
            }
        }
    }

  vgpWriteLine (vf, '&', vf->object, NULL); // number of objects in file = length of index

  if (vf->groupType > 0 && vf->group > 0)
    { ((I64 *) vf->lineInfo['*']->buffer)[vf->group] = vf->object;
      vgpWriteLine (vf, '*', vf->group+1, NULL); // number of groups in file + 1 = length of index
    }

  fprintf (vf->f, "^\n");

  if (fwrite (&footOff, sizeof(off_t), 1, vf->f) != 1)
    die ("failed writing footer offset");
}

void vgpFinalizeCounts(VgpFile *vf)
{ int       i, j, n, k, len;
  LineInfo *li, *ln;

  if (vf->share < 0)
    die ("Cannot call vgpFileClose on a slave VgpFile");

  vf->isFinal = TRUE;

  if (vf->share == 0)
    { updateGroupCount(vf,FALSE);
      return;
    }

  len = vf->share;
  
  //  Close current groups at the end of each part (if any)

  if (vf->groupType > 0)
    for (i = 'A'; i <= 'Z'; i++)
      if (vf->lineInfo[i] != NULL)
        for (j = 0; j < len; j++)
          if (vf[j].inGroup)
            { I64 oc, ot;

              ot = oc = 0;
              for (k = j+1; k < len; k++)
                if (vf[k].inGroup)
                  { oc += vf[k].lineInfo[i]->oCount;
                    ot += vf[k].lineInfo[i]->oTotal; 
                    break;
                  }
                else
                  { oc += vf[k].lineInfo[i]->accum.count;
                    ot += vf[k].lineInfo[i]->accum.total;
                  }

              li = vf[j].lineInfo[i];
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
    { ln = vf->lineInfo[i];
      for (j = 1; j < len; j++)
        { li = (vf+j)->lineInfo[i];
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
        { li = vf[j].lineInfo['*'];
          gi = (I64 *) (li->buffer);
          for (i = 0; i < vf[j].group; i++)
            gb[ns++] = gi[i] + off;
          off += vf[j].object;
        }
      gb[ns] = off;
      li = vf->lineInfo['*'];
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
      { li = vf[j].lineInfo['&'];
        gi = (I64 *) (li->buffer);
        for (i = 0; i < vf[j].object; i++)
          gb[ns++] = gi[i] + off;
        off += ftello(vf[j].f);
      }

    li = vf->lineInfo['&'];
    free(li->buffer);
    li->buffer  = gb;
    li->bufSize = ns;
    vf->object  = ns;
  }
}

// automatically rewrites header if allowed when writing

void vgpFileClose (VgpFile *vf)
{
  if (vf->share < 0)
    die ("Cannot call vgpFileClose on a slave VgpFile");

  if ( ! vf->isFinal)
    vgpFinalizeCounts (vf);

  if (vf->isWrite)
    { if (vf->share > 0)
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
                  die ("write error cat'ing thread bits (vgpFileClose)");
              if (unlink(name) < 0)
                die ("could not delete thread file %s", name);
            }
          free(buf);
        }
      fputc ('\n', vf->f);  // end of file if ascii, end of data marker if binary
      if (vf->isBinary) // write the footer
        vgpWriteFooter (vf);
    }
  vgpFileDestroy (vf);
}
  
/***********************************************************************************
 *
 *    UTILITIES: memory allocation, file opening, timer
 *
 **********************************************************************************/

void die(char *format, ...)
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

void *myalloc(size_t size)
{ void *p;

  p = malloc(size);
  if (p == NULL)
    die("myalloc failure requesting %d bytes", size);
  nAlloc     += 1;
  totalAlloc += size;
  return (p);
}

void *mycalloc(size_t number, size_t size)
{ void *p;

  p = calloc(number,size);
  if (p == NULL)
    die("mycalloc failure requesting %d objects of size %d", number, size);
  nAlloc     += 1;
  totalAlloc += size*number;
  return p;
}

FILE *fzopen(const char *path, const char *mode)
{
#ifdef WITH_ZLIB   // very cool from https://stackoverflow.com/users/3306211/fernando-mut
  gzFile zfp;      // fernando said *zfp - makes me worry....

  zfp = gzopen(path,mode);     // try gzopen
  if (zfp == NULL)
    return fopen(path,mode);

  return funopen(zfp,                                 // open file pointer
                 (int(*)(void*,char*,int))gzread,
                 (int(*)(void*,const char*,int))gzwrite,
                 (fpos_t(*)(void*,fpos_t,int))gzseek,
                 (int(*)(void*))gzclose);
#else
  return fopen(path,mode);
#endif
}

FILE *fopenTag (char* root, char* tag, char* mode)
{ char *fileName;
  FILE *f;

  if (strlen(tag) > 30)
    die("tag %s in fopenTag too long - should be < 30 chars", tag);
  fileName = new (strlen(root) + 32, char);
  sprintf (fileName, "%s.%s", root, tag);
  f = fzopen (fileName, mode);
  free (fileName);
  return (f);
}

  // NB this is not threadsafe - do not call inside threads

static struct rusage rLast, rFirst;

void timeUpdate(FILE *f)
{ static BOOL   isFirst = 1;
  struct rusage rNew;
  int           secs, mics;

  getrusage(RUSAGE_SELF, &rNew);

  if (isFirst)
    { rFirst  = rNew;
      isFirst = FALSE;
    }
  else
    { secs = rNew.ru_utime.tv_sec  - rLast.ru_utime.tv_sec;
      mics = rNew.ru_utime.tv_usec - rLast.ru_utime.tv_usec;
      if (mics < 0)
        { mics += 1000000; secs -= 1 ; }
      fprintf (f, "user\t%d.%06d", secs, mics);

      secs = rNew.ru_stime.tv_sec  - rLast.ru_stime.tv_sec;
      mics = rNew.ru_stime.tv_usec - rLast.ru_stime.tv_usec;
      if (mics < 0)
        { mics += 1000000; secs -= 1 ; }
      fprintf (f, "\tsystem\t%d.%06d", secs, mics);

      fprintf (f, "\tmax_RSS\t%ld", rNew.ru_maxrss - rLast.ru_maxrss);
      fprintf (f, "\tnnew\t%lld", nAlloc);   
      fprintf (f, "\ttotnew\t%lld", totalAlloc);   
      fprintf (f, "\n");
    }

  rLast = rNew;
}

void timeTotal(FILE *f)
{ rLast = rFirst;
  timeUpdate(f);
}

/********************* end of file ***********************/
