/*  File: vgpvalidate.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jul 22 13:25 2019 (rd109)
 * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "vgprd.h"

#include <assert.h>

#include <string.h>		/* strcmp etc. */
#include <stdlib.h>		/* for exit() */

static char *commandLine (int argc, char **argv)
{
  if (!argc) die ("commandLine needs at least one arg") ;
  int i, totLen = 0 ;
  for (i = 0 ; i < argc ; ++i) totLen += 1 + strlen(argv[i]) ;
  char *buf = new (totLen, char) ;
  strcpy (buf, argv[0]) ;
  for (i = 1 ; i < argc ; ++i) { strcat (buf, " ") ; strcat (buf, argv[i]) ; }
  return buf ;
}

typedef struct ObjectListStruct {
  I64 i0, iN ;
  struct ObjectListStruct *next ;
} ObjectList ;

static ObjectList *parseObjectList (char *s)
{
  ObjectList *ol, *ol0 = ol = new0 (1, ObjectList) ;
  while (*s)
    { while (*s >= '0' && *s <= '9') ol->i0 = ol->i0*10 + (*s++ - '0') ;
      if (*s == '-')
	{ ++s ; while (*s >= '0' && *s <= '9') ol->iN = ol->iN*10 + (*s++ - '0') ;
	  if (ol->iN <= ol->i0) die ("end index %lld <= start index %lld", ol->iN, ol->i0) ;
	}
      else
	ol->iN = ol->i0 + 1 ;
      if (*s == ',') { ol->next = new0 (1, ObjectList) ; ol = ol->next ; ++s ; }
      else if (*s) die ("unrecognised character %c at %s in object list\n", *s, s) ;
    }
  return ol0 ; 
}

int main (int argc, char **argv)
{
  I64 i ;
  FileType fileType = 0 ;
  char *outFileName = "-" ;
  BOOL isHeader = FALSE, isHeaderOnly = FALSE, isBinary = FALSE, isUsage = FALSE ;
  ObjectList *objList = 0 ;
  
  timeUpdate (0) ;

  char *command = commandLine (argc, argv) ;
  --argc ; ++argv ;		/* drop the program name */

  if (!argc)
    { fprintf (stderr, "vgpview [options] vgpfile\n") ;
      fprintf (stderr, "  -t --type <abc>           file type, e.g. seq, aln - required if no header\n") ;
      fprintf (stderr, "  -h --header               include the header in ascii output\n") ;
      fprintf (stderr, "  -H --headerOnly           only write the header (in ascii)\n") ;
      fprintf (stderr, "  -b --binary               write in binary (default is ascii)\n") ;
      fprintf (stderr, "  -o --output <filename>    output file name (default stdout)\n") ;
      fprintf (stderr, "  -i --index x[-y](,x[-y])* write specified objects\n") ;
      fprintf (stderr, "  -u --usage                byte usage per line type; no other output\n") ;
      fprintf (stderr, "index works only for binary files; '-i 0-10' outputs first 10 objects\n") ;
      exit (0) ;
    }
  
  while (argc && **argv == '-')
    if (argc > 1 && (!strcmp (*argv, "-t") || !strcmp (*argv, "--type")))
      { for (i = SEQ ; i < MAX_FILE ; ++i)
	  if (!strcmp (argv[1], fileTypeName[i])) fileType = i ;
	if (!fileType) die ("unknown file type %s requested", argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-h") || !strcmp (*argv, "--header"))
      { isHeader = TRUE ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-H") || !strcmp (*argv, "--headerOnly"))
      { isHeaderOnly = TRUE ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-b") || !strcmp (*argv, "--binary"))
      { isBinary = TRUE ; --argc ; ++argv ; }
    else if (argc > 1 && (!strcmp (*argv, "-o") || !strcmp (*argv, "--output")))
      { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (argc > 1 && (!strcmp (*argv, "-i") || !strcmp (*argv, "--index")))
      { objList = parseObjectList (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-u") || !strcmp (*argv, "--usage"))
      { isUsage = TRUE ; --argc ; ++argv ; }
    else die ("unknown option %s - run without arguments to see options", *argv) ;

  if (isBinary) isHeader = TRUE ;
  if (isHeaderOnly) { isHeaderOnly = TRUE ; isBinary = FALSE ; }
  
  if (argc != 1) die ("can currently only take one input file") ;

  VgpFile *vfIn = vgpFileOpenRead (*argv, fileType) ; /* reads the header */
  if (!vfIn) die ("failed to open vgp file %s", *argv) ;
  
  if (isUsage)
    { I64 usage[128] ; memset (usage, 0, 128*sizeof(I64)) ; 
      off_t u, uLast = ftello (vfIn->f) ;
#define UPDATE_U { u = ftello (vfIn->f) ; usage[vfIn->lineType] += u-uLast ; uLast = u ; }
      
      if (objList)
	{ while (objList)
	    { if (!vgpGotoObject (vfIn, objList->i0)) die ("bad seek to %lld", objList->i0 ) ;
	      if (!vgpReadLine (vfIn)) die ("can't read object %lld", objList->i0) ;
	      UPDATE_U ;
	      for (i = objList->i0 ; i < objList->iN ; ) // conditional increment end of loop
		{ if (!vgpReadLine (vfIn)) break ;
		  UPDATE_U ;
		  if (vfIn->lineType == vfIn->objectType) ++i ;
		}
	      objList = objList->next ;
	    }
	}
      else
	while (vgpReadLine (vfIn))
	  UPDATE_U ;

      UPDATE_U ;
      for (i = 'A' ; i < 128 ; ++i)
	if (usage[i]) printf ("usage line type %c bytes %lld\n", (char)i, usage[i]) ;
     }
  else
    { VgpFile *vfOut = vgpFileOpenWriteFrom (outFileName, vfIn, isBinary) ;
      if (!vfOut) die ("failed to open output file %s", outFileName) ;

      if (isHeaderOnly)
	{ vgpWriteHeader (vfOut) ; fputc ('\n', vfOut->f) ; }
      else
	{ vgpAddProvenance (vfOut, "vgpview", "0.0", command, 0) ;
	  if (isHeader) vgpWriteHeader (vfOut) ;

	  static size_t fieldSize[128] ;
	  for (i = 0 ; i < 128 ; ++i)
	    if (vfIn->lineInfo[i]) fieldSize[i] = vfIn->lineInfo[i]->nField*sizeof(Field) ;

	  if (objList)
	    { vfOut->isLastLineBinary = TRUE ; /* prevents leading '\n' */
	      while (objList)
		{ if (!vgpGotoObject (vfIn, objList->i0)) die ("bad seek to %lld", objList->i0 ) ;
		  if (!vgpReadLine (vfIn)) die ("can't read object %lld", objList->i0) ;
		  for (i = objList->i0 ; i < objList->iN ; ) // conditional increment end of loop
		    { memcpy (vfOut->field, vfIn->field, fieldSize[vfIn->lineType]) ;
		      vgpWriteLine (vfOut, vfIn->lineType, vfIn->lineInfo[vfIn->lineType]->buffer) ;
		      if (!vgpReadLine (vfIn)) break ;
		      if (vfIn->lineType == vfIn->objectType) ++i ;
		    }
		  objList = objList->next ;
		}
	    }
	  else
	    while (vgpReadLine (vfIn))
	      { memcpy (vfOut->field, vfIn->field, fieldSize[vfIn->lineType]) ;
		vgpWriteLine (vfOut, vfIn->lineType, vfIn->lineInfo[vfIn->lineType]->buffer) ;
	      }
	}
      
      vgpClose (vfOut) ;
    }

  timeTotal (stderr) ;
}

/********************* end of file ***********************/
