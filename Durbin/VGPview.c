/*  File: VGPview.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Feb  4 01:05 2020 (rd109)
 * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "VGPlib.h"

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

typedef struct IndexListStruct {
  I64 i0, iN ;
  struct IndexListStruct *next ;
} IndexList ;

static IndexList *parseIndexList (char *s)
{
  IndexList *ol, *ol0 = ol = new0 (1, IndexList) ;
  while (*s)
    { while (*s >= '0' && *s <= '9') ol->i0 = ol->i0*10 + (*s++ - '0') ;
      if (*s == '-')
	{ ++s ; while (*s >= '0' && *s <= '9') ol->iN = ol->iN*10 + (*s++ - '0') ;
	  if (ol->iN <= ol->i0) die ("end index %lld <= start index %lld", ol->iN, ol->i0) ;
	}
      else
	ol->iN = ol->i0 + 1 ;
      if (*s == ',') { ol->next = new0 (1, IndexList) ; ol = ol->next ; ++s ; }
      else if (*s) die ("unrecognised character %c at %s in object list\n", *s, s) ;
    }
  return ol0 ; 
}

int main (int argc, char **argv)
{
  I64 i ;
  FileType fileType = 0 ;
  char *outFileName = "-" ;
  BOOL isNoHeader = FALSE, isHeaderOnly = FALSE, isBinary = FALSE, isUsage = FALSE ;
  IndexList *objList = 0, *groupList = 0 ;
  
  timeUpdate (0) ;

  char *command = commandLine (argc, argv) ;
  --argc ; ++argv ;		/* drop the program name */

  if (!argc)
    { fprintf (stderr, "VGPview [options] vgpfile\n") ;
      fprintf (stderr, "  -t --type <abc>           file type, e.g. seq, aln - required if no header\n") ;
      fprintf (stderr, "  -h --noHeader             skip the header in ascii output\n") ;
      fprintf (stderr, "  -H --headerOnly           only write the header (in ascii)\n") ;
      fprintf (stderr, "  -b --binary               write in binary (default is ascii)\n") ;
      fprintf (stderr, "  -o --output <filename>    output file name (default stdout)\n") ;
      fprintf (stderr, "  -i --index x[-y](,x[-y])* write specified objects\n") ;
      fprintf (stderr, "  -g --group x[-y](,x[-y])* write specified groups\n") ;
      fprintf (stderr, "index and group only work for binary files; '-i 0-10' outputs first 10 objects\n") ;
      exit (0) ;
    }
  
  while (argc > 1)
    if (!strcmp (*argv, "-t") || !strcmp (*argv, "--type"))
      { for (i = SEQ ; i < MAX_FILE ; ++i)
	  if (!strcmp (argv[1], fileTypeName[i])) fileType = i ;
	if (!fileType) die ("unknown file type %s requested", argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-h") || !strcmp (*argv, "--header"))
      { isNoHeader = TRUE ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-H") || !strcmp (*argv, "--headerOnly"))
      { isHeaderOnly = TRUE ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-b") || !strcmp (*argv, "--binary"))
      { isBinary = TRUE ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-o") || !strcmp (*argv, "--output"))
      { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-i") || !strcmp (*argv, "--index"))
      { objList = parseIndexList (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-g") || !strcmp (*argv, "--group"))
      { objList = parseIndexList (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else die ("unknown option %s - run without arguments to see options", *argv) ;

  if (isBinary) isNoHeader = FALSE ;
  if (isHeaderOnly) { isHeaderOnly = TRUE ; isBinary = FALSE ; }
  
  if (argc != 1) die ("can currently only take one input file") ;

  VgpFile *vfIn = vgpFileOpenRead (*argv, fileType, 1) ; /* reads the header */
  if (!vfIn) die ("failed to open vgp file %s", *argv) ;
  
    { VgpFile *vfOut = vgpFileOpenWriteFrom (outFileName, vfIn, FALSE, isBinary, 1) ;
      if (!vfOut) die ("failed to open output file %s", outFileName) ;

      if (isHeaderOnly)
	vgpWriteHeader (vfOut) ;
      else
	{ vgpAddProvenance (vfOut, "vgpview", "0.0", command, 0) ;
	  if (!isNoHeader) vgpWriteHeader (vfOut) ;

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
		      vgpWriteLine (vfOut, vfIn->lineType, vgpString(vfIn)) ;
		      if (!vgpReadLine (vfIn)) break ;
		      if (vfIn->lineType == vfIn->objectType) ++i ;
		    }
		  objList = objList->next ;
		}
	    }
	  else
	    while (vgpReadLine (vfIn))
	      { memcpy (vfOut->field, vfIn->field, fieldSize[vfIn->lineType]) ;
		vgpWriteLine (vfOut, vfIn->lineType, vgpString(vfIn)) ;	      }
	}
      
      vgpFileClose (vfOut) ;
    }

  timeTotal (stderr) ;
}

/********************* end of file ***********************/
