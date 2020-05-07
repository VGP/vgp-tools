/*  File: VGPview.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Apr 16 20:48 2020 (rd109)
 * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "VGPlib.h"

#include <assert.h>
#include <stdbool.h>  /* bool, true, false */
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

static void transferLine (VgpFile *vfIn, VgpFile *vfOut, size_t *fieldSize)
{ memcpy (vfOut->field, vfIn->field, fieldSize[(int)vfIn->lineType]) ;
  vgpWriteLine (vfOut, vfIn->lineType, vgpLen(vfIn), vgpString(vfIn)) ;
  char *s = vgpReadComment (vfIn) ; if (s) vgpWriteComment (vfOut, s) ;
}

int main (int argc, char **argv)
{
  I64 i ;
  FileType fileType = 0 ;
  char *outFileName = "-" ;
  bool isNoHeader = false, isHeaderOnly = false, isBinary = false ;
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
      { isNoHeader = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-H") || !strcmp (*argv, "--headerOnly"))
      { isHeaderOnly = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-b") || !strcmp (*argv, "--binary"))
      { isBinary = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-o") || !strcmp (*argv, "--output"))
      { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-i") || !strcmp (*argv, "--index"))
      { objList = parseIndexList (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-g") || !strcmp (*argv, "--group"))
      { groupList = parseIndexList (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else die ("unknown option %s - run without arguments to see options", *argv) ;

  if (isBinary) isNoHeader = false ;
  if (isHeaderOnly) isBinary = false ;
  
  if (argc != 1) die ("can currently only take one input file") ;

  VgpFile *vfIn = vgpFileOpenRead (*argv, fileType, 1) ; /* reads the header */
  if (!vfIn) die ("failed to open vgp file %s", *argv) ;

  if ((objList || groupList) && !vfIn->isBinary)
    die ("%s is ascii - you can only access objects and groups by index in binary files", *argv) ;
  
  VgpFile *vfOut = vgpFileOpenWriteFrom (outFileName, vfIn, false, isBinary, 1) ;
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
	{ while (objList)
	    { if (!vgpGotoObject (vfIn, objList->i0))
		die ("can't locate to object %lld", objList->i0 ) ;
	      if (!vgpReadLine (vfIn))
		die ("can't read object %lld", objList->i0) ;
	      while (objList->i0 < objList->iN)
		{ transferLine (vfIn, vfOut, fieldSize) ;
		  if (!vgpReadLine (vfIn)) break ;
		  if (vfIn->lineType == vfIn->objectType) ++objList->i0 ;
		}
	      objList = objList->next ;
	    }
	}
      else if (groupList)
	{ while (groupList)
	    { if (!vgpGotoGroup (vfIn, groupList->i0))
		die ("can't locate to group %lld", groupList->i0 ) ;
	      if (!vgpReadLine (vfIn))
		die ("can't read group %lld", groupList->i0) ;
	      while (groupList->i0 < groupList->iN)
		{ transferLine (vfIn, vfOut, fieldSize) ;
		  if (!vgpReadLine (vfIn)) break ;
		  if (vfIn->lineType == vfIn->groupType) ++groupList->i0 ;
		}
	      groupList = groupList->next ;
	    }
	}
      else
	while (vgpReadLine (vfIn))
	  transferLine (vfIn, vfOut, fieldSize) ;
    }
  
  vgpFileClose (vfIn) ;
  vgpFileClose (vfOut) ;

  free (command) ;
  timeTotal (stderr) ;
}

/********************* end of file ***********************/
