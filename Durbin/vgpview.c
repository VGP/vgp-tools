/*  File: vgpvalidate.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jul  7 15:00 2019 (rd109)
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

int main (int argc, char **argv)
{
  int i ;
  FileType fileType = 0 ;
  char *outFileName = "-" ;
  BOOL isHeader = FALSE, isHeaderOnly = FALSE, isBinary = FALSE ;
  
  timeUpdate (0) ;

  char *command = commandLine (argc, argv) ;
  --argc ; ++argv ;		/* drop the program name */

  if (!argc)
    { fprintf (stderr, "vgpview [options] vgpfile\n") ;
      fprintf (stderr, "  -t --type <abc>         file type, e.g. seq, aln - required if no header\n") ;
      fprintf (stderr, "  -h --header             include the header in ascii output\n") ;
      fprintf (stderr, "  -H --headerOnly         only write the header (in ascii)\n") ;
      fprintf (stderr, "  -b --binary             write in binary (default is ascii)\n") ;
      fprintf (stderr, "  -o --output <filename>  output file name (default stdout)\n") ;
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
    else die ("unknown option %s - run without arguments to see options", *argv) ;

  if (isBinary) isHeader = TRUE ;
  if (isHeaderOnly) { isHeaderOnly = TRUE ; isBinary = FALSE ; }
  
  if (argc != 1) die ("can currently only take one input file") ;
  
  VgpFile *vfIn = vgpFileOpenRead (*argv, fileType) ; /* reads the header */
  if (!vfIn) die ("failed to open vgp file %s", *argv) ;

  VgpFile *vfOut = vgpFileOpenWriteFrom (outFileName, vfIn, isBinary) ;
  if (!vfOut) die ("failed to open output file %s", outFileName) ;

  if (isHeaderOnly)
    { vgpWriteHeader (vfOut) ; fputc ('\n', vfOut->f) ;
    }
  else
    { vgpAddProvenance (vfOut, "vgpview", "0.0", command, 0) ;
      if (isHeader) vgpWriteHeader (vfOut) ;

      static size_t fieldSize[128] ;
      for (i = 0 ; i < 128 ; ++i)
	if (vfIn->spec->line[i]) fieldSize[i] = vfIn->spec->line[i]->nField*sizeof(Field) ;
      
      while (vgpReadLine (vfIn))
	{ memcpy (vfOut->field, vfIn->field, fieldSize[vfIn->lineType]) ;
	  vgpWriteLine (vfOut, vfIn->lineType, vfIn->buffer[vfIn->lineType]) ;
	}
    }

  vgpClose (vfOut) ;

  timeTotal (stderr) ;
}

/********************* end of file ***********************/
