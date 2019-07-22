/*  File: vgpvalidate.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jul 20 17:11 2019 (rd109)
 * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "vgprd.h"

#include <assert.h>

#include <string.h>		/* strcmp etc. */
#include <stdlib.h>		/* for exit() */

int main (int argc, char **argv)
{
  int i ;
  FileType fileType = 0 ;
  char *outFileName = "-" ;
  BOOL isHeader = FALSE, isWrite = FALSE ;
  
  timeUpdate (0) ;
  
  --argc ; ++argv ;

  if (!argc)
    { fprintf (stderr, "vgpvalidate [options] vgpfile\n") ;
      fprintf (stderr, "  -t --type <abc>         file type, e.g. seq, aln - required if no header\n") ;
      fprintf (stderr, "  -h --header             write out the header\n") ;
      fprintf (stderr, "  -w --write              rewrite in ascii\n") ;
      fprintf (stderr, "  -o --output <filename>  output file name\n") ;
      fprintf (stderr, "vgpvalidate aborts on a syntactic parse error with a (hopefully informative) message.\n") ;
      fprintf (stderr, "Otherwise information is written to stderr about any inconsistencies between the header and\n") ;
      fprintf (stderr, "the data in the body of the file.\n") ;
      fprintf (stderr, "You can print out a correct header with -h (default to stdout if not set with -o),\n") ;
      fprintf (stderr, "and also rewrite the entire file with the header using '-w -h'.\n") ;
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
    else if (!strcmp (*argv, "-w") || !strcmp (*argv, "--write"))
      { isWrite = TRUE ; --argc ; ++argv ; }
    else if (argc > 1 && (!strcmp (*argv, "-o") || !strcmp (*argv, "--output")))
      { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }
    else die ("unknown option %s - run without arguments to see options", *argv) ;
  
  if (argc != 1) die ("can currently only take one input file") ;
  
  VgpFile *vfIn = vgpFileOpenRead (*argv, fileType) ; /* reads the header */
  if (!vfIn) die ("failed to open vgp file %s", *argv) ;
  vfIn->isCheckString = TRUE ;

  if (vfIn->line == 1) fprintf (stderr, "header missing\n") ;
  else fprintf (stderr, "read %lld header lines\n", vfIn->line) ;

  I64 lastGroupObject = 0, lastGroupSize, lastGroupLine = 0 ;
  char lastGroupChar ;
  while (vgpReadLine (vfIn))
    if (vfIn->lineType == vfIn->groupType) /* a group line */
      { if (lastGroupLine && vfIn->object - lastGroupObject != lastGroupSize)
	  fprintf (stderr, "group size mismatch: group %c line %lld said %lld objects, but found %lld\n",
		   vfIn->lineType, lastGroupLine, lastGroupSize, vfIn->object-lastGroupObject) ;
	lastGroupLine = vfIn->line ;
	lastGroupSize = vgpInt(vfIn,0) ;
	lastGroupChar = vfIn->lineType ;
	lastGroupObject = vfIn->object ;
      }
  if (lastGroupLine && vfIn->object - lastGroupObject != lastGroupSize)
    fprintf (stderr, "last group size mismatch: group %c line %lld said %lld objects, but found %lld\n",
	     lastGroupChar, lastGroupLine, lastGroupSize, vfIn->object-lastGroupObject) ;

  fprintf (stderr, "read %lld objects in %lld lines from VGP file %s type %s\n",
	   vfIn->object, vfIn->line, *argv, fileTypeName[vfIn->fileType]) ;

  {
#define CHECK(X,Y,Z) if (li->X && li->X != li->Y)		\
      { fprintf (stderr, "header mismatch %s %c: header %lld data %lld\n", Z, i, li->X, li->Y) ; ++nBad ; } \
    else if (li->Y && !li->X) { fprintf (stderr, "header %s line missing for %c, value is %lld\n", Z, i, li->Y) ; ++nMissing ; } \
      if (li->Y) ++nTotal

    int i, nTotal = 0, nBad = 0, nMissing = 0 ;
    for (i = 0 ; i < 128 ; ++i)
      if (((i >= 'A' && i <= 'Z') || i == vfIn->groupType) && vfIn->lineInfo[i])
	{ LineInfo *li = vfIn->lineInfo[i] ;
	  CHECK(expectCount, count, "count") ;
	  CHECK(expectMax, max, "max") ;
	  CHECK(expectTotal, total, "total") ;
	  CHECK(expectGroupCount, groupCount, "group count") ;
	  CHECK(expectGroupTotal, groupTotal, "group total") ;
	}
    fprintf (stderr, "total %d header content lines expected, of which %d bad and %d missing\n",
	     nTotal, nBad, nMissing) ;
  }

  if (isHeader || isWrite)
    { VgpFile *vfOut = vgpFileOpenWriteFrom (outFileName, vfIn, FALSE) ;
      if (!vfOut) die ("failed to open output file %s", outFileName) ;
  
      if (isHeader)
	{ for (i = 0 ; i < 128 ; ++i)
	    if (vfIn->lineInfo[i] && vfIn->lineInfo[i]->count)
	      { LineInfo *liIn = vfIn->lineInfo[i], *liOut = vfOut->lineInfo[i] ;
		liOut->expectCount = liIn->count ;
		liOut->expectMax = liIn->max ;
		liOut->expectTotal = liIn->total ;
		liOut->expectGroupCount = liIn->groupCount ;
		liOut->expectGroupTotal = liIn->groupTotal ;
	      }
	  vgpWriteHeader (vfOut) ; fputc ('\n', vfOut->f) ; fflush (vfOut->f) ;
	}

      if (isWrite)
	{ vgpClose (vfIn) ;
	  vfIn = vgpFileOpenRead (*argv, fileType) ; /* reread the header */
	  if (!vfIn) die ("failed to reopen vgp file %s", *argv) ;
	  int bufsize = 2 << 24, nRead ;
	  char *buf = new (bufsize, char) ;
	  while ((nRead = fread (buf, 1, bufsize, vfIn->f)))  /* copy remainder of vfIn to vfOut */
	    if (fwrite (buf, 1, nRead, vfOut->f) != nRead) die ("failed to write to outFile") ;
	  free (buf) ;
	}
      vgpClose (vfOut) ;
    }

  vgpClose (vfIn) ;

  timeTotal (stderr) ;
}

/********************* end of file ***********************/
