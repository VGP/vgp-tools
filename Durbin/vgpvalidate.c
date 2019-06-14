/*  File: vgpvalidate.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jun 14 10:46 2019 (rd109)
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
  FILE *outFile = stdout ;
  BOOL isHeader = FALSE, isWrite = FALSE ;
  
  timeUpdate (0) ;
  
  --argc ; ++argv ;

  if (!argc)
    { fprintf (stderr, "vgpvalidate [options] vgpfile\n") ;
      fprintf (stderr, "  -t --type <abc>         file type, e.g. seq, aln - required if no header\n") ;
      fprintf (stderr, "  -h --header             write out the header\n") ;
      fprintf (stderr, "  -w --write              rewrite the full file\n") ;
      fprintf (stderr, "  -o --output <filename>  output file name\n") ;
      fprintf (stderr, "vgpvalidate aborts with a syntactic parse error with a (hopefully informative) message.\n") ;
      fprintf (stderr, "Then information is written to stderr about any inconsistencies between the header and\n") ;
      fprintf (stderr, "the data in the body of the file.\n") ;
      fprintf (stderr, "You can print out a correct header (default to stdout if not set with -o), and also\n") ;
      fprintf (stderr, "rewrite the entire file with the header using '-w -h', which takes a second pass.\n") ;
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
      { isWrite = TRUE ; --argc ; ++argv ;
	fprintf (stderr, "Sorry, rewriting option is not present yet.  Write out a header and cat the body of the data onto it as a temporary measure.\n") ;
      }
    else if (argc > 1 && (!strcmp (*argv, "-o") || !strcmp (*argv, "--output")))
      { if (!(outFile = fopen (argv[1], "w"))) die ("failed to open output file %s", argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else die ("unknown option %s - run without arguments to see options", *argv) ;

  if (argc != 1) die ("can currently only take one input file") ;
  
  VgpFile *vf = vgpFileOpenRead (*argv, fileType) ; /* reads the header */
  if (!vf) die ("failed to open vgp file %s", *argv) ;

  if (vf->line == 1) fprintf (stderr, "header missing\n") ;
  else fprintf (stderr, "read %lld header lines\n", vf->line) ;

  while (vgpReadLine (vf))
    {
    }
  printf ("read %lld objects in %lld lines from VGP file %s type %s\n",
	  vf->count[vf->spec->objectType], vf->line, *argv, fileTypeName[vf->type]) ;

  {
#define CHECK(X,Y,Z) if (vf->X[i] && vf->X[i] != vf->Y[i])		\
      { fprintf (stderr, "header mismatch %s %c: header %lld data %lld\n", Z, i, vf->X[i], vf->Y[i]) ; ++nBad ; } \
    else if (vf->Y[i] && !vf->X[i]) { fprintf (stderr, "header %s line missing for %c, value is %lld\n", Z, i, vf->Y[i]) ; ++nMissing ; } \
      if (vf->Y[i]) ++nTotal

    int i, nTotal = 0, nBad = 0, nMissing = 0 ;
    for (i = 0 ; i < 128 ; ++i)
      if ((i >= 'A' && i <= 'Z') || (i >= 'a' && i <= 'z'))
	{ CHECK(expectCount, count, "count") ;
	  CHECK(expectMax, max, "max") ;
	  CHECK(expectTotal, total, "total") ;
	  CHECK(expectGroupCount, groupCount, "group count") ;
	  CHECK(expectGroupTotal, groupTotal, "group total") ;
	}
    fprintf (stderr, "total %d header content lines expected, of which %d bad and %d missing\n",
	     nTotal, nBad, nMissing) ;
  }
  
  if (isHeader) { vgpWriteHeader (vf, outFile) ; fprintf (outFile, "\n") ; }

  vgpFileClose (vf) ;

  timeTotal (stderr) ;
}

/********************* end of file ***********************/
