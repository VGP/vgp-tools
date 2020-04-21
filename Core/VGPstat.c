/*****************************************************************************************
 *
 *  File: VGPstat.c
 *    vgp format validator and header generator
 *
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *
 * HISTORY:
 * Last edited: Apr 16 20:45 2020 (rd109)
 *   * Dec 27 09:20 2019 (gene): style edits
 *   * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *
 ****************************************************************************************/

#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "VGPlib.h"

int main (int argc, char **argv)
{ int      i ;
  FileType fileType = 0 ;
  char    *outFileName = "-" ;
  BOOL     isHeader = FALSE, isUsage = FALSE ;
  VgpFile *vfIn ;
  
  timeUpdate (0) ;

  //  Process command line arguments

  --argc ; ++argv ;		/* drop the program name */
  if (argc == 0)
    { fprintf (stderr, "VGPstat [options] vgpfile\n") ;
      fprintf (stderr, "  -t --type <abc>         file type, e.g. seq - required if no header\n") ;
      fprintf (stderr, "  -H --header             output header accumulated from data\n") ;
      fprintf (stderr, "  -o --output <filename>  output to filename\n") ;
      fprintf (stderr, "  -u --usage              byte usage per line type; no other output\n") ;
      fprintf (stderr, "VGPstat aborts on a syntactic parse error with a message.\n") ;
      fprintf (stderr, "Otherwise information is written to stderr about any inconsistencies\n") ;
      fprintf (stderr, "between the header and the data in the body of the file.\n") ;
      fprintf (stderr, "Output is to stdout by default, use -o to overide\n");
      exit (0) ;
    }
  
  fileType = 0;
  while (argc && **argv == '-')
    if (!strcmp (*argv, "-H") || !strcmp (*argv, "--header"))
      { isHeader = TRUE ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-u") || !strcmp (*argv, "--usage"))
      { isUsage = TRUE ; --argc ; ++argv ; }
    else if (argc > 1 && (!strcmp (*argv, "-t") || !strcmp (*argv, "--type")))
      { for (i = SEQ; i < MAX_FILE; i++)
	  if (!strcmp (argv[1], fileTypeName[i]))
	    fileType = i;
	if (fileType == 0)
	  die ("unknown file type %s requested", argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else if (argc > 1 && (!strcmp (*argv, "-o") || !strcmp (*argv, "--output")))
      { outFileName = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else die ("unknown option %s - run without arguments to see options", *argv) ;
  
  if (argc != 1)
    die ("can currently only take one input file") ;


  //  Open subject file for reading and read header (if present)
  
  vfIn = vgpFileOpenRead (*argv, fileType, 1) ;
  if (vfIn == NULL)
    die ("failed to open vgp file %s", *argv) ;
  vfIn->isCheckString = TRUE ;

  if (vfIn->line == 1)
    fprintf (stderr, "header missing\n") ;
  else
    fprintf (stderr, "read %lld header lines\n", vfIn->line) ;

  // if requesting usage, then 

  if (isUsage)
    { I64 usage[128] ; memset (usage, 0, 128*sizeof(I64)) ; 
      off_t u, uLast = ftello (vfIn->f) ;

      while (vgpReadLine (vfIn))
	{ u = ftello (vfIn->f) ; usage[(int)vfIn->lineType] += u-uLast ; uLast = u ; }
      u = ftello (vfIn->f) ; usage[(int)vfIn->lineType] += u-uLast ; uLast = u ;

      FILE *f ;
      if (strcmp (outFileName, "-") && !(f = fopen (outFileName, "w")))
	die ("failed to open output file %s", outFileName) ;
      else
	f = stdout ;
      
      for (i = 'A' ; i < 128 ; ++i)
	if (usage[i]) fprintf (f, "usage line type %c bytes %lld\n", (char)i, usage[i]) ;

      if (f != stdout) fclose (f) ;
     }

  else
    {
      //  Read data portion of file checking syntax and group sizes (if present)

      { I64 lastObj = 0, lastSize = 0, lastLine = 0 ;
	
	lastObj = lastSize = lastLine = 0 ;
	while (vgpReadLine (vfIn))
	  if (vfIn->lineType == vfIn->groupType)
	    { if (lastLine > 0 && vfIn->object - lastObj != lastSize)
		{ fprintf (stderr, "group size mismatch: group %c at line %lld asserted %lld objects",
			   vfIn->groupType, lastLine, lastSize);
		  fprintf (stderr, " but found %lld\n", vfIn->object-lastObj) ;
		}
	      lastLine = vfIn->line ;
	      lastSize = vgpInt(vfIn,0) ;
	      lastObj  = vfIn->object ;
	    }
	if (lastLine && vfIn->object - lastObj != lastSize)
	  { fprintf (stderr, "group size mismatch: group %c at line %lld asserted %lld objects",
		     vfIn->groupType, lastLine, lastSize) ;
	    fprintf (stderr, " but found %lld\n", vfIn->object-lastObj) ;
	  }
      }
      
      fprintf (stderr, "read %lld objects in %lld lines from VGP file %s type %s\n",
	       vfIn->object, vfIn->line, *argv, fileTypeName[vfIn->fileType]) ;

      vgpFinalizeCounts (vfIn) ;
    
      //  Check count statistics for each line type versus those in header (if was present)

      { I64 nTotal = 0, nBad = 0, nMissing = 0 ;
	  
#define CHECK(X,Y,Z)							                       \
  if (li->X > 0 && li->X != li->Y)				                               \
    { fprintf (stderr, "header mismatch %s %c: header %lld data %lld\n", Z, i, li->X, li->Y) ; \
      nBad += 1 ;							                       \
   } 											       \
 else if (li->Y > 0 && li->X == 0)							       \
   { fprintf (stderr, "header %s line missing for %c, value is %lld\n", Z, i, li->Y) ;	       \
     nMissing += 1 ;									       \
   } 											       \
 if (li->Y > 0)										       \
   nTotal += 1 ;

	for (i = 0; i < 128; i++)
	  if (((i >= 'A' && i <= 'Z') || i == vfIn->groupType) && vfIn->lineInfo[i] != NULL)
	    { LineInfo *li = vfIn->lineInfo[i] ;
	      CHECK(given.count, accum.count, "count") ;
	      CHECK(given.max, accum.max, "max") ;
	      CHECK(given.total, accum.total, "total") ;
	      CHECK(given.groupCount, accum.groupCount, "group count") ;
	      CHECK(given.groupTotal, accum.groupTotal, "group total") ;
	  }
	fprintf (stderr, "expected %lld header content lines, of which %lld bad and %lld missing\n",
		 nTotal, nBad, nMissing) ;
      }

  //  Write header if requested

      if (isHeader)
	{ VgpFile *vfOut = vgpFileOpenWriteFrom (outFileName, vfIn, TRUE, FALSE, 1) ;
	  if (vfOut == NULL)
	    die ("failed to open output file %s", outFileName) ;
  
	  vgpWriteHeader (vfOut) ;
	  fflush (vfOut->f) ;
	  vgpFileClose (vfOut) ;
	}
    }

  vgpFileClose (vfIn) ;

  timeTotal (stderr) ;

  exit (0) ;
}

/********************* end of file ***********************/
