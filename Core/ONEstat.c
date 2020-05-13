/*****************************************************************************************
 *
 *  File: ONEstat.c
 *    one format validator and header generator
 *
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *
 * HISTORY:
 * Last edited: May 12 19:36 2020 (rd109)
 *   * Dec 27 09:20 2019 (gene): style edits
 *   * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *
 ****************************************************************************************/

#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "utils.h"
#include "ONElib.h"

extern void oneFinalizeCounts (OneFile *vf) ; // secret connection into ONElib.c for checking

int main (int argc, char **argv)
{ int        i ;
  char      *fileType = 0 ;
  char      *outFileName = "-" ;
  bool       isHeader = false, isUsage = false, isVerbose = false ;
  char      *schemaFileName = 0 ;
  char      *checkText = 0 ;
  
  timeUpdate (0) ;

  //  Process command line arguments

  --argc ; ++argv ;		/* drop the program name */
  if (argc == 0)
    { fprintf (stderr, "ONEstat [options] onefile\n") ;
      fprintf (stderr, "  -t --type <abc>          file type, e.g. seq - required if no header\n") ;
      fprintf (stderr, "  -S --schema <schema>     schema file - required if not in file\n") ;
      fprintf (stderr, "  -C --check 'schematext'  check for a limited set of features\n") ;
      fprintf (stderr, "  -H --header              output header accumulated from data\n") ;
      fprintf (stderr, "  -o --output <filename>   output to filename\n") ;
      fprintf (stderr, "  -u --usage               byte usage per line type; no other output\n") ;
      fprintf (stderr, "  -v --verbose             else only errors and requested output\n") ;
      fprintf (stderr, "ONEstat aborts on a syntactic parse error with a message.\n") ;
      fprintf (stderr, "Otherwise information is written to stderr about any inconsistencies\n") ;
      fprintf (stderr, "between the header and the data in the body of the file.\n") ;
      fprintf (stderr, "Output is to stdout by default, use -o to overide\n");
      exit (0) ;
    }
  
  while (argc && **argv == '-')
    if (!strcmp (*argv, "-H") || !strcmp (*argv, "--header"))
      { isHeader = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-u") || !strcmp (*argv, "--usage"))
      { isUsage = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-v") || !strcmp (*argv, "--verbose"))
      { isVerbose = true ; --argc ; ++argv ; }
    else if (argc > 1 && (!strcmp (*argv, "-t") || !strcmp (*argv, "--type")))
      { fileType = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (argc > 1 && (!strcmp (*argv, "-S") || !strcmp (*argv, "--schema")))
      { schemaFileName = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (argc > 1 && (!strcmp (*argv, "-C") || !strcmp (*argv, "--check")))
      { checkText = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (argc > 1 && (!strcmp (*argv, "-o") || !strcmp (*argv, "--output")))
      { outFileName = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else die ("unknown option %s - run without arguments to see options", *argv) ;
  
  if (argc != 1)
    die ("need to give a single data file as argument") ;

  //  Open subject file for reading and read header (if present)
  OneSchema *vs = 0 ;
  if (schemaFileName)
    { vs = oneSchemaCreateFromFile (schemaFileName) ;
      if (!vs) die ("failed to read schema file %s", schemaFileName) ;
    }
  OneFile *vf = oneFileOpenRead (argv[0], vs, fileType, 1) ;
  if (!vf) die ("failed to open OneFile %s", argv[0]) ;
  oneSchemaDestroy (vs) ; // no longer needed

  if (isVerbose)
    { if (vf->line == 1)
	fprintf (stderr, "header missing\n") ;
      else
	fprintf (stderr, "read %" PRId64 " header lines\n", vf->line) ;
    }

  if (checkText)
    oneFileCheckSchema (vf, checkText) ;

  vf->isCheckString = true ;

  // if requesting usage, then 

  if (isUsage)
    { I64 usage[128] ; memset (usage, 0, 128*sizeof(I64)) ; 
      off_t u, uLast = ftello (vf->f) ;

      while (oneReadLine (vf))
	{ u = ftello (vf->f) ; usage[(int)vf->lineType] += u-uLast ; uLast = u ; }
      u = ftello (vf->f) ; usage[(int)vf->lineType] += u-uLast ; uLast = u ;

      FILE *f ;
      if (strcmp (outFileName, "-") && !(f = fopen (outFileName, "w")))
	die ("failed to open output file %s", outFileName) ;
      else
	f = stdout ;
      
      for (i = 'A' ; i < 128 ; ++i)
	if (usage[i]) fprintf (f, "usage line type %c bytes %" PRId64 "\n", (char)i,  usage[i]) ;

      if (f != stdout) fclose (f) ;
     }

  else
    {
      //  Read data portion of file checking syntax and group sizes (if present)

      { I64 lastObj = 0, lastSize = 0, lastLine = 0 ;
	
	lastObj = lastSize = lastLine = 0 ;
	while (oneReadLine (vf))
	  if (vf->lineType == vf->groupType)
	    { if (lastLine > 0 && vf->object - lastObj != lastSize)
		{ fprintf (stderr, "group size mismatch: group %c at line %" PRId64 " asserted %" PRId64 " objects",
			   vf->groupType,  lastLine,  lastSize);
		  fprintf (stderr, " but found %" PRId64 "\n", vf->object-lastObj) ;
		}
	      lastLine = vf->line ;
	      lastSize = oneInt(vf,0) ;
	      lastObj  = vf->object ;
	    }
	if (lastLine && vf->object - lastObj != lastSize)
	  { fprintf (stderr, "group size mismatch: group %c at line %" PRId64 " asserted %" PRId64 " objects",
		     vf->groupType,  lastLine,  lastSize) ;
	    fprintf (stderr, " but found %" PRId64 "\n",  vf->object-lastObj) ;
	  }
      }

      if (isVerbose)
	fprintf (stderr, "read %" PRId64 " objects in %" PRId64 " lines from OneFile %s type %s\n",
		  vf->object,  vf->line, argv[0], vf->fileType) ;

      oneFinalizeCounts (vf) ;
    
      //  Check count statistics for each line type versus those in header (if was present)

      { I64 nTotal = 0, nBad = 0, nMissing = 0 ;
	  
#define CHECK(X,Y,Z)							                       \
  if (li->X > 0 && li->X != li->Y)				                               \
    { fprintf (stderr, "header mismatch %s %c: header %" PRId64 " data %" PRId64 "\n", Z, i,  li->X,  li->Y) ; \
      nBad += 1 ;							                       \
   } 											       \
 else if (li->Y > 0 && li->X == 0)							       \
   { fprintf (stderr, "header %s line missing for %c, value is %" PRId64 "\n", Z, i,  li->Y) ;	       \
     nMissing += 1 ;									       \
   } 											       \
 if (li->Y > 0)										       \
   nTotal += 1 ;

	for (i = 0; i < 128; i++)
	  if (((i >= 'A' && i <= 'Z') || i == vf->groupType) && vf->info[i] != NULL)
	    { OneInfo *li = vf->info[i] ;
	      CHECK(given.count, accum.count, "count") ;
	      CHECK(given.max, accum.max, "max") ;
	      CHECK(given.total, accum.total, "total") ;
	      CHECK(given.groupCount, accum.groupCount, "group count") ;
	      CHECK(given.groupTotal, accum.groupTotal, "group total") ;
	  }
	if (isVerbose || nBad || nMissing)
	  fprintf (stderr, "expected %" PRId64 " header content lines, of which %" PRId64 " bad and %" PRId64 " missing\n",
		    nTotal,  nBad,  nMissing) ;
      }

  //  Write header if requested

      if (isHeader)
	{ OneFile *vfOut = oneFileOpenWriteFrom (outFileName, vf, false, 1) ;
	  if (vfOut == NULL)
	    die ("failed to open output file %s", outFileName) ;

	  for (i = 0 ; i < 128 ; i++) // transfer accumulated counts for vfIn to given for vfOut
	    if (vfOut->info[i])
	      vfOut->info[i]->given = vf->info[i]->accum ;
  
	  oneWriteHeader (vfOut) ;
	  fflush (vfOut->f) ;
	  oneFileClose (vfOut) ;
	}
    }

  oneFileClose (vf) ;

  if (isVerbose) timeTotal (stderr) ;

  exit (0) ;
}

/********************* end of file ***********************/
