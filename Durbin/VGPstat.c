/*****************************************************************************************
 *
 *  File: vgpvalidate.c
 *    vgp format validator and header generator
 *
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *
 * HISTORY:
 * Last edited: Dec 27 09:46 2019 (gene)
 *   * Dec 27 09:20 2019 (gene): style edits
 *   * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *
 ****************************************************************************************/

#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "vgprd.h"

int main (int argc, char **argv)
{ int      i;
  FileType fileType;
  char    *outFile;
  BOOL     isHeader, isWrite;
  VgpFile *vfIn;
  
  timeUpdate (0);


  //  Process command line arguments
  
  argc -= 1;
  argv += 1;

  if (argc == 0)
    { fprintf (stderr, "vgpvalidate [options] vgpfile\n") ;
      fprintf (stderr, "  -t --type <abc>         file type, e.g. seq - required if no header\n") ;
      fprintf (stderr, "  -h --header             write out the header\n") ;
      fprintf (stderr, "  -w --write              rewrite in ascii\n") ;
      fprintf (stderr, "  -o --output <filename>  output to filename\n") ;
      fprintf (stderr, "vgpvalidate aborts on a syntactic parse error with a message.\n") ;
      fprintf (stderr, "Otherwise information is written to stderr about any inconsistencies\n") ;
      fprintf (stderr, "between the header and the data in the body of the file.\n") ;
      fprintf (stderr, "You can print out a correct header with -h followed by data section if\n");
      fprintf (stderr, "-w is also set.  Output is to stdout by default, use -o to overide\n");
      exit (0) ;
    }
  
  outFile  = "-";
  isHeader = FALSE;
  isWrite  = FALSE;
  fileType = 0;
  while (argc && **argv == '-')
    { if (!strcmp (*argv, "-h") || !strcmp (*argv, "--header"))
        isHeader = TRUE;
      else if (!strcmp (*argv, "-w") || !strcmp (*argv, "--write"))
        isWrite = TRUE;
      else
        { if (argc > 1 && (!strcmp (*argv, "-t") || !strcmp (*argv, "--type")))
            { for (i = SEQ; i < MAX_FILE; i++)
                if (!strcmp (argv[1], fileTypeName[i]))
                  fileType = i;
              if (fileType == 0)
                die ("unknown file type %s requested", argv[1]) ;
            }
          else if (argc > 1 && (!strcmp (*argv, "-o") || !strcmp (*argv, "--output")))
            outFile = argv[1];
          else
            die ("unknown option %s - run without arguments to see options", *argv);
          argc -= 1;
          argv += 1;
        }
      argc -= 1;
      argv += 1;
    }
  
  if (argc != 1)
    die ("can currently only take one input file");


  //  Open subject file for reading and read header (if present)
  
  vfIn = vgpFileOpenRead (*argv, fileType, 1);
  if (vfIn == NULL)
    die ("failed to open vgp file %s", *argv);
  vfIn->isCheckString = TRUE;

  if (vfIn->line == 1)
    fprintf (stderr, "header missing\n");
  else
    fprintf (stderr, "read %lld header lines\n", vfIn->line);


  //  Read data portion of file checking syntax and group sizes (if present)

  { I64 lastObj, lastSize, lastLine;

    lastObj = lastSize = lastLine = 0;
    while (vgpReadLine (vfIn))
      if (vfIn->lineType == vfIn->groupType)
        { if (lastLine > 0 && vfIn->object - lastObj != lastSize)
            { fprintf (stderr, "group size mismatch: group %c at line %lld asserted %lld objects",
                               vfIn->groupType, lastLine, lastSize);
              fprintf (stderr, " but found %lld\n", vfIn->object-lastObj);
            }
          lastLine = vfIn->line;
          lastSize = vgpInt(vfIn,0);
          lastObj  = vfIn->object;
      }
    if (lastLine && vfIn->object - lastObj != lastSize)
      { fprintf (stderr, "group size mismatch: group %c at line %lld asserted %lld objects",
                         vfIn->groupType, lastLine, lastSize);
        fprintf (stderr, " but found %lld\n", vfIn->object-lastObj);
      }

    fprintf (stderr, "read %lld objects in %lld lines from VGP file %s type %s\n",
             vfIn->object, vfIn->line, *argv, fileTypeName[vfIn->fileType]);

    vgpFinalizeCounts (vfIn);
  }


  //  Check count statistics for each line type versus those in header (if was present)

#define CHECK(X,Y,Z)										\
 if (li->X > 0 && li->X != li->Y)								\
   { fprintf (stderr, "header mismatch %s %c: header %lld data %lld\n", Z, i, li->X, li->Y);	\
     nBad += 1;											\
   } 												\
 else if (li->Y > 0 && li->X == 0)								\
   { fprintf (stderr, "header %s line missing for %c, value is %lld\n", Z, i, li->Y);		\
     nMissing += 1;										\
   } 												\
 if (li->Y > 0)											\
   nTotal += 1;

  { int i, nTotal, nBad, nMissing;

    nTotal = nBad = nMissing = 0;
    for (i = 0; i < 128; i++)
      if (((i >= 'A' && i <= 'Z') || i == vfIn->groupType) && vfIn->lineInfo[i] != NULL)
        { LineInfo *li = vfIn->lineInfo[i];
          CHECK(given.count, accum.count, "count");
          CHECK(given.max, accum.max, "max");
          CHECK(given.total, accum.total, "total");
          CHECK(given.groupCount, accum.groupCount, "group count");
          CHECK(given.groupTotal, accum.groupTotal, "group total");
        }
    fprintf (stderr, "total %d header content lines expected, of which %d bad and %d missing\n",
                     nTotal, nBad, nMissing);
  }


  //  Write header and data as directed to outFile

  if (isHeader || isWrite)
    { VgpFile *vfOut;

      vfOut = vgpFileOpenWriteFrom (outFile, vfIn, TRUE, FALSE, 1);
      if (vfOut == NULL)
        die ("failed to open output file %s", outFile);
  
      if (isHeader)
        { vgpWriteHeader (vfOut);
          fputc ('\n', vfOut->f);
          fflush (vfOut->f);
        }

      if (isWrite)
        { vgpFileClose (vfIn);
          vfIn = vgpFileOpenRead (*argv, fileType, 1);
          if (vfIn == NULL)
            die ("failed to reopen vgp file %s", *argv);

          { int   bufsize, nRead;
            char *buf;

            bufsize = 0x2000000;
            buf     = new (bufsize, char);
            while ((nRead = fread (buf, 1, bufsize, vfIn->f)) > 0)
              if (fwrite (buf, 1, nRead, vfOut->f) != nRead)
                die ("failed to write to outFile");
            free (buf);
          }
        }
      vgpFileClose (vfOut);
    }

  vgpFileClose (vfIn);

  timeTotal (stderr);

  exit (0);
}

/********************* end of file ***********************/
