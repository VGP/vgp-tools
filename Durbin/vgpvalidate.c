/*  File: vgpvalidate.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Apr 28 17:49 2019 (rd109)
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
  
  --argc ; ++argv ;

  while (argc && **argv == '-')
    if (argc > 1 && !strcmp (*argv, "-t"))
      { for (i = SEQ ; i < MAX_FILE ; ++i)
	  if (!strcmp (argv[1], fileTypeName[i])) fileType = i ;
	if (!fileType) die ("unknown file type %s requested with -t option", argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else die ("unknown option %s", *argv) ;

  if (argc != 1) die ("can currently only take one input file") ;
  
  VgpFile *vf = vgpFileOpenRead (*argv, fileType) ;
  if (!vf) die ("failed to open vgp file *s", *argv) ;
  while (vgpReadLine (vf)) ;
  printf ("read %lld lines from VGP file %s type %s\n", vf->line, *argv, fileTypeName[vf->type]) ;
  vgpFileClose (vf) ;
  
  exit (0) ;
}

/********************* end of file ***********************/
