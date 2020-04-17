/*  Last edited: Apr 16 22:08 2020 (rd109) */
/*******************************************************************************************
 *
 *  VGPpipe: Convert VGP files from ASCII to binary or vice versa
 *
 *  Author:  Gene Myers
 *  Date  :  Jan. 1, 2020
 *
 ********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include <pthread.h>
#include <fcntl.h>
#include <sys/stat.h>

#include "gene_core.h"
#include "../include/VGPlib.h"

static char *Usage = "<in >out";


  //  Main

int main(int argc, char* argv[])
{
  //  Process command line arguments

  { int   i, j, k;
    int   flags[128];
    
    ARG_INIT("VGPpipe")
    
    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("")
            break;
        }   
      else
        argv[j++] = argv[i];
    argc = j;
    
    if (argc != 1)
      { fprintf(stderr,"Usage: %s %s\n",Prog_Name,Usage);
        exit (1);
      } 
  } 

  { VgpFile  *inf, *ouf;
    LineInfo *li;
    int       i, t;

    inf = vgpFileOpenRead("-",0,1);
    ouf = vgpFileOpenWriteFrom("-",inf,FALSE,!inf->isBinary,1);
    vgpWriteHeader(ouf);
    while ((t = vgpReadLine(inf)) > 0)
      { li = inf->lineInfo[t];
        for (i = 0; i < li->nField; i++)
          ouf->field[i] = inf->field[i];
        vgpWriteLine(ouf,t,vgpLen(inf),li->buffer);
      }
    vgpFileClose(ouf);
    vgpFileClose(inf);
  }

  exit (0);
}
