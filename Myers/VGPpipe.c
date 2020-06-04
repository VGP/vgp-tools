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
#include <stdbool.h>

#include "gene_core.h"
#include "../Core/ONElib.h"

#include "VGPschema.h"

static char *Usage = "<in >out";


  //  Main

int main(int argc, char* argv[])
{ OneSchema *schema;

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

    schema = oneSchemaCreateFromText(vgpSchemaText);
  } 

  { OneFile  *inf, *ouf;
    OneInfo  *li;
    int       i, t;

    inf = oneFileOpenRead("-",schema,NULL,1);
    ouf = oneFileOpenWriteFrom("-",inf,!inf->isBinary,1);
    oneWriteHeader(ouf);
    while ((t = oneReadLine(inf)) > 0)
      { li = inf->info[t];
        for (i = 0; i < li->nField; i++)
          ouf->field[i] = inf->field[i];
        oneWriteLine(ouf,t,oneLen(inf),li->buffer);
      }
    oneFileClose(ouf);
    oneFileClose(inf);
  }

  oneSchemaDestroy(schema);

  exit (0);
}
