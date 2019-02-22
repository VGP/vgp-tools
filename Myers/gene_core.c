#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <dirent.h>

#include "gene_core.h"

/*******************************************************************************************
 *
 *  GENERAL UTILITIES
 *
 ********************************************************************************************/

char *Prog_Name;

void *Malloc(int64 size, char *mesg)
{ void *p;

  if ((p = malloc(size)) == NULL)
    { if (mesg == NULL)
        fprintf(stderr,"%s: Out of memory\n",Prog_Name);
      else
        fprintf(stderr,"%s: Out of memory (%s)\n",Prog_Name,mesg);
    }
  return (p);
}

void *Realloc(void *p, int64 size, char *mesg)
{ if (size <= 0)
    size = 1;
  if ((p = realloc(p,size)) == NULL)
    { if (mesg == NULL)
        fprintf(stderr,"%s: Out of memory\n",Prog_Name);
      else
        fprintf(stderr,"%s: Out of memory (%s)\n",Prog_Name,mesg);
    }
  return (p);
}

char *Strdup(char *name, char *mesg)
{ char *s;

  if (name == NULL)
    return (NULL);
  if ((s = strdup(name)) == NULL)
    { if (mesg == NULL)
        fprintf(stderr,"%s: Out of memory\n",Prog_Name);
      else
        fprintf(stderr,"%s: Out of memory (%s)\n",Prog_Name,mesg);
    }
  return (s);
}

FILE *Fopen(char *name, char *mode)
{ FILE *f;

  if (name == NULL || mode == NULL)
    return (NULL);
  if ((f = fopen(name,mode)) == NULL)
    fprintf(stderr,"%s: Cannot open %s for '%s'\n",Prog_Name,name,mode);
  return (f);
}

char *PathTo(char *name)
{ char *path, *find;

  if (name == NULL)
    return (NULL);
  if ((find = rindex(name,'/')) != NULL)
    { *find = '\0';
      path = Strdup(name,"Extracting path from");
      *find = '/';
    }
  else
    path = Strdup(".","Allocating default path");
  return (path);
}

char *Root(char *name, char *suffix)
{ char *path, *find, *dot;
  int   epos;

  if (name == NULL)
    return (NULL);
  find = rindex(name,'/');
  if (find == NULL)
    find = name;
  else
    find += 1;
  if (suffix == NULL)
    { dot = strchr(find,'.');
      if (dot != NULL)
        *dot = '\0';
      path = Strdup(find,"Extracting root from");
      if (dot != NULL)
        *dot = '.';
    }
  else
    { epos  = strlen(find);
      epos -= strlen(suffix);
      if (epos > 0 && strcasecmp(find+epos,suffix) == 0)
        { find[epos] = '\0';
          path = Strdup(find,"Extracting root from");
          find[epos] = suffix[0];
        }
      else
        path = Strdup(find,"Allocating root");
    }
  return (path);
}

char *Catenate(char *path, char *sep, char *root, char *suffix)
{ static char *cat = NULL;
  static int   max = -1;
  int len;

  if (path == NULL || root == NULL || sep == NULL || suffix == NULL)
    return (NULL);
  len =  strlen(path);
  len += strlen(sep);
  len += strlen(root);
  len += strlen(suffix);
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      if ((cat = (char *) realloc(cat,max+1)) == NULL)
        { fprintf(stderr,"%s: Out of memory (Making path name for %s)\n",Prog_Name,root);
          return (NULL);
        }
    }
  sprintf(cat,"%s%s%s%s",path,sep,root,suffix);
  return (cat);
}

char *Numbered_Suffix(char *left, int num, char *right)
{ static char *suffix = NULL;
  static int   max = -1;
  int len;

  if (left == NULL || right == NULL)
    return (NULL);
  len =  strlen(left);
  len += strlen(right) + 40;
  if (len > max)
    { max = ((int) (1.2*len)) + 100;
      if ((suffix = (char *) realloc(suffix,max+1)) == NULL)
        { fprintf(stderr,"%s: Out of memory (Making number suffix for %d)\n",Prog_Name,num);
          return (NULL);
        }
    }
  sprintf(suffix,"%s%d%s",left,num,right);
  return (suffix);
}


#define  COMMA  ','

//  Print big integers with commas/periods for better readability

void Print_Number(int64 num, int width, FILE *out)
{ if (width == 0)
    { if (num < 1000ll)
        fprintf(out,"%lld",num);
      else if (num < 1000000ll)
        fprintf(out,"%lld%c%03lld",num/1000ll,COMMA,num%1000ll);
      else if (num < 1000000000ll)
        fprintf(out,"%lld%c%03lld%c%03lld",num/1000000ll,
                                           COMMA,(num%1000000ll)/1000ll,COMMA,num%1000ll);
      else
        fprintf(out,"%lld%c%03lld%c%03lld%c%03lld",num/1000000000ll,
                                                   COMMA,(num%1000000000ll)/1000000ll,
                                                   COMMA,(num%1000000ll)/1000ll,COMMA,num%1000ll);
    }
  else
    { if (num < 1000ll)
        fprintf(out,"%*lld",width,num);
      else if (num < 1000000ll)
        { if (width <= 4)
            fprintf(out,"%lld%c%03lld",num/1000ll,COMMA,num%1000ll);
          else
            fprintf(out,"%*lld%c%03lld",width-4,num/1000ll,COMMA,num%1000ll);
        }
      else if (num < 1000000000ll)
        { if (width <= 8)
            fprintf(out,"%lld%c%03lld%c%03lld",num/1000000ll,COMMA,(num%1000000ll)/1000ll,
                                               COMMA,num%1000ll);
          else
            fprintf(out,"%*lld%c%03lld%c%03lld",width-8,num/1000000ll,COMMA,(num%1000000ll)/1000ll,
                                                COMMA,num%1000ll);
        }
      else
        { if (width <= 12)
            fprintf(out,"%lld%c%03lld%c%03lld%c%03lld",num/1000000000ll,COMMA,
                                                       (num%1000000000ll)/1000000ll,COMMA,
                                                       (num%1000000ll)/1000ll,COMMA,num%1000ll);
          else
            fprintf(out,"%*lld%c%03lld%c%03lld%c%03lld",width-12,num/1000000000ll,COMMA,
                                                        (num%1000000000ll)/1000000ll,COMMA,
                                                        (num%1000000ll)/1000ll,COMMA,num%1000ll);
        }
    }
}

//  Return the number of digits, base 10, of num

int  Number_Digits(int64 num)
{ int digit;

  digit = 0;
  while (num >= 1)
    { num /= 10;
      digit += 1;
    }
  return (digit);
}
