#ifndef _CORE

#define _CORE

#include <stdio.h>

typedef unsigned char      uint8;
typedef unsigned short     uint16;
typedef unsigned int       uint32;
typedef unsigned long long uint64;
typedef signed char        int8;
typedef signed short       int16;
typedef signed int         int32;
typedef signed long long   int64;
typedef float              float32;
typedef double             float64;

extern char *Prog_Name;   //  Name of program

#define ARG_INIT(name)                  \
  Prog_Name = Strdup(name,"");          \
  for (i = 0; i < 128; i++)             \
    flags[i] = 0;

#define ARG_FLAGS(set)                                                                  \
  for (k = 1; argv[i][k] != '\0'; k++)                                                  \
    { if (index(set,argv[i][k]) == NULL)                                                \
        { fprintf(stderr,"%s: -%c is an illegal option\n",Prog_Name,argv[i][k]);        \
          exit (1);                                                                     \
        }                                                                               \
      flags[(int) argv[i][k]] = 1;                                                      \
    }

#define ARG_POSITIVE(var,name)                                                          \
  var = strtol(argv[i]+2,&eptr,10);                                                     \
  if (*eptr != '\0' || argv[i][2] == '\0')                                              \
    { fprintf(stderr,"%s: -%c '%s' argument is not an integer\n",			\
                     Prog_Name,argv[i][1],argv[i]+2);      				\
      exit (1);                                                                         \
    }                                                                                   \
  if (var <= 0)                                                                         \
    { fprintf(stderr,"%s: %s must be positive (%d)\n",Prog_Name,name,var);              \
      exit (1);                                                                         \
    }

#define ARG_NON_NEGATIVE(var,name)                                                      \
  var = strtol(argv[i]+2,&eptr,10);                                                     \
  if (*eptr != '\0' || argv[i][2] == '\0')                                              \
    { fprintf(stderr,"%s: -%c '%s' argument is not an integer\n",			\
                     Prog_Name,argv[i][1],argv[i]+2);      				\
      exit (1);                                                                         \
    }                                                                                   \
  if (var < 0)	                                                                        \
    { fprintf(stderr,"%s: %s must be non-negative (%d)\n",Prog_Name,name,var);          \
      exit (1);                                                                         \
    }

#define ARG_REAL(var)                                                                   \
  var = strtod(argv[i]+2,&eptr);                                                        \
  if (*eptr != '\0' || argv[i][2] == '\0')                                              \
    { fprintf(stderr,"%s: -%c '%s' argument is not a real number\n",			\
                     Prog_Name,argv[i][1],argv[i]+2);      				\
      exit (1);                                                                         \
    }

//  The following general utilities return NULL if any of their input pointers are NULL, or if they
//    could not perform their function (in which case they also print an error to stderr).

void *Malloc(int64 size, char *mesg);                    //  Guarded versions of malloc, realloc
void *Realloc(void *object, int64 size, char *mesg);     //  and strdup, that output "mesg" to
char *Strdup(char *string, char *mesg);                  //  stderr if out of memory

FILE *Fzopen(char *path, char *mode);    // Open file path for "mode", will directly read .gz
char *PathTo(char *path);                // Return path portion of file name "path"
char *Root(char *path, char *suffix);    // Return the root name, excluding suffix, of "path"

#define OPEN(arg,pwd,root,input,suffix,nsuf)                    \
  for (i = 0; i < nsuf; i++)                                    \
    { root  = Root(arg,suffix[i]);                              \
      input = Fzopen(Catenate(pwd,"/",root,suffix[i]),"r");     \
      if (input != NULL) break;                                 \
      free(root);                                               \
    }

// Catenate returns concatenation of path.sep.root.suffix in a *temporary* buffer
// Numbered_Suffix returns concatenation of left.<num>.right in a *temporary* buffer

char *Catenate(char *path, char *sep, char *root, char *suffix);
char *Numbered_Suffix(char *left, int num, char *right);

void Print_Number(int64 num, int width, FILE *out);   //  Print readable big integer
int  Number_Digits(int64 num);                        //  Return # of digits in printed number

#define COMPRESSED_LEN(len)  (((len)+3) >> 2)

void   Compress_Read(int len, char *s);   //  Compress read in-place into 2-bit form
void Uncompress_Read(int len, char *s);   //  Uncompress read in-place into numeric form
void      Print_Read(char *s, int width);

void Lower_Read(char *s);     //  Convert read from numbers to lowercase letters (0-3 to acgt)
void Upper_Read(char *s);     //  Convert read from numbers to uppercase letters (0-3 to ACGT)
void Number_Read(char *s);    //  Convert read from letters to numbers
void Change_Read(char *s);    //  Convert read from one case to the other

void Letter_Arrow(char *s);   //  Convert arrow pw's from numbers to uppercase letters (0-3 to 1234)
void Number_Arrow(char *s);   //  Convert arrow pw string from letters to numbers

#endif // _CORE
