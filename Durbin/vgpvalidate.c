/*  File: vgpvalidate.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Feb 22 14:52 2019 (rd109)
 * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"

typedef enum { NULL_FILE = 0, SEQ, LIS, IRP, PBR, X10, RMD, SXS, RXR, SCF, MAX_FILE } FileType ;

char *fileTypeString[] = { 0, "seq", "lis", "irp", "pbr", "10x", "rmd", "sxs", "rxr", "scf" } ;

typedef enum { NULL_FIELD = 0, INT, REAL, CHAR, STRING, INT_LIST, REAL_LIST, STRING_LIST } FieldType ;

typedef struct {
  FieldType field[6] ;		/* 0-terminated list of fields */
  BOOL isList ;
} LineInfo ;

typedef struct {
  char *name ;
  LineInfo *line[256] ;
} FileInfo ;

U64 *count, *max, *total ;
FileInfo *fileInfo ;
LineInfo **header ;
FileType fileType = NULL_FILE ;
U64 major = -1, minor = -1 ;

void defineFormat (void)
{
  int i, j ;
  
  header = new0 (128, LineInfo*) ;
  header['1'] = &(LineInfo) { { STRING, INT, INT, 0, 0, 0 }, FALSE } ;
  header['#'] = &(LineInfo) { { CHAR, INT, 0, 0, 0, 0 }, FALSE } ;
  header['@'] = &(LineInfo) { { CHAR, INT, 0, 0, 0, 0 }, FALSE } ;
  header['+'] = &(LineInfo) { { CHAR, INT, 0, 0, 0, 0 }, FALSE } ;
  header['%'] = &(LineInfo) { { CHAR, CHAR, CHAR, INT, 0, 0 }, FALSE } ;
  header['!'] = &(LineInfo) { { STRING, STRING, STRING, STRING, 0, 0 }, FALSE } ;
  header['<'] = &(LineInfo) { { STRING, CHAR, INT, 0, 0, 0 }, FALSE } ;
  header['>'] = &(LineInfo) { { STRING, 0, 0, 0, 0, 0 }, FALSE } ;
  /* strictly !, < and > lines contain strings so isList should be TRUE for them */

  fileInfo = new0 (MAX_FILE, FileInfo) ;
  for (i = 0 ; i < MAX_FILE ; ++i)
    { fileInfo[i].name = fileTypeString[i] ;
      for (j = 0 ; j < 128 ; ++j) if (header[j]) fileInfo[i].line[j] = header[j] ;
    }

  fileInfo[SEQ].line['S'] = &(LineInfo) { { STRING, 0, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[SEQ].line['I'] = &(LineInfo) { { STRING, 0, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[SEQ].line['Q'] = &(LineInfo) { { STRING, 0, 0, 0, 0, 0 }, TRUE } ;

  fileInfo[LIS].line['L'] = &(LineInfo) { { INT_LIST, 0, 0, 0, 0, 0 }, TRUE } ;

  fileInfo[IRP].line['g'] = &(LineInfo) { { INT, STRING, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[IRP].line['P'] = &(LineInfo) { { 0, 0, 0, 0, 0, 0 }, FALSE } ;
  fileInfo[IRP].line['S'] = &(LineInfo) { { STRING, 0, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[IRP].line['Q'] = &(LineInfo) { { STRING, 0, 0, 0, 0, 0 }, TRUE } ;

  fileInfo[X10].line['c'] = &(LineInfo) { { INT, STRING, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[X10].line['P'] = &(LineInfo) { { 0, 0, 0, 0, 0, 0 }, FALSE } ;
  fileInfo[X10].line['S'] = &(LineInfo) { { STRING, 0, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[X10].line['Q'] = &(LineInfo) { { STRING, 0, 0, 0, 0, 0 }, TRUE } ;
  /* NB syntactically .10x files are the same as .irp files */

  fileInfo[PBR].line['g'] = &(LineInfo) { { INT, STRING, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[PBR].line['S'] = &(LineInfo) { { STRING, 0, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[PBR].line['W'] = &(LineInfo) { { INT, INT, INT, REAL, 0, 0 }, FALSE } ;
  fileInfo[PBR].line['N'] = &(LineInfo) { { REAL, REAL, REAL, REAL, 0, 0 }, TRUE } ;
  fileInfo[PBR].line['A'] = &(LineInfo) { { STRING, 0, 0, 0, 0, 0 }, FALSE } ;

  fileInfo[RMD].line['r'] = &(LineInfo) { { INT, STRING_LIST, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[RMD].line['M'] = &(LineInfo) { { INT, INT_LIST, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[RMD].line['R'] = &(LineInfo) { { INT_LIST, 0, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[RMD].line['I'] = &(LineInfo) { { REAL_LIST, 0, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[RMD].line['N'] = &(LineInfo) { { REAL_LIST, 0, 0, 0, 0, 0 }, TRUE } ;

  fileInfo[SXS].line['A'] = &(LineInfo) { { INT, INT, 0, 0, 0, 0 }, FALSE } ;
  fileInfo[SXS].line['I'] = &(LineInfo) { { INT, INT, INT, INT, 0, 0 }, FALSE } ;
  fileInfo[SXS].line['Q'] = &(LineInfo) { { INT, 0, 0, 0, 0, 0 }, FALSE } ;
  fileInfo[SXS].line['C'] = &(LineInfo) { { STRING, 0, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[SXS].line['T'] = &(LineInfo) { { INT_LIST, 0, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[SXS].line['M'] = &(LineInfo) { { INT, 0, 0, 0, 0, 0 }, FALSE } ;
  fileInfo[SXS].line['D'] = &(LineInfo) { { INT, 0, 0, 0, 0, 0 }, FALSE } ;

  fileInfo[RXR].line['A'] = &(LineInfo) { { INT, INT, 0, 0, 0, 0 }, FALSE } ;
  fileInfo[RXR].line['I'] = &(LineInfo) { { INT, INT, INT, INT, 0, 0 }, FALSE } ;
  fileInfo[RXR].line['Q'] = &(LineInfo) { { INT, 0, 0, 0, 0, 0 }, FALSE } ;
  fileInfo[RXR].line['M'] = &(LineInfo) { { INT, 0, 0, 0, 0, 0 }, FALSE } ;
  fileInfo[RXR].line['U'] = &(LineInfo) { { INT_LIST, 0, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[RXR].line['V'] = &(LineInfo) { { INT_LIST, 0, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[RXR].line['C'] = &(LineInfo) { { STRING, 0, 0, 0, 0, 0 }, TRUE } ;
  fileInfo[RXR].line['D'] = &(LineInfo) { { INT, 0, 0, 0, 0, 0 }, FALSE } ;
  
  fileInfo[SCF].line['B'] = &(LineInfo) { { INT, INT, 0, 0, 0, 0 }, FALSE } ;
  fileInfo[SCF].line['Q'] = &(LineInfo) { { INT, INT, INT, INT, 0, 0 }, FALSE } ;
  fileInfo[SCF].line['J'] = &(LineInfo) { { INT, 0, 0, 0, 0, 0 }, FALSE } ;
  fileInfo[SCF].line['G'] = &(LineInfo) { { INT, 0, 0, 0, 0, 0 }, FALSE } ;
}

static inline char readChar (FILE *f) ;
static inline U64 readInt (FILE *f) ;
static inline double readReal (FILE *f) ;
static inline char* readString (FILE *f) ;
static inline void readFlush (FILE *f) ; /* reads to the end of the line */

int main (int argc, char **argv)
{
  int i, line = 0 ;
  FILE *in = stdin ;
  FILE *out = stdout ;
  
  --argc ; ++argv ;
  defineFormat () ;

  while (argc && **argv == '-')
    if (argc > 1 && !strcmp (*argv, "-t"))
      { for (i = SEQ ; i < MAX_FILE ; ++i)
	  if (!strcmp (argv[1], fileTypeString[i])) fileType = i ;
	if (!fileType) die ("unknown file type %s requested with -t option", argv[1]) ;
	argc -= 2 ; argv += 2 ;
      }
    else die ("unknown option %s", *argv) ;
  if (argc > 1) die ("can currently only take one input file") ;
  
  if (argc && !(in = fzopen (*argv, "r"))) die ("failed to open input file %s", *argv) ;

  BOOL isFirst = TRUE ;
  BOOL isHeader = FALSE ;
  FileInfo *finfo ;
  while (!feof (in))
    { char x = readChar(in) ;
      if (isFirst)
	{ isFirst = FALSE ;
	  if (x == '.')
	    { char *suffix = readString (in) ;
	      for (i = SEQ ; i < MAX_FILE ; ++i)
		if (!strcmp (suffix, fileTypeString[i])) fileType = i ;
	      major = readInt (in) ; minor = readInt (in) ;
	    }
	  if (!fileType) die ("file type not given on first line; must specify with -t option") ;
	  finfo = &fileInfo[fileType] ;
	  finfo->line['1'] = 0 ; /* don't allow to parse first line any more */
	}
      readFlush (in) ;
    }
    
  fclose (in) ;
  fclose (out) ;
  
  exit (0) ;
}

static inline char readChar (FILE *f)
{ char x = getc (f) ; char y = getc (f) ; if (y == '\n') ungetc (y, f) ; return x ; }

static inline U64 readInt (FILE *f) { return (0); }
static inline double readReal (FILE *f) { return (0.); }
static inline char* readString (FILE *f) { return (NULL); }

static inline void readFlush (FILE *f) /* reads to the end of the line */
{ char x ; while ((x = getc (f)) && x != '\n') if (x == EOF) die ("premature end of file") ; }
