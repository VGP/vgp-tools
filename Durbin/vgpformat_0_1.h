/*  File: vgpformat_0_1.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description: this contains the VGP format specification 0.1
 *   the idea is to include this file into vgprd.c, so it is compiled into the executable
 *   it is a separate file so as to make it easy to change the format and version
 * Exported functions:
 * HISTORY:
 * Last edited: May  1 00:05 2019 (rd109)
 * Created: Sun Feb 24 14:48:21 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "vgprd.h"

static LineSpecification *vgpDefineLine (FieldType f0, FieldType f1, FieldType f2, FieldType f3)
{
  LineSpecification *ls = new0 (1, LineSpecification) ;
  ls->field[0] = f0 ; ls->field[1] = f1 ; ls->field[2] = f2 ; ls->field[3] = f3 ;
  return ls ;
}

static FileSpecification *vgpDefineFormat (void)
{
  static int MajorVersion = 1 ;
  static int MinorVersion = 0 ;
  
  int i, j, k ;
  FileSpecification *fileSpec = new0 (MAX_FILE, FileSpecification) ;
  LineSpecification **header = new0 (128, LineSpecification*) ;

  /* header line type 1 is treated separately below */
  header['2'] = vgpDefineLine (STRING, 0, 0, 0) ;
  header['#'] = vgpDefineLine (CHAR, INT, 0, 0) ;
  header['@'] = vgpDefineLine (CHAR, INT, 0, 0) ;
  header['+'] = vgpDefineLine (CHAR, INT, 0, 0) ;
  header['%'] = vgpDefineLine (CHAR, CHAR, INT, 0) ;
  header['!'] = vgpDefineLine (STRING_LIST, 0, 0, 0) ;
  header['<'] = vgpDefineLine (STRING, INT, 0, 0) ;
  header['>'] = vgpDefineLine (STRING, 0, 0, 0) ;
  /* NB strictly !, < and > lines contain strings so isList should be TRUE for them */

  for (i = 0 ; i < MAX_FILE ; ++i)
    { fileSpec[i].major = MajorVersion ; fileSpec[i].minor = MinorVersion ;
      if (i > 0) for (j = 0 ; j < 128 ; ++j) if (header[j]) fileSpec[i].line[j] = header[j] ;
    }

  fileSpec[0].line['1'] = vgpDefineLine (STRING, INT, INT, 0) ;

  fileSpec[SEQ].objectType = 'S' ;
  fileSpec[SEQ].line['S'] = vgpDefineLine (STRING, 0, 0, 0) ;
  fileSpec[SEQ].line['Q'] = vgpDefineLine (STRING, 0, 0, 0) ;
  fileSpec[SEQ].line['P'] = vgpDefineLine (0, 0, 0, 0) ;
  fileSpec[SEQ].line['W'] = vgpDefineLine (INT, INT, INT, REAL) ;
  fileSpec[SEQ].line['N'] = vgpDefineLine (REAL, REAL, REAL, REAL) ;
  fileSpec[SEQ].line['A'] = vgpDefineLine (STRING, 0, 0, 0) ;
  fileSpec[SEQ].line['g'] = vgpDefineLine (INT, STRING, 0, 0) ;

  fileSpec[RMD].objectType = 'R' ;
  fileSpec[RMD].line['R'] = vgpDefineLine (INT, INT_LIST, 0, 0) ;
  fileSpec[RMD].line['E'] = vgpDefineLine (INT_LIST, 0, 0, 0) ;
  fileSpec[RMD].line['I'] = vgpDefineLine (REAL_LIST, 0, 0, 0) ;
  fileSpec[RMD].line['N'] = vgpDefineLine (REAL_LIST, 0, 0, 0) ;
  fileSpec[RMD].line['r'] = vgpDefineLine (INT, STRING_LIST, 0, 0) ;

  fileSpec[SXS].objectType = 'X' ;
  fileSpec[SXS].line['A'] = vgpDefineLine (INT, INT, 0, 0) ;
  fileSpec[SXS].line['I'] = vgpDefineLine (INT, INT, INT, INT) ;
  fileSpec[SXS].line['Q'] = vgpDefineLine (INT, 0, 0, 0) ;
  fileSpec[SXS].line['M'] = vgpDefineLine (INT, 0, 0, 0) ;
  fileSpec[SXS].line['D'] = vgpDefineLine (INT, 0, 0, 0) ;
  fileSpec[SXS].line['C'] = vgpDefineLine (STRING, 0, 0, 0) ;
  fileSpec[SXS].line['T'] = vgpDefineLine (INT_LIST, 0, 0, 0) ;
  fileSpec[SXS].line['U'] = vgpDefineLine (INT_LIST, 0, 0, 0) ;
  fileSpec[SXS].line['V'] = vgpDefineLine (INT_LIST, 0, 0, 0) ;

  fileSpec[JNS].objectType = 'J' ;
  fileSpec[JNS].line['J'] = vgpDefineLine (INT, INT, CHAR, CHAR) ;
  fileSpec[JNS].line['G'] = vgpDefineLine (INT, 0, 0, 0) ;
  fileSpec[JNS].line['Q'] = vgpDefineLine (INT, 0, 0, 0) ;
  fileSpec[JNS].line['X'] = vgpDefineLine (INT_LIST, 0, 0, 0) ;
  
  fileSpec[BRK].objectType = 'B' ;
  fileSpec[BRK].line['B'] = vgpDefineLine (INT, INT, INT, 0) ;
  fileSpec[BRK].line['Q'] = vgpDefineLine (INT, 0, 0, 0) ;
  fileSpec[BRK].line['X'] = vgpDefineLine (INT_LIST, 0, 0, 0) ;
  
  fileSpec[JNS].objectType = 'F' ;
  fileSpec[JNS].line['F'] = vgpDefineLine (INT, CHAR, INT_LIST, 0) ;
  fileSpec[JNS].line['N'] = vgpDefineLine (INT, INT_LIST, 0, 0) ;

  fileSpec[LIS].line['L'] = vgpDefineLine (INT_LIST, 0, 0, 0) ;

  /* fill in the listByteSize information */
  int listSize[8] = { 0, 0, 0, 0, 1, sizeof(I64), sizeof(double), 1 } ;
  for (i = 0 ; i < MAX_FILE ; ++i)
    for (j = 0 ; j < 128 ; ++j)
      if (fileSpec[i].line[j])
	for (k = 0 ; k < MAX_FIELD ; k++)
	  { if (listSize[fileSpec[i].line[j]->field[k]] && fileSpec[i].line[j]->listByteSize)
	      die ("VGP format %s spec error: two list types in record %c", fileTypeName[i], j) ;
	    else
	      fileSpec[i].line[j]->listByteSize = listSize[fileSpec[i].line[j]->field[k]] ;
	  }
  printf ("format defined\n") ;
  return fileSpec ;
}

/********* end of file ***********/
