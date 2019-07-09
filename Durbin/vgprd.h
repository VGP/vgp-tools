/*  File: vgprd.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description: header for VGP file reading and writing
 * Exported functions:
 * HISTORY:
 * Last edited: Jul  7 22:50 2019 (rd109)
 * Created: Sat Feb 23 10:12:43 2019 (rd109)
 *-------------------------------------------------------------------
 */

#ifndef VGPRD_DEFINED
#define VGPRD_DEFINED

#include "vgptypes_1_0.h"	/* must match vgpformat_x_y.h */ 
#include "../Myers/compression.h"

/*** basic types ***/

#include <stdio.h>		/* for FILE etc. */
#include <stdint.h>		/* for standard size int types */
#include <limits.h>		/* for INT_MAX etc. */

typedef char BOOL ;
#define TRUE 1
#define FALSE 0

typedef int64_t I64 ;
const static I64 I64MAX = 0x7fffffffffffffff ;

typedef enum { INT = 1, REAL, CHAR, STRING, INT_LIST, REAL_LIST, STRING_LIST } FieldType ;
typedef union { I64 i ; double r ; char c ; I64 len ; } Field ; /* len for lists */
typedef struct { char *program, *version, *command, *date ; } Provenance ;
typedef struct { char *filename ; I64 count ; } Reference ;

typedef struct {
  FieldType field[MAX_FIELD] ;
  int listByteSize ;
  int nField ;			/* number of fields in this line type */
  int listField ;		/* 1 + field number for length of list (so 0 is no list) */
  //  char *description ;    	/* should really add this */
} LineSpecification ;

typedef struct {
  I64 major, minor ;
  LineSpecification *line[128] ;
  char binaryTypeUnpack[256], binaryTypePack[128] ;
  char objectType ;		/* line designation character for primary objects */
  char groupType ;		/* line designation character for groups (optional) */
} FileSpecification ;

/*** the main VGP file type ***/

typedef struct {
  /* these fields may be read by user - but don't change them! */
  FileType fileType ;
  SubType subType ;
  FileSpecification *spec ;	/* contains the file type and line specifications */
  I64 major, minor ;		/* actual major and minor versions of this file */
  char lineType ;		/* current lineType */
  I64 line ;			/* current line number */
  I64 object ;			/* current object - incremented when object line read */
  I64 group ; 			/* current group - incremented when group line read */
  I64 count[128], max[128], total[128] ; /* amounts read or written so far */
  I64 groupCount[128], groupTotal[128] ; /* assumes only one group type per file; true for now */
  I64 expectCount[128], expectMax[128], expectTotal[128] ; /* values read from header */
  I64 expectGroupCount[128], expectGroupTotal[128] ; /* NB the group* are maxes per group */
  Provenance *provenance ;      /* if non-zero then count['!'] entries */
  Reference *reference ;	/* if non-zero then count['<'] entries */
  Reference *deferred ; 	/* if non-zero then count['>'] entries */
  Field field[MAX_FIELD] ;	/* used to hold the current line - accessed by macros */
  /* fields below here are private to the package */
  FILE *f ;
  BOOL isWrite ;
  BOOL isHeader ;
  BOOL isBinary ;		/* true if writing a binary file */
  BOOL inGroup ;       		/* set once inside a group */
  BOOL isLastLineBinary ;	/* needed to deal with newlines on ascii files */
  void *buffer[128] ;
  I64 bufSize[128] ;
  /* buffer is used for list contents when reading - only nonzero if lineSpec->listType
     buffer can be owned by VgpFile, in which case size is expectMax | max
     or owned by user */
  BOOL isUserBuf[128] ;		/* flag for whether buffer is owned by user */
  I64 gCount[128], gTotal[128] ; /* used internally to calculate groupCount and groupTotal */
  char lineBuf[128], numberBuf[32], *codecBuf ;
  I64 linePos ;
  VGPcodec *fieldCodec[128], *listCodec[128] ;
  BOOL isUseFieldCodec[128], isUseListCodec[128] ;
  /* need the isUse*Codec because we accumulate in the codec object before we can use it */
} VgpFile ;

/*** function definitions for reading and writing VgpFiles ***/

VgpFile *vgpFileOpenRead (const char *path, FileType type) ;
/* opens file and reads header if it exists
   if there is no header then type must be given
   if there is a header and type is non-zero then it must match
   will read ascii or vgp binary files transparently
*/
void vgpClose (VgpFile *vf) ;
/* call this for files opened either read or write
   when writing automatically writes footer, if there was a $ directive in the header
*/
BOOL vgpReadLine (VgpFile *vf) ;
/* this reads the next line and returns FALSE at end of file or on error
   the line is parsed according to its linetype and contents accessed by macros that follow
*/
#define vgpInt(vf,x) ((vf)->field[x].i)
#define vgpReal(vf,i) ((vf)->field[i].r)
#define vgpChar(vf,i) ((vf)->field[i].c)
#define vgpLen(vf,i) ((vf)->field[i].len)
#define vgpString(vf) (char*)((vf)->buffer[(vf)->lineType])
#define vgpIntList(vf) (I64*)((vf)->buffer[(vf)->lineType])
#define vgpRealList(vf) (double*)((vf)->buffer[(vf)->lineType])
void vgpUserBuffer (VgpFile *vf, char lineType, void* buffer) ;
/* this lets the user reassign the buffer that lists are read into
   prior to this function being called, or if buffer==0, a default is provided
   in either case, the buffer will be overwritten each line 
   this can be called repeatedly, so the location can be changed, e.g. for each line
   NB the user must allocate enough memory for user buffers (package buffers are always safe)
*/

VgpFile *vgpFileOpenWriteFrom (const char *path, VgpFile *vfIn, BOOL isBinary) ;
VgpFile *vgpFileOpenWriteNew (const char *path, FileType type, SubType sub, BOOL isBinary) ;
void vgpWriteLine (VgpFile *vf, char lineType, void *buf) ;
/* process is to fill fields by assigning to macros, then call */
/* list contents are in buf - string is 0-terminated, though len must also be set */
/* if lineType is stringList then buf is concatenation of null-terminated strings */
/* NB adds '\n' before writing line not after, so user fprintf() can add extra material */
BOOL vgpInheritProvenance (VgpFile *vf, VgpFile *source) ;
/* add all provenance lines in source - fails after WriteLine called */
BOOL vgpAddProvenance (VgpFile *vf, char *prog, char *version, char *command, char *dateTime) ;
/* if dateTime is 0 then the current date-time is inserted - fails after WriteLine called */
BOOL vgpInheritReference (VgpFile *vf, VgpFile *source) ; /* as for provenance */
BOOL vgpAddReference (VgpFile *vf, char *filename, I64 count) ;
BOOL vgpInheritDeferred (VgpFile *vf, VgpFile *source) ;
BOOL vgpAddDeferred (VgpFile *vf, char *filename) ;
void vgpWriteHeader (VgpFile *vf) ;

/*** utility declarations that RD uses in his parsing code and so users get for free ***/

void die (char *format, ...) ;
void *myalloc (size_t size) ;
void *mycalloc (size_t number, size_t size) ;
#define	new(n,type)	(type*)myalloc((n)*sizeof(type))
#define	new0(n,type)	(type*)mycalloc((n),sizeof(type))
FILE *fzopen (const char* path, const char* mode) ; /* will open gzip files silently */
FILE *fopenTag (char* root, char* tag, char* mode) ; /* uses fzopen, silently handling .gz */
void timeUpdate (FILE *f) ;	/* print time usage since last call to file f */
void timeTotal (FILE *f) ;	/* print full time usage since first call of timeUpdate */

#endif	/* VGPRD_DEFINED */

/******************* end of file **************/
