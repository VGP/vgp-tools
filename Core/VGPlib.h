/******************************************************************************************
 *
 *  File: VGPlib.h
 *    Header for VGP file reading and writing
 *
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *
 * HISTORY:
 * Last edited: Apr 16 20:42 2020 (rd109)
 *   * Dec 27 09:46 2019 (gene): style edits
 *   * Created: Sat Feb 23 10:12:43 2019 (rd109)
 *
 *****************************************************************************************/

#ifndef VGPRD_DEFINED
#define VGPRD_DEFINED

#include <stdio.h>    // for FILE etc.
#include <stdint.h>   // for standard size int types
#include <stdbool.h>  // for true, false, bool
#include <limits.h>   // for INT_MAX etc.
#include <pthread.h>  // for mutex locks

#include "VGPtypes_1_0.h"

/***********************************************************************************
 *
 *    DATA TYPES
 *
 **********************************************************************************/

#ifndef UTILS_DEFINED
  // Basic Types

typedef int64_t       I64;
typedef unsigned char U8;


static const I64 I64MAX = 0x7fffffffffffffffll;
#endif // UTILS_DEFINED

typedef enum { INT = 1, REAL, CHAR, STRING, INT_LIST, REAL_LIST, STRING_LIST } FieldType;

typedef union
  { I64    i;
    double r;
    char   c;
    I64    len; // For lists : top 8 bits encode excess bytes, low 56 length
  } Field;

typedef struct
  { char *program;
    char *version;
    char *command;
    char *date;
  } Provenance;

typedef struct
  { char *filename; 
    I64   count;
  } Reference;

typedef struct
  { I64 count;
    I64 max;
    I64 total;
    I64 groupCount;
    I64 groupTotal;
  } Counts;

  // VGPcodecs are a private package for binary vgp file compression

typedef void VGPcodec; // forward declaration of opaque type for compression codecs

  // DNAcodec is a special pre-existing compressor one should use for DNA.
  // It compresses every base to 2-bits, where any non-ACGT letter is
  // effectively converted to an A.  Compression is case insensitive,
  // but decompression always delivers lower-case.

extern  VGPcodec *DNAcodec;

  // Record for a particular line type.  There is at most one list element.

typedef struct
  { Counts     accum;         // counts read or written to this moment
    Counts     given;         // counts read from header
    I64        gCount;        // used internally to calculate groupCount and groupTotal
    I64        gTotal;
    I64        oCount;        // # of objects in prefix before first group (if any)
    I64        oTotal;        // + of objects in prefix (these 2 needed for thread parallel apps)

    int        nField;        // number of fields
    FieldType  fieldSpec[MAX_FIELD];   //  type of each datum, 0 for [nField..MAX_FIELD-1]
    int        listByteSize;  // size of list field elements (if present)
    int        listField;     // field index of list + 1, 0 if no list

    bool       isUserBuf;     // flag for whether buffer is owned by user
    I64        bufSize;       // system buffer and size if not user supplied
    void      *buffer;

    VGPcodec *fieldCodec;       // compression codecs and flags
    VGPcodec *listCodec;
    bool      isUseFieldCodec;  // on once enough data collected to train associated codec
    bool      isUseListCodec;
    bool      isIntListDiff;    // diff int lists before compressing with codec
    char      binaryTypePack;   // binary code for line type, bit 8 set.
                                //     bit 0: fields compressed
                                //     bit 1: list compressed
    I64       fieldTack;        // accumulated training data for this threads fieldCodec (master)
    I64       listTack;         // accumulated training data for this threads codeCodec (master)
  } LineInfo;

  // The Main VGP file type

typedef struct
  {
    // this field may be set by the user

    bool        isCheckString;       // set if want to validate string char by char, else fread()

    // these fields may be read by user - but don't change them!

    FileType    fileType;
    SubType     subType;
    I64         major;               // actual major and minor versions of this file
    I64         minor;
    char        lineType;            // current lineType
    char        objectType;          // line designation character for primary objects
    char        groupType;           // line designation character for groups (optional)
    I64         line;                // current line number
    I64         byte;                // current byte position when writing binary
    I64         object;              // current object - incremented when object line read
    I64         group;               // current group - incremented when group line read
    Provenance *provenance;          // if non-zero then count['!'] entries
    Reference  *reference;           // if non-zero then count['<'] entries
    Reference  *deferred;            // if non-zero then count['>'] entries
    Field       field[MAX_FIELD];    // used to hold the current line - accessed by macros
    LineInfo   *lineInfo[128];       // all the per-linetype information
    I64         codecTrainingSize;   // amount of data to see before building Huffman code table

    // fields below here are private to the package

    FILE *f;
    bool  isWrite;                // true if open for writing
    bool  isHeaderOut;            // true if header already written
    bool  isBinary;               // true if writing a binary file
    bool  inGroup;                // set once inside a group
    bool  isLastLineBinary;       // needed to deal with newlines on ascii files
    bool  isIndexIn;              // index read in
    bool  isBig;                  // are we on a big-endian machine?
    char  lineBuf[128];           // working buffers
    char  numberBuf[32];
    char *codecBuf;
    I64   codecBufSize;
    I64   linePos;                // current line position

    char  binaryTypeUnpack[256];  // invert binary line code to ASCII line character.
    int   share;                  // index if slave of threaded write, +nthreads > 0 if master
    int   isFinal;                // vgpFinalizeCounts has been called on file
    pthread_mutex_t fieldLock;    // Mutexs to protect training accumumulation stats when threadded
    pthread_mutex_t listLock;
  } VgpFile;                      //   the footer will be in the concatenated result.


/***********************************************************************************
 *
 *    ROUTINES FOR READING & WRITING VGP FILES IN BOTH ASCII & BINARY (TRANSPARENTLY)
 *
 **********************************************************************************/

//  READING VGP FILES:

VgpFile *vgpFileOpenRead (const char *path, FileType type, int nthreads);

  // Open VGP file 'path', either binary or ascii encoded, for reading.
  //   If the file doesn't have a header, then 'type' must be specified,
  //   otherwise, if 'type' is non-zero it must match the header type.
  //   All header information (if present) is read.
  // If nthreads > 1 then nthreadds VgpFiles are generated as an array and the pointer
  //   to the first, called the master, is returned.  The other nthreads-1 files are
  //   called slaves.  The package routines are aware of when a VgpFile argument is a
  //   slave or master in a parallel group.  The master recieves provenance, counts, etc.
  //   The slaves only read data and have the virture of sharing indices and codecs with
  //   the master if relevant.

char vgpReadLine (VgpFile *vf);

  // Read the next VGP formatted line returning the line type of the line, or 0
  //   if at the end of the data section.  The content macros immediately below can be
  //   used to access the information of the line just read.

#define vgpInt(vf,x)        ((vf)->field[x].i)
#define vgpReal(vf,x)       ((vf)->field[x].r)
#define vgpChar(vf,x)       ((vf)->field[x].c)
#define _LF(vf)             ((vf)->lineInfo[(int)(vf)->lineType]->listField)
#define vgpLen(vf)          ((vf)->field[_LF(vf)?_LF(vf)-1:0].len & 0xffffffffffffffll)
#define vgpString(vf)       (char *) ((vf)->lineInfo[(int) (vf)->lineType]->buffer)
#define vgpIntList(vf)      (I64 *) ((vf)->lineInfo[(int) (vf)->lineType]->buffer)
#define vgpRealList(vf)     (double *) ((vf)->lineInfo[(int) (vf)->lineType]->buffer)
#define vgpNextString(vf,s) (s + strlen(s) + 1)

  // Access field information.  The index x of a list object is not required as there is
  //   only one list per line, stored in ->buffer.
  //   A "string list" is implicitly supported, get the first string with vgpString, and
  //   subsequent strings sequentially with vgpNextString, e.g.:
  //
  //       char *s = vgpString(vf);
  //       for (i = 0; i < vgpLen(vf); i++)
  //         { // do something with i'th string
  //           s = vgpNextString(vf,s);
  //         }

char *vgpReadComment (VgpFile *vf);

  // Can be called after vgpReadLine() to read any optional comment text after the fixed fields.
  // Returns NULL if there is no comment.

//  WRITING VGP FILES:

VgpFile *vgpFileOpenWriteNew(const char *path, FileType type, SubType sub,
                             bool isBinary, int nthreads);
VgpFile *vgpFileOpenWriteFrom(const char *path, VgpFile *vfIn, bool useAccum,
                              bool isBinary, int nthreads);

  // Create a new vgpFile that will be written to 'path'.  For the 'New' variant supply
  //   the file type, subtype (if non-zero), and whether it should be binary or ASCII.
  //   For the 'From' variant, specify binary or ASCII, all other header information is
  //   inherited from 'vfIn', where the count stats are from vfIn's accumulation (assumes
  //   vfIn has been fully read or written) if useAccum is true, and from vfIn's header
  //   otherwise.
  // If nthreads > 1 then nthreads VgpFiles are generated as an array and the pointer
  //   to the first, called the master, is returned.  The other nthreads-1 files are
  //   called slaves.  The package routines are aware of when a VgpFile argument is a
  //   slave or master in a parallel group.  The slaves are expected to only write data
  //   lines, with the master adding provenance, producing the header, and then some
  //   segment of the initial data lines.  Upon close the final result is effectively
  //   the concatenation of the master, followed by the output of each slave in sequence.

bool vgpInheritProvenance (VgpFile *vf, VgpFile *source);
bool vgpInheritReference  (VgpFile *vf, VgpFile *source);
bool vgpInheritDeferred   (VgpFile *vf, VgpFile *source);

  // Add all provenance/reference/deferred entries in source to header of vf.  Must be
  //   called before call to vgpWriteHeader.

bool vgpAddProvenance (VgpFile *vf, char *prog, char *version, char *command, char *dateTime);
bool vgpAddReference  (VgpFile *vf, char *filename, I64 count);
bool vgpAddDeferred   (VgpFile *vf, char *filename);

  // Append provenance/reference/deferred to header information.  Must be called before
  //   call to vgpWriteHeader.  Current data & time filled in if dateTime == NULL.

void vgpWriteHeader (VgpFile *vf);

  // Write out the header for file.  For ASCII output, if you want the header to contain
  //   count information then you must create and fill the relevant Counts objects before
  //   calling this. For binary output, the counts will be accumulated and output in a
  //   footer upon vgpClose.

void vgpWriteLine (VgpFile *vf, char lineType, I64 listLen, void *listBuf);

  // Set up a line for output just as it would be returned by vgpReadLine and then call
  //   this routine to output the line (ASCII or binary).
  // Use the macros above on the l.h.s. of assignments to fill fields (e.g. vgpInt(vf,2) = 3).
  // For lists, give the length in the listLen argument, and either place the list data in your
  //   own buffer and give it as listBuf, or put in the line's buffer and set listBuf == NULL.

void vgpWriteComment (VgpFile *vf, char *comment);

  // Adds a comment to the current line. Need to use this not fprintf() so as to keep the
  // index correct in binary mode.

// CLOSING FILES (FOR BOTH READ & WRITE)

void vgpFinalizeCounts (VgpFile *vf);

  // After all input has been read, or all data has been written, this routine will finish
  //   accumulating counts/statistics for the file and merge thread stats into those for
  //   the master file (if a parallel VgpFile).

void vgpFileClose (VgpFile *vf);

  // Close vf (opened either for reading or writing).  Finalizes counts if not explicitly
  //   requested, merges theaded files, and writes footer if binary.   Frees all non-user
  //   memory associated with vf.

//  GOTO & BUFFER MANAGEMENT

void vgpUserBuffer (VgpFile *vf, char lineType, void *buffer);

  // A buffer is used to capture the list element of each line type that has one.
  //   This routine allows you to reassign the buffer to one you've allocated, or
  //   to revert to a default system buffer if 'buffer' = NULL.  The previous buffer
  //   (if any) is freed.  The user must ensure that a buffer they supply is large
  //   enough.  BTW, this buffer is overwritten with each new line read of the given type.

bool vgpGotoObject (VgpFile *vf, I64 i);

  // Goto i'th object in the file.  This only works on binary files, which have an index.

I64  vgpGotoGroup  (VgpFile *vf, I64 i);

  // Goto the first object in group i.  Return the size (in objects) of the group, or 0
  //   if an error (i out of range or vf has not group type).  Only works for binary files.


#ifndef UTILS_DEFINED
//  UTILITIES: RD utilities used in package for your use if desired

void  die(char *format, ...);                  //  print message to stderr and exit -1
void *myalloc(size_t size);                    //  allocate block, NULL if not available
void *mycalloc(size_t number, size_t size);    //  allocate & zero # objects of size, NULL if NA

#define new(n,type)  (type *) myalloc((n)*sizeof(type))
#define new0(n,type)  (type *) mycalloc((n),sizeof(type))

FILE *fzopen(const char* path, const char* mode);   // will open gzip files silently
FILE *fopenTag(char* root, char* tag, char* mode);  // uses fzopen, silently handling .gz
void  timeUpdate(FILE *f);                          // print time and memory usage since last call
void  timeTotal(FILE *f);                           // print full usage since first call
#endif

/***********************************************************************************
 *
 *    A BIT ABOUT THE FORMAT OF BINARY FILES
 *
 **********************************************************************************/

 //   <bin file> <- <ASCII Prolog> <$-line> <binary data> <footer> <^-line> <footer-size:int64>
 //
 // '$'-line flags file is binary and gives endian
 // The data block ends with a blank line consisting of '\n'
 //
 // EWM: Removed '-' line, simply write off_t to footer start
 //
 //   <ASCII Prolog> <- <'1'-line> [<'2'-line>] ( <'!'-line> | <'<'-line> | <'>'-line> )*
 //
 // The ASCII prolog contains the type, subtype, provenance, reference, and deferred lines
 //   in the ASCII format.  The VGP count statistic lines for each data line type are found
 //   in the footer along with binary '\01' and '\02' lines that encode their compressors as
 //   needed.  The footer also contains binary '&' and '*' lines that encode the object index
 //   and group indices, respectively.
 //
 //   <Binary line> <- <Binary line code + tags> <fields> [<list data>]
 //
 // Line codes are >= 128 for binary encoded lines.  The low two order bits of these are flags,
 //   so each binary-encoded line type has 4 codes and a table maps these to the ASCII code.
 //   Bit 0 indicates if the fields of the line type are compressed, and Bit 1 indicates if
 //   the list data (if present) is compressed.
 //
 // If a field is a list, then the field array element for that field is the list's length
 //   where the low 56 bits encode length, and the high 8 bits encode the # of high-order
 //   0-bytes in every list element if an INT_LIST (0 otherwise).

#endif  // VGPRD_DEFINED

/******************* end of file **************/
