# The ONE-Code C-library Interface

### Authors:  Gene Myers, Richard Durbin, and the Vertebrate Genome Project Assembly Group
### Last Update: April 13, 2020


The interface is defined in `ONElib.h` (in subdirectory Core/).  There are 19 functions and 8 macros, with one primary type `OneFile` which maintains information about the current line.

As a brief synopsis, the following reads a sequence file, prints out some simple stats, and writes a binary file containing the reverse-complemented sequences.

```
 { OneFile *in = oneFileOpenRead(inFile, "seq", 1);
	if (!in)
	  { fprintf(stderr,"Can't open sequence file %s to read\n",inFile);
	    exit (1);
	  }

	OneFile *out = oneFileOpenWriteFrom(outFile, in, FALSE, TRUE, 1);
	if (!out)
	  { fprintf(stderr,"Can't open sequence file %s to write\n",outFile);
	  	 exit (1);
	  }
	oneAddProvenance(out,"revcomp","1.0","revcomp inFile outFile",0);
	oneWriteHeader(out);
	
	int totLen = 0, totCount = 0;
	while (oneReadLine(in))
	  if (in->lineType == 'S')
	    { totLen += oneLen(in);
	      reverseComplement(oneString(in), oneLen(in)); // user-provided, assume acts in place
	      oneWriteLine(out, 'S', onepLen(in), oneString(in));
	    }
	  else if (in->lineType == 'C')
	    totCount += oneInt(in,0);
	printf("total sequence length %d and counts %d\n", totLen, totCount);
	printf("auto-accumulated length %d should be the same\n", in->lineInfo['S']->accum.total);
	oneFileClose(in);
	oneFileClose(out);
 }
```

For now, more details are provided in `ONElib.h`.


```
#define TRUE  1
#define FALSE 0

typedef char          BOOL;
typedef int64_t       I64;
typedef unsigned char U8;

static const I64 I64MAX = 0x7fffffffffffffffll;
```

```
typedef enum { vINT = 1, vREAL, vCHAR, vSTRING, vINT_LIST, vREAL_LIST, vSTRING_LIST } OneType;

typedef union
  { I64    i;
    double r;
    char   c;
    I64    len;
  } OneField;

typedef struct
  { char *program;
    char *version;
    char *command;
    char *date;
  } OneProvenance;

typedef struct
  { char *filename; 
    I64   count;
  } OneReference;

typedef struct
  { I64 count;
    I64 max;
    I64 total;
    I64 groupCount;
    I64 groupTotal;
  } OneCounts;

  // OneCodecs are a private package for binary one file compression

typedef void OneCodec; // forward declaration of opaque type for compression codecs

  // DNAcodec is a special pre-existing compressor one should use for DNA.
  // It compresses every base to 2-bits, where any non-ACGT letter is
  // effectively converted to an A.  Compression is case insensitive,
  // but decompression always delivers lower-case.

extern  OneCodec *DNAcodec;

  // Record for a particular line type.  There is at most one list element.

typedef struct
  { OneCounts  accum;           // counts read or written to this moment
    OneCounts  given;           // counts read from header
    I64        gCount;          // used internally to calculate groupCount and groupTotal
    I64        gTotal;
    I64        oCount;          // # of objects in prefix before first group (if any)
    I64        oTotal;          // + of objects in prefix (these 2 are for thread parallel apps)

    int        nField;          // number of fields
    OneType   *fieldType;       // type of each field
    int        listEltSize;     // size of list field elements (if present, else 0)
    int        listField;       // field index of list
    BOOL       isIntListDiff;   // diff int lists before compressing with codec

    BOOL       isUserBuf;       // flag for whether buffer is owned by user
    I64        bufSize;         // system buffer and size if not user supplied
    void      *buffer;

    OneCodec *fieldCodec;       // compression codecs and flags
    OneCodec *listCodec;
    BOOL      isUseFieldCodec;  // on once enough data collected to train associated codec
    BOOL      isUseListCodec;
    char      binaryTypePack;   // binary code for line type, bit 8 set.
                                //     bit 0: fields compressed
                                //     bit 1: list compressed
    I64       fieldTack;        // accumulated training data for this threads fieldCodec (master)
    I64       listTack;         // accumulated training data for this threads codeCodec (master)
  } OneInfo;

  // the schema type - the first record is the header spec, then a linked list of primary classes

typedef struct OneSchema
  { int               major, minor;
    char              primary[4];
    int               nSecondary;
    char            **secondary;
    OneInfo          *info[128];
    int               nFieldMax;
    char              objectType;
    char              groupType;
    int               nBinary;  // number of line types which allow binary encoding
    int               nBinaryHeader ; // need to start counting nBinary from here
    struct OneSchema *nxt;
  } OneSchema ;

  // The main OneFile type - this is the primary handle used by the end user

typedef struct
  {
    // this field may be set by the user

    BOOL           isCheckString;      // set if want to validate string char by char

    // these fields may be read by user - but don't change them!

    char           fileType[4];
    char           subType[4];
    I64            major;              // actual major and minor versions of this file
    I64            minor;
    char           lineType;           // current lineType
    char           objectType;         // line designation character for primary objects
    char           groupType;          // line designation character for groups (optional)
    I64            line;               // current line number
    I64            byte;               // current byte position when writing binary
    I64            object;             // current object - incremented when object line read
    I64            group;              // current group - incremented when group line read
    OneProvenance *provenance;         // if non-zero then count['!'] entries
    OneReference  *reference;          // if non-zero then count['<'] entries
    OneReference  *deferred;           // if non-zero then count['>'] entries
    OneField      *field;              // used to hold the current line - accessed by macros
    OneInfo       *info[128];          // all the per-linetype information
    I64            codecTrainingSize;  // amount of data to see before building codec

    // fields below here are private to the package

    FILE *f;
    BOOL  isWrite;                // true if open for writing
    BOOL  isHeaderOut;            // true if header already written
    BOOL  isBinary;               // true if writing a binary file
    BOOL  inGroup;                // set once inside a group
    BOOL  isLastLineBinary;       // needed to deal with newlines on ascii files
    BOOL  isIndexIn;              // index read in
    BOOL  isBig;                  // are we on a big-endian machine?
    char  lineBuf[128];           // working buffers
    char  numberBuf[32];
    int   nFieldMax;
    I64   codecBufSize;
    char *codecBuf;
    I64   linePos;                // current line position

    char  binaryTypeUnpack[256];  // invert binary line code to ASCII line character.
    int   share;                  // index if slave of threaded write, +nthreads > 0 if master
    int   isFinal;                // oneFinalizeCounts has been called on file
    pthread_mutex_t fieldLock;    // Mutexs to protect training accumumulation stats when threadded
    pthread_mutex_t listLock;
  } OneFile;                      //   the footer will be in the concatenated result.
```


# THE ROUTINES

## Creating and Destroying Schemas

```
OneSchema *oneSchemaCreateFromFile (char *path);
void oneSchemaDestroy (OneSchema *vs);
```

## Reading ONE Files

```
OneFile *oneFileOpenRead (const char *path, OneSchema *vs, char *type, int nthreads);
```

Open ONE file 'path', either binary or ascii encoded, for reading.
If the file doesn't have a header, then 'type' must be specified,
otherwise, if 'type' is non-zero it must match the header type.
All header information (if present) is read.
If nthreads > 1 then nthreadds OneFiles are generated as an array and the pointer
to the first, called the master, is returned.  The other nthreads-1 files are
called slaves.  The package routines are aware of when a OneFile argument is a
slave or master in a parallel group.  The master recieves provenance, counts, etc.
The slaves only read data and have the virture of sharing indices and codecs with
the master if relevant.

```
char oneReadLine (OneFile *vf);
```

Read the next ONE formatted line returning the line type of the line, or 0
if at the end of the data section.  The content macros immediately below can be
used to access the information of the line just read.

```
#define oneInt(vf,x)        ((vf)->field[x].i)
#define oneReal(vf,x)       ((vf)->field[x].r)
#define oneChar(vf,x)       ((vf)->field[x].c)
#define _LF(vf)             ((vf)->info[(int)(vf)->lineType]->listField)
#define oneLen(vf)          ((vf)->field[_LF(vf)].len & 0xffffffffffffffll)
#define oneString(vf)       (char *) ((vf)->info[(int) (vf)->lineType]->buffer)
#define oneIntList(vf)      (I64 *) ((vf)->info[(int) (vf)->lineType]->buffer)
#define oneRealList(vf)     (double *) ((vf)->info[(int) (vf)->lineType]->buffer)
#define oneNextString(vf,s) (s + strlen(s) + 1)
```

Access field information.  The index x of a list object is not required as there is
only one list per line, stored in ->buffer.
A "string list" is implicitly supported, get the first string with oneString, and
subsequent strings sequentially with oneNextString, e.g.:

```
        char *s = oneString(vf);
        for (i = 0; i < oneLen(vf); i++)
          { // do something with i'th string
            s = oneNextString(vf,s);
          }
```

```
char *oneReadComment (OneFile *vf);
```

Can be called after oneReadLine() to read any optional comment text after the fixed fields.
Returns NULL if there is no comment.

## Writing One Files

```
OneFile *oneFileOpenWriteNew (const char *path, OneSchema *vs, char *type,
			      BOOL isBinary, int nthreads);
OneFile *oneFileOpenWriteFrom (const char *path, OneSchema *vs, OneFile *vfIn,
			       BOOL useAccum, BOOL isBinary, int nthreads);
```		

Create a new oneFile that will be written to 'path'.  For the 'New' variant supply
the file type, subtype (if non-zero), and whether it should be binary or ASCII.
For the 'From' variant, specify binary or ASCII, all other header information is
inherited from 'vfIn', where the count stats are from vfIn's accumulation (assumes
vfIn has been fully read or written) if useAccum is true, and from vfIn's header
otherwise.

If nthreads > 1 then nthreads OneFiles are generated as an array and the pointer
to the first, called the master, is returned.  The other nthreads-1 files are
called slaves.  The package routines are aware of when a OneFile argument is a
slave or master in a parallel group.  The slaves are expected to only write data
lines, with the master adding provenance, producing the header, and then some
segment of the initial data lines.  Upon close the final result is effectively
the concatenation of the master, followed by the output of each slave in sequence.

```
BOOL oneInheritProvenance (OneFile *vf, OneFile *source);
BOOL oneInheritReference  (OneFile *vf, OneFile *source);
BOOL oneInheritDeferred   (OneFile *vf, OneFile *source);
```

Add all provenance/reference/deferred entries in source to header of vf.  Must be
called before call to oneWriteHeader.

```
BOOL oneAddProvenance (OneFile *vf, char *prog, char *version, char *command, char *dateTime);
BOOL oneAddReference  (OneFile *vf, char *filename, I64 count);
BOOL oneAddDeferred   (OneFile *vf, char *filename);
```

Append provenance/reference/deferred to header information.  Must be called before
call to oneWriteHeader.  Current data & time filled in if dateTime == NULL.

```
void oneWriteHeader (OneFile *vf);
```

Write out the header for file.  For ASCII output, if you want the header to contain
count information then you must create and fill the relevant OneCounts objects before
calling this. For binary output, the counts will be accumulated and output in a
footer upon oneClose.

```
void oneWriteLine (OneFile *vf, char lineType, I64 listLen, void *listBuf);
```

Set up a line for output just as it would be returned by oneReadLine and then call
this routine to output the line (ASCII or binary).
Use the macros above on the l.h.s. of assignments to fill fields (e.g. oneInt(vf,2) = 3).
For lists, give the length in the listLen argument, and either place the list data in your
own buffer and give it as listBuf, or put in the line's buffer and set listBuf == NULL.

```
void oneWriteComment (OneFile *vf, char *comment);
```

Adds a comment to the current line. Need to use this not fprintf() so as to keep the
index correct in binary mode.

## Closing Files (for both read & write)

```
void oneFinalizeCounts (OneFile *vf);
```

After all input has been read, or all data has been written, this routine will finish
accumulating counts/statistics for the file and merge thread stats into those for
the master file (if a parallel OneFile).

```
void oneFileClose (OneFile *vf);
```

Close vf (opened either for reading or writing). Finalizes counts if not explicitly
requested, merges theaded files, and writes footer if binary. Frees all non-user
memory associated with vf.

## Indexing & Buffer Management

```
void oneUserBuffer (OneFile *vf, char lineType, void *buffer);
```

A buffer is used to capture the list element of each line type that has one.
This routine allows you to reassign the buffer to one you've allocated, or
to revert to a default system buffer if 'buffer' = NULL.  The previous buffer
(if any) is freed.  The user must ensure that a buffer they supply is large
enough. BTW, this buffer is overwritten with each new line read of the given type.

```
BOOL oneGotoObject (OneFile *vf, I64 i);
```

Goto i'th object in the file. This only works on binary files, which have an index.

```
I64  oneGotoGroup  (OneFile *vf, I64 i);
```
Goto the first object in group i. Return the size (in objects) of the group, or 0
if an error (i out of range or vf has not group type). Only works for binary files.

# A BIT ABOUT THE FORMAT OF BINARY FILES

```
<bin file> <- <ASCII Prolog> <$-line> <binary data> <footer> <^-line> <footer-size:int64>
```

'$'-line flags file is binary and gives endian
The data block ends with a blank line consisting of '\n'

```
<ASCII Prolog> <- <'1'-line> [<'2'-line>] ( <'!'-line> | <'<'-line> | <'>'-line> )*
```

The ASCII prolog contains the type, subtype, provenance, reference, and deferred lines
in the ASCII format.  The ONE count statistic lines for each data line type are found
in the footer along with binary '\01' and '\02' lines that encode their compressors as
needed.  The footer also contains binary '&' and '*' lines that encode the object index
and group indices, respectively.

```
<Binary line> <- <Binary line code + tags> <fields> [<list data>]
```

Line codes are >= 128 for binary encoded lines.  The low two order bits of these are flags,
so each binary-encoded line type has 4 codes and a table maps these to the ASCII code.
Bit 0 indicates if the fields of the line type are compressed, and Bit 1 indicates if
the list data (if present) is compressed.

If a field is a list, then the field array element for that field is the list's length
where the low 56 bits encode length, and the high 8 bits encode the # of high-order
0-bytes in every list element if an INT_LIST (0 otherwise).
