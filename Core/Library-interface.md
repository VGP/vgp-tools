# The One Tools C library interface

The interface is defined in `Onelib.h`.  There are 23 functions and 9 macros, with one primary
type `OneFile` which maintains information about the file being read or written, including the current line.

## Synopsis

As a brief synopsis, the following reads a sequence file, prints out some simple stats, and writes a binary file containing the reverse-complemented sequences.

```
{  int totLen   = 0;
   int totCount = 0;
   
   OneFile *in = oneFileOpenRead(inFile, 0, "seq", 1); // 0 for read schema from file, 1 for single thread
   if (in == NULL) 
     die("can't open sequence file %s to read", inFile);
   
   OneFile *out = oneFileOpenWriteFrom(outFile, in, false, true, 1); // false for don't use counts from in's header, true for binary, 1 for single thread
   if (out == NULL)
     die("can't open VGP sequence file %s to write", outFile);
     
   oneAddProvenance(out,"revcomp","1.0","revcomp inFile outFile",0);
   oneWriteHeader(out);
   while (oneReadLine(in))
     if (in->lineType == 'S')
       { totLen += oneLen(in);
         reverseComplement(oneString(in), oneLen(in)); // user-provided, assume acts in place
         oneWriteLine(out, 'S', oneLen(in), oneString(in));
       }
     else if (in->lineType == 'C')
       totCount += oneInt(in,0);
	    
   printf("total sequence length %d and counts %d\n", totLen, totCount);
   printf("auto-accumulated length %d should be the same\n", in->lineInfo['S']->accum.total);
   oneFileClose(in);
   oneFileClose(out); // NB this writes out the footer as well as closing the file - don't omit!
}
```
In the above, there is no check that the schema of the file fits the expectations in the code below.  It would have been possible to carry out such a check using

```
  if (! oneFileCheckSchema (in, "D S 1 3 DNA\nD C 1 3 INT\n")) die ("schema mismatch") ;
```
which confirms that there are S lines with a single field encoding DNA, and C lines with a single field encoding an integer.  Alternatively, one could define the schema ahead of opening the file as in
```
  OneSchema *schema = oneSchemaCreateFromText ("P 3 seq\nD S 1 3 DNA\nD C 1 3 INT\n") ;
  OneFile *in = oneFileOpenRead (inFile, schema, "seq", 1) ;
  oneSchemaDestroy (schema) ;
```
Note that in this case it is necessary to define the file type "seq" in the schema, since a general schema can specify multiple file types.  The schema can also be read from file using oneSchemacreateFromFile().  

Also there is a more subtle difference, in that the first version checks that the S and C lines are present and specified as required while allowing additional unspecified line types, while the second version requires that the file only contain S and C lines.  i.e. for oneFileCheckSchema() all defined lines must be in the file, and for a schema given as an argument to oneFileOpenRead all lines the in file must be in the schema.


# INTERFACE

The following is derived from the file `ONElib.h` which provides the entire interface.  First we provide the subroutine interface, then the data types (reversing the order in a normal C header file).

## Routines for manipulating ONE files

### Creating and destroying schemas

```
OneSchema *oneSchemaCreateFromFile (char *path) ;
OneSchema *oneSchemaCreateFromText (char *text) ;
```

These functions create a schema handle that can be used to open One-code data files 
for reading and writing.  A schema file is itself a One-code file, consisting of
a set of objects, one per primary file type.  Valid lines in this file are:

```
   P <primary file type>   // a string of length 3
   S <secondary file type> // a string of length 3 - any number of these
   D <char> <field_list>   // definition of line with uncompressed fields
   C <char> <field_list>   // definition of line with compressed fields
```

`<char>` must be a lower or upper case letter.  The first upper case letter definition determines 
the objects in this file type. A maximum of one lower case letter determines the group type. 
`<field_list>` is a list of field types from: `CHAR, INT, REAL, STRING, INT_LIST, REAL_LIST, STRING_LIST, DNA`.
By convention comments on each line explain the definition.  
Example, with lists and strings preceded by their length in OneCode style:

```
   P 3 seq                            this is a sequence file
   D S 1 3 DNA                        the DNA sequence - each S line starts an object
   D Q 1 6 STRING                     the phred encoded quality score + ASCII 33
   C N 4 4 REAL 4 REAL 4 REAL 4 REAL  signal to noise ratio in A, C, G, T channels
   D g 2 3 INT 6 STRING               group designator: number of objects, name
```

The `oneSchemaCreateFromText()` alternative writes the text to a temp file and reads it with 
`oneSchemaCreateFromFile()`. This allows code to set the schema.

```
void oneSchemaDestroy (OneSchema *schema);
```
Recovers the memory allocated for the schema object.

### Reading ONE files

```
OneFile *oneFileOpenRead (const char *path, OneSchema *schema, char *type, int nthreads) ;
```
Open ONE file 'path', either binary or ascii encoded, for reading.
If the file doesn't have a header, then 'type' must be specified,
otherwise, if 'type' is non-zero it must match the header type.
All header information (if present) is read.
'schema' is also optional.  If it is NULL then the file must contain its own schema.  
If 'schema' is present then it must support 'type', and if the file contains its 
own schema, then that must be a subset of the one for this type in 'schema'.

If nthreads > 1 then nthreads OneFiles are generated as an array and the pointer
to the first, called the master, is returned.  The other nthreads-1 files are
called slaves.  The package routines are aware of when a OneFile argument is a
slave or master in a parallel group.  The master recieves provenance, counts, etc.
The slaves only read data and have the virture of sharing indices and codecs with
the master if relevant.

```
BOOL oneFileCheckSchema (OneFile *vf, char *textSchema) ; // EXPERIMENTAL
```
Checks if file schema is consistent with text schema.  Mismatches are reported to stderr.
Filetype and all linetypes in text must match.  File schema can contain additional linetypes.
e.g. if (! oneFileCheckSchema (vf, "P 3 seq\nD S 1 3 DNA\nD Q 1 6 STRING\nD P 0\n")) die () ;
This is provided to enable a program to ensure that its assumptions about data layout
are satisfied.

```
char oneReadLine (OneFile *vf);
```
Read the next ONE formatted line returning the line type of the line, or 0
if at the end of the data section.  

The content macros immediately below are
used to access the information of the line most recently read.

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

### Writing ONE files

```
OneFile *oneFileOpenWriteNew (const char *path, OneSchema *schema, char *type,
			      BOOL isBinary, int nthreads);
OneFile *oneFileOpenWriteFrom (const char *path, OneFile *vfIn,
			       BOOL isBinary, int nthreads);
```
Create a new oneFile that will be written to 'path'.  For the 'New' variant supply
the file type, subtype (if non-zero), and whether it should be binary or ASCII.
For the 'From' variant, specify binary or ASCII, schema and all other header 
information is inherited from 'vfIn'.

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
call to oneWriteHeader.  Current data & time are filled in if 'dateTime' == NULL.

```
void oneWriteHeader (OneFile *vf);
```
Write out the header for file.  For ASCII output, if you want the header to contain
count information then you must create and fill the relevant OneCounts objects before
calling this. For binary output, the counts will be accumulated and output in a
footer upon oneFileClose().

```
void oneWriteLine (OneFile *vf, char lineType, I64 listLen, void *listBuf);
```
Set up a line for output just as it would be returned by oneReadLine and then call
this routine to output the line (ASCII or binary).
Use the macros above on the l.h.s. of assignments to fill fields (e.g. oneInt(vf,2) = 3).
For lists, give the length in the listLen argument, and either place the list data in your
own buffer and give it as listBuf, or put it in the line's buffer and set listBuf == NULL.

```
void oneWriteComment (OneFile *vf, char *comment);
```
Adds a comment to the current line. Need to use this not fprintf() so as to keep the
index correct in binary mode.

### Closing files (for both read and write)

```
void oneFileClose (OneFile *vf);
```
Close vf (opened either for reading or writing).  Merges theaded
files, and writes footer if binary. Frees all non-user memory
associated with vf. 

## Goto and buffer management

```
BOOL oneGotoObject (OneFile *vf, I64 i);
```
Goto i'th object in the file. This only works on binary files, which have an index.

```
I64  oneGotoGroup  (OneFile *vf, I64 i);
```
Goto the first object in group i. Return the size (in objects) of the group, or 0
if an error (i out of range or vf has not group type). Only works for binary files.

```
void oneUserBuffer (OneFile *vf, char lineType, void *buffer);
```
A buffer is used to capture the list element of each line type that has one.
This routine allows you to reassign the buffer to one you've allocated, or
to revert to a default system buffer if 'buffer' = NULL.  The previous buffer
(if any) is freed.  The user must ensure that a buffer they supply is large
enough. BTW, this buffer is overwritten with each new line read of the given type.

# DATA TYPES

```
typedef int64_t       I64;
typedef unsigned char U8;

typedef enum { oneINT = 1, oneREAL, oneCHAR, oneSTRING, oneINT_LIST, oneREAL_LIST, oneSTRING_LIST, oneDNA } OneType;
static char* oneTypeString[] = { 0, "INT", "REAL", "CHAR", "STRING", "INT_LIST", "REAL_LIST", "STRING_LIST", "DNA" };
```
Basic data types.

```
typedef union
  { I64    i;
    double r;
    char   c;
    I64    len; // for lists: top 8 bits encode excess bytes, low 56 bits encode length
  } OneField;
```
Encoding of a data value.

```
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
```
Natural data structures for programmatic access to information in the header. Note that all of these should be used read only.

```
typedef void OneCodec; // forward declaration of opaque type for compression codecs
extern  OneCodec *DNAcodec;
```
OneCodecs are a private package for binary one file compression. DNAcodec is a special pre-existing compressor one should use for DNA. It compresses every base to 2-bits, where any non-acgt letter is effectively converted to an a.  DNA compression is case insensitive, but decompression always delivers lower-case.

```
typedef struct
  { OneCounts accum;            // counts read or written to this moment
    OneCounts given;            // counts read from header
    int       nField;           // number of fields
    OneType  *fieldType;        // type of each field
    char     *comment;          // the comment on the definition line in the schema
	// plus private fields
} OneInfo;
```
Record for a particular line type.  There is at most one list element per line type.  Again, all read only.

```
typedef struct OneSchema {} OneSchema ;
```
The schema type, all private to the package.  Internally a schema is stored as a linked list of OneSchema objects, with the first holding the (hard-coded) schema for the header and footer, and the remainder each holding the schema definition data for one primary file type.

And, finally, the main OneFile type - this is the primary handle used by the end user.
```
typedef struct
  {
    // these fielde may be set by the user

    BOOL           isCheckString;      // set if want to validate strings char by char - slows down reading
    I64            codecTrainingSize;  // number of bytes to see before building codec - default 100k - can set before writing

    // these fields may be read by user - but don't change them!

    char           fileType[4];
    char           subType[4];
    char           objectType;         // line designation character for primary objects
    char           groupType;          // line designation character for groups (optional)

    I64            line;               // current line number
    char           lineType;           // current lineType
    I64            object;             // current object - incremented when object line read
    I64            group;              // current group - incremented when group line read
  
    OneProvenance *provenance;         // if non-zero then count['!'] entries
    OneReference  *reference;          // if non-zero then count['<'] entries
    OneReference  *deferred;           // if non-zero then count['>'] entries
 
    OneField      *field;              // used to hold the current line - accessed by macros
    OneInfo       *info[128];          // all the per-linetype information
 
    // the remainder is private to the package
    
  } OneFile;                      //   the footer will be in the concatenated result.
```

# A BIT ABOUT THE FORMAT OF BINARY FILES

```
<bin file> <- <ASCII Prolog> <$-line> <binary data> <footer> <^-line> <footer-size:int64>
```
'$'-line flags file is binary and gives endian. The data block ends with a blank line consisting of '\n' only.

```
<ASCII Prolog> <- <'1'-line> [<'2'-line>] ( <'!'-line> | <'<'-line> | <'>'-line> )*
```
The ASCII prolog contains the type, subtype, provenance, reference, and deferred lines and schema
in the ASCII format.  The ONE count statistic lines for each data line type are found
in the footer along with binary ';' and ':' lines that encode their compressors as
needed.  The footer also contains binary '&' and '*' lines that encode the object index
and group indices, respectively.

```
<Binary line> <- <Binary line code + tags> <fields> [<list data>]
```
Binary line codes are >= 128.  The low two order bits of these are flags,
so each binary-encoded line type has 4 codes and a table maps these to the ASCII code.
Bit 0 indicates if the fields of the line type are compressed, and Bit 1 indicates if
the list data (if present) is compressed.

If a field is a list, then the field array element for that field is the list's length
where the low 56 bits encode length, and the high 8 bits encode the # of high-order
0-bytes in every list element if an INT_LIST (0 otherwise).
