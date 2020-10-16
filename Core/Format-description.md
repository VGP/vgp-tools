# VGP-tools format description and schema definition

### Authors:  Gene Myers, Richard Durbin, and the Vertebrate Genome Project Assembly Group
### Last Update: April 13, 2020

## 0. Introduction

This document describes the generic VGP-Tools file structure.
For other VGP tools documentation see the top level [README.md](https://github.com/VGP/vgp-tools/blob/master/README.md).

VGP files are container for objects of a single *type*.  There are *primary* types, each of
which has its own schema, and *secondary* types which are specialisations of a primary type for
particular applications, requiring particular properties or fields to be present.
Every primary and secondary file type is assigned a 3-letter lower case string, that can be used as a filename suffix, e.g. dataset1.seq
or dataset2.pbr.  For example, from the [standard  VGP sequence sequence schema](https://github.com/VGP/vgp-tools/blob/master/VGP-sequence-schema.md),

- **.seq**   primary file type for sequence objects
	- **.irp** for Illumina read pairs
	- **.pbr** for PacBio long reads and relevant meta-data
	- **.10x** for 10X Genomics read clouds with their extracted barcodes
	- **.ctg** for contigs from an assembly
	- **.kmr** for sets of kmers

We describe here the ASCII version of the format, which gives the specification.  There is also a binary version encoding the same information that is more compact and efficient.  

Every VGP file begins with a header segment followed by a data segment.  Headers, relevant for all
file types, are described in Chapter 1 following this introduction, and the encoding of the data is documented in Chapter 2.
The utility program **VGPstat** validates the format of any given file against the schema, and also
can optionally reconstruct a file's header segment given only the data part.

Each line in either the header or data section starts with a single character ("1-code") which determines the encoding of the remainder of the line. Header lines are specified by *non-alphabetic* characters, and data lines by *alphabetic* letters.  Consider as an example the following Illumina read pair (.irp) file:

```
     1 3 seq 1 0            header start: primary file type, major and minor version
     2 3 irp                optional secondary subtype
     ! 7 VGPpair 3 0.1 ...  provenance line saying how file arose
     # P 3                  number of read pairs in file
     # S 6                  number of sequences in file
     @ S 5                  maximum number of bp’s in a read
     + S 26                 total number of bp's in file
     P                      data start: separator for read pairs - helps human interpretability 
     S 5 acgta              forward sequence of pair 1
     S 3 gtt                reverse sequence of pair 1
     P                      pair 2
     S 4 gcta               sequence 3
     S 5 ggtac              sequence 4
     P                      pair 3
     S 4 atta               sequence 5
     S 5 cctac              sequence 6
```

The design of VGP formats is based on the following principles:

1.	The format should be trivial to parse as input.  All the burden of encoding is placed
on the software that produces the formatted file.

2.	The length of a list should always precede the list (e.g. 3 xxx).  Bracket (.eg [xxx]) or
terminator (e.g. xxx0) constructions are not permitted.

3.	The total number of items in the file, and the maximum size or total size of variable-length items,
should be specified in the header of a file so that memory allocation can be performed
once at the start by the reader if desired.

4.	There should be a VGP tool that translates the ASCII format
to a terse, indexable binary form and another that inverts this conversion, as well as a code-level interface to the binary encoding.

5.	The ASCII form should not be overly verbose. It is not primarily for human consumption, but must be easily interpretable by a human being when and if they need to look at the data.

6.	Complex identifiers or symbolic names for objects are not used.  The nth occurrence of an
object of a given type is referred to by the number n in any future context within the file.

7.	A reference to an object is not permitted before its definition.

8.	The provenance of a file must be given as the sequence of tools applied to obtain it,
where each tool application is described by its name, version number, command line applied,
and date & time performed. 

Syntactically, the goal is that an ASCII VGP file should be readable by a simple parser that can
(a) read the next character, (b) read the next integer, (c) read the next real number, (d) read the next n symbols, and
(e) skip to the next line, where tokens are separated by a single white space character.
Moreover, any program reading in VGP data should never require dynamic memory
allocation to store objects – the size of any object or object collection should be given
before the first object of that type is described.  We realized this by designing VGP formats
so that every file begins with a header section that gives the number of objects, maximum size
of an object, and total size of all objects, along with the provenance information of Principle 8.

The "1-code" design supports Principle 5.  The sequence of items on the remainder of the line is determined by the code.  Tokens are separated by a single space character.  Variable length lists,
including strings, are preceded by their length in keeping with Principle 2.  Thus, the 1-code and
subsequent tokens determine when the encoding of information on a line is at an end. 
Any additional information to the end of the line, for example the comments in
the example above, is optional, providing an extensibility
mechanism where one can add auxiliary information if desired.

The first header line must always be a 1-line confirming the primary file type and specifying
the major and minor version numbers separated by whitespace.
This can optionally be followed immediately by a file subtype line of type ```2``` that gives the
secondary file type.
Additional header lines give information about the number of items
in the file (#-lines), the maximum length of lists (@-lines), and the total number of items in
a given list class (+-lines).  Finally, provenance lines (!-lines) inform one about how the
particular file came to be.  The data segment in this simple example then consists of pairs
of reads encoded by P-lines indicating the beginning of a pair, followed by two S-lines
giving the two sequence reads in the pair.

Conceptually VGP files are *immutable*, meaning that we do not expect their contents to change.
This means that subsequently in the same file, or more often in future files in a pipeline, we
can refer to objects by their ordinal position in a file, 1...n, not requiring named identifiers
(Principle 6).  

Together, these design decisions simplify many things.


## 1. File Headers

Considerable effort is invested on headers in VGP-formats in keeping with Principles 3 and 8.
We will specify VGP format syntax with "casual" context-free grammar rules.
Tokens are between angle-braces ```<>``` and may optionally have a descriptive tag following
a colon after the token name, e.g. ```<string:file_name>```.  The token name can be a predefined
type, such as ```int```, ```real```, or ```char```, in which case the syntax of the token is
the usual lexical syntax for such literals.  The right hand side of a rule can also be a
regular expression, where ```*``` denotes 0-or-more, ```+``` 1-or-more, ```^n``` exactly n
repetitions, ```|``` separates two alternatives, an item between square brackets ```[]``` is
optional, and parentheses may be used to disambiguate precedence.  For example, a string
is defined by the rule:
```
    <string> = <int:n> <char>^n
```

The first line of any header must be the primary file type declaration that has
the synatx:

```
    <version_header> = 1 <string:file_type> <int:major> <int:minor>
```
where the initial ```1``` indicates that this is a "1-code" file (as well as this being line 1
&#x1F609;)
and ```<file_type>``` is one of the five 3-letter file suffixes given at the start.

The initial header line can be followed by an optional subtype line

```
    <subtype_header> = 2 <string:file_subtype>
```
where now file type is one of the suffixes above for a secondary file type.

Next there are three header line types - #, +, and @ - that allow one to specify the number, total size, and maximum size of objects across a file.  These all have the syntax:

```
    <size_header> = [#+@] <symbol:S> <int>
```
\#-lines tell you the number of lines in the file of type S.  For line types that encode a
list of items, such as ```<string>``` (a list of characters) or say a list of restriction map
sites, a +-line tells you the total number of items in all the lists in the file,
e.g. "<code>+ S 26</code>" in the example above indicates that altogether the sequences in the
file total 26 bases.  Similarly, an @-line indicates the length of the largest list that
occurs in any given line of the specified type.
These "limit" lines are only present for line types that occur in the body and not lines
that occur in the header itself.  In almost all cases there is only one list object
per line type, however, a notable exception is the restriction site line (see 2.2), which has
a variable number of strings in a line.  In such cases, ```+``` is the sum of all of the list
element lengths over all lines of the given type, and <code>@</code> is the maximum of the sum of the
list lengths of each line.

Often the objects in a file are naturally partitioned into groups, e.g. all the read pairs in a
flow-cell lane, or all the read pairs in a "cloud".   VGP formats support this with the concept
of "group" data lines that by convention have lower case symbols as the line type indicator.
For a group, one would like to know the maximum number of items of a given type in any group
and the maximum size needed to contain all the list objects occurring in any group.  So in the
header, a ```%``` designator indicates the maximum number of objects or maximum total size of
list objects within any given group type.  The syntax for these lines is:  

```
    <group_header> = % <symbol:g> [#+] <symbol:S> <int>
```
where ```g``` is the group line designator and ```S``` is the line type in question.
(The ```@``` reduction is not necessary as it is the same as the global maximum.)

Another important header line type indicates that references in this file are made to objects
in another file. This has the syntax:

```
    <reference_header> = < <string:file_name> <int:nx>
```
All the objects (all of the same type) in the specified file are available and
and ```nx``` indicates the number of these items in the file (hence
the range of reference indices is ```[1,nx]```).  For example a ```.sxs``` alignment file refers
to sequence objects designated by S-lines in another file.

A related concept is to refer to another file upon which the objects in the current file depend.
We denote these with a '>'-line that has the opposite direction to the '<' of the reference line
above.

```
    <forward_header> = > <string:file_name>
```
In this case there is no need to indicate the number of objects in the file, since the current
file will not refer to them.

The final header line type is the provenance or !-line that records a processing step that
was involved in producing the current file.  Each line contains four strings giving (a) the
program name, (b) the version of that program as a string, (c) the command line that was executed,
and (d) the date and time it was run.

```
    <provenance_step> = ! <string:name> <string:version> <string:command> <string:date>
```

In summary, every VGP formatted file begins with a header.  Every header starts with a version
line optionally followed by a subtype line.  Then ensue a number of size lines for every relevant
data line of the primary file type.  And finally, at the end, any relevant reference-, forward-
and provenance lines.  In a rule:

```
    <header> = <version_header> [<subtype_header>] (<size_header>|<group_header>)+
                   (<reference_header>|<forward_header>|<provenance_step>)+
```

## 2. File Data and Schema Definition

Every data line has a type designated by a letter. The types of data
lines permitted depend on the schema associated to the primary file type.
At most one data line type can be used to start a new *group* of objects.





