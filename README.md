# VGP Tools: The VGP tool kit and formats

## Authors:  Gene Myers, Richard Durbin and the VGP assembly group
## _Date:   February 15, 2019_

VGP-Tools is a collection of tools that operate on DNA sequencing data encoded in a collection
of file formats called collectively VGP-Formats.  The encodings include descriptions of source
data, process intermediates, and the ultimate reconstructed genome assemblies.  There is a file
type, indicated in the file name by a 3-letter extension, for each desired data object collection
within which a specific, very simple ASCII format is used to encode the data.  For example, the
VGP currently collects 4 core data sets for a given genome project which can be encoded in the
following file types:

**.pbr**	A collection of PacBio long reads and relevant meta-data about each are
encoded in such a file (e.g. myfish.pbr)

**.brm**	An encoding of the Bionano restriction maps for a collection of molecules.

**.irp**	An encoding of Illumina read pairs.  We generate such data using a Hi-C protocol, but this applies more generally to any read pair protocol.

**.10x**	An encoding of 10X Genomics read clouds. 

The VGP tool library contains programs that takes the files produced by each instrument supplier
and converts it into the relevant VGP format and file type above.  Then, more importantly, we support
the following downstream file types that allow one to encode assemblies, scaffolds, and important
intermediates therein:

**.sxs** a collection of sequence to sequence matches.

**.sxr** a collection of sequence to restriction map matches.

**.sml** subsets of matches that define a coherent object, e.g. an assembly contig, a haplotype branch, etc.

**.seq** a simple collection of sequences, used to represent contigs (effectively a basic type)

**.scf** a collection of scaffolding suggestsions representing putative joins and breaks between contigs

**.sfl** subsets of scaffold suggestions, which can be used for example to generate a new set of contigs,
         or a set of putative chromosomes (two different files for the different haplotypes).

The VGP format encoding of each file type is documented in a separate chapter of this document following the overview
given here.  The formal definition of all currently valid specs is provided in the code for the utility **vgpvalidate**
which also provides utilities to reconstruct headers.  The design of VGP formats is based on the following principles:

1.	The format should be trivial to parse as input.  All the burden of encoding is placed on the software that produces the formatted file.

2.	The length of a list should always precede the list (e.g. 3 xxx).  Bracket (.eg [xxx]) or terminator (e.g. xxx0) constructions are not permitted.

3.	The maximum size or total size of items to follow should be specified at the start of a file so that any desired memory allocation can be performed once at the start by the reader.

4.	For every format file type, there is a VGP tool that translates the ASCII format to a terse binary form and another that inverts this conversion.

5.	The ASCII form should not be overly verbose, it is not for human consumption, but simply must be easily interpretable by a human being when and if they need to look at the data.

6.	Complex identifiers or symbolic names for objects are not used.  The nth occurrence of an object of a given type is referred to by the number n in any future context within the file.

7.	A reference to an object is not permitted before its definition.

8.	The provenance of a file must be given as the sequence of tools applied to obtain it, where each tool application is described by its name, version number, command line applied, and date & time performed. 

From another point of view, a VGP formatted file should be readable by a simple parser that can
(a) read the next integer, (b) read the next real number, (c) read the next n symbols, and
(d) skip to the next line, where each read can begin with white space that is skipped.  The
parser should further never require dynamic memory allocation to store objects – the size of
any object or object collection should be given before the first object of that type is described.

Following these principles, VGP formats uses a very simple “1-code” schema in which the data is
encoded in an ASCII file in a sequence of lines where the first character of each line determines
what kind of information the line encodes.  Consider as a working example,
the following Illumina Read Pair (.irp) file:

```
	1 irp 1.0          file type, major and minor version
	# ! 1                  number of provenance lines in the file
	# P 2                  number of read pairs in file
	# S 4                  number of sequences in file
	@ S 5                  maximum number of bp’s in a read
	+ S 17                 total number of bp's in file
	! 7 VGPpair 3 0.1 ...  provenance line saying how file arose
	P                      separator for read pairs - helps human interpretability 
	S 5 acgta              forward sequence of pair 1
	S 3 gtt                reverse sequence of pair 1
	P
	S 4 gcta
	S 5 ggtac 
```
Non-alphabetic "1-code" symbols indicate header lines, which must precede data lines that always begin with an
alphabetic "1-code".  Tokens on a line are separated by a single whitespace character (space or tab).  Variable length lists, including strings, are preceded by their length in keeping with principle 2. The
1-code and subsequent tokens determine when the encoding of expected information on a line is at an end. 
Therefore optional additional information to the end of the line provides an extensibility
mechanism where one can add auxiliary information if desired.  In the example above, the whitespace and 
comments following the encoded information on each line are ignored by a VGP parser.

### Header lines

Considerable effort is invested on headers in VGP-formats in keeping with principles 3 and 8.
The first header line must always be a version line confirming the file type and specifying the major and minor
version numbers separated by a ```.```.  Additional header lines give information about the number of items
in the file (#-lines), the maximum length of lists (@-lines), and the total number of items in a given list class
(+-lines).  In addition, provenance lines (!-lines) inform one about how the particular file came to be.

We now introduce how we formally describe VGP formats, defining the header lines to illustrate. Using a "casual" context free grammar rule the inital file type declaration has the syntax:

```
    <version_header> = 1 <file_suffix> <major>.<minor>
```
where the initial ```1``` indicates that this is a a "1-code" file (as well as this being line 1) and ```<file_suffix>``` is one of the ten 3-letter file suffixes above.

There are three header line types -- #, +, and @ -- that allow one to specify the number, total size, and maximum size of objects across a file.  These all have the syntax:

```
    <size_header> = [#+@] <symbol:S> <int>
```
```#```-lines tell you the number of lines in the file of type S.  For line types that encode a list of items, such as
string (a list of characters) or say a list of digest sites, a ```+```-line tells you the total number of items
in all the lists in the file, e.g. "<code>+ S 17</code>" in the example above indicates that altogether the
sequences in the file total 17 bases.  Similarly, an ```@```-line indicates the length of the largest list that
occurs in any given line of the specified type.  

Often the objects in a file are naturally partitioned into groups, e.g. all the read pairs in a flow-cell lane, and the VGP formats support this with lower case symbols as data line indicators.  For these the ```%``` designator indicates the maximum number of objects or total size of list objects within any given group type.  The syntax for these lines is:  

```
    <group_header> = % <symbol:g> [#+] <symbol:S> <int>
```
where ```<g>``` is the group line designator (always a lower case letter by convention) and ```<S>``` is the line type in question.

Another important header line type indicates that a given file should be included and has the following syntax

```
    <include_header> = '<' <string:file_name> <symbol:S> <int:nS>
```
The effect is as if the data portion of the named file had been inserted at this point and any size information in its header
had been appropriately factored into the header of the current file.  Note carefully, that the indexing/numbering of
lines/objects is according to their order after inclusion as per principle 6. In this example the S-objects from the included file are the ones that are to be referred to subsequently, and the number of them is given here, so that actual inclusion is not required.

In some cases there is an another file type that is closely associated with the current file, and which depends on objects defined in the current file.  It is possible to indicate such files with a forward reference which looks similar to the include line, but with a '>' symbol rather than a '<' symbol.

```
	<forward_header> = '>' <string:file_name>
```
The semantics of this can be taken to be that the file referred to will be included *after* the end of the current file, rather than at the current location, i.e. deferred inclusion.  In this case there is no need to indicate the primary objects in the file and their number, since the current file will not refer to them.

The final header line type is provenance or '!'-lines that record a processing step that was involved
in producing the current file.  Each line contains 3 strings and 2 integers giving (a) the program name, (b) the major
and minor version numbers separated by a '.', (c) the command line that was executed, and (d) the date and time it
was run.  Each \<string\>, as per principle 2, is an integer followed, possibly after white space, by that number
of characters given by the integer.
```
    <step>   = <string:name> <string:version> <string:command> <string:date>
    <string> = <int:n> <char> ^n
```
Since there can be multiple provenance lines in a header, one expects to see "<code># ! <int></code>" before any
of these lines in order to know how many there are.

Below we summarise the proposed file types, dropping the generic header
line types already introduced above.  There are a few additional header
line types that we will introduce as they become necessary.

## Data lines defined according to file type

The types of data lines permitted depend on the file type.  We therefore specify here each file type in turn.  We remind the reader that while we aim to keep this README up to date, the formal versioned definition of each VGP file type is determined by the validator **vgpvalidate** and documented in its code.

### Illumina Read Pair files, .irp

The data lines are as follows:

```
   g <int:r> <string:group_name>    read group with r read pairs
   P              a pair follows immediately
   S <string>     forward or reverse sequence
   Q <string>     QV scores in single character phred encoding (ASCII 33+q) 
```
The data portion of a .irp file consists of a sequence of line bundles defining a read pair with group lines
interspersed.  Each read pair is defined by a P-line followed by two S-lines given the forward and reverse
sequences, and optionally two Q-lines giving the QV's for the sequences.  By convention the first S/Q-line is
for the forward read and the second for the reverse read.  Between these bundles one may find interspersed
g-lines that specify a user defined group of reads where the group consists off the next r-pairs of reads, where
r is given by the first integer in the line.  This is followed by a string giving a user-selected
name for the group.
In terms of a regular expression:
```
    (<g-line> | <P-line><S-line>[<Q-line>]<S-line>[<Q-line>]) *
```

The S strings are over the alphabet A, C, G, T.  Each symbol of a Q string is an
ASCII character as in FASTQ.  So Q-strings are over the 94 printable ASCII
characters [!-~] and correspond to the range 0 to 93.  IUPAC ambiguity codes, except for N,
are not supported, any unsupported code is turned into an N by VGP tools.  DNA strings may be
upper or lower case.

If group lines are present they do not need to form a partitioning of the data (i.e. every pair belongs
to one and only one group), although this property is often true of the programs that produce groups
in .irp's.

**VGPpair** is a VGP tool that takes two .fastq files containing corresponding forward and reverse
reads (the usual bioinformatic convention for such) and produces a .irp from the inputs.  Reads
within the same lane are placed in a given group. 

## 10X Read Cloud files, .10x

The data lines are as follows:

```
   c <int:n> <string:bar>    the next n read (pairs) are in the cloud with the given barcode

   P                         a pair follows immediately
   S <string>                forward and reverse sequence
   Q <string>                QV scores
```
The P, S and Q lines are exactly as described for an .irp file previously.  The c-lines group read
pair lines into clouds which are all the read pairs
in which the forward read contained the indicated barcode sequence for the cloud (which is removed
for the .10x file).  In this case the c-lines are guaranteed to partition the read pairs (i.e. every
read pair is guaranteed to be in one and only one cloud).

In terms of a regular expression:
```
      (<C-line> (<S-line>[<Q-line>]<R-line>[<Q-line>]) ^n ) *
```

**VGPcloud** is a VGP tool that takes an .irp-file containing the bar-coded Illumina read pairs
and re-organizes the data into clouds and removes the bar codes.  It understands how barcodes
are encoded in the forward read of each pair, and takes into account potential errors in the
barcode sequence.

## PacBio Long Reads, .pbr

The data lines are as follows:

```
g <int:r> <string:SMRT_header>				SMRT cell with r reads

L <int:well> <int:1st.pulse> <int:last.pulse> <real:score>	Read well and pulse range
S <string>				                	sequence
N <real:A> <real:C> <real:G> <real:T>   	SNR in each base channel for read
A <string>				                	capped pulse widths
```
The data is grouped into SMRT cells by g-lines that begin each group and declare how many reads
were produced in that cell.  Each read is described by an L- and S-line and optionally an
N- and A-line.  The L-line gives the well and pulse range from which the sequence in the S-line
was extracted.  If present, the N-line give the average SNR in the channels for each base and
the A-line gives the capped pulse width for each base as the character 1, 2, 3 or 4.  Pulse
widths larger than 4 are clipped to 4 as Arrow does not treat widths larger than 4 differently.
Basically the very short pulses are indications of potential error but any pulse over 4 units
long is almost certainly a good call.  The A and S strings for a given read have the same
length, therefore the S header information applies also to the A strings.

In terms of a regular expression:
```
   (<g-line> (<L-line><F-line>[<N-line><A-line>])r)c
```

**VGPextract** is a VGP tool that takes one or more Pacbio subreads.bam or .sam files as input
and extracts the information to make a .pbr file.

# Bionano Restriction Maps, .brm

This file encodes the primary Bionano data that is required for downstream map building or alignment to sequence.
It only allows a single restriction enzyme per data set, consistent with current usage of the BioNano instrument, although the BioNano .bnx format in principle supports multiple enzymes.  

```
= <string>                     number of channels and their recognition sites

M <read:len> <int:n>           molecule, length and number of sites
D <int:n> <real1> ... <realn>         site locations in digest
N <int:n> <real1> ... <realn>         signal to oise R of each site in digest
A <int:n> <real1> ... <realn>         Average intensity of sites in digest
```
The =-line declares how many distinct probes were used (typically only 1) and their recognition
sequence.  The other header lines tell one how many molecules were mapped and for each
probe/channel, the total number of sites in the file, and the maximum number of sites that
occurred in any molecule.  The data portion of the file consists of a sequence of molecule
descriptions, where the first M-line tells one how many sites there are in that molecule for
each probe.  The ensuing lines then give for each channel, the location of the sites and
length of the molecule (D-line) to within the accuracy of the machine, the SNR at each site
(N-line), and the average intensity of the spot at each site (A-line).  So for example, if
a file encodes maps for 2 channels, then each molecule is described by 7 lines.  The lists
of site information are presumed to be in order across each molecule, so for example, the
positions of the D-line should always be in increasing order ending with the length of the
molecule.  Recognition sites are assumed to be simple sequences.  The D,N, and A lines can
occur in any order within a molecule description.

In terms of a regular expression:

```
(<M-line> (<D-line><N-line><A-line>)c )m
```

**VGPdigest** is a VGP tool that takes a Bionano .bnx-file as input and extracts the information to make a .brm-file.

## Sequence match files, .sxs

```
< <string:filename_a> <S> <int:na>
< <string:filename_b> <S> <int:nb>

A <int:a> <int:b>			           	indexes of aligned sequences
I <int:as> <int:ae> <int:bs> <int:be>	start and end in a and in b
Q <int>									mapping quality in phred units
C <string>								cigar string
T <int:t> <int>t						list of trace points
M <int>									number of matching aligned bases
D <int>									number of differences = substitutions + indel bases
```
There must be two included files, each containing some number of S lines designating sequences.  It is these sequences that the ```a``` and ```b``` index fields on the ```A``` alignment lines refer to.

## Sequence Match List file, .sml


## Sequence file, .seq


## Scaffolding information file, .scf

```
< <string:contig_file_name> <int:nseqs>
B <int:seq> <int:start> <int:end>
Q <int:phred_confidence>	additional free information after phred score
J <int:seq_a> <int:seq_b> <[s|e]:end_a> <[s|e]:end_b>  
G <int:gap>
```

## Scaffolding List file, .slf

## Alignment/assembly files, .ala

We view the assembly task as choice of a subset of read alignments that
result in a consistent layout.  These can then be used to generate
consensus contig sequences, which are the subject of the next file
type, which represents contigs and scaffolding.

We also support subsets of alignments which are believed to represent
consistent groups of connections between .

The syntax and meaning of each line is as follows, where <string> denotes a length followed
after white space by a sequence of characters of that length.
```
   <  <string> <int>   filename of source reads for alignment, and number of reads in it
   + L <int>      total number of alignments (links) in the file
   + A <int>     total number of assembly components in the file
   + H <int>     total number of homology groups
   L <int:a> <int:b> <int:sa> <int:ea> <int:sb> <int:eb> link between  read a interval [sa,ea) and read b interval [sb,eb)
   D <int>          number of differences
   T <int:t> <int>t    trace points
   C <string>     cigar string
   A <int:n> <int>n  list of links in this assembly component
```
The F, R and Q lines are exactly as described for an .irp file previously.  Provenance lines may
occur in the header.  The C-lines group read pair lines into clouds which are all the read pairs
in which the forward read contained the indicated barcode sequence for the cloud (which is removed
for the .10x file)

The + and @ lines are header lines and must proceed the data lines.  The data portion consists of
a C-line declaring how many read pairs are in this cloud, followed by that many 
F- and R-lines optionally followed by a corresponding pair of Q-lines.  In terms of a regular
expression:
```
(<+-line>|<@-line>)+ (<C-line> (<F-line><R-line>[<Q-line>2])n )c
```
Provenance lines beginning with ! as defined in the prolog may also occur in the header section
of the file.

# VGP Tool Manuals

```
1. VGPpair [-vsg] <forward:fastq> <reverse:fastq>
```

VGPpair reads two correlated fastq files and outputs the pairs in .irp format to the
standard output.

The input file names must have a .fastq suffix, e.g. somewhere/foo.fastq, however
the suffix need not be given on the command line as VGPpair will automatically add the
appropriate suffix if it is not present.  That is, "<code>VGPpair foo1 foo2</code>" will
operate on <code>foo1.fastq</code> and <code>foo2.fastq</code>.

The -v option asks VGPpair to output information on its progress to the standard error.
The -s option asks VGPpair to *not* output the quality values or Q-lines, but just the
forward and reverse sequences in F- and R-lines.
The -g option asks VGPpair to group the data into lanes.  In this case the files must
have been produced by standard Illumina software from their more basic .bcl files, and
therefore the .fastq headers encode the instrument, flow cell, lane, etc. in fields
between :'s where the data is in order of flow cell and lane.  VGPpair uses this
information to group reads into lanes.

VGPpair checks the syntax of the .fastq files but does not verify that the DNA and QV
strings are over the appropriate symbols.

```
2. VGPpacbio [-vaU] <data:subreads.[bam|sam]> ...
```

VGPpacbio reads a sequence of Pacbio .bam or .sam files and outputs a .brp file to
standard output.

The input file names must have a .bam or .sam suffix, e.g. somewhere/foo.bam, however
the suffix need not be given on the command line as VGPpacbio will automatically add the
appropriate suffix if it is not present.  That is, "<code>VGPpacbio foo1 foo2</code>" will
operate on <code>foo1.bam</code> and <code>foo2.sam</code> if those are the files present.

The -v option asks VGPpacbio to output information on its progress to the standard error.
The -a option asks VGPpacbio to  output the arrow information in N- and A-lines per bundle,
the default is to not output this information.  The -U option requests that the DNA seuqences
are in upper case (the default is lower case).

Each file is forms a separate group in the output file.
