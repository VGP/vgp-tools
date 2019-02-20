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

**.irp**	An encoding of Illumina read pairs (linked by a Hi-C protocol, but applies
more generally to any read pair protocol).

**.10x**	An encoding of lllumina read 10X clouds. 

The VGP tool library contains a program that takes the files produced by each instrument supplier
and converts it into the relevant VGP format and file type above.  Then we support downstream assembly-related file types:

**.sxs** a collection of sequence to sequence matches.

**.sxl** containing lists of matches that be used to represent, for example, the matches that should contribute to an assembly

**.seq** a simple collection of sequences, used to represent contigs

**.scf** a collection of scaffolding suggestsions representing putative joins and breaks between contigs

**.scl** contains lists of scaffold suggestions, which can be used for example to generate a new set of contigs, or a set of putative chromosomes (two different files for the different haplotypes).

The VGP format encoding of
each file type is documented in a separate chapter of this document following the overview
given here.  The formal definition of all currently valid specs is provided in the code for the utility **VGPformatCheck** which also provides utilities to reconstruct headers.  The design of VGP formats is based on the following principles:

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
what kind of object the line encodes.  Consider as a working example,
the following simple Illumina Read Pair (.irp) file:

```
	. VGPirp 1 0 file type, major and minor version
	# P 2        number of read pairs in file
	# S 4        number of sequences in file
	@ S 5        maximum number of bp’s in a read
	+ S 17       total number of bp's in file
	P            separator for read pairs - helps 
	S 5 acgta    forward sequence of pair 1
	S 3 gtt      reverse sequence of pair 1
	P
	S 4 gcta
	S 5 ggtac 
```

Non-alphabetic "1-code" symbols indicate header lines, which must precede data lines denoted by alphabetic symbols. The first line must be a '.' header line confirming the file type. Subsequent header symbols used here are # for the number of lines of a particular type, @ for the maximum
number of list items (bases in this case) on a line of the given type, and + for the total
number of list items on lines of this type. In this file P data lines indicate read pairs, and S lines sequence
objects.  Tokens on a line are separated by a single whitespace character (space or tab).  The
1-code and subsequent tokens determine when the encoding of information on a line is at an end,
therefore optional additional information to the end of the line provides an extensibility
mechanism where one can add auxiliary information if desired.  In the example above, the
comments following the information on each line are ignored by a VGP parser.

In line with principle 8, every VGP formatted file should have a provenance line in its header
that documents how the file came to be.  We illustrate how we formally describe VGP formats
using this multi-line type as an example.  The first line begins with a ! and is followed by
the number of ensuing lines that will document a step in the processing that led to the current
file.  We use basically a “casual” context free rewrite rule where (a) line breaks are indicated
by actually using lines in the rule,  (b) tokens in angle brackets can define their value,
e.g. \<int:n\> denotes that n holds the value of the integer, and (c) token names can be indexed
to indicate count, e.g. step.3:  ** RICHARD COMMENT: I WANT TO CHANGE THIS **

```
    <!-line> = ! <int:n>
                   <step.1>
                   . . .
                   <step.n>
```
Each step-line contains 4 strings giving the program name, version #, command line that was
executed and date and time it was run.  Each \<string\>, as per principle 2, is an integer followed,
possibly after white space, by that number of characters given by the integer.  The style of our
rewrite rules for this are as follows:

```
    <step>   = <string:name> <string:version> <string:command> <string:date>
    <string> = <int:n> <char>n
```
Provenance lines may be in the header section of any VGP-Format file and for brevity will not
be explicitly described in each of the format descriptions.

Below we summarise proposed file types, dropping the generic header
line types already introduced.  There are a few additional header
line types that we will introduce as they become necessary.

## Illumina Read Pair files, .irp

```
   G <string:lane_header> <string:file> <int:r>  read group with r read pairs
   P              a pair follows immediately
   S <string>     forward and reverse sequence
   Q <string>     QV scores in single character phred encoding (ASCII 33+q) 
```
The G-line contains a header string for the group which is the standard Illumina read identifier through to the
lane number, and there is a group for each lane.  The G-line also gives the file name or names the group
was read in from (separated by a |) and the number of read pairs that follow.
The P string simply indicates the next two S-lines are a forward/reverse pair is present primarily
for human interpretability.
The S strings are over the alphabet A, C, G, T, and N.  Each symbol of a Q string is the
ASCII character as in FASTQ.  So Q-strings are over the 94 printable ASCII
characters [!-~] and correspond to the range 0 to 93.  IUPAC ambiguity codes, except for N,
are not supported, any unsupported code is turned into an N by VGP tools.  DNA strings may be
upper or lower case.

In terms of a regular
expression:

```
    (<G-line> (<P-line><S-line>[<Q-line>]<S-line>[<Q-line>])r)*
```

**VGPpair** is a VGP tool that takes two .fastq files containing corresponding forward and reverse
reads (the usual bioinformatic convention for such) and produces a .irp from the inputs.  Reads
with the same header prefix are placed in a given group. 

## 10X Read Cloud files, .10x

The syntax and meaning of each line is as follows, where <string> denotes a length followed
after white space by a sequence of characters of that length.
```
   + S <int>		total bases of sequence data in file
   + C <int:c>		total number of clouds in file
   + P <int>		total number of read (pairs) in the file
   @ P <int>		maximum number of read (pairs) in any cloud
   @ S <int>		maximum number of read bases in any cloud
   @ F <int>		maximum forward read length
   @ R <int> 		maximum reverse read length

   C <int:n> <string:bar>    the next n read (pairs) are in the cloud with the given barcode
   P                         boundary marker for read pair
   S <string>                forward sequence
   Q <string>                QV scores
   S <string>	             reverse sequence
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
Provenance lines beginning with ! as defined in the prolog are also expected to occur in the
header section of the file.

**VGPcloud** is a VGP tool that takes an .irp-file containing the bar-coded Illumina read pairs
and re-organizes the data into clouds and removes the bar codes.  It understands how barcodes
are encoded in the forward read of each pair, and takes into account potential errors in the
barcode sequence.

## PacBio Long Reads, .pbr

The syntax and meaning of each line is as follows, where <string> denotes a length followed
after white space by a sequence of characters of that length.  One should be particularly
careful to observe the subscripted meta-values.

```
% C <int>					max. number of bases in a SMRT cell
@ R <int>					max. number of reads in a SMRT cell
@ S <int>					max. number of bases in a read

G <string:SMRT_header> <string:file_name> <int:r>	        SMRT cell with r reads
L <int:well> <int:1st.pulse> <int:last.pulse> <real:score>	Read well and pulse range
F <string>				                	sequence
N <real:A> <real:C> <real:G> <real:T>                   	SNR in each base channel for read
A <string>				                	capped pulse widths
```
The data is grouped into SMRT cells by G-lines that begin each group and declare how many reads
were produced in that cell.  Each read is described by an L- and S-line and optionally an
N- and A-line.  The L-line gives the well and pulse range from which the sequence in the F-line
was extracted.  If present, the N-line give the average SNR in the channels for each base and
the A-line gives the capped pulse width for each base as the character 1, 2, 3 or 4.  Pulse
widths larger than 4 are clipped to 4 as Arrow does not treat widths larger than 4 differently.
Basically the very short pulses are indications of potential error but any pulse over 4 units
long is almost certainly a good call.  The A and S strings for a given read have the same
length, therefore the S header information applies also to the A strings.

The + and @ lines are header lines and must proceed the data lines.  The data portion consists of
a C-line declaring how many read pairs are in the ensuing read group, followed by that many 
l- and F-lines optionally followed by an N- and A-line.  In terms of a regular
expression:
```
(<+-line>|@=-line)+ (<C-line> (<L-line><F-line>[<N-line><A-line>])r)c
```
Provenance lines beginning with ! as defined in the prolog are also expected to occur in the
header section of the file.
**VGPextract** is a VGP tool that takes one or more Pacbio subreads.bam or .sam files as input
and extracts the information to make a .pbr file.

# Bionano Restriction Maps, .brm

This file encodes the primary Bionano data that is required for downstream map building or alignment to sequence. It refects the current BioNano .bnx format in supporting multiple enzymes used concurrently, although we don't currently see this in use.  As a consequence one should be particularly careful to observe the subscripted meta-values.  **Richard: are we sure we really want to support this?  See below for an alternative**

```
= <int:c> <string1> .. <stringc>     number of channels and their recognition sites
+ D <int:t∈[1,c]> <int>	             total sites for channel t
+ M <int:m>                          total number of molecules in file
@ D <int:t∈[1,c]> <int>              max. number of sites for channel t in any molecule

M <real:len> <int:n1> .. <int:nc>    molecule length, number of sites in each channel
D <int:t∈[1,c]> <int:n> <real>n      length and locations for site n
N <int:t∈[1,c]> <real1> .. <realnt>  signal:noise ratio of each site
A <int:t∈[1,c]> <real1> .. <realnt>  average intensity of sites 
```
with regular expression structure

```
(<M-line> (<D-line><N-line><A-line>)c )m
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

**VGPdigest** is a VGP tool that takes a Bionano .bnx-file as input and extracts the information to make a .brm-file.

A simpler alternative proposal requiring just a single restriction site per file would be

```
= <string>                           site recognition sequence                         
M <real:len> <int>                   molecule length, number of sites 
D <int:n> <real1> .. <realn>         length and locations for site n
N <int:n> <real1> .. <realn>         signal:noise ratio of each site
A <int:n> <real1> .. <realn>         average intensity of sites 
```

## Sequence match files, .sxs


```
+ R <int>					total number of reads in file
+ C <int:c>					total number of SMRT cells in file
+ S <int>					total bases of sequence in file
@ C <int>					max. number of bases in a SMRT cell
@ R <int>					max. number of reads in a SMRT cell
@ S <int>					max. number of bases in a read

G <string:SMRT_header> <string:file_name> <int:r>	        SMRT cell with r reads
L <int:well> <int:1st.pulse> <int:last.pulse> <real:score>	Read well and pulse range
F <string>				                	sequence
N <real:A> <real:C> <real:G> <real:T>                   	SNR in each base channel for read
A <string>				                	capped pulse widths
```


## Sequence Match List file, .sml


## Sequence file, .seq


## Scaffolding information file, .scf

```
< <string:contig_file_name> <int:nseqs>
B <int:seq> <int:start> <int:end>
Q <int:phred_confidence>	additional free information after phred score
J <int:seq_a> <int:seq_b> <[s|e]:end_a> <[s|e]:end_b>  
G <int:gap>
Q <int:phred_confidence>
```

## Scaffolding List file, .slf





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
1. VGPpair [-vs] <forward:fastq> <reverse:fastq>
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

VGPpair checks the syntax of the .fastq files but does not verify that the DNA and QV
strings are over the appropriate symbols.
