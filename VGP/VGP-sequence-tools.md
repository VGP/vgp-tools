# VGP Tools: The VGP command line tools for sequence data and genome assembly

### Authors:  Gene Myers, Richard Durbin, and the Vertebrate Genome Project Assembly Group
### Last Update: April 13, 2020

# 0. Introduction

This document describes the VGP tools command line tools developed to
interconvert with and operate on the
[sequencing data and genome assembly VGP tools formats](file://VGP-assembly-schema.md)
developed for the
[Vertebrate Genomes Project](http://www.vertebrategenomes.org).

Apart from genome assembly, we believe that these formats and some of
the early step tools will be potentially useful for other high
throughput DNA sequencing operations, and hope that the core elements
of the schema will become standardised to facilitate broad modular
data processing in large scale sequence data analysis.

For descriptions of the basic principles and properties of VGP tools
and formats, and other VGP tools documentation, see the top level
[README.md](file://README.md). Included in these is a an illustrative hypothetical [work flow](file://vgp-sequence-workflow.md)
using some of the tools listed here.

Many of the tools in the VGP repertoire, in particular data import
tools such as **VGPseq**, can perform parallel threaded processing on
VGPzip'd files, and thus operate much more efficiently.

# 1. List of tools

#### <code>1.1. VGPseq [-vsg] [-T\<int(4)\>] \<name:.fast[aq][.gz]> ...</code>

VGPseq reads one or more, possibly compressed, fasta or fastq files and outputs a single file in .seq format to the standard output.  If more than one file is given then they must
all be either .fasta or .fastq files, a mix is not allowed.
The output has Q-lines if .fastq is input and the -s option is *not* set, and does not otherwise.  It is preferrable that the inputs if compressed, were compressed
with VGPzip, as otherwise a temporary uncompressed version of each input must be created
in order to take advantage of parallel threads.

The file names given to VGPseq do not need to have a complete suffix designation, the
program will find the appropriate extension.  That is, if a user wishes to refer to a
file ```foo.fastq.gz``` then simply saying ```foo``` or ```foo.fastq``` on the command
line will suffice.

VGPseq is threaded 4 ways by default, but the number of threads can be explicitly controlled
with the -T parameter.

The -v option asks VGPseq to output information on its progress to the standard error output.
The -s option asks VGPseq to *not* output the quality values or Q-lines, if present, but just the
sequences in S-lines.
The -g option asks VGPseq to group the data into lanes.  In this case the files must
have been produced by standard Illumina software from their more basic .bcl files, and
therefore the .fastq headers encode the instrument, flow cell, lane, etc. in fields
between :'s where the data is in order of flow cell and lane.  VGPseq uses this
information to group reads into lanes.

#### <code>1.2. VGPpair [-v] [-T\<int(4)\>] \<forward:.seq> \<reverse:.seq></code>

VGPpair reads two, presumably paired .seq files and outputs to stdout a compressed binary
.irp file in which the sequences with the same indices are paired together, with the forward sequence (and any qualifying lines, e.g. 'Q', 'W', etc) immediately preceding the reverse sequence (and its modulating lines if any).  The only condition is that the two files have
the same number of sequences.  The group structure, if any, is taken from the forward file.

#### <code>1.3. VGPpacbio [-vaq] [-T\<int(4)\>] [-e<expr(ln>=500 && rq>=750)>] \<data:.subreads.[bam|sam]> ...</code>

VGPpacbio reads a sequence of Pacbio .subread.bam or subread.sam files and outputs a compressed
binary VGP .pbr file to
standard output.  As is typical for Myers' tools, the suffixes may be dropped on the command
line and the tool will find and add them.  The program is threaded, by default with 4 threads,
but the number can be explicitly set with the -T option.  For example, with six threads the program
runs about 5.5 times faster than with only one.

The -v option asks VGPpacbio to output information on its progress to the standard error output.
The -a option asks VGPpacbio to  output the arrow information in N- and A-lines per read bundle,
the default is to not output this information.  The -q option asks VGPpacbio to out Phred quality string in Q-lines, the default is to not output this information.  Please note that -q only really makes sense for HiFi data.  The reads are grouped into SMRT cells where each
input file is assumed to contain the data produced by a single cell.

#### <code>1.4. VGPcloud [-v] [-P\<dir(/tmp)>] [-T\<int(4)>] \<forward:.seq> \<reverse:.seq></code>

*Still under development but operational*

VGPcloud takes pair of .seq files containing paired Illumina data produced with the
10x Genomics read cloud technology.  As such the first 23bp of the forward reads are assumed
to consist of a 16bp barcode followed by a 7bp linker.  VGPcloud examines the bar codes and
considers valid any that occur in sufficient number and do not contain a low quality base or
'N'.  It then corrects any bar code that has one difference from a unique valid code.  Those
read pairs with valid or corrected bar codes are then sortedd
according to their barcodes in the subdirectory given by the -P option (/tmp
by default).  It then groups the pairs by barcode and outputs them as an .10x file to the
standard output after trimming off the first 23bp of the forward read.

The -v option asks VGPcloud to output information on its progress to the standard error output.
The sorts of VGPcloud are threaded and the -T option specifies how many threads to use (4 by
default).

#### <code>1.5. VGPbionano [-v] \<source:.bnx></code>

*Not yet.*

#### <code>1.6. Dazz2pbr [-vagu] [-T\<int(4)\>] \<dazzler:.db\></code>

Dazz2pbr takes a Dazzler database of a Pacbio long read data set and outputs
a compressed, binary encoding of a VGP .pbr file to the standard output whose contents
are modulated by the option flags.
As usual, suffixes are added as necessary to the command line arguments.   The -T argument specifies the number of threads to use.

The -v option asks Dazz2pbr to output information on its progress to the standard error output.
The -a option asks Dazz2pbr to output the arrow information in N- and A-lines per read bundle.
The -g option asks Dazz2pbr to organize the reads into SMRT cell groups.
The -u option asks Dazz2pbr to output the entire data set rather than the default choice
of the trimmed data set.  The Dazzler suite encodes the set of all reads output from the
instrument for a given project, but then allows a user to effectively remove from further
consideration all reads less than a given threshold length, and optionally to take only
one (the longest) read from a given well.  This is the trimmed data set that is by default
output by Dazz2pbr.

#### <code>1.7. Dazz2sxs [-vidtg] [-T\<int(4)\>] \<src1:.pbr> [\<src2:.pbr>] \<dazzler:.las\> ...</code>

Dazz2sxs takes one or more Dazzler .las file encoding a collection of local alignments found by daligner
and outputs a single compressed, binary VGP .sxs file to the standard output.  To do so, it also needs .pbr
files containing the A- and B-collections of reads that were compared, presumably produced
by Dazz2pbr above.  If the comparison is symmetric then only one .pbr file need be given.
As usual all suffix extensions are auto-completed by the program as necessary. Moreover,
the Dazzler '@' notation in the .las file names is supported.

The -v option asks Dazz2sxs to output information on its progress to the standard error output.
The -i option asks Dazz2sxs to output I-lines.
The -d option asks Dazz2sxs to output D-lines.
The -t option asks Dazz2sxs to output alignment information using T-, R- and Q-lines.
The -g option asks Dazz2sxs to output the alignments grouped according to the A-read
(facilitating the Dazzler concept of a read pile).
