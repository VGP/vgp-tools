
# VGP Tools: Command line tools for <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; sequence data and genome assembly

### Authors:  Gene Myers, Richard Durbin, and the Vertebrate Genome Project Assembly Group
### Last Update: April 27, 2020

# Introduction

This document describes the VGP command line tools developed to
import aand operate on the
sequencing and assembly data in the [VGP tools format](file://VGP-assembly-schema.md)
that has been developed for the
[Vertebrate Genomes Project](http://www.vertebrategenomes.org).

All of the tools of the library output compressed binary One-Code files and accept as input either binary or ascii One-Code files (when the input is such).  Furthermore,
every tool takes advantage of threading to obtain speedups dependent on the number
of threads requested.  In otherwords, the tools here are designed for large-scale
bioinformatics applications and are competitive or superior to similar tools in terms of compute times and the storage consumed by intermediate files in a workflow.

# List of Tools

### <code>1. VGPzip [-x] [-C\<int(6)\>] [-T\<int(4)\>] \<file\></code>

This tool is a replacement for gzip.  It offers two big advantages:

* it is threaded and uses the libdeflate library as opposed to zlib, thus
providing much faster compression.  For example compressing a 35GB file
with ```gzip``` takes over an hour, but VGPzip takes just 2.25 minutes with the 6 cores
on a 2019 MacBook Pro (at the default compression levels).

* it automatically produces a small .vzi index file that allows applications
(like VGPseq below) to read and decompress a file with parallel threads.

Moreover, it produces a .gz file that can be used exactly as a file produced with gzip,
e.g., gunzip works on the file.  So we highly recommend that if you are going to use
the VGP tools, that you compress your fasta or fastq data with VGPzip.

VGPzip compresses the given file ```<file>``` into a blocked gzip file with the name ```<file>.gz``` and
produces an associated index in ```<file>.vzi``` if the -x option is *not* set.  The ```.gz``` file can be treated just like any other gzip'd file.  The compression is parallized with threads,
either 4 by default or the number requested by the ```-T``` option.  The default compression level is 6 but may
be set with the -C option to a level between 1 and 12, where lower numbers are less time and compression, and higher
number are more time and better compresssion.

### <code>2. VGPseq [-viqp] [-g#x] [-T\<int(4)\>] \<name:cram|[bs]am|f{ast}[aq][.gz]> ...</code>

VGPseq reads one or more cram, bam, sam, fastq, or fasta files and outputs a single file in .seq format to the standard output.
If more than one file is given then they must all be of the same type, a mix is currently not allowed.  In addition, fasta or
fastq files can optionally be compressed with either VGPzip or gzip.  VGPzip'd files are decompressed on the fly, whereas gzip'd
files must be less-efficiently expanded into a temporary file.

The file names given to VGPseq do not need to have a complete suffix designation, the
program will find the appropriate extension.  That is, if a user wishes to refer to a
file ```foo.fastq.gz``` then simply saying ```foo``` or ```foo.fastq``` on the command
line will suffice.

VGPseq is threaded 4 ways by default, but the number of threads can be explicitly
controlled with the -T parameter.  With the -v option set, VGPseq verbosely outputs to
stderr information about its progress and operation.

By default the output is a series of S-lines containing the sequence reads.  Additional
information is output as directed by the following flags:

* ```-q```: Q-lines encoding the QV scores of each base are output, *if* the information is in the files.

* ```-i```:  I-lines giving the read identifier are output.  For fasta and fastq the identifier is considered
to be the header line (that begins with > or @) up to the first white space character or new-line.

* ```-p```: Read pairs are output indicated by a P-line proceeding each pair.  This
option is only if the data in the file is *pairable*.  Basically, this amounts to the
entries or records being in order so that each read pair is consecutive in the file.
For sam, bam, and cram, the "flags" field of the underlying bam records are used to
identify the forward and reverse reads.  For fasta and fastq, the identifiers of a
pair are expected to be identical.  Pairability is checked and the program exits
prematurely with an error message if it cannot parse all the input files into pairs.

* ```-g#x```: The -g option asks VGPseq to group the data (using g-lines, see the specification) according to
a prefix of a sequence's identifier.  Specifically, the prefix ends at the ```#```<sup>th</sup> occurence of the
symbol ```x```.  For example, to group Illumina data into lanes use ```-g3:```, or to group Pacbio data
into cells use ```-g1/```.

### <code>3. VGPpair [-v] [-T\<int(4)\>] \<forward:.seq> \<reverse:.seq></code>

VGPpair reads two, presumably paired .seq files and outputs to stdout a compressed binary
.irp file in which the sequences with the same indices are paired together, with the forward sequence (and any qualifying lines, e.g. 'Q', 'W', etc) immediately preceding the reverse sequence (and its modulating lines if any).  The only condition is that the two files have
the same number of sequences.  The group structure, if any, is taken from the forward file.

### <code>4. VGPpacbio [-vaq] [-T\<int(4)\>] [-e\<expr(ln>=500 && rq>=750)>]</code> <br><code>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; \<data:.subreads.[bam|sam]> ...</code>

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

### <code>5. VGPcloud [-vH] [-t\<int(100)>] [-T\<int(4)>] \<pairs:.irp></code>

VGPcloud takes an .irp sequence file containing paired Illumina data produced with the
10x Genomics read cloud technology.  As such the first 23bp of the forward reads consist
of a 16bp barcode followed by a 7bp linker.  VGPcloud examines the bar codes and
considers **valid** any that occur 100 or more times.  This threshold can be reset by -t and to
help you choose a threshold, the -H option outputs a histogram of barcode frequencies.
Any bar code that has one difference from a unique valid code is corrected to that valid code.  Those
read pairs with valid bar codes are then sorted
and output into a compressed, binary .10x sequence file, grouped according to
their barcordes which have been removed from the forward sequences along with the linker.
If the input .irp has name \<X>.irp, then VGPcloud creates and places its output in a
file called \<X>.10x.  The -v option asks VGPcloud to output information on its progress to the standard error output.  All the steps of VGPcloud are threaded and the -T option specifies how many threads
to use (4 by default).

### <code>6. Dazz2pbr [-vagu] [-T\<int(4)\>] \<dazzler:.db\></code>

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

### <code>7. Dazz2sxs [-vidtg] [-T\<int(4)\>] \<src1:.pbr> [\<src2:.pbr>] \<dazzler:.las\> ...</code>

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
