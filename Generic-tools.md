# VGP Tools: The generic command line tools

Authors:  Gene Myers, Richard Durbin, and the Vertebrate Genome Project Assembly Group

Last Update: April 13, 2020


This document describes generic command line tools for interacting with VGP-Tools files.
For other VGP tools documentation see the top level [README.md](https://github.com/VGP/vgp-tools/blob/master/README.md).

#### <code>1.0. VGPstat [-Hu] [-o \<name>] [-t <3-code>] \<input:VGP-file></code>

VGPstat provides information about a VGP file.  Without arguments it validates an ascii file, including reporting any missing header information, and states how many objects, groups, and lines it contains.  Details of how many lines of each type are present are available in the count '@' header lines output by the -H option.

The -H option calculates a full and correct header from the contents of the file, and writes it out in ascii.  A header can be added to an ascii file that is lacking one using VGPview as explained below in the next section.

The -u option outputs the number of bytes used by each line type.

The -o option redirects the output to the named file. The default is stdout.

The -t option specifies the file type, and is required if the inspected file is an ascii file without a header, but is not needed for a binary file or an ascii file with a proper header.

#### <code>1.1. VGPview [-bhH] [-o \<filename>] [-t <3-code>] [-i \<ranges>] [-g \<ranges>] \<input:VGP-file></code>
	
VGPview is the standard utility to extract data from VGP files and convert between ascii and binary forms of the format.

The -b option outputs the file in binary.  The default is ascii.  Note that the binary form is compressed and indexed, and should be the standard form for programmatic access.

The -h option drops the header lines from ascii output.  It has no effect when writing binary because all binary files automatically generate a full header (much of which is actually written as a footer, since it can only be created once all the data is processed).

The -H option just prints out the header, in ascii.

The -o option redirects the output to the named file. The default is stdout.

The -t option specifies the file type, and is required if the inspected file is an ascii file without a header, but is not needed for a binary file or an ascii file with a proper header.

The -i and -g options make use of the binary file indices to allow random access to arbitrary sets of ojects or groups.  Legal range arguments include "0-10" which outputs the first 10 items, "7" which outputs the eighth item (remember numbering starts at 0), or compound ranges such as "3,5,9,20-33,4" which returns the requested items in the specified order.

It is possible to stream from a binary file to ascii and back from ascii to binary, so a standard pattern is 
```
   VGPview -h <binary-file> | <script operating on ascii> | VGPview -b -t <type> - > <new-binary-file>
```
Note that the script can ignore the header information, which is reconstructed by the second call to VGPview.

Regrettably it is not possible to stream a binary file into VGPview, since reading a binary file requires a seek operation to the end of the file to read footer information before reading the rest of the file.  An intermediate binary file must therefore be created to pass from ascii to binary and back again, as follows:
```
   VGPview -b -t <type> <ascii-file> > <binary-file>
   VGPview <binary-file> > <new-ascii-file>
```
This pattern has the effect of standardising an ascii file, and is the recommended way to add a header to an ascii VGP file that lacks a header.  Although some format consistency checks will be performed, if you want to fully validate a VGP file then use VGPstat.

#### <code>1.2. VGPzip [-x] [-T\<int(4)\>] \<file\></code>

This utility is not required for VGP files, given the generic compressed binary format now available
for VGP format data sets. However we retain it here for more general use, including for compressing
files in other formats prior to using VGP import tools such as VGPseq.

VGPzip compresseses the given file into a blocked gzip file with the name ```<file>.gz``` and
produces an associated index in ```<file>.vzi``` if the -x option is *not* set.  The ```.gz``` file can be decompressed with ```gunzip``` just like any other gzip file.
The file is compressed in 10MB blocks save for the last.   This program is like ```bgzip``` except that it parallelizes the compression with threads,
either 4 by default or the number requested by the ```-T``` option, and employes the much speedier LIBDEFLATE library as opposed to zlib.  This saves a great deal of
waiting around: compressing a 35GB file with ```gzip``` takes over an hour, but VGPzip takes
just 2.25 minutes with the 6 cores on a 2019 MacBook Pro (at the default compression levels).

The other distinguishing feature of ```VGPzip``` is the very large block size (```bgzip``` uses
64KB blocks).  Our goal is to have a compressed file that can be operated on in several
large pieces by parallel threads.  For this purpose, small blocks are
just a nuisance and would make the associated ```.vzi``` decompression
index excessively large.
