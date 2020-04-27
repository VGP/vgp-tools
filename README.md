# The One-Code Data Framework

One-Code is a data representation framework with a growing collection of associated software,
initially designed in the context of the
[Vertebrate Genomes Project](http://vertebrategenomesproject.org) (VGP) to operate on all the forms of 
data involved in a DNA sequencing and assembly project.   Data is represented in a very simple ASCII file format that is easy for both humans and
programs to read and interpret.  Moreover, there is a corresponding compressed and indexed binary
version for each ASCII file so that production tools built with the One-Code library are very efficient in time and the
size of the data files they manipulate.  A specific collection of data types constituting a *schema* is
specified in a schema file, that is itself a One-Code data file.  A generic converter allows one to move between the ASCII and
binary representations of data, and anoter core tool is provided to validate files against a 
given schema.

There are currently two packagess.  **Core** contains the general One-Code command-line utilities
alluded to above, and a C-library supporting the developement of One-Code applications.
**VGP** contains the VGP-formats schema designed specifically for DNA sequencing and assembly app's, and a growwing number of VGP-tools for importing and operating on data in this schema.
To make all libraries and command line tools just type ```make``` in this top
level directory.  Both packages have no dependencies on other software.

Further documentation is available in the following files:

- [Core/Format-description.md](https://github.com/VGP/vgp-tools/blob/master/Core/Format-description.md) definition of One-Code schemas and data representation
- [Core/Generic-tools.md](https://github.com/VGP/vgp-tools/blob/master/Core/Generic-tools.md) generic One-Code command line tools
- [Core/Library-interface.md](https://github.com/VGP/vgp-tools/blob/master/Core/Library-interface.md) the standard C library interface for One-Code developers
- [VGP/VGP-sequence-schema.md](https://github.com/VGP/vgp-tools/blob/master/VGP/VGP-sequence-schema.md) the VGP-formats schema for sequence data and genome assembly
- [VGP/VGP-sequence-tools.md](https://github.com/VGP/vgp-tools/blob/master/VGP/VGP-sequence-tools.md) command line tools for importing and operating on data in VGP-formats.  Many of these are useful for general file conversions and other
  basic operations on arbitrary DNA sequence data sets
- [VGP/VGP-assembly-workflow.md](https://github.com/VGP/vgp-tools/blob/master/VGP/VGP-assembly-workflow.md) a hypothetical work flow for genome assembly to illustrate
  potential use of VGP formats and tools.  Note that this is *not* all implemented - it is just for illustrative purposes.

in addition to specific technical documentation within individual files.

Authors:  Gene Myers, Richard Durbin, and the Vertebrate Genome Project Assembly Group

Date: February 2019
