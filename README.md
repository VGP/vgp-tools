# The VGP-tools formats and tool kit

VGP-tools is a data representation framework with a growing collection of associated software,
initially designed in the context of the
[Vertebrate Genomes Project](http://vertebrategenomesproject.org) (VGP) to operate on all the forms of 
data involved in a DNA sequencing and assembly project.   The data have a representation in a very simple ASCII file format that is easy for both humans and
programs to read and interpret.  Moreover, there is a corresponding compressed and indexed binary
version for each ASCII file so that production VGP tools are very efficient in time and the
size of the data files they manipulate.  A simple converter allows one to move between the ASCII and
binary representations of data, and a tool is provided to validate files against a schema file. A standard library with no dependencies supports reading and writing arbitrary VGP tools files.

There are currently two packages, Core which contains the library and
general utilities, and VGP, which contains the schema and tools developed for
the Vertebrate Genomes Project.

To make the libraries and executables just type ```make``` in this top
level directory.

Further documentation is available in the following files:

- [Core/Format-description.md](https://github.com/VGP/vgp-tools/blob/master/Core/Format-description.md) format and schema definition
- [Core/Generic-tools.md](https://github.com/VGP/vgp-tools/blob/master/Core/Generic-tools.md) generic command line tools
- [Core/Library-interface.md](https://github.com/VGP/vgp-tools/blob/master/Core/Library-interface.md) the standard C library interface
- [VGP/VGP-sequence-schema.md](https://github.com/VGP/vgp-tools/blob/master/VGP/VGP-sequence-schema.md) the schema developed by the Vertebrate Genomes Project for sequence data and genome assembly
- [VGP/VGP-sequence-tools.md](https://github.com/VGP/vgp-tools/blob/master/VGP/VGP-sequence-tools.md) command line tools developed for the
  Vertebrate Genomes Project to go with the sequence schema, many of which are useful for general file conversions and other
  basic operations on arbitrary DNA sequence data sets
- [VGP/VGP-assembly-workflow.md](https://github.com/VGP/vgp-tools/blob/master/VGP/VGP-assembly-workflow.md) a hypothetical work flow for genome assembly to illustrate
  potential use of the format and tools.  Note that this is all implemented - it is just for illustrative purposes.

in addition to specific technical documentation within individual files.

Authors:  Gene Myers, Richard Durbin, and the Vertebrate Genome Project Assembly Group

Date: February 2019
