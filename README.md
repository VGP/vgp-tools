# The VGP-tools formats and tool kit

VGP-tools is a data representation framework with a growing collection of associated software,
initially designed in the context of the
[Vertebrate Genomes Project](http://vertebrategenomesproject.org) to operate on all the forms of 
data involved in a DNA sequencing and assembly project.   The data have a representation in a very simple ASCII file format that is easy for both humans and
programs to read and interpret.  Moreover, there is a corresponding compressed and indexed binary
version for each ASCII file so that production VGP tools are very efficient in time and the
size of the data files they manipulate.  A simple converter allows one to move between the ASCII and
binary representations of data, and a tool is provided to validate files against a schema file. A standard library with no dependencies supports reading and writing arbitrary VGP tools files.

Further documentation is available in the following files which should be in the same directory as this README

- [Format-description.md](https://github.com/VGP/vgp-tools/blob/master/Format-description.md) format and schema definition
- [Generic-tools.md](https://github.com/VGP/vgp-tools/blob/master/Generic-tools.md) generic command line tools
- [Library-interface.md](https://github.com/VGP/vgp-tools/blob/master/Library-interface.md) the standard C library interface
- [VGP-sequence-schema.md](https://github.com/VGP/vgp-tools/blob/master/VGP-sequence-schema.md) the schema developed by the Vertebrate Genomes Project for sequence data and genome assembly
- [VGP-sequence-tools.md](https://github.com/VGP/vgp-tools/blob/master/VGP-sequence-tools.md) command line tools developed for the
  Vertebrate Genomes Project to go with the sequence schema, many of which are useful for general file conversions and other
  basic operations on arbitrary DNA sequence data sets
- [VGP-assembly-workflow.md](https://github.com/VGP/vgp-tools/blob/master/VGP-assembly-workflow.md) a hypothetical work flow for genome assembly to illustrate
  potential use of the format and tools.  Note that this is all implemented - it is just for illustrative purposes.

in addition to specific technical documentation within individual files.

Authors:  Gene Myers, Richard Durbin, and the Vertebrate Genome Project Assembly Group

Date: February 2019
