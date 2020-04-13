# VGP Tools: The VGP tool formats and tool kit

VGP-Tools is a data representation framework with a growing collection of associated software,
initially designed in the context of the
[Vertebrate Genomes Project](http://vertebrategenomesproject.org) to operate on all the forms of 
data involved in a DNA sequencing and assembly project.   The specification of the content of the
several types of data files involved has a very simple ASCII format that is easy for both humans and
programs to read and interpret.  Moreover, there is a corresponding compressed and indexed binary
representation for each ASCII datum so that production VGP tools are very efficient in time and the
size of the data files they manipulate.  A simple converter allows one to move between the ASCII and
binary representations of data. A standard library with no dependencies supports reading, writing
and validation of VGP tools files.

Further overview documentation includes the following files which should be in the same directory as this README

- [VGP-file-structure.md](https://github.com/VGP/vgp-tools/blob/master/VGP-file-structure.md) format structure and schema definition
- [VGP-generic-tools.md](https://github.com/VGP/vgp-tools/blob/master/VGP-generic-tools.md) generic command line tools
- [VGP-library-interface.md](https://github.com/VGP/vgp-tools/blob/master/library-interface.md) the standard C library interface
- [VGP-sequence-schema.md](https://github.com/VGP/vgp-tools/blob/master/sequence-schema.md) the schema developed by the Vertebrate Genomes Project for sequence data and genome assembly
- [VGP-sequence-tools.md](https://github.com/VGP/vgp-tools/blob/master/VGP-sequence-tools.md) command line tools developed for the
  Vertebrate Genomes Project to go with the sequence schema, many of which are useful for general file conversions and other
  basic operations on arbitrary DNA sequence data sets
- [VGP-assembly-workflow.md](https://github.com/VGP/vgp-tools/blob/master/VGP-assembly-workflow.md) a hypothetical work flow for genome assembly to illustrate
  potential use of the format and tools.  Note that this is not up to date or complete - it is just for illustrative purposes.

in addition to specific technical documentation within individual files.

Authors:  Gene Myers, Richard Durbin, and the Vertebrate Genome Project Assembly Group

Date: February 2019
