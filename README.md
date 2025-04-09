# The One-Code Data Framework

One-Code is a data representation framework with a growing collection of associated software,
initially designed in the context of the
[Vertebrate Genomes Project](http://vertebrategenomesproject.org) (VGP) to operate on all the forms of 
data involved in a DNA sequencing and assembly project. It developed into a general purpose schema-driven data representation framework, providing interconvertible binary and ascii record support, with (for the binary form) compression and indexing.

**The version of OneCode in this repository is frozen and out of date.  The latest version can be found at [https:github.com/thegenemyers/ONEcode](https:github.com/thegenemyers/ONEcode).**

Data is represented in a very simple ASCII file format that is easy for both humans and
programs to read and interpret.  Moreover, there is a corresponding compressed and indexed binary
version for each ASCII file so that production tools built with the One-Code library are very efficient in time and the
size of the data files they manipulate.  All fields are strongly typed and a specific collection of data types constituting a *schema* is
specified in a schema file, that is itself a One-Code data file.  A generic converter allows one to move between the ASCII and
binary representations of data, and another core tool is provided to validate files against a 
given schema.

There are currently two packages.  **Core** contains the general One-Code command-line utilities
alluded to above, and a C-library supporting the development of One-Code applications.
**VGP** contains the VGP-formats schema designed specifically for DNA sequencing and assembly applications, and a growing number of VGP-tools for importing and operating on data in this schema.

To make all libraries and command line tools just type ```make``` in this top
level directory.  Both packages have no dependencies on other software.

Further documentation is available as follows:

**CORE**

- [Framework Description](https://github.com/VGP/vgp-tools/blob/master/Core/Format-description.md): definition of One-Code schemas and data representation
- [Generic Tools](https://github.com/VGP/vgp-tools/blob/master/Core/Generic-tools.md): generic One-Code command line tools
- [Library Interface](https://github.com/VGP/vgp-tools/blob/master/Core/Library-interface.md): the standard C library interface for One-Code developers

**VGP**

- [The VGP Schema](https://github.com/VGP/vgp-tools/blob/master/VGP/docs/VGP-sequence-schema.md): the VGP formats schema for sequence data and genome assembly
- [The VGP Tools](https://github.com/VGP/vgp-tools/blob/master/VGP/docs/VGP-sequence-tools.md): command line tools for importing and operating on data in VGP-formats.  Many of these are useful for general file conversions and other
  basic operations on arbitrary DNA sequence data sets
- [Assembly Workflow Example](https://github.com/VGP/vgp-tools/blob/master/VGP/docs/VGP-assembly-workflow.md): a hypothetical work flow for genome assembly to illustrate
  potential use of VGP formats and tools.  Note that this is *not* all implemented - it is just for illustrative purposes.

In addition, specific technical documentation can be found within individual source files.

Authors:  Gene Myers, Richard Durbin, and the Vertebrate Genome Project Assembly Group

Date: February 2019
