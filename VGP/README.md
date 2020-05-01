# VGP-Tools: The VGP formats and tools for genome <br>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; sequence assembly and related activities

### Authors:  Gene Myers, Richard Durbin, and the Vertebrate Genome Project Assembly Group
### Last Update: April 13, 2020

This subdirectory contains the documents and tools designed within the One-Code
framework to support the [Vertebrate Genomes Project](http://www.vertebrategenomes.org)
with data representations for all the source data, process intermediates, and resulting
reconstructed genome assemblies involved in a large-scale DNA sequencing project.
Along with these formats there is a growing set of documented
[command line tools](https://github.com/VGP/vgp-tools/blob/master/VGP/docs/VGP-sequence-tools.md)
for carrying out file conversions and other operations on these objects. There is also
available an illustrative hypothetical
[work flow](https://github.com/VGP/vgp-tools/blob/master/VGP/docs/VGP-assembly-workflow.md)
using these formats and tools (including some tools not yet written).

Apart from genome assembly, we believe that these formats and some of
the current tools will be potentially useful for other high
throughput DNA sequencing operations, and hope that the core elements
of the schema will become standardised to facilitate broad modular
data processing in large scale sequence data analysis.

As mentioned initially, the VGP tools data schema is a specific instance of a "One-Code" schema
which is a general framework for representing and manipulating data whose
primary advantages are a very simple ASCII format for each data object type
that is easy for both humans and programs to read and interpret,
with a corresponding compressed and indexed binary version so that production
software can be used for large scale data intensive applications.  Furthermore,
adding and modifying data types is simply a matter of editing the simple
specification of the schema [VGP_1\_0.def](https://github.com/VGP/vgp-tools/blob/master/VGP/VGP_1_0.def).
For descriptions of the basic principles and properties of the One-Code
framework, we highly recommend you begin by reading the
[top level document](https://github.com/VGP/vgp-tools/blob/master/README.md).  The current documentation specifically for the VGP formats and tools is as follows:

- [The VGP Formats](https://github.com/VGP/vgp-tools/blob/master/VGP/docs/VGP-sequence-schema.md): the VGP formats schema for sequence data and genome assembly
- [The VGP Tools](https://github.com/VGP/vgp-tools/blob/master/VGP/docs/VGP-sequence-tools.md): command line tools for importing and operating on data in VGP-formats.  Many of these are useful for general file conversions and other
  basic operations on arbitrary DNA sequence data sets
- [Assembly Workflow Example](https://github.com/VGP/vgp-tools/blob/master/VGP/docs/VGP-assembly-workflow.md): a hypothetical work flow for genome assembly to illustrate
  potential use of VGP formats and tools.  Note that this is *not* all implemented - it is just for illustrative purposes.

In addition, specific technical documentation can be found within individual source files.

