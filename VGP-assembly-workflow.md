# VGP Tools: The VGP tool kit and formats

### Authors:  Gene Myers, Richard Durbin, and the Vertebrate Genome Project Assembly Group
### Date: Spring, 2019

This document presents an illustrative hypothetical application workflow for genome
assembly.  It was written relatively early in the VGP tools
development process, and has never been fully implemented, so it is
both out of date and non-functional.  However we have left it as a
document to give a sense of the thinking behind VGP tools design.

For other VGP tools documentation see the top level [README.md](https://github.com/VGP/vgp-tools/blob/master/README.md).

# A Hypothetical Assembly Work Flow

We view the assembly task as choice of a subset of read alignments that
result in a consistent layout.  These can then be used to generate
consensus contig sequences, with possible joins at branch points.  While a read may overlap multiple contigs, our convention is that contigs do not overlap, but rather abut at junctions.  Explicitly representing the subset of alignments to be used allows different tools to be used to generate alignments, select the spanning subset, then generate the graph and consensus contig sequences.

Following this different scaffolding programs may generate lists of possible breaks and joins, which can be applied in some order to generate longer contigs

Below is pseudocode for a somewhat fanciful standard assembly pipeline. In practice we do not envision working the whole time with VGP formats.  Certain series of operations will naturally be carried out on data in a native binary format.  However, we propose to require that all of our tools can convert back and forth into VGP formats to permit interchange of information.

```
// edit and assemble PacBio data
VGPpacbio pb1.bam pb2.bam pb3.bam > pb-raw.pbr
long_read_compare pb-raw.pbr pb-raw.pbr > pb-raw.sxs
pb_edit pb-raw.sxs > (pb-clean.seq pb-clean.map)
long_read_compare pb-clean.pbr pb-clean.pbr > pb-clean.sxs
transitive_reduction pb-clean.sxs > pb-clean.lis
build_consensus pb-clean.lis > (contigs1.ctg contigs1.map)

// short range scaffold with 10x data
VGPpair tenx.R1.fastq.gz tenx.R2.fastq.gz > tenx.irp
VGPcloud tenx.irp > tenx.10x
read_map contigs1.ctg tenx.10x > tenx.sxs
scaff_10x contigs1.ctg tenx.sxs > (tenx.brk tenx.jns)
apply_breaks contigs1.ctg tenx.brk > (contigs2.ctg contigs2.map)
apply_joins contigs2.ctg tenx.jns > (contigs3.ctg contigs3.map)

// long range scaffold with BioNano data
bnx2rmm bionan.bnx > bionan.rmm
rmm_assemble bionan.rmm > (bionan.rma, bionan.map)
digest_compare contigs3.ctg bionan.rma > bionan.sxr
scaff_rmc contigs3.ctg bionan.sxr > (bionan.brk bionan.jns)
apply_breaks contigs3.ctg bionan.brk > (contigs4.ctg, contigs4.map)
apply_joins contigs4.ctg bionan.jns > (contigs5.ctg contigs5.map)

// scaffold to chromosomes with HiC data
VGPpair hic.R1.fastq.gz hic2.R2.fastq.gz > hic.irp
read_map contigs5.ctg hic.irp > hic.sxs
scaff_hic contigs5.ctg hic.sxs > (hic.brk hic.jns)
apply_breaks contigs5.ctg hic.brk > (contigs6.ctg contigs6.map) 
apply_joins contigs6.ctg hic.jns > (contigs7.ctg contigs7.map) 

// build scaffold objects
scaffolds contigs7.ctg hic.jns > scaffolds7.scf

// and polish, first with PacBio data
compose_maps contigs7.map contigs6.map contigs5.map contigs4.map contigs3.map contigs2.map contigs1.map pb-clean.map > contigs7-pbraw.map
pb_polish contigs7.ctg scaffolds1.scf contigs7-pbraw.map pb-raw.pbr > (contigs8.ctg contigs8.map scaffolds8.scf)

// and then with 10x Illumina data
read_map contigs8.ctg tenx.10x > contigs8-tenx.sxs
illumina_polish contigs8.ctg scaffolds8.scf contigs8-tenx.sxs > (contigs9.ctg contigs9.map scaffolds9.scf)

// finally export the scaffolds into fasta for everyone to use
export_fasta contigs9.ctg > assembly.fa
export_agp scaffolds9.scf > assembly.agp
```

An alternative, enabled by having proposed scaffolding operations in VGP formats, would be to first make all the .jns and .brk files then select a consensus of them, and apply this to generate the final contig and scaffold sets.  This process could be repeated once, to allow scaffolding operations to act on revised contig sets.
