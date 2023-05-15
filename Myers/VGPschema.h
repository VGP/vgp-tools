/*  File: VGPschema.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2020
 *-------------------------------------------------------------------
 * Description: schema to include in VGP ONEcode applications
 * Exported functions:
 * HISTORY:
 * Last edited: May 14 03:11 2020 (rd109)
 * Created: Wed May 13 23:34:04 2020 (rd109)
 *-------------------------------------------------------------------
 */

static char *vgpSchemaText =
  "1 3 def 1 0  OneCode schema for VGP genome assembly pipeline and related purposes\n"
  "! 4 14 Richard Durbin 1 0 1 0 7 by hand 24 Wed Apr 22 20:00:00 2020\n"
  ".\n"
  ". these are blank lines, which can also be used for arbitrary comments, anywhere in the file\n"
  ".\n"
  ". below P lines give primary file types, and define the objects in this file\n"
  ". then S is for secondary file types\n"
  ". note we can also add arbitrary comments/descriptions after the data on each line\n"
  ".\n"
  "P 3 seq SEQUENCE\n"
  "S 3 irp   read pairs\n"
  "S 3 pbr   pacbio reads\n"
  "S 3 10x   10X Genomics data\n"
  "S 3 ctg   contigs from an assembly\n"
  "S 3 kmr   kmers\n"
  "G g 2 3 INT 6 STRING               group: count, name (e.g. use for flow cell/lane grouping)\n"
  "O S 1 3 DNA                        sequence: the DNA string\n"
  "D I 1 6 STRING                     id: sequence identifier\n"
  "D Q 1 6 STRING                     quality: Q values (ascii string = q+33)\n"
  "D P 0                              marks start of a readpair\n"
  "D W 4 3 INT 3 INT 3 INT 4 REAL     PacBio metadata: well, pulse start, pulse end, score\n"
  "D N 4 4 REAL 4 REAL 4 REAL 4 REAL  read Signal to Noise Ratio: values in A,C,G,T channels  \n"
  "D A 1 6 STRING                     PacBio capped pulse widths: values in [1,4] inclusive\n"
  "D C 1 3 INT                        count (for kmers)    \n"
  ".\n"
  "P 3 rmp RESTRICTION MAP\n"
  "S 3 rmm   maps of single molecules, e.g. BioNano primary data\n"
  "S 3 rms   maps from sequence\n"
  "S 3 rma   maps from assembly of molecule level data\n"
  "G r 3 3 INT 3 INT 11 STRING_LIST   group: count, number of 'enzymes', corresponding patterns\n"
  "O R 2 3 INT 8 INT_LIST             map: length, site locations (both in base pairs)\n"
  "D E 1 8 INT_LIST                   enzyme: for each site the enzyme index into list in line r\n"
  "D I 1 9 REAL_LIST                  intensities: intensity of signal for each site in map\n"
  "D N 1 9 REAL_LIST                  noise: signal to noise ratio at each site in map\n"
  "D O 1 3 INT                        index of the sequence object from which map is calculated\n"
  ".\n"
  "P 3 aln ALIGNMENT\n"
  "S 3 sxs   sequence to sequence alignment\n"
  "S 3 rxr   restriction map to restriction map alignment\n"
  "S 3 sxr   sequence to restriction map alignment\n"
  "S 3 map   relationship between two versions of an object, e.g. contigs before/after polishing\n"
  "G g 2 3 INT 6 STRING                       group: count, name\n"
  "O A 2 3 INT 3 INT                          alignment: a index in its file, b index in its file\n"
  "D I 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT  start, end, len for A, start, end, len for B\n"
  "D Q 1 3 INT                                quality: alignment confidence in phred units\n"
  "D M 1 3 INT                                match: number of matching bases\n"
  "D D 1 3 INT                                differences: number of diffs = substitions + indels\n"
  "D C 1 6 STRING                             cigar string: encodes precise alignment\n"
  "D U 1 8 INT_LIST                           trace points in a\n"
  "D V 1 8 INT_LIST                           trace points in b\n"
  "D T 1 3 INT                                trace point spacing in a - global until reset\n"
  "D W 1 8 INT_LIST                           trace point separations in b\n"
  "D X 1 8 INT_LIST                           number of differences in alignment per trace panel\n"
  ".\n"
  "P 3 hit HIT_LIST\n"
  "S 3 s2k   sequence to kmer hit list\n"
  "S 3 k2s   kmer to sequence hit list\n"
  "O H 2 3 INT 8 INT_LIST     hits: indices of query a in its file, targets b in their file\n"
  "D O 1 8 INT_LIST           offsets: in query a of each target b\n"
  "D P 1 8 INT_LIST           positions: in each target b of query a\n"
  ".\n"
  "P 3 jns JOIN\n"
  "O J 6 3 INT 3 INT 4 CHAR 3 INT 3 INT 4 CHAR join: a, a-pos, a-dir s or e, b, b-pos, b-dir\n"
  "D G 2 3 INT 3 INT                           gap: size estimate, standard dev. estimate (bp)\n"
  "D Q 1 3 INT                                 quality: confidence in phred units\n"
  "D E 1 8 INT_LIST                            evidence: list of alignments supporting the join\n"
  ".\n"
  "P 3 brk BREAK\n"
  "O B 3 3 INT 3 INT 3 INT    break: object, start, end - material in [start,end] uncertain\n"
  "D Q 1 3 INT                quality: confidence in phred units\n"
  "D E 1 8 INT_LIST           evidence: list of alignments supporting the break\n"
  ".\n"
  "P 3 lis LIST    \n"
  "S 3 lyo   layout for assembly: selection of alignments generating a contig\n"
  "S 3 scf   scaffold: list over joins\n"
  "O L 1 8 INT_LIST           list: indexes of objects in list in reference file\n"
  "D N 1 6 STRING             name: optional name for list\n"
  "D S 1 3 INT                seed: optional seed sequence for scaffold\n"
  ".\n"
   ". end of file\n" ;
