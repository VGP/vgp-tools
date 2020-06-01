static char *VGP_SPEC = 

"1 3 def 1 0  OneCode schema for VGP genome assembly pipeline and related purposes \n\
! 14 Richard Durbin 1 0 1 0 7 by hand 24 Wed Apr 22 20:00:00 2020 \n\
. \n\
. these are blank lines, which can also be used for arbitrary comments \n\
. \n\
. below P lines give primary file types, and define the objects in this file \n\
. then S is for secondary file types \n\
. and L, C, F and B define the legal record types, all taking a CHAR, then STRING_LIST \n\
.     L for standard linetypes saved with no compression in binary mode \n\
.     C means compress the list - perhaps all lists should be compressed by default? \n\
.     F means compress the fields - not default because it slows IO down for small gain \n\
.     B means compress both \n\
. note we can, and do, also add arbitrary comments/descriptions after the data on each line \n\
. \n\
V 1 0  versions: major, minor.  This is a global that carries forward until changed. \n\
. \n\
P 3 seq SEQUENCE \n\
S 3 irp   read pairs \n\
S 3 pbr   pacbio reads \n\
S 3 10x   10X Genomics data \n\
S 3 ctg   contigs from an assembly \n\
S 3 kmr   kmers \n\
L g 2 3 INT 6 STRING               group: count, name (e.g. use for flow cell/lane grouping) \n\
L S 1 3 DNA                        sequence: the DNA string \n\
C I 1 6 STRING                     id: (optional) sequence identifier \n\
C Q 1 6 STRING                     quality: Q values (ascii string = q+33) \n\
L P 0                              marks start of a readpair \n\
F W 4 3 INT 3 INT 3 INT 4 REAL     PacBio metadata: well, pulse start, pulse end, score \n\
F N 4 4 REAL 4 REAL 4 REAL 4 REAL  read Signal to Noise Ratio: values in A,C,G,T channels  \n\
C A 1 6 STRING                     PacBio capped pulse widths: values between 1 and 4 inclusive \n\
L C 1 3 INT                        count (for kmers)    \n\
. \n\
P 3 rmp RESTRICTION MAP \n\
S 3 rmm   maps of single molecules, e.g. BioNano primary data \n\
S 3 rms   maps from sequence \n\
S 3 rma   maps from assembly of molecule level data \n\
L r 3 3 INT 3 INT 11 STRING_LIST   group: count, number of 'enzymes', corresponding patterns \n\
C R 2 3 INT 8 INT_LIST             map: length, site locations (both in base pairs) \n\
C E 1 8 INT_LIST                   enzyme: for each site the enzyme index from the list in line r \n\
C I 1 9 REAL_LIST                  intensities: intensity of signal for each site in map \n\
C N 1 9 REAL_LIST                  noise: signal to noise ratio at each site in map \n\
L O 1 3 INT                        index of the sequence object from which map is calculated \n\
. \n\
P 3 aln ALIGNMENT \n\
S 3 sxs   sequence to sequence alignment \n\
S 3 rxr   restriction map to restriction map alignment \n\
S 3 sxr   sequence to restriction map alignment \n\
S 3 map   relationship between two versions of an object, e.g. contigs before/after polishing \n\
L g 2 3 INT 6 STRING                       group: count, name \n\
F A 2 3 INT 3 INT                          alignment: a index in its file, b index in its file \n\
F I 6 3 INT 3 INT 3 INT 3 INT 3 INT 3 INT  indices: aStart, aEnd, aLength, bStart, bEnd, bLength \n\
F Q 1 3 INT                                quality: alignment confidence in phred units \n\
F M 1 3 INT                                match: number of matching bases \n\
F D 1 3 INT                                differences: number of diffs = substitions + indels \n\
C C 1 6 STRING                             cigar string: encodes precise alignment \n\
C U 1 8 INT_LIST                           trace points in a \n\
C V 1 8 INT_LIST                           trace points in b \n\
L T 1 3 INT                                trace point spacing in a - global until reset \n\
C W 1 8 INT_LIST                           trace point separations in b \n\
C X 1 8 INT_LIST                           number of differences in alignment per trace interval \n\
. \n\
P 3 hit HIT_LIST \n\
S 3 s2k   sequence to kmer hit list \n\
S 3 k2s   kmer to sequence hit list \n\
C H 2 3 INT 8 INT_LIST     hits: indices of query a in its file, targets b in their file \n\
C O 1 8 INT_LIST           offsets: in query a of each target b \n\
C P 1 8 INT_LIST           positions: in each target b of query a \n\
. \n\
P 3 jns JOIN \n\
F J 6 3 INT 3 INT 4 CHAR 3 INT 3 INT 4 CHAR join: a, pos in a, dirxn s or e, b, pos in b, dirxn  \n\
F G 2 3 INT 3 INT                           gap: size estimate, standard deviation estimate (bp) \n\
L Q 1 3 INT                                 quality: confidence in phred units \n\
C E 1 8 INT_LIST                            evidence: list of alignments supporting the join \n\
. \n\
P 3 brk BREAK \n\
L B 3 3 INT 3 INT 3 INT    break: object, start, end - material between start and end uncertain \n\
L Q 1 3 INT                quality: confidence in phred units \n\
C E 1 8 INT_LIST           evidence: list of alignments supporting the break \n\
. \n\
P 3 lis LIST    \n\
S 3 lyo   layout for assembly: selection of alignments generating a contig \n\
S 3 scf   scaffold: list over joins \n\
C L 1 8 INT_LIST           list: indexes of objects in list in reference file   \n\
C N 1 6 STRING             name: optional name for list \n\
L S 1 3 INT                seed: optional seed sequence for scaffold \n\
. \n\
. end of file \n";
