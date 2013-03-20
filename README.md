#Overview

SCUBAT (Scaffolding Contigs Using BLAT And Transcripts) uses any set of transcripts to identify cases where a transcript is split over multiple genome fragments and attempts to use this information to scaffold the genome. 

#Software Requirements:

[CAP3](http://www.ncbi.nlm.nih.gov/pubmed/10508846) - to assemble scaffolds 

[GNU parallel](http://www.gnu.org/software/parallel/) - to run the CAP3 assemblies in parallel 

#Data requirments

1. A genome in FASTA format
2. A set  of transcripts in FASTA format
3. A BLAT psl file aligning the transcripts to the genome

#Quick guide

`scubat.pl -t genome.fa -q transcripts.fa -p blat_output.psl`

To see all options just run `scubat.pl`

#Details

SCUBAT procedes via 6 steps:

###1)  Align the transcripts to the genome (not part of the script):

  - Externally - Use BLAT to align the transcripts to the genome with psl output, e.g. `blat -t=dna genome.fa -q=dna transcripts.fa -noHead out.psl`

###2)  Identify informative split transcripts:

  - All transcripts that mapped to more than one contig are flagged. Each mapping is then analysed to identify those that map consecutive non-overlapping sections of the transcript on separate contigs allowing for an overlap buffer zone of up to 20bp.

###3) Create scaffolds:

- Each transcript-contig complex is assembled by orientating the contigs based on the
BLAT information and adding 10 N’s in-between the contigs. 

###4) Cluster scaffolds into groups and assemble:

- The contigs used in each transcript-contig complex are then cross-referenced and any complexes sharing a contig are grouped. Groups are then assembled using CAP3 with default parameters.

###5) Filter the assemblies:

- The assemblies are parsed to identify successful assembly events. Failed assemblies are identified and the largest complex kept, the remainder removed.

###6) Create new contig set:

- A new contig set is generated using the successful CAP3 assemblies, the remaining complexes and any contig not involved.
