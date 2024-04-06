# My Project 

This is my project investigating the hybrid lineage of *Cryptococcus neoformans* called VNIII Serotype AD Hybrid. All analysis conducted locally or on the University of Exeter's ISCA Server

**Data Collection**
- 73 _Cryptococcus_ isolates have come from the Farrer Lab
- 39 Samples gathered from the SRA Database

**Genome Alignment and Variant Calling**
- Genome Alignment using the GATK_no_scattergather.sh -r reference.fasta -b folder with paired FASTQ -s sample name
- Done through the use of a mamba environment:
  -   mamba create -n SNP
  -   mamba install BWA
  -   mamba install Picard

