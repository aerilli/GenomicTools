# GRCoverage

This script:

  1. takes genomic reads input files.fastq
  2. maps them to the reference genome.fasta, via bwa if Illumina paired-end reads files are given or via minimap2 if Oxford Nanopore/PacBio reads file is given
  3. assesses genome coverage, using bedtools
  4. outputs a .genomecoverage file



## Requirements

  - **Python 3.X**\
  
  - **bwa** to align paired-end reads to the reference genome

  - **minimap2**  to align long reads

  - **samtools**  to convert sam files to bam files and to sort bam files

  - **bedtools**   to retrieve genome coverage




## Installation

```
git clone https://github.com/aerilli/GenomicTools.git
```

The script will be in the GenomeCoverage subdirectory




## Options:

   **-i**    Index genome before alignment\
   **-I**    Paired-end reads input: 2nd and 3rd positional arguments\
   **-L**    Long reads input: 2nd positional argument 
   
Note: Reference genome as first positional argument



## Quick Guide

If genome coverage of Illumina paired-end reads is to be assessed, use:

```
python3 GRCoverage.py -I refgenome.fasta read1.fastq read2.fastq
```


If genome coverage of OxfordNanopore/PacBio long reads is to be assessed, use:

```
python3 GRCoverage.py -L genome_target.fasta reads_query.fastq
```

