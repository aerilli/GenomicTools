# genomestats

This script takes one or more genome fasta files as input and returns several stats in a csv file.

Stats retrieved:
n° of seqs, total size in nt, min size, max size, average size, min size over xbp (default=100), N50/L50, N90/L90.
Stats about gaps can be retrieved as well. Automatically detected default gap sizes (>20% entire set of gap values), can be excluded from gap stats calculation.


## Requirements

  - **Python 3.X**
  
  - **Biopython** module to process the input file
  
  - **pandas**  module to send data to a dataframe and convert it to the output file

  - *N_GAPS.py* script is required for gaps data retrieval


## Installation

```
git clone https://github.com/aerilli/GenomicTools.git
```

The scripts (genomestats.py and N_GAPS.py) will be in the Genomestats subdirectory




## Options:

   **-i**     Input file(s)\
   **-g**     Retrieve stats about gaps\
   **-D**     Include default gap values in output file (To be run with -g), otherwise redirected to stderr\
   **-f**     Output format - either csv (default) or tsv\
   **-o**     Output name (default stats)\
   **-m**     min size - accepts more than one argument\
   **-v**     Version\
   **-h**     Help



## Quick Guide


```
python3 genomestats.py -i genome1.fasta
```
To simply collect stats from the genome1.fasta file and output them in a stats.csv file.


```
python3 genomestats.py -gi genome1.fasta genome2.fasta -o genome1genome2 -f tsv
```
In this command:\
genome1 and genome2 inputs are given, gap stats are assessed, and a genome1genome2.tsv tsv file is the output.

