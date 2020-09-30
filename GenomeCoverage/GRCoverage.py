#!/usr/bin/env python3

'''
GRCoverage
takes genomic reads input files.fastq, 
maps them to the reference genome.fasta (via bwa if paired end reads are given, via minimap2 if long reads file is given),
assesses genome coverage with bedtools and
outputs a .genomecoverage file

Instructions to run the script:
GRCoverage.py [-h --help] ILLUMINAPAIREDEND[-I --IPE] LONGREADS[-L --LongReads]  [-v -version]

First positional input: reference genome
Second: Read1 (IPE), Read (LongRead)
Third: Read2 (IPE)

Version 0.0.2 (Beta)
'''

#libraries
import subprocess
import os
import argparse

##Setting flags
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script takes genomic reads input files.fastq, maps them to the reference genome.fasta (via bwa if paired end reads are given, via minimap2 if one long reads file is given), assesses genome coverage with bedtools and outputs a .genomecoverage file')
    parser.add_argument('input_genome', type=str, help='input: reference genome')
    parser.add_argument('input_reads', nargs='*', type=str, help="input read(s) - read1 read2, as Illumina PE reads")
    parser.add_argument('-i','--index', help='indexing genome first', action='store_true', default=False)
    parser.add_argument('-I','--IPE', help='Illumina PE reads as reads input (read1 read2)', action='store_true', default=False)
    parser.add_argument('-L','--LongReads', help='Long reads as query input', action='store_true', default=False)
    parser.add_argument('-v','--version', action='version', version='0.0.1 (Beta)')
    Args = parser.parse_args()


##Setting functions
def short_indexing(g): #indexing genome function for short reads
    print("Indexing genome first..."+'\n\n')
    ind = subprocess.run(['bwa', 'index', g]) #bwa index cmd
    print(ind.args)

#function for short PE reads processing
def shortPE(read1, read2, genome):
    if '_R1' in read1: #setting name variable to name downstream files 
        name = os.path.splitext(read1)[0].replace('_R1','')
    elif 'R1' in read1:
        name = os.path.splitext(read1)[0].replace('R1','')
    else:
        name = os.path.splitext(read1)[0]

    if Args.index == True: #if -i is selected genome is indexed first
        short_indexing(genome)
 
    print('Mapping reads...'+'\n') #bwa mem for mapping reads
    mem = subprocess.run('bwa mem %s %s %s' % (genome, read1, read2), stdout=subprocess.PIPE, shell=True)
    print(mem.args+'\n')
    #samtools view for converting sam to bam file
    sam = subprocess.run('samtools view -F 4 -Sb -o %s.bam - ' % name, input=mem.stdout, shell=True)
    print(sam.args+'\n')
    print('Sorting BAM file...'+'\n') #samtools sort for sorting bam file
    sort = subprocess.run('samtools sort -@30 -o %s_sorted.bam %s.bam' % (name,name), shell=True)
    print(sort.args+'\n')
    print('Creating genomecoverage file...'+'\n') #bedtools genomecov to retrieve genome coverage data
    bed = subprocess.run('bedtools genomecov -ibam %s_sorted.bam > %s.genomecoverage' % (name,name), shell=True)
    print(bed.args+'\n')

def long(query, target):
    t = os.path.splitext(target)[0] #setting variables to name downstream file
    q = os.path.splitext(query)[0]
    if Args.index == True: #if -i, genome is indexed first
            ind = subprocess.run('minimap2 -d %s.mmi %s' % (t, target), shell=True)
            print(ind.args+'\n')
    mini =subprocess.run('minimap2 -a %s -o %s.sam %s' % (target, q, query), shell=True) #minimap2 for mapping long reads
    #samtools view for converting sam to bam file
    sam = subprocess.run('samtools view -F 4 -Sb -o %s.sam %s.bam' % (q, q), shell=True) 
    print(mini.args+'\n')
    print('Sorting BAM file...'+'\n') #samtools sort for sorting bam file
    sam = subprocess.run('samtools sort -@30 -o %s_sorted.bam %s.bam' % (q,q), shell=True)
    print(sam.args+'\n')
    print('Creating genomecoverage file...'+'\n') #bedtools genomecov to retrieve genome coverage data
    bed = subprocess.run('bedtools genomecov -ibam %s_sorted.bam > %s.genomecoverage' % (q,q), shell=True)
    print(bed.args+'\n')




#If -I is selected => PE reads are processed
if Args.IPE == True and len(Args.input_reads)==2:
    shortPE(read1=Args.input_reads[0], read2=Args.input_reads[1], genome=Args.input_genome)

#If -L is selected => long reads are processed instead
elif Args.LongReads == True and len(Args.LongReads)==1:
    long(target=Args.input_genome, query=Args.input_reads[0])

