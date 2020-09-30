#!/usr/bin/env python3

'''

genomestats takes one or more genome fasta input files and returns several stats in a csv file:
n° of seqs, total size in nt, min size, max size, average size, min size over xbp (default=100, for x value see MINSIZEBP), N50/L50, N90/L90.
Stats about gaps can be retrieved (see GAPSTATS). 
Automatically detected default gap sizes (>20% entire set of gap values), can be excluded from gap stats calculation (see DEFVAL(GAPSTATS)).

Instructions to run the script:
genomestats.py [-h --help] INPUTFILE[--infile -i] MINSIZEBP[--minsize -m] OUTPUTNAME[--name -n] OUTPUTFORMAT[--output -o] GAPSTATS[--gapstats -g] DEFVALUES(GAPSTATS)[--DEFVAL -D] [-v --version]
INPUTFILE is required

Version 0.0.2 (Beta)

'''

#importing libraries
import argparse
from os import path
from os import stat
from Bio import SeqIO
import pandas as pd
import N_GAPS

#SETTING FLAGS
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script takes a fasta input file (-i: required!) and returns several stats in a csv file: n° of seqs, total size in nt, min_size, max_size, average_size, min-size_xbp(for x value, see MINSIZE), N50/L50, N90/L90; Stats about N gaps can be retrieved (see GAPSTAT) and automatically detected default gap sizes (>20% entire set of gap values), can be excluded from gap stats calculation (see DEFVAL)\nusage: creepyfasta.py [-h] INPUTFILE[--infile -i] MINSIZEBP[--minsize -m] OUTPUTNAME[--name -n] OUTPUTFORMAT[--output -o] GAPSTATS[--gapstats -g] DEFVALUES(GAPSTATS)[--DEFVAL -D] [-v]')
    parser.add_argument('-i','--infile', nargs='*', type=str, help="input fasta file(s)", required=True)
    parser.add_argument('-g', '--gapstats', help='Retrieve stats about gaps (N)', action='store_true')
    parser.add_argument('-D','--DEFVAL', help='To be run with -g - Consider default gap size values, otherwise removed (stderr)', action='store_true', default=False)
    parser.add_argument("-f","--format", help='output format - either csv or tsv', type=str, default='csv')
    parser.add_argument("-o","--output", help='output name', type=str, default='stats')
    parser.add_argument("-m","--minsize", nargs='*', help='min size(bp) - accepts more than one argument', type=int, default=[100])
    parser.add_argument('-v','--version', action='version', version='0.0.1 (Beta)')
    args = parser.parse_args()



#type of output file
if args.output == 'csv':
    output = ','
elif args.output == 'tsv':
    output = '\t'
else:
    raise Exception('Invalid output format')


##FUNCTIONS

#Counting SeqN
def count_seq(d):
    return len(d)+1

#Calculating the total nucleotide size 
def count_nt(d):
    d_list = list(d.values())
    tot_nt = 0
    for dl in d_list:
        tot_nt += len(dl) #incrementing by the length of each Seq
    return tot_nt

#sorting the values to retrieve... 
def sort_seq(d):
    d_list = list(d.values())
    len_list = []
    for dl in d_list:
        len_list.append(len(dl))
    len_list.sort()
    return len_list

#...min size
def min_seq(d):
    return sort_seq(d)[0]

#...max size
def max_seq(d):
    return sort_seq(d)[-1]

#...and avg size
def avg_seq(d):
    average = sum(sort_seq(d))/len(sort_seq(d))
    return round(average,2)

#min size sample n
def min_size_Xbp(n,d):
    above_n = [seq for seq in sort_seq(d) if seq >= n] #sorted Seq lengths above n threshold
    Xseq = len(above_n) 
    sum_above_n = sum(above_n) #total summation of seq length above n
    percent_set = sum_above_n/count_nt(d)*100 #percentage of the set containing min size of n
    return n, Xseq, round(percent_set,2)


#Define NX/LX, with X = 50 OR 90
def NXoverLX(x, d):
    x/=100
    S = sort_seq(d)
    S.reverse() #sorted and reverted Seq lengths
    sum_seq= 0
    counter = 0
    for seq in S: #starting from bigger length, summation until threshold is reached
        sum_seq += seq
        counter += 1
        if sum_seq > x:
            break
        else:
            continue
    return sum_seq, counter #N50 and L50


#Creating dataframe
row_list_init = {'Stat':
    ['Seq-N',
    'Tot_Size',
    'Min_Size',
    'Max_Size',
    'Avg_Size',
    'N50',
    'L50',
    'N90',
    'L90']}



q = 5
for a in args.minsize:
    row_list_init['Stat'].insert(q,'min_size_%ibp_Xseq' % a) 
    q+=1
    row_list_init['Stat'].insert(q,'min_size_%ibp_SetPercentage' % a)
    q+=1

df = pd.DataFrame(row_list_init)

for fastafile in args.infile:

#Does the file exist?
    if path.exists(fastafile): 
        print('I am opening your file...')
    elif not path.exists(fastafile):
        print("Cannot find %s file\n" % fastafile)
        quit()
    elif stat(fastafile).st_size == 0: #OR is it empty?
        print('File %s empty'% fastafile)
        quit()

    print('\t\t...%s' % fastafile)

##Step1: Parse the input

#creating a dictionary containing ids as keys and seqs as values
    D={}
    with open(fastafile, 'r') as handle:
        parse = SeqIO.parse(handle, 'fasta')
        if any(parse): #But... Is the file .fasta?
            for seq_record in parse:
                    D[seq_record.id] = seq_record.seq
        else: #If not kill the program!
            print('The input file is not .fasta\n')
            quit()

    
    print('Retrieving stats...')


#Printing stats in terminal
    print('Seq-N is %i' % count_seq(D))
    print('Total size is %int' % count_nt(D))
    print('min is %ibp' % min_seq(D))
    print('max is %ibp' % max_seq(D))
    print('avg is %fbp' % avg_seq(D))
    for a in args.minsize:
        print('min size %ibp: Xseq = %ibp, setpercentage = %.2fbp' % min_size_Xbp(a,D))
    print('L50 is %i, N50 is %i' % NXoverLX(50, D))
    print('L90 is %i, N90 is %i' % NXoverLX(90, D))


#Updating dataframe with a new column
    row_list_fastafile = [count_seq(D),
                          count_nt(D),
                          min_seq(D),
                          max_seq(D),
                          avg_seq(D),
                          NXoverLX(50, D)[0],
                          NXoverLX(50, D)[1],
                          NXoverLX(90, D)[0],
                          NXoverLX(90, D)[1]
                          ]
    q = 5
    for a in args.minsize:
        row_list_fastafile.insert(q,min_size_Xbp(a, D)[1]) 
        q+=1
        row_list_fastafile.insert(q,min_size_Xbp(a, D)[2])
        q+=1

    df[fastafile[2::]] = row_list_fastafile

    if args.gapstats == True: #-g flag to retrieve gap stats as well
        df = N_GAPS.GAPS(D, df, args.DEFVAL,fastafile[2::])
    
    df.set_index('Stat',inplace=True)

#Sending df to csv file
df.to_csv ('./{}.csv'.format(args.name), index = True, float_format='%.16g', header=True, sep = output)


#checking if the output is ok...
if path.exists('./{}.csv'.format(args.name)):  
    print('\nstats file succesfully created')
elif not path.exists('./{}.csv'.format(args.name)) or stat('./{}.csv'.format(args.name)).st_size == 0:
    print("\nstats file not created")
    raise Exception ('Output problem: file either empty or not existing')
