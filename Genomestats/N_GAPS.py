#Stats about gaps 'N' nucleotides
#To be imported by genomestats.py

#libraries
from Bio.Seq import Seq
import pandas as pd


def GAPS(d, dafr, defval,ff):
        
    totalgapsize = 0
    Lgapsizes=[]
    numberd=0
    for v in list(d.values()):
        numberd+=1 #which seq is being analyzed
        totalN= v.count('N') #total Ns in a seq
        stpoint = v.find('N') # starting point for assessing gap length
        N1=0
        while N1 != -1:
            N1 = v.find('N',start=stpoint) #first N found starting from stpoint, if no Ns are found -1 is returned
            posL=[]
            for nt in ['A','C','G','T', 'U', 'a', 'c', 'g', 't', 'u']:
                pos = v.find(nt, start=N1) #from the first N found, when can we find a nt
                if pos != -1:
                    posL.append(pos)
            gap_length = min(posL)- N1 #the length of the retrieved gap is the position of the first nt found after the first found N
            stpoint = min(posL)
            if gap_length !=0: #and if the gap is different from 0, append the size to a list
                Lgapsizes.append(gap_length)
        totalgapsize += totalN #increment this value of the total Ns found in a given seq

    #Gap length default values
    if defval == False:
        SLgapsizes = set(Lgapsizes)
        Lgapsizes_unique = list(SLgapsizes)
        for lgs in Lgapsizes_unique:
            percent_lgs = Lgapsizes.count(lgs)/len(Lgapsizes)
            if percent_lgs > 0.2:
                list(filter(lambda a: a == lgs, Lgapsizes))
                import sys
                sys.stderr.write('Gaps of default gap size of %i were removed to calculate the stats' % lgs)
        

    #Defining stats
    totLgapsize = sum(Lgapsizes)
    lenLgapsizes = len(Lgapsizes)
    meanLgapsize = round(totalgapsize/lenLgapsizes,2)
    mingap = min(Lgapsizes)
    maxgap = max(Lgapsizes)
    gaps_majorthan100 = len ([x for x in Lgapsizes if x>100])
    gaps_majorthan1000 = len ([x for x in Lgapsizes if x>1000])

    #Putting values in a list
    list_gapstats = [totLgapsize,lenLgapsizes,mingap,maxgap,meanLgapsize,gaps_majorthan100,gaps_majorthan1000]

    #appending gapstats to genomestats df
    df_gap = pd.DataFrame({'Stat':['total gapsize','number of gaps','min gap size', 'max gapsize', 'mean gap size', 'gaps >100nt', 'gaps <1000nt'], ff:list_gapstats})
    result = pd.concat([dafr,df_gap]) 
    print('\nGAPSTATS: total gap size is %i, total n of gaps is %i, min gap size is %i, max gap size is %i, mean gap size is %.2f, n of gaps major than 100 nt is %i, n of gaps major than 1000 nt is %i'%(totalgapsize, lenLgapsizes, mingap, maxgap, meanLgapsize, gaps_majorthan100, gaps_majorthan1000))
    return result
