import pandas as pd
import argparse
from itertools import cycle
from collections import Counter
import glob
import sys
import os
import shutil

out='/Users/luke/genomes/genomes/hg19/AncestralRef/RunAll/'
names=pd.read_table('/Users/luke/bin/smaller_mut_spectrum_pipeline/Name.BigPop.Pop.txt', 
                    sep=' ', names=['Name','BigPop','Pop'])
names = names[names.Pop != 'NAG']
Pops=list(names['Pop'].unique())
BigPops=list(names['BigPop'].unique())

print("AlleleFrequency")
for pop in Pops:
    list_of_files = glob.glob(out+'AlleleFrequency_Chr*'+pop+'_HeatMap.txt')
    AF=pd.DataFrame()
    for file_name in list_of_files:
        AF=AF.append(pd.read_table(file_name, sep=',', index_col=0), ignore_index=True)
    AF=AF.groupby([list(AF)[0], list(AF)[1]]).sum().reset_index()
    fname=out+'AlleleFrequency_GenomeWide_'+pop+'_HeatMap.txt'
    AF.to_csv(fname)
print("Using a total of "+str(len(list_of_files))+" Chromosomes")    
print("Genotypes")    
HomHet=['HomRef','Hetero','HomAlt']
for h in range(len(HomHet)):
    allele=HomHet[h]
    for BP in BigPops:
        SmallPops = names[names['BigPop'] == BP]['Pop'].unique()
        for spop in SmallPops:
            list_of_files = glob.glob(
                out+'Chr*_'+BP+'_'+spop+'_'+allele+'.MutTypeCountPerSample.txt')
            Genotype=pd.DataFrame()
            for file_name in list_of_files:
                Genotype=Genotype.append(pd.read_table(file_name, sep='\t'))
            Genotype = Genotype.groupby('Mut').sum().reset_index()
            fname=out+'MutationSpectrum_GenomeWide_'+spop+'.txt'
            Genotype.to_csv(fname, sep="\t")
print("Using a total of "+str(len(list_of_files))+" Chromosomes") 