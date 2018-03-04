import pandas as pd
import argparse
from itertools import cycle
import sys
import pickle
from collections import Counter

def main(args):
    Genotype = []
    AlleleFreq = Counter()
    #Genotype : #CHROM POS HG00096 HG00097 ...
    Chunk = pd.read_table(args.gt, compression = 'gzip', sep='\t', chunksize=10000)
    days = pd.read_table('~/bin/smaller_mut_spectrum_pipeline/Name.Pop.Day.txt', sep=' ',
                         names=['Name','Pop_Day']).set_index('Name').to_dict()
    print("Reading Genotypes for Chr : "+str(args.chrom))
    for GT in Chunk:
        #print("## GenotypeTable :")
        #print(GT.iloc[0:3,0:7])
        name = GT.iloc[:,6:].columns
        #groupby MutType then count number of muts per category per individual
        gb = GT.groupby(['CONTEXT','ALT']).agg(dict(zip(name,cycle(['value_counts'])))).fillna(0)
        Genotype.append(gb)
        
        #Allele Frequency
        GT.rename(columns=days['Pop_Day'],inplace=True)
        af = GT.iloc[:,6:].transpose().groupby(level=0).sum()
        tuples = [tuple(x) for x in af.transpose().values]
        AlleleFreq = AlleleFreq + Counter(map(tuple, tuples))
    
    #Allele Frequency
    print("AlleleFreq for Chr : "+str(args.chrom))
    AFname = args.out + 'Chr' + args.chrom + '_AlleleFrequencyPerPopDay.txt'
    with open(AFname, 'wb') as outputfile:
        pickle.dump(AlleleFreq, outputfile)
    
    print("Formatting table for Chr : "+str(args.chrom))
    #from each chunk, add each category together
    Genotype=pd.concat(Genotype).sort_index()
    Genotype = Genotype.groupby(level=[0,1,2]).agg(dict(zip(name,cycle(['sum']))))    
    #do some formatting, flatten multi index
    Genotype=Genotype.reset_index()
    Genotype['Mut'] = Genotype[['level_0','level_1']].apply(lambda x: '_'.join(x),axis=1)
    cols = Genotype.columns.tolist()
    ncols = [cols[-1]] + cols[:-1]
    Genotype = Genotype[ncols]    

    
    print("Done AF and GT format for Chr : "+str(args.chrom))
    #read namePop file
    NamePop = pd.read_csv('~/bin/smaller_mut_spectrum_pipeline/Name.BigPop.Pop.txt', 
                          names=['Name','BigPop','Pop'], sep=' ')
    NamePop=NamePop[NamePop.Pop != 'NAG'] #get rid of NAGs
    
    print("Writing output for each pop to tables for Chr : "+str(args.chrom))
    #make a file per BigPop
    BigPops = list(NamePop['BigPop'].unique())
    for pop in BigPops:
        names=list(['level_0', 'level_1', 'level_2', 'Mut'])
        names.extend(NamePop[NamePop['BigPop'] == pop]['Name'])
        sub=Genotype[names]
        HomRef=sub[sub.level_2 == 0].drop(['level_0','level_1','level_2'], axis=1)
        Het=sub[sub.level_2 == 1].drop(['level_0','level_1','level_2'], axis=1)
        HomAlt=sub[sub.level_2 == 2].drop(['level_0','level_1','level_2'], axis=1)    

        #write table
        fileName1= args.out + 'Chr' + args.chrom + '_' + pop + '_HomAlt.MutTypeCountPerSample.txt'
        fileName2= args.out + 'Chr' + args.chrom + '_' + pop + '_Hetero.MutTypeCountPerSample.txt'
        fileName3= args.out + 'Chr' + args.chrom + '_' + pop + '_HomRef.MutTypeCountPerSample.txt'

        print("OutFile : " + args.out + 'Chr' + args.chrom + '_' + pop)

        HomRef.to_csv(fileName1, sep="\t", index=False)
        Het.to_csv(fileName2, sep="\t", index=False)
        HomAlt.to_csv(fileName3, sep="\t", index=False)
        
        SmallPops = NamePop[NamePop['BigPop'] == pop]['Pop']
        for spop in SmallPops:
            names=list(['level_0', 'level_1', 'level_2', 'Mut'])
            names.extend(NamePop[NamePop['Pop'] == spop]['Name'])
            single=sub[names]
            HomRef=single[single.level_2 == 0].drop(['level_0','level_1','level_2'], axis=1)
            Het=single[single.level_2 == 1].drop(['level_0','level_1','level_2'], axis=1)
            HomAlt=single[single.level_2 == 2].drop(['level_0','level_1','level_2'], axis=1)    

            #write table
            fileName1= args.out+'Chr'+args.chrom+'_'+pop+'_'+spop+'_HomAlt.MutTypeCountPerSample.txt'
            fileName2= args.out+'Chr'+args.chrom+'_'+pop+'_'+spop+'_Hetero.MutTypeCountPerSample.txt'
            fileName3= args.out+'Chr'+args.chrom+'_'+pop+'_'+spop+'_HomRef.MutTypeCountPerSample.txt'

            HomRef.to_csv(fileName1, sep="\t", index=False)
            Het.to_csv(fileName2, sep="\t", index=False)
            HomAlt.to_csv(fileName3, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Count Genotype per MutSpect')
    parser.add_argument('-gt', help = 'input genotype')
    parser.add_argument('-chrom', help = 'chromosome')
    parser.add_argument('-out', help= 'output directory')
    args = parser.parse_args()
    main(args)
