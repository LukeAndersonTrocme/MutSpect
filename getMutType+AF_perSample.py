import pandas as pd
import argparse
from itertools import cycle
import sys
import pickle
from collections import Counter
import os

def main(args):
    Genotype = []
    AlleleFreq = Counter()
    names=pd.read_table('/Users/luke/bin/smaller_mut_spectrum_pipeline/Name.Pop.Day.txt', 
                        sep=' ', names=['Name','Pop_Day'])
    Pops=list(names['Pop_Day'].unique())

    days = pd.read_table('~/bin/smaller_mut_spectrum_pipeline/Name.Pop.Day.txt', sep=' ',
                             names=['Name','Pop_Day']).set_index('Name').to_dict()
    #list_of_files = glob.glob(args.i + 'Chr2?.FlippedGenotypes.Context_10000.txt.gz')
    if args.only is None:
        print('####Get MutSpect and AlleleFrequency')
        getGT=True
        getAF=True
    if args.only == 'AF':
        print('####Get ONLY AlleleFrequency')
        getGT=False
        getAF=True        
    if args.only == 'GT':
        print('####Get ONLY MutSpect')
        getGT=True
        getAF=False
        
    #for file_name in list_of_files:
    file_name = args.i + 'Chr' + args.chrom + '.FlippedGenotypes.Context.txt.gz'
    Chunk = pd.read_table(file_name, compression = 'gzip', sep='\t', chunksize=10000)
    print("##Reading Genotypes for Chr : "+ args.chrom, end=" ", flush=True)
    for GT in Chunk:
        if getGT:
            name = GT.iloc[:,6:].columns
            #groupby MutType then count number of muts per category per individual
            gb = GT.groupby(['CONTEXT','ALT']).agg(dict(zip(
                name,cycle(['value_counts'])))).fillna(0)
            Genotype.append(gb)

        if getAF:
            GT.rename(columns=days['Pop_Day'],inplace=True) #rename to groupby
            af = GT.iloc[:,6:].transpose().groupby(level=0).sum() #get AF per group
            tuples = [tuple(x) for x in af.transpose().values] #turn it into tuple
            AlleleFreq = AlleleFreq + Counter(map(tuple, tuples))
    print("[DONE]")       
    if getAF:        
        print("##AF Processing for Chr : "+ args.chrom)    
        D=pd.DataFrame.from_dict(AlleleFreq,orient='index').reset_index() #dict to DF
        D[Pops] = D['index'].apply(pd.Series) #split tuple to columns
        D=D.drop(['index'], axis=1) #drop tuple
        D.columns.values[0] = 'Count' #rename
        for i in range(0,52,2):
            print(Pops[i][:3], end=" ", flush=True)
            Sub=D.loc[:,['Count',Pops[i], Pops[i+1]]] #subset per pop
            Sub=Sub.groupby([Pops[i], Pops[i+1]]).sum().reset_index() #get count of AF
            fname=args.out+'AlleleFrequency_'+ 'Chr' + args.chrom +Pops[i][:3]+'_HeatMap.txt'
            # if file does not exist write header 
            if not os.path.isfile(fname):
                Sub.to_csv(fname, mode = 'w')
            else: # else it exists so append without writing the header
                Sub.to_csv(fname, mode = 'a',header=False)
        print("[DONE]")             
    if getGT:
        print("##Formatting table for Chr : "+ args.chrom, end=" ", flush=True)
        #from each chunk, add each category together
        Genotype=pd.concat(Genotype).sort_index()
        Genotype = Genotype.groupby(level=[0,1,2]).agg(dict(zip(name,cycle(['sum']))))    
        #do some formatting, flatten multi index
        Genotype=Genotype.reset_index()
        Genotype['Mut'] = Genotype[['level_0','level_1']].apply(lambda x: '_'.join(x),axis=1)
        cols = Genotype.columns.tolist()
        ncols = [cols[-1]] + cols[:-1]
        Genotype = Genotype[ncols]    

        print("[DONE]")
        #read namePop file
        NamePop = pd.read_csv('~/bin/smaller_mut_spectrum_pipeline/Name.BigPop.Pop.txt', 
                              names=['Name','BigPop','Pop'], sep=' ')
        NamePop=NamePop[NamePop.Pop != 'NAG'] #get rid of NAGs

        print("##Writing output for each pop to tables")
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
            fileName1= args.out + 'Chr' + args.chrom + pop + '_HomAlt.MutTypeCountPerSample.txt'
            fileName2= args.out + 'Chr' + args.chrom + pop + '_Hetero.MutTypeCountPerSample.txt'
            fileName3= args.out + 'Chr' + args.chrom + pop + '_HomRef.MutTypeCountPerSample.txt'

            print(pop, end=" ", flush=True)

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
                fileName1= args.out+ 'Chr' + args.chrom +pop+'_'+spop+'_HomAlt.MutTypeCountPerSample.txt'
                fileName2= args.out+ 'Chr' + args.chrom +pop+'_'+spop+'_Hetero .MutTypeCountPerSample.txt'
                fileName3= args.out+ 'Chr' + args.chrom +pop+'_'+spop+'_HomRef.MutTypeCountPerSample.txt'

                HomRef.to_csv(fileName1, sep="\t", index=False)
                Het.to_csv(fileName2, sep="\t", index=False)
                HomAlt.to_csv(fileName3, sep="\t", index=False)
        print("[DONE]")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Count Genotype per MutSpect')
    parser.add_argument('-i', help = 'input Directory')
    parser.add_argument('-chrom', help = 'Chrom')
    parser.add_argument('-only', help = 'AF or GT')
    parser.add_argument('-out', help= 'output directory')
    args = parser.parse_args()
    main(args)
