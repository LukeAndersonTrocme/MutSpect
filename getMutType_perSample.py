import pandas as pd
import argparse
from itertools import cycle
import sys

def main(args):
    #context : Chrom Pos Ref Alt Derived Allele Context
    context=pd.read_table(args.context, sep='\t')
    #Genotype : #CHROM POS HG00096 HG00097 ...
    GT=pd.read_table(args.gt, compression = 'gzip', sep='\t')
    print("## GenotypeTable :")
    print(GT.iloc[0:2,0:5])

    #get sites shared by both files
    iPos = set(GT['POS']) & set(context['Pos'])
    #subset each table to get just shared sites
    iGT = GT[GT['POS'].isin(list(iPos))]
    icontext = context[context['Pos'].isin(list(iPos))]
    #join the two tables together based on position
    joined=icontext.merge(iGT, left_on='Pos', right_on='POS')
    print("## Tables Merged")
    print(joined.iloc[0:2,0:6])
    print('sys.getsizeof(joined) : '+str(sys.getsizeof(joined)))
    
    v=(joined['DerivedAllele'] == 0)
    joined.loc[v,joined.columns[8:]]=abs(joined.loc[v,joined.columns[8:]] - 2)

    #Norm=joined[joined['DerivedAllele'] == 1] #get HumRef
    #print('sys.getsizeof(Norm) : '+str(sys.getsizeof(Norm)))
    #Der=joined[joined['DerivedAllele'] == 0] #get ChimpRef
    #print('sys.getsizeof(Der) : '+str(sys.getsizeof(Der)))
    #incorrect_to_correct = {0:2, 2:0} #swap HomRef with HomAlt
    #D=Der.iloc[:,8:].replace(incorrect_to_correct) #swap
    #Fixed=pd.concat([Der.iloc[:,0:8], D], axis=1) #column bind
    #FixedMerged=pd.concat([Fixed,Norm], axis=0).sort_values(by=['Pos']) #rowbind
    print("## Fixed Merged")
    ##grouped by the Alt -> Context
    ##each table of MutType is then summed per individual
    name = FixedMerged.iloc[:,8:].columns
    #gb = FixedMerged.groupby(['Alt','Context']).agg(dict(zip(name,cycle(['value_counts'])))).fillna(0)
    gb = joined.groupby(['Alt','Context'])
    gb = gb.agg(dict(zip(name,cycle(['value_counts']))))
    gb = gb.fillna(0)

    print(gb.iloc[0:5,0:6])
    #write table
    fileName= args.out + 'Chr' + args.chrom + 'MutTypeCountPerSample.txt'
    print("FileName : " + fileName)
    gb.to_csv(fileName, index = False, sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Count Genotype per MutSpect')
    parser.add_argument('-context', help = 'input context')
    parser.add_argument('-gt', help = 'input genotype')
    parser.add_argument('-chrom', help = 'chromosome')
    parser.add_argument('-out', help= 'output directory')
    args = parser.parse_args()
    main(args)
