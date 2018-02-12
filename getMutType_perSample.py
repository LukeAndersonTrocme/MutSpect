import pandas as pd
import argparse

#total count of mutations per individual for each type
def getCounts(table):
    n=table.iloc[:,8:].apply(pd.value_counts)
    n.insert(0, 'Type', str(table.iloc[0,3])+"->"+str(table.iloc[0,5]))
    return n

parser = argparse.ArgumentParser(description = 'Count Genotype per MutSpect')
parser.add_argument('-context', help = 'input context')
parser.add_argument('-gt', help = 'input genotype')
parser.add_argument('-chrom', help = 'chromosome')
parser.add_argument('-out', help= 'output directory')
args = parser.parse_args()

#context : Chrom Pos Ref Alt Derived Allele Context
context=pd.read_table(args.context, sep='\t')
#Genotype : #CHROM POS HG00096 HG00097 ...
GT=pd.read_table(args.gt, compression = 'gzip', sep='\t')
print("## GenotypeTable :")
print(GT.iloc[0,0:5])
#get sites shared by both files
iPos = set(GT['POS']) & set(context['Pos'])
#subset each table to get just shared sites
iGT = GT[GT['POS'].isin(list(iPos))]
icontext = context[context['Pos'].isin(list(iPos))]
#join the two tables together based on position
joined=icontext.merge(iGT, left_on='Pos', right_on='POS', how='inner')
print("## Tables Merged")
print(joined.iloc[0,0:6])

##fancy GroupBy :
##this splits the long table into a list of tables,
##each table is grouped by the Alt -> Context
gb = joined.groupby(['Alt','Context'])
t = [gb.get_group(x) for x in gb.groups]

print(t[1].iloc[0,0:6])
print("## Tables Split by MutType")

#concatenate each of these back into a single table
#each row is a count of mutation type per individual
f=pd.concat([getCounts(t[x]) for x in range(1,len(t))])
print("## GetCounts Complete")
#write table
fileName= args.out + 'Chr' + args.chrom + 'MutTypeCountPerSample.txt'
print("FileName : " + fileName)
f.to_csv(fileName, index = False, sep="\t")
