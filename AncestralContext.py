import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = 'get Ancestral Context')
parser.add_argument('-input', help = 'input file')
parser.add_argument('-chrom', help = 'chromosome')
parser.add_argument('-out', help= 'output directory')
args = parser.parse_args()
#main(args)

#read input file
vcfPos=pd.read_table(args.input,
                     names=["Chrom","Pos", "Ref", "Alt"])
#only keep SNPs
DNA=['A','C','G','T']
vcfPos= vcfPos[vcfPos['Ref'].isin(DNA)]
vcfPos= vcfPos[vcfPos['Alt'].isin(DNA)]
#set Derived Allele to 1
vcfPos['DerivedAllele'] = 1

#read chimp human alignment
chimp=pd.read_table('/Users/luke/bin/smaller_mut_spectrum_pipeline/\
hg19_chimp_align/human_chimp_diffs_chr' \
+ args.chrom + '.txt', sep=" ")
chimp= chimp[chimp['SNP/Indel'] == 'SNP']

#read oneLine fasta ref file
ref=open('/Users/luke/bin/smaller_mut_spectrum_pipeline/\
hg19_reference/chr' + args.chrom + '_oneline.txt')
refseq=ref.read()

#get context
vcfPos['Context'] = vcfPos.apply(\
lambda row : refseq[row.Pos - 2 : row.Pos + 1], axis=1)

#get sites shared by both files
iPos = set(chimp['Pos']) & set(vcfPos['Pos'])
iChimp = chimp[chimp['Pos'].isin(list(iPos))]
iVCF = vcfPos[vcfPos['Pos'].isin(list(iPos))]
#find sites where Chimp is same as Human Alt
ChimpDiff = iVCF[iVCF['Alt'] == list(iChimp['Chimp'])]['Pos']

#For sites that have ChimpDiff, change DerivedAllele and Context
vcfPos.loc[vcfPos['Pos'].isin(ChimpDiff), 'DerivedAllele'] = 0
vcfPos.loc[vcfPos['Pos'].isin(ChimpDiff), 'Context'] = \
vcfPos['Context'].str[0] \
+ vcfPos['Alt'].str[0] \
+ vcfPos['Context'].str[2]

#write table
vcfPos.to_csv(args.out + 'Chr' + args.chrom + \
'AncestralContext.txt', index = False, sep="\t")
