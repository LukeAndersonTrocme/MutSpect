#!/bin/bash
myVcf=$1
name=$2

PathToData='/Users/luke/Documents/MutSpect'

gzcat $myVcf | vcf2bed | cut -f1-7 - \
> $PathToData/$name+TestVariants.bed

bedops --everything \
--range 1 $PathToData/$name+TestVariants.bed \
> $PathToData/$name+Context_TestVariants.bed

bedtools getfasta \
-fi $PathToData/chr22.fa \
-bed $PathToData/$name+Context_TestVariants.bed \
> $PathToData/$name+Context_TestVariants.fasta

grep -v '^>' $PathToData/$name+Context_TestVariants.fasta \
> $PathToData/$name+Context_TestVariants.nt

awk '{print $1,$2,$6,$7}' $PathToData/$name+TestVariants.bed \
| paste -d" " - $PathToData/$name+Context_TestVariants.nt \
> $PathToData/$name+Context_TestVariants_position.txt
