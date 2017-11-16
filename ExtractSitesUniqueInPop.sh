#!/bin/bash


for file in /Users/luke/Documents/MutSpect/data/SitesIn_NAG_Tohoku.txt \
/Users/luke/Documents/MutSpect/data/SitesIn_NAG.unique.txt \
/Users/luke/Documents/MutSpect/data/SitesIn_ThouGenome_NAG.txt \
/Users/luke/Documents/MutSpect/data/SitesIn_ThouGenome.unique.txt \
/Users/luke/Documents/MutSpect/data/SitesIn_Tohoku_ThouGenome.txt \
/Users/luke/Documents/MutSpect/data/SitesIn_Tohoku.unique.txt
do echo $file
vcftools \
--vcf $pathToData/1000Genome+Japan_merged_filtered.vcf \
--positions $file \
--recode --recode-INFO-all \
--out  $pathToData/1000Genome+Japan$(basename "$file" .txt).vcf
 done
