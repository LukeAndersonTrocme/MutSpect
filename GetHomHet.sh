#!/bin/bash
myVcf=$1
chrom=$2
outputDIR=$3

cd /Users/luke/bin/snpEff_latest_core/snpEff
if [ -f $outputDIR/Chr$chrom.Genotypes.txt ]; then
  echo "File $outputDIR/Chr$chrom.Genotypes.txt exists, no need to get GT"
else
  java -jar SnpSift.jar extractFields \
  $myVcf POS "GEN[*].GT" \
  | sed 's/\//|/g ; s/0|0/0/g ; s/0|1/1/g ; s/1|0/1/g ; s/1|1/2/g ; /|/d' \
  | tail -n +2 > $outputDIR/Chr$chrom.Genotypes.txt
fi

#check that the Genotype file doesn't have a header
line=$(head -n 1 $outputDIR/Chr$chrom.Genotypes.txt)
if [[ $line = *"POS"* ]]; \
then echo "Have to remove fist line of Genotypes"
tail -n +2 $outputDIR/Chr$chrom.Genotypes.txt > \
tail -n +2 $outputDIR/Temp.Genotypes.txt && \
mv $outputDIR/Temp.Genotypes.txt $outputDIR/Chr$chrom.Genotypes.txt
fi

echo "#### Attempting Join"
join -1 2 -2 1 -t $'\t' \
$outputDIR/Chr$chrom.AncestralContext.txt \
$outputDIR/Chr$chrom.Genotypes.txt \
| awk 'BEGIN{ OFS = "\t" } { if ($5 == 0) {for (i = 7; i <= NF; i++) $i == $i - 2} print $0}' \
| sed 's/-//g' \
| sed "1s/.*/$(head -n 1 $outputDIR/header.txt )/" \
| gzip > $outputDIR/Chr$chrom.FlippedGenotypes.Context.txt.gz

gzcat $outputDIR/Chr$chrom.FlippedGenotypes.Context.txt.gz | wc -l

echo "#### Attempted Join"
#remove this file because it's big
#rm $outputDIR/Chr$chrom.Genotypes.txt
