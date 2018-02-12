#!/bin/bash
myVcf=$1
chrom=$2
outputDIR=$3

cd /Users/luke/bin/snpEff_latest_core/snpEff

# gzcat $myVcf | head -500 \
# | grep "^#CHROM" \
# | cut -f1,2,10- \
# > $outputDIR/header.txt

java -jar SnpSift.jar extractFields \
$myVcf CHROM POS "GEN[*].GT" \
| sed 's/\//|/g' \
| sed 's/0|0/0/g;s/0|1/1/g;s/1|0/1/g;s/1|1/2/g' \
| sed '/|/d' \
| sed "1s/.*/$(head -n 1 $outputDIR/header.txt )/" \
| gzip -f > $outputDIR/Chr$chrom.Filtered_Genotypes.txt.gz
