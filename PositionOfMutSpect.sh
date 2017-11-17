#!/bin/bash
set -e
myVcf=$1
name=$2
chrom=$3
outputDIR=$4

#inspired by https://www.biostars.org/p/129502/

#convert vcf to bed
gzcat $myVcf | vcf2bed | cut -f1-7 - \
> $outputDIR/$name+chr$chrom+TestVariants.bed
#extend bed range to +1 on either end
bedops --everything \
--range 1 $outputDIR/$name+chr$chrom+TestVariants.bed \
> $outputDIR/$name+chr$chrom+Context_TestVariants.bed
#get sequences of those positions
bedtools getfasta \
-fi /Users/luke/genomes/fastaRef/chr$chrom.fa \
-bed $outputDIR/$name+chr$chrom+Context_TestVariants.bed \
> $outputDIR/$name+chr$chrom+Context_TestVariants.fasta
#get each position
grep -v '^>' $outputDIR/$name+chr$chrom+Context_TestVariants.fasta \
> $outputDIR/$name+chr$chrom+Context_TestVariants.nt
#put it all together
awk '{print $1,$2,$6,$7}' $outputDIR/$name+chr$chrom+TestVariants.bed \
| paste -d" " - $outputDIR/$name+chr$chrom+Context_TestVariants.nt \
> $outputDIR/$name+chr$chrom+Context_TestVariants_position.txt

awk '{print $1,$2,$3,$6,$7}' $outputDIR/$name+chr$chrom+TestVariants.bed \
| paste -d"\t" - $outputDIR/$name+chr$chrom+Context_TestVariants.nt  \
| awk '{print $1,$2,$3,$4"#"$5"#"$6}' \
| awk -v OFS="\t" '$1=$1'  \
> $outputDIR/$name+chr$chrom+Context_TestVariants_position.bed

bedtools unionbedg \
-i $outputDIR/$name+chr$chrom+Context_TestVariants_position.bed \
/Users/luke/genomes/genomes/NAG/misc/NAG14231.callable.bed \
-filler N/A \
> $outputDIR/$name+chr$chrom+Context_TestVariants_position_ANNOTATED.bed

grep -v "N/A" $outputDIR/$name+chr$chrom+Context_TestVariants_position_ANNOTATED.bed \
| sed 's/#/'$'\t/g' | awk '{print toupper($0)}'\
> $outputDIR/$name+Context_FINAL_position_ANNOTATED.bed
