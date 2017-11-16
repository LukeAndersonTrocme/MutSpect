#!/bin/bash
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

# ls /Users/luke/Documents/MutSpect/MutSpectPosition/BED/ | gshuf -n 5
# vcftools \
# --gzvcf $inputVCF \
# --keep /Users/luke/Documents/MutSpect/data/1000genomes_phase3_sample_IDs_JUST_JPT.txt \
# --freq2 \
# --out $outputDIR/1000Genome_filtered_JUST_JPT_freq
#
# vcftools \
# --gzvcf $inputVCF \
# --keep /Users/luke/Documents/MutSpect/data/1000genomes_phase3_sample_IDs_JUST_NAG.txt \
# --freq2  \
# --out $outputDIR/NAGJapan_filtered_freq

# for f in /Users/luke/Documents/MutSpect/MutSpectPosition/BED/*.bed
# do /usr/local/bin/bcftools-1.6/bcftools filter \
#  --regions-file $f \
# /Users/luke/Documents/MutSpect/FilteredData/2017-10-30_1kGenome_NAG_filtered_chr22.vcf.gz \
# > /Users/luke/Documents/MutSpect/FilteredData/2017-10-30_1kGenome_NAG_filtered_chr22_CACtoCCC.vcf

# Rscript FreqSpectPlot.R \
# $outputDIR/1000Genome_filtered_JUST_JPT_freq.frq \
# $outputDIR/NAGJapan_filtered_freq.frq \
# $outputDIR/plots/$name
