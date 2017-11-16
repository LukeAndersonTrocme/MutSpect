#!/bin/bash

##Merged three bed files using bedtools unionbedg
##see evernote for documentation
TimeStamp=$(date +%Y-%m-%d)
inputVCF=$1
name=$2
outputDIR=$3

pathToBed="/Users/luke/genomes/BED_MASKS"
mkdir /Users/luke/Documents/MutSpect/FilteredData/$TimeStamp
outputDIR="/Users/luke/Documents/MutSpect/FilteredData/$TimeStamp"

cp /Users/luke/bin/smaller_mut_spectrum_pipeline/FiltrationPipeLine2.0.sh \
$outputDIR/FiltrationPipeLine2.0_$name.sh

#this is to time the whole thing
start=`date +%s`
#filter the input VCF
# 3 bed masks : nestedRepeats, phastCons100way and strict_mask
# Remove : low Qual, missing > 1%, max AF 98%
# Keep only biallelic snps
echo "########## Filtering $name "
echo '/usr/local/bin/bcftools-1.6/bcftools filter \
--regions-file $pathToBed/Strict.Cons100way.Repeats.bed \
$inputVCF \
| /usr/local/bin/bcftools-1.6/bcftools filter \
-e "QUAL < 10 || F_MISSING > 0.01 || MAF > 0.98" \
| /usr/local/bin/bcftools-1.6/bcftools view \
-m2 -M2 -v snps --output-type z \
--output-file \
$outputDIR/$name.3bed_filtered.vcf.gz \'

/usr/local/bin/bcftools-1.6/bcftools filter \
--regions-file $pathToBed/Strict.Cons100way.Repeats.bed \
$inputVCF \
| /usr/local/bin/bcftools-1.6/bcftools filter \
-e "QUAL < 10 || F_MISSING > 0.01 || MAF > 0.98" \
| /usr/local/bin/bcftools-1.6/bcftools view \
-m2 -M2 -v snps --output-type z \
--output-file \
$outputDIR/$name.3bed_filtered.vcf.gz

tabix -f -p vcf $outputDIR/$name.3bed_filtered.vcf.gz

# end=`date +%s`
# runtime=$((end-start))
# echo "Time (seconds) to run FiltrationPipeLine2.0 on is : "
# echo $runtime
