#!/bin/bash
##set python environment to python2
source activate python2
set -e
inputVCF=$1
chrom=$2
name=$3
outputDirectory=$4


repos="/Users/luke/bin/smaller_mut_spectrum_pipeline/"
mkdir $outputDirectory/$name.MutSpect
mkdir $outputDirectory/$name.MutSpect/plots
mkdir $outputDirectory/$name.MutSpect/files

outputDIR="$outputDirectory/$name.MutSpect"
##change to the respository
cd $repos
##copy all the files to the directory just to keep track of any changes to files.
cp get_finescale_mut_spectra_pep8.py $outputDirectory/scripts/get_finescale_mut_spectra_pep8_$name.py
cp make_heatmap_ASN_filtered.py $outputDirectory/scripts/make_heatmap_ASN_filtered_$name.py
cp Make_PCA.R $outputDirectory/scripts/Make_PCA_$name.R
cp FreqSpectPlot.R $outputDirectory/scripts/FreqSpectPlot_$name.R
cp MutSpect_PipeLine.sh $outputDirectory/scripts/MutSpect_PipeLine_$name.sh


 echo "########## FIRST STEP : \
 python get_finescale_mut_spectra_pep8.py \
 -chrom $chrom \
 -vcf $inputVCF \
 -id $repos/1000genomes_phase3_sample_IDs_NAG_SGDP.txt \
 -repos $repos \
 -out $outputDirectory/$name.MutSpect/files
 "
##convert VCF to popMutSpect and indMutSpect
python get_finescale_mut_spectra_pep8.py \
-chrom $chrom \
-vcf $inputVCF \
-id $repos/1000genomes_phase3_sample_IDs_NAG_SGDP.txt \
-repos $repos \
-out $outputDirectory/$name.MutSpect/files/

echo "########## SECOND STEP : \
python make_heatmap_BIGPOP_filtered.py \
-i $outputDIR/files/ \
-out $outputDIR/plots/$name"
##make heatmaps
python make_heatmap_BIGPOP_filtered.py \
-i $outputDIR/files/ \
-out $outputDIR/plots/$name \
-chrom $chrom

echo "########## SECOND STEP : \
python make_heatmap_ALLPOP_filtered.py \
-i $outputDIR/files/ \
-out $outputDIR/plots/$name"
##make heatmaps
python make_heatmap_ALLPOP_filtered.py \
-i $outputDIR/files/ \
-out $outputDIR/plots/$name \
-chrom $chrom

 ##make R plots
echo "########## THIRD STEP : \
Rscript Make_PCA.R \
$outputDIR/files/ \
$chrom \
$repos \
$outputDIR/plots/$name"

Rscript Make_PCA.R \
$outputDIR/files/ \
$chrom \
$repos \
$outputDIR/plots/$name
