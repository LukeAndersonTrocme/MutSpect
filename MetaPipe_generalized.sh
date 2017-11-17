#!/bin/bash
set -e
echo "\
                    __ __ __ __ __
                   /__/__/__/__/__/|
                  /__/__/__/__/__/|/
                  |__'__'__'__'__|/
"
#check if I haven't saved on git
# if git diff-index --quiet HEAD --; then
#     echo 'no changes'
# else
#     echo 'changes'
# fi
#test for emptiness
chrom=$1
TimeStamp=$(date +%Y-%m-%d) #TimeStamp used to make the output folder

cd /Users/luke/bin/smaller_mut_spectrum_pipeline

PathToGenome="/Users/luke/genomes/genomes" #this is where the data is
#PathToGenome="/Users/luke/genomes/genomes/test22" #this is where the data is
pathToBed="/Users/luke/genomes/BED_MASKS" #this is where the bed masks are
#check if folder exists
mkdir -p /Users/luke/Documents/MutSpect/FilteredData/$TimeStamp
#make a new directory for hour/minutes
outputDIR="/Users/luke/Documents/MutSpect/FilteredData/$TimeStamp/$chrom" #set outputDIR
mkdir -p $outputDIR
mkdir -p $outputDIR/scripts
#git rev-parse HEAD > $outputDIR/git_info.txt #(pointing to script directory)
#git checkout ^ this number to get back to that state

#copy all scripts for reproducibility
cp MetaPipe_generalized.sh $outputDIR/scripts/
cp FiltrationPipeLine2.0.sh $outputDIR/scripts/
cp MutSpect_PipeLine.sh $outputDIR/scripts/
cp PositionOfMutSpect.sh $outputDIR/scripts/
cp FreqSpectPlot.R $outputDIR/scripts/
cp PlotFilter.R $outputDIR/scripts/

#if [ -f $outputDIR/NAG_chr$chrom.3bed_filtered.vcf.gz.tbi ]; then
  echo "########### 1 : Nag Filtration ###########"
  bash FiltrationPipeLine2.0.sh \
  $PathToGenome/NAG/allSamples.$chrom.genotyped.vcf.gz \
  NAG_chr$chrom \
  $outputDIR
  # /usr/local/bin/bcftools-1.6/bcftools filter \
  # --regions-file $pathToBed/Strict.Cons100way.Repeats.bed \
  # $PathToGenome/NAG/allSamples.$chrom.genotyped.vcf.gz \
  # | /usr/local/bin/bcftools-1.6/bcftools filter \
  # -e "QUAL < 10 || F_MISSING > 0.01 || MAF > 0.98" \
  # | /usr/local/bin/bcftools-1.6/bcftools view \
  # -m2 -M2 -v snps --output-type v \
  # --output-file \
  # $outputDIR/NAG_chr$chrom.3bed_filtered.vcf \
  # && bgzip -f $outputDIR/NAG_chr$chrom.3bed_filtered.vcf \
  # && tabix -f -p vcf $outputDIR/NAG_chr$chrom.3bed_filtered.vcf.gz #\
  #& echo "  #   #   #    NAG_chr$chrom.3bed_filtered.vcf.gz Kill PID : $!"
#fi
echo "hey"

#if [ -f $outputDIR/1kG_chr$chrom.3bed_filtered.vcf.gz.tib ]; then
  echo "########### 2 : 1000Genome Filtration ###########"
  bash FiltrationPipeLine2.0.sh \
  $PathToGenome/hg19/phase3/ALL.chr$chrom.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
  1kG_chr$chrom \
  $outputDIR

  # /usr/local/bin/bcftools-1.6/bcftools filter \
  # --regions-file $pathToBed/Strict.Cons100way.Repeats.bed \
  # $PathToGenome/hg19/phase3/ALL.chr$chrom.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
  # | /usr/local/bin/bcftools-1.6/bcftools filter \
  # -e "QUAL < 10 || F_MISSING > 0.01 || MAF > 0.98" \
  # | /usr/local/bin/bcftools-1.6/bcftools view \
  # -m2 -M2 -v snps --output-type v \
  # --output-file \
  # $outputDIR/1kG_chr$chrom.3bed_filtered.vcf \
  # && bgzip -f $outputDIR/1kG_chr$chrom.3bed_filtered.vcf \
  # && tabix -f -p vcf $outputDIR/1kG_chr$chrom.3bed_filtered.vcf.gz #\
  #& echo "  #   #   #    1kG_chr$chrom.3bed_filtered.vcf.gz Kill PID : $!"
#fi

##wait until each step has completed
wait

fileTest="$outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf.gz"
if [ -f $fileTest ] && echo "$fileTest Found" || echo "$fileTest Not found"; then
  echo "########### 3 : Merge Nag and 1000Genome ###########"
  /usr/local/bin/bcftools-1.6/bcftools merge \
  -0 $outputDIR/1kG_chr$chrom.3bed_filtered.vcf.gz \
  $outputDIR/NAG_chr$chrom.3bed_filtered.vcf.gz \
  | /usr/local/bin/bcftools-1.6/bcftools view \
  --output-type z \
  > $outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf.gz

  /usr/local/bin/bcftools-1.6/bcftools index \
  -t --output-file $outputDIR/$name.3bed_filtered.vcf.gz.tbi \
  $outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf.gz
  # bgzip -f $outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf
  # tabix -f -p vcf $outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf.gz
fi

echo "########### MutSpect PipeLine ###########"
bash MutSpect_PipeLine.sh \
$outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf.gz \
$chrom \
$TimeStamp.chr$chrom \
$outputDIR \
& echo "  #   #   #  MutSpect_PipeLine  Kill PID : $!"

echo "########### Position Of MutSpect (unfiltered)###########"
bash PositionOfMutSpect.sh \
$outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf.gz \
unfiltered \
$chrom \
$outputDIR

fileTest="$outputDIR/1000Genome_filtered_JUST_JPT_freq.frq"
if [ -f $fileTest ] && echo "$fileTest Found" || echo "$fileTest Not found"; then
  echo "########### SITE FREQUENCY SPECTRUM PLOT ###########"
  vcftools \
  --gzvcf $inputVCF \
  --keep /Users/luke/Documents/MutSpect/data/1000genomes_phase3_sample_IDs_JUST_JPT.txt \
  --freq2 \
  --out $outputDIR/1000Genome_filtered_JUST_JPT_freq \
  & echo "  #   #   #  vcftools JPT  Kill PID : $!"

  vcftools \
  --gzvcf $inputVCF \
  --keep /Users/luke/Documents/MutSpect/data/1000genomes_phase3_sample_IDs_JUST_NAG.txt \
  --freq2  \
  --out $outputDIR/NAGJapan_filtered_freq \
  & echo "  #   #   #  vcftools NAG  Kill PID : $!"

  wait

  echo "Rscript FreqSpectPlot.R \
  $outputDIR/1000Genome_filtered_JUST_JPT_freq.frq \
  $outputDIR/NAGJapan_filtered_freq.frq \
  $outputDIR/plots/$name
  "
  Rscript FreqSpectPlot.R \
  $outputDIR/1000Genome_filtered_JUST_JPT_freq.frq \
  $outputDIR/NAGJapan_filtered_freq.frq \
  $outputDIR/plots/$name

  Rscript Get_Fixed_Sites.R \
  $outputDIR/1000Genome_filtered_JUST_JPT_freq.frq \
  $outputDIR/NAGJapan_filtered_freq.frq \
  $outputDIR
fi

fileTest="$outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01.recode.vcf.gz"
if [ -f $fileTest ] && echo "$fileTest Found" || echo "$fileTest Not found"; then
  echo "########### Exclude 0.01 Sites ###########"
  vcftools \
  --gzvcf $outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf.gz \
  --chr $chrom \
  --exclude-positions $outputDIR/RemoveSites_0.01.txt \
  --recode --recode-INFO-all \
  --out $outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01

  bgzip -f $outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01.recode.vcf
  tabix -f -p vcf $outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01.recode.vcf.gz
fi

echo "########### Position Of MutSpect (filtered)###########"
bash PositionOfMutSpect.sh \
$outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01.recode.vcf.gz \
filtered \
$chrom \
$outputDIR

echo "########### Plot Filter ###########"
Rscript PlotFilter.R \
$outputDIR/unfiltered+chr$chrom+Context_TestVariants_position.txt \
$outputDIR/filtered+chr$chrom+Context_TestVariants_position.txt \
$outputDIR

echo "########### MutSpect PipeLine (filtered)###########"
bash MutSpect_PipeLine.sh \
$outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01.recode.vcf.gz \
$chrom \
$TimeStamp.chr$chrom.RemoveSites_0.01 \
$outputDIR
