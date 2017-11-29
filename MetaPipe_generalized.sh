#!/bin/bash
set -e
echo "\
                    __ __ __ __ __
                   /__/__/__/__/__/|
                  /__/__/__/__/__/|/
                  |__'__'__'__'__|/
"

chrom=$1 #this is the chromosome we're dealing with
if [ $2 = "test" ]; then
	TestMode=true
	echo "#%$%@%#$# TEST  MODE #%$%@%#$#"
fi

if [ -z "$3" ] #check if specified TimeStamp is supplied
  then
    echo "No TimeStamp argument supplied. Using Today."
    TimeStamp=$(date +%Y-%m-%d) #TimeStamp used to make the output folder
else
  echo "Using TimeStamp input argument"
  TimeStamp=$3
fi

if [ -z "$4" ]
  then
    echo "No HourStamp argument supplied. Using Now."
    HourStamp=$(date +%H-%M)
else
  echo "Using HourStamp input argument"
  HourStamp=$4
fi

cd /Users/luke/bin/smaller_mut_spectrum_pipeline

PathToGenome="/Users/luke/genomes/genomes" #this is where the data is
#PathToGenome="/Users/luke/genomes/genomes/test22" #this is where the data is
pathToBed="/Users/luke/genomes/BED_MASKS" #this is where the bed masks are
#check if folder exists
mkdir -p /Users/luke/Documents/MutSpect/FilteredData/$TimeStamp
#make a new directory for hour/minutes
mkdir -p "/Users/luke/Documents/MutSpect/FilteredData/$TimeStamp/$chrom"
#set outputDIR
outputDIR="/Users/luke/Documents/MutSpect/FilteredData/$TimeStamp/$chrom/$HourStamp"
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

echo "# # # # # # Time since last Git Commit ###########"
echo "MetaPipe_generalized.sh"
git log -1 --format=%cd MetaPipe_generalized.sh
echo "FiltrationPipeLine2.0.sh"
git log -1 --format=%cd FiltrationPipeLine2.0.sh
echo "MutSpect_PipeLine.sh"
git log -1 --format=%cd MutSpect_PipeLine.sh
echo "PositionOfMutSpect.sh"
git log -1 --format=%cd PositionOfMutSpect.sh
echo "FreqSpectPlot.R"
git log -1 --format=%cd FreqSpectPlot.R
echo "PlotFilter.R"
git log -1 --format=%cd PlotFilter.R

start=`date +%s` #timer
startStep=`date +%s` #timer
#Take inputVCF filter low qual sites

if [ "$TestMode" = true ]; then
	echo "#%$%@%#$# TEST  MODE inside if statement #%$%@%#$#" 
	mkdir -p TestMode/{Nag,hg19/phase3}
	if [ ! -f TestMode/Nag/allSamples.$chrom.genotyped.vcf.gz.tbi ]; then

		gzcat $PathToGenome/NAG/allSamples.$chrom.genotyped.vcf.gz \
		| head -100000 | bgzip -f -c \
		> TestMode/Nag/allSamples.$chrom.genotyped.vcf.gz

		tabix -f TestMode/Nag/allSamples.$chrom.genotyped.vcf.gz
		echo "Got first 100000 lines in NAG"
	fi

	if [ ! -f TestMode/hg19/phase3/\
ALL.chr$chrom.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi ]; then
        	gzcat $PathToGenome/hg19/phase3/\
ALL.chr$chrom.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
        	| head -100000 | bgzip -f -c \
        	> TestMode/hg19/phase3/\
ALL.chr$chrom.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 
	
	 	tabix -f TestMode/hg19/phase3/\
ALL.chr$chrom.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
		echo "Got first 100000 lines in 1000G"	 	
	fi
	PathToGenome="TestMode"
fi
	
if [ ! -f $outputDIR/NAG_chr$chrom.3bed_filtered.vcf.gz ]; then
  echo "########### 1 : Nag Filtration ###########"
  bash FiltrationPipeLine2.0.sh \
  $PathToGenome/NAG/allSamples.$chrom.genotyped.vcf.gz \
  NAG_chr$chrom \
  $outputDIR \
  & echo "  #   #   #    NAG_chr$chrom.3bed_filtered.vcf.gz Kill PID : $!"

else echo "File Exists : $outputDIR/NAG_chr$chrom.3bed_filtered.vcf.gz "
fi

if [ ! -f $outputDIR/1kG_chr$chrom.3bed_filtered.vcf.gz ]; then
  echo "########### 2 : 1000Genome Filtration ###########"
  sleep 10
  bash FiltrationPipeLine2.0.sh \
  $PathToGenome/hg19/phase3/ALL.chr$chrom.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
  1kG_chr$chrom \
  $outputDIR \
  & echo "  #   #   #    1kG_chr$chrom.3bed_filtered.vcf.gz Kill PID : $!"
else echo "File Exists : $outputDIR/1kG_chr$chrom.3bed_filtered.vcf.gz "
fi

wait
echo "# # # # # # Time to run Filtration pipeline on chr $chrom is : \
$((($(date +%s)-$startStep)/60)) minutes or $((($(date +%s)-$startStep)/60/60)) hours"
startStep=`date +%s` #timer

#merge both files
if [ ! -f $outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf.gz ]; then
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
else echo "File Exists : $outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf.gz"
fi

echo "# # # # # # Time to run merging step on chr $chrom is : \
$((($(date +%s)-$startStep)/60)) minutes or $((($(date +%s)-$startStep)/60/60)) hours"


#run MutSpect PipeLine on first sample
echo "########### MutSpect PipeLine ###########"
bash MutSpect_PipeLine.sh \
$outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf.gz \
$chrom \
$TimeStamp.chr$chrom \
$outputDIR \
& echo "  #   #   #  MutSpect_PipeLine  Kill PID : $!"

if [ ! -f $outputDIR/1000Genome_filtered_JUST_JPT_freq.frq ]; then
  startStep=`date +%s` #timer
  echo "########### Get Fixed Sites ###########"
  vcftools \
  --gzvcf $outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf.gz \
  --keep /Users/luke/Documents/MutSpect/data/1000genomes_phase3_sample_IDs_JUST_JPT.txt \
  --freq2 \
  --out $outputDIR/1000Genome_filtered_JUST_JPT_freq \
  & echo "  #   #   #  vcftools JPT  Kill PID : $!"

  sleep 10

  vcftools \
  --gzvcf $outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf.gz \
  --keep /Users/luke/Documents/MutSpect/data/1000genomes_phase3_sample_IDs_JUST_NAG.txt \
  --freq2  \
  --out $outputDIR/NAGJapan_filtered_freq \
  & echo "  #   #   #  vcftools NAG  Kill PID : $!"

  wait
  echo "Rscript Get_Fixed_Sites.R \
  $outputDIR/1000Genome_filtered_JUST_JPT_freq.frq \
  $outputDIR/NAGJapan_filtered_freq.frq \
  $outputDIR
  "
  Rscript Get_Fixed_Sites.R \
  $outputDIR/1000Genome_filtered_JUST_JPT_freq.frq \
  $outputDIR/NAGJapan_filtered_freq.frq \
  $outputDIR
  echo "# # # # # # Time to run Get_Fixed_Sites step on chr $chrom is : \
  $((($(date +%s)-$startStep)/60)) minutes or $((($(date +%s)-$startStep)/60/60)) hours"

else echo "File Exists : $outputDIR/1000Genome_filtered_JUST_JPT_freq.frq "
fi

if [ ! -f $outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01.recode.vcf.gz ]; then
  startStep=`date +%s` #timer
  echo "########### Exclude 0.01 Sites ###########"
  vcftools \
  --gzvcf $outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf.gz \
  --chr $chrom \
  --exclude-positions $outputDIR/RemoveSites_0.01.txt \
  --recode --recode-INFO-all \
  --out $outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01 \
  && bgzip -f $outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01.recode.vcf \
  && tabix -f -p vcf $outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01.recode.vcf.gz #\
  #& echo "  #   #   #  Exclude 0.01 Sites  Kill PID : $!"
  echo "# # # # # # Time to run Exclude 0.01 Sites step on chr $chrom is : \
  $((($(date +%s)-$startStep)/60)) minutes or $((($(date +%s)-$startStep)/60/60)) hours"
else echo "File Exists : $outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01.recode.vcf.gz"
fi

# pre=($(zgrep -Ec "$" /$outputDIR/1kGenome_NAG_filtered_chr$chrom.recode.vcf.gz))
# post=($(zgrep -Ec "$" $outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01.recode.vcf.gz))
# diff=$((pre-post))
# prop=$((diff / pre))
# if [[ ($pre-$post)/$pre > 0.15 ]]; then
#     echo "## ## ## ## ## Checked diff, more than 15%"
#     exit
#
# else echo "## ## ## ## ## Checked diff, less than 15%"
# fi

echo "########### MutSpect PipeLine (filtered 0.01)###########"
bash MutSpect_PipeLine.sh \
$outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01.recode.vcf.gz \
$chrom \
$TimeStamp.chr$chrom.RemoveSites_0.01 \
$outputDIR \
& echo "  #   #   #  MutSpect PipeLine (filtered 0.01)  Kill PID : $!"

if [ ! -f $outputDIR/1000Genome_filtered0.01_JUST_JPT_freq.frq ]; then
  startStep=`date +%s` #timer
  echo "########### SITE FREQUENCY SPECTRUM PLOT (filtered 0.01) ###########"

  vcftools \
  --gzvcf $outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01.recode.vcf.gz \
  --keep /Users/luke/Documents/MutSpect/data/1000genomes_phase3_sample_IDs_JUST_JPT.txt \
  --freq2 \
  --out $outputDIR/1000Genome_filtered0.01_JUST_JPT_freq \
  & echo "  #   #   #  vcftools JPT  Kill PID : $!"

  sleep 10

  vcftools \
  --gzvcf $outputDIR/1kGenome_NAG_filtered_chr$chrom.RemoveSites_0.01.recode.vcf.gz \
  --keep /Users/luke/Documents/MutSpect/data/1000genomes_phase3_sample_IDs_JUST_NAG.txt \
  --freq2  \
  --out $outputDIR/NAGJapan_filtered0.01_freq \
  & echo "  #   #   #  vcftools NAG  Kill PID : $!"
else echo "File Exists : $outputDIR/1000Genome_filtered0.01_JUST_JPT_freq.frq"
fi

wait
echo "# # # # # # Time to run SITE FREQUENCY SPECTRUM PLOT (filtered 0.01) step on chr $chrom is : \
$((($(date +%s)-$startStep)/60)) minutes or $((($(date +%s)-$startStep)/60/60)) hours"
echo "########### Position Of MutSpect (unfiltered)###########"
bash PositionOfMutSpect.sh \
$outputDIR/1kGenome_NAG_filtered_chr$chrom.vcf.gz \
unfiltered \
$chrom \
$outputDIR

echo "########### Position Of MutSpect (filtered 0.01) ###########"
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

echo "Rscript FreqSpectPlot.R \
$outputDIR/1000Genome_filtered_JUST_JPT_freq.frq \
$outputDIR/NAGJapan_filtered_freq.frq \
$outputDIR/filtered_0.01
"
Rscript FreqSpectPlot.R \
$outputDIR/1000Genome_filtered0.01_JUST_JPT_freq.frq \
$outputDIR/NAGJapan_filtered0.01_freq.frq \
$outputDIR/filtered_0.01

echo "Rscript FreqSpectPlot.R \
$outputDIR/1000Genome_filtered_JUST_JPT_freq.frq \
$outputDIR/NAGJapan_filtered_freq.frq \
$outputDIR/filtered\
"
Rscript FreqSpectPlot.R \
$outputDIR/1000Genome_filtered_JUST_JPT_freq.frq \
$outputDIR/NAGJapan_filtered_freq.frq \
$outputDIR/filtered

#if [ TestMode = true ]; then
#	rm TestMode/Nag/allSamples.$chrom.genotyped.vcf.gz
#	rm TestMode/hg19/phase3/\
#        ALL.chr$chrom.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
#fi


echo "# # # # # # Time to run pipeline on chr $chrom is : $((($(date +%s)-$start)/60)) minutes or $((($(date +%s)-$start)/60/60)) hours"
