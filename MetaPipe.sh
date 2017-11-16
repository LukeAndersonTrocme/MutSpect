#!/bin/bash

#metaPipe
cd /Users/luke/bin/smaller_mut_spectrum_pipeline
PathToGenome="/Users/luke/genomes/genomes"
PathToData="/Users/luke/Documents/MutSpect/FilteredData"

bash FiltrationPipeLine2.0.sh \
$PathToGenome/NAG/allSamples.20.genotyped.vcf.gz \
NAG_chr22

bash FiltrationPipeLine2.0.sh \
$PathToGenome/hg19/phase3/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
1kG_chr22

/usr/local/bin/bcftools-1.6/bcftools merge \
-0 $PathToData/2017-10-30+1kG_chr22/1kG_chr22.3bed_filtered.vcf.gz \
$PathToData/2017-10-30+NAG_chr22/NAG_chr22.3bed_filtered.vcf.gz \
| /usr/local/bin/bcftools-1.6/bcftools view \
--output-type z \
> $PathToData/2017-10-30_1kGenome_NAG_filtered_chr22.vcf

/usr/local/bin/bcftools-1.6/bcftools filter \
--regions-file $pathToBed/Illu_PE.exclusion_complement.bed \
--output-type z $PathToData/2017-10-30_1kGenome_NAG_filtered_chr22.vcf.gz \
| /usr/local/bin/bcftools-1.6/bcftools view \
--output-file $PathToData/2017-10-30_1kGenome_NAG_filtered_chr22_callable.vcf.gz \
| /usr/local/bin/bcftools-1.6/bcftools index -t \
--output-file $PathToData/2017-10-30_1kGenome_NAG_filtered_chr22_callable.vcf.gz.tbi

bash MutSpect_PipeLine.sh $PathToData/2017-10-30_1kGenome_NAG_filtered_chr22.vcf.gz 22 Oct30_Filter2.0
