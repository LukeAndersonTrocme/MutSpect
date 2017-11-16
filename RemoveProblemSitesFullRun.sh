#!/bin/bash

##this script will run the pipeline with certain sites removed.

#awk '$3 == "POOR_MAPPING_QUALITY" {print $1,$2}' \
#/Users/luke/Documents/MutSpect/1kG_problemSites/breakdownOfCallability.txt \
#> /Users/luke/Documents/MutSpect/1kG_problemSites/breakdownOfCallability_POOR_MAPPING_QUALITY.txt

vcftools --vcf /Users/luke/Documents/MutSpect/data/1000Genome+JapanSitesIn_ThouGenome.unique.vcf.recode.vcf \
--exclude-positions /Users/luke/Documents/MutSpect/1kG_problemSites/breakdownOfCallability_POOR_MAPPING_QUALITY.txt \
--recode --recode-INFO-all \
--out /Users/luke/Documents/MutSpect/1kG_problemSites/POOR_MAPPING_QUALITY

bgzip -f /Users/luke/Documents/MutSpect/1kG_problemSites/POOR_MAPPING_QUALITY.recode.vcf
tabix -f -p vcf /Users/luke/Documents/MutSpect/1kG_problemSites/POOR_MAPPING_QUALITY.recode.vcf.gz

bash MutSpect_PipeLine.sh \
/Users/luke/Documents/MutSpect/1kG_problemSites/POOR_MAPPING_QUALITY.recode.vcf.gz \
Oct4_POOR_MAPPING_QUALITY_

#awk '$3 == “LOW_COVERAGE” {print $1,$2}' \
#/Users/luke/Documents/MutSpect/1kG_problemSites/breakdownOfCallability.txt \
#> /Users/luke/Documents/MutSpect/1kG_problemSites/breakdownOfCallability_LOW_COVERAGE.txt

vcftools --vcf /Users/luke/Documents/MutSpect/data/1000Genome+JapanSitesIn_ThouGenome.unique.vcf.recode.vcf \
--exclude-positions /Users/luke/Documents/MutSpect/1kG_problemSites/breakdownOfCallability_LOW_COVERAGE.txt \
--recode --recode-INFO-all \
--out /Users/luke/Documents/MutSpect/1kG_problemSites/LOW_COVERAGE

bgzip -f /Users/luke/Documents/MutSpect/1kG_problemSites/LOW_COVERAGE.recode.vcf
tabix -f -p vcf /Users/luke/Documents/MutSpect/1kG_problemSites/LOW_COVERAGE.recode.vcf.gz

bash MutSpect_PipeLine.sh \
/Users/luke/Documents/MutSpect/1kG_problemSites/LOW_COVERAGE.recode.vcf.gz \
Oct4_LOW_COVERAGE_
