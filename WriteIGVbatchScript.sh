#!/bin/bash

# gshuf -n 50 /Users/luke/Documents/MutSpect/data/Dec4_SitesIn_ThouGenome.unique.txt \
# | awk -v s=500 '{print $1, $2-s"-"$2+s}' \
# | sed 's/22 /chr22:/g' \
# > /Users/luke/Documents/MutSpect/data/Dec4_SitesIn_ThouGenome.unique_random50.txt


echo "new
genome hg19
load  /Users/luke/genomes/alignments/hg19/NA18939.mapped.ILLUMINA.bwa.JPT.low_coverage.20130415.bam
load /Users/luke/genomes/alignments/NAG14231.matefixed.sorted.bam
snapshotDirectory /Users/luke/igv/MutSpect/" \
> /Users/luke/Documents/MutSpect/IGV_ThouGenome_siteCheck_nag_1kg_NOV15.txt

while read p; do
  echo "goto $p
  sort position
  expand
  snapshot" >> /Users/luke/Documents/MutSpect/IGV_ThouGenome_siteCheck_nag_1kg_Dec4.txt
done < /Users/luke/Documents/MutSpect/data/Dec4_SitesIn_ThouGenome.unique_random50.txt
#/Users/luke/Documents/MutSpect/OutputFiles/Nov13_NAG_1KG_CHR1/MissingInNAG_random50.txt


# new
# genome hg19
# load  /Users/luke/genomes/alignments/hg19/NA18939.mapped.ILLUMINA.bwa.JPT.low_coverage.20130415.bam
# snapshotDirectory mySnapshotDirectory
# goto chr1:65,289,335-65,309,335
# sort position
# collapse
# snapshot
# goto chr1:113,144,120-113,164,120
# sort base
# collapse
# snapshot
# goto chr4:68,457,006-68,467,006
# sort strand
# collapse
# snapshot
