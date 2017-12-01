from copy import deepcopy
import sys
import gzip
import argparse
from tqdm import tqdm
import os
import logging
import re

def eprint(*args, **kwargs):
    print '*args, file=sys.stderr, **kwargs'


def main(args):
    ##Read Sample ID file : Columns of file :
    ##SAMPLE_ID FAMILY_ID Population_ID Populaion_description
    ##File was made using AddSampleIDs.R saved in pipeline folder
    sampleIDs = open(args.sampleIDs)
    IDs = sampleIDs.readlines()
    sampleIDs.close()

    ##Read Human Reference
    infile = open(args.repos
                + '/hg19_reference/chr'
                + args.chrom
                + '_oneline.txt')

    refseq = infile.read()
    infile.close()
    ##Read Chimp Reference
    infile=open(args.repos
                + '/hg19_chimp_align/human_chimp_diffs_chr'
                + args.chrom
                + '.txt')
    anc_lines = infile.readlines()
    #print(len(anc_lines))
    infile.close()

    ##Make Population Dictionary
    big_populs = ['EAS', 'SAS', 'EUR', 'AMR', 'AFR']
    populs = dict({})
    #removed CHD and SAG
    populs['EAS'] = ['CHB', 'JPT', 'CHS', 'CDX', 'KHV', 'NAG', 'EAS']
    populs['EUR'] = ['CEU', 'TSI', 'GBR', 'FIN', 'IBS', 'EUR']
    populs['AFR'] = ['YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ACB', 'ASW', 'AFR']
    populs['SAS'] = ['GIH', 'PJL', 'BEB', 'STU', 'ITU', 'SAS']
    populs['AMR'] = ['CLM', 'MXL', 'PUR', 'PEL', 'AMR']
    group = dict({})
    allpops = []
    for bigpop in big_populs:
        for pop in populs[bigpop]:
            allpops.append(pop)
            group[pop] = bigpop

    pop_thisID = dict({})
    for line_IDs in IDs:
        s_IDs = line_IDs.split('\t')
        if len(s_IDs) > 3 :
            pop_thisID[s_IDs[0]] = s_IDs[2]

    muts=[]
    for b1 in 'ACGT':
        for b2 in 'ACGT':
            for b3 in 'ACGT':
                for b4 in 'ACGT':
                    if not b2 == b4:
                        muts.append((b1 + b2 + b3, b4))
    #print muts
    ##Open VCF file
    #print 'opening file'
    vcf = gzip.open(args.vcf)
    #print 'file open : ' + args.vcf
    ##Get Line number of VCF (used for progress bar)
    lineString = os.popen("gzcat "
                        + args.vcf
                        + " | grep -v '^#' | wc -l ").read()
    #print 'this is the lineNumbers  :' + lineString
    lineNo = int(lineString.split()[0])

    line_vcf = vcf.readline()
    while not line_vcf.startswith('#CHROM'):
        line_vcf = vcf.readline()

    #print 'fast forwarded through file'
    s_vcf = line_vcf.strip('\n').split('\t')
    num_lineages = 2 * (len(s_vcf) - 9)

    mut_count_derived = dict({})
    for hap_ind in range(num_lineages):
        for mut in muts:
            mut_count_derived[(mut, hap_ind)] = 0

    output = 'Mut_type'

    for i in range(9, len(s_vcf)):
        output += ' ' + pop_thisID[s_vcf[i]]
    output += '\n'

    popul = dict({})
    indices = dict({})

    for pop in allpops:
        indices[pop] = []

    #print "count indices"
    for i in range(9, len(s_vcf)):
        popul[i] = pop_thisID[s_vcf[i]]
        indices[popul[i]].append(i)

    count=dict({})

    #print "ancestral lines"
    anc_lines.pop(0)
    anc_ind = 0
    countS=0
    while anc_ind < len(anc_lines):
	s = anc_lines[anc_ind].strip('\n').split(' ')
        #print('s is :  ' +str(s))
        try :
            if s[1] == 'SNP':
                anc_lines[anc_ind] = deepcopy(s)
                anc_ind += 1
            else:
                anc_lines.pop(anc_ind)
        except IndexError:
            print "indexError.. Check inputVCF"
            countS+=1
	    if countS >= 100:
		print "Something is wrong"
		sys.exit()
	    continue
    anc_ind = 0

    #print "mutations per pop"
    output = dict({})
    mut_count = dict({})
    for pop in allpops:
        output[pop] = 'Mut'
        for i in range(1, 2 * len(indices[pop]) + 1):
            output[pop] += ' ' + str(i)
            for mut in muts:
                mut_count[(mut,pop,i)] = 0
        output[pop] += '\n'

    ##DERIVED
    output_derived = 'Mut_type'
    for i in range(9, len(s)):
        output_derived += ' ' + pop_thisID[s[i]]
        output_derived += '\n'

    #print 'get context'

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    ##SANITY CHECKS :
    ReverseAnc = 0
    NonReverseAnc = 0
    MissCount = 0
    Ref_count = 0
    Alt_count = 0
    for line in tqdm(vcf, total = lineNo):
        s=line.strip('\n').split()
        pos=int(s[1])
        context=refseq[pos - 2 : pos + 1]

        if s[3] in ['A','C','G','T']\
                and s[4] in ['A','C','G','T']\
                and s[3]!= s[4]\
                and len(s[3] + s[4]) == 2\
                and not 'N' in context:
            #go through the chimp file until it matches the position
            while anc_ind < len(anc_lines) - 1\
                    and int(anc_lines[anc_ind][0]) <   pos:
                anc_ind += 1
            #if chimps differ than humans for this position
            #then we use the chimp as reference
            if int(anc_lines[anc_ind][0]) == pos\
                    and s[4] == anc_lines[anc_ind][3]:
                reverse=True
                der_allele = '0'
                this_mut = (context[0] + s[4] + context[2], s[3])
                ReverseAnc += 1
            #if the position doesn't exist in the chimp file
            #humans and chimps are the same for this position
            #use the human ref as the ref.
            else:
                reverse = False
                der_allele = '1'
                this_mut = (context, s[4])
                NonReverseAnc += 1

            i = 9
            der_observed = 0
            for pop in allpops:
                count[pop] = 0

            while i < len(s):
                for j in [0,2]:
                    if s[i][j] == der_allele:
                        mut_count_derived[(this_mut,
                                           2 * (i - 9) + j / 2)] += 1
                        count[popul[i]] += 1
                        der_observed += 1
                        Alt_count += 1
                    elif s[i][j] == ".":
                        MissCount += 1
                    else:
                        Ref_count += 1

                i += 1

            for pop in allpops:
                if count[pop] > 0:
                    mut_count[(this_mut, pop, count[pop])] += 1

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
    # print "ReverseAnc : " + str(ReverseAnc)
    # print "NonReverseAnc : " + str(NonReverseAnc)
    # print "skippedLines : " + str(lineNo - NonReverseAnc - ReverseAnc)
    # print "Ref_count : " + str(Ref_count)
    # print "Alt_count : " + str(Alt_count)
    # print "MissCount : " + str(MissCount)
    # print "total Genotypes per line : " +str((Ref_count + Alt_count + MissCount)/(ReverseAnc + NonReverseAnc))
    # print "number of Ind : " + str(len(s))
    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

    for pop in allpops:
        for mut in muts:
            output[pop] += mut[0] + '_' + mut[1]
            for i in range(1, 2 * len(indices[pop]) + 1):
                output[pop] += ' '\
                            + str(mut_count[(mut, pop, i)])
            output[pop] += '\n'

        outfile=open(args.outputroot\
                    + 'mut_type_v_allele_freq_'\
                    + pop\
                    + '_chr'\
                    + args.chrom\
                    + '_nosingle.txt','w')

        outfile.write(output[pop])
        outfile.close()

    print 'finished chrom ', args.chrom
    print 'output files are in : ', args.outputroot

    for mut in muts:
        output_derived += mut[0] + '_' + mut[1]
        for i in range(num_lineages):
            output_derived += ' '\
                            + str(mut_count_derived[(mut,i )])
        output_derived += '\n'

    outfile=open(args.outputroot\
                + '_derived_each_lineage_chr'\
                + args.chrom\
                + '_nosingle.txt','w')

    outfile.write(output_derived)
    outfile.close()

    ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

if __name__ == "__main__":
        parser = argparse.ArgumentParser(description = 'get finescale mutation spectra')
        parser.add_argument('-chrom',
                            help = 'chromosome')
        parser.add_argument('-vcf',
                            help = 'input vcf')
        parser.add_argument('-id',
                            dest = 'sampleIDs',
                            help='sample_IDs.txt \
                            ( /Users/luke/bin/smaller_mut_spectrum_pipeline/\
                            1000genomes_phase3_sample_IDs_NAG.txt)')
        parser.add_argument('-repos',
                            help = 'directory containing scripts and input files \
                            (/Users/luke/bin/smaller_mut_spectrum_pipeline/)')
        parser.add_argument('-out',
                            dest = 'outputroot', help='output directory')
        args = parser.parse_args()
        main(args)
