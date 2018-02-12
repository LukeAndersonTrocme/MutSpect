import multiprocessing as mp
import sys, os

def Bash_cmd(i):
    print('Working on Chrom : '+str(i))
    #Get CHROM POS REF ALT from vcf
    #os.system("gzcat {1}/phase3/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | awk '/^#CHROM/{{y=1;next}}y' | cut -f1,2,4,5 > {1}/AncestralRef/Chr{0}.txt".format(i,path))

    print("#### GetHomHet")

    os.system("bash GetHomHet.sh \
    {1}/phase3/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
    {0} \
    {1}/AncestralRef/".format(i,path))

    print("####AncestralContext")

    os.system("python AncestralContext.py \
    -input {1}/AncestralRef/Chr{0}.txt \
    -chrom {0} -out {1}/AncestralRef/".format(i,path))

    print("#### getMutType_perSample")

    os.system("python getMutType_perSample.py \
    -context {1}/AncestralRef/Chr{0}AncestralContext.txt \
    -gt {1}/AncestralRef/Chr{0}.Filtered_Genotypes.txt.gz \
    -chrom {0} \
    -out {1}/AncestralRef/".format(i,path))

if __name__ == '__main__':
    path = '/Users/luke/genomes/genomes/hg19'
    ListOfChrom=list(range(21,23))
    pool_size=7 #mp.cpu_count()
    pool=mp.Pool(processes=pool_size)
    pool_outputs= pool.map(Bash_cmd, ListOfChrom)
    pool.close()
    pool.join()
