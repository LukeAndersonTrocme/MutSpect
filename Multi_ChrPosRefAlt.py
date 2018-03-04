import multiprocessing as mp
import sys, os

def Bash_cmd(i):
    print('Working on Chrom : '+str(i))

    myVcf = "{1}/phase3/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz".format(i,path)

    my_path = "{0}/AncestralRef/header.txt".format(path)
    if os.path.isfile(my_path) == False :
        print("#### Get header from vcf".format(i))
        os.system("gzcat myVcf | head -500 | grep '^#CHROM'| cut -f10- | awk '{print 'POS\tCHROM\tREF\tALT\tCHIMP.FLIP\tCONTEXT\t' $0}' > {1}/AncestralRef/header.txt")

    my_path1 = "{1}/AncestralRef/Chr{0}.txt".format(i,path)
    if os.path.isfile(my_path1) == False :
        print("#### Get CHROM POS REF ALT from vcf Chrom : {0}".format(i))
        os.system("gzcat myVcf | awk '/^#CHROM/{{y=1;next}}y' | cut -f1,2,4,5 > {1}/AncestralRef/Chr{0}.txt".format(i,path))

    my_path2 ="{1}/AncestralRef/Chr{0}.AncestralContext.txt".format(i,path)
    #print(my_path2)
    if (os.path.isfile(my_path2) == False) or ((os.path.getsize(my_path2) < 100) == True):
        print("#### AncestralContext Chrom : {0}".format(i))

        os.system("python AncestralContext.py \
        -input {1}/AncestralRef/Chr{0}.txt \
        -chrom {0} -out {1}/AncestralRef/".format(i,path))

    my_path3 = "{1}/AncestralRef/Chr{0}.FlippedGenotypes.Context.txt.gz".format(i,path)
    #print((os.path.getsize(my_path3) < 100))
    if (os.path.isfile(my_path3) == False) or ((os.path.getsize(my_path3) < 100) == True):
        print("#### GetHomHet Chrom : {0}".format(i))

        os.system("bash GetHomHet.sh {2} {0} \
        {1}/AncestralRef/".format(i,path,myVcf))

    print("#### getMutType_perSample Chrom : {0}".format(i))

    os.system("python getMutType_perSample.py \
    -gt {1}/AncestralRef/Chr{0}.FlippedGenotypes.Context.txt.gz \
    -chrom {0} \
    -out {1}/AncestralRef/RunAll/".format(i,path))

if __name__ == '__main__':
    path = '/Users/luke/genomes/genomes/hg19'
    ListOfChrom=list(range(1,5))
    pool_size=6 #mp.cpu_count()
    pool=mp.Pool(processes=pool_size)
    pool_outputs= pool.map(Bash_cmd, ListOfChrom)
    pool.close()
    pool.join()
