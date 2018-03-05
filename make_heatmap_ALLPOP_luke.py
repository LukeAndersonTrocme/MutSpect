import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import chi2_contingency
import math
import argparse
import os

if __name__ == "__main__":
        parser = argparse.ArgumentParser(description='make heatmaps for mutation spectra')
        parser.add_argument('-i', dest='inputroot', help='should be same as outputroot from get_finescale_mut_spectra_nosingle_argparse.py')
        parser.add_argument('-out', dest='outputroot', help='output directory')
        #parser.add_argument('-chrom', dest='chrom')

        args = parser.parse_args()

comp=dict({})
comp['A'],comp['C'],comp['G'],comp['T']='T','G','C','A'
ypos, ylabel=[],[]

inv_mut_index=dict({})
mut_index=dict({})
row, col = 0,0

for (b2,d) in [('A','T'),('A','C'),('A','G'),('C','T'),('C','G'),('C','A')]:
    for b1 in 'ACGT':
        col=0
        ypos.append(row+0.5)
        if b1=='T' and b2=='C' and d=='A':
            ylabel.append('5\'-'+b1)
        elif b1=='C':
            ylabel.append(b2+r'$\to$'+d+r'  '+b1)
        else:
            ylabel.append(b1)
#        ylabel.append('5\'-'+b1+b2+r'$\to$'+d)
        for b3 in 'ACGT':
            mut_index[(b1+b2+b3,d)]=(row,col)
            inv_mut_index[(row,col)]=b1+b2+b3+'_'+d
            mut_index[(comp[b3]+comp[b2]+comp[b1],comp[d])]=(row,col)
            col+=1
        row+=1

def frequency_breakdown(allele,pop,spop,start_chr, end_chr):
    count_array=np.full((row,col),0.0000001)
    chromsIncluded=[]
    for chrom in range(start_chr,end_chr):
        fName=args.inputroot+'Chr'+str(chrom)+'_'+pop+'_'+spop+'_'+allele+'.MutTypeCountPerSample.txt'
        if os.path.exists(fName):
            chromsIncluded.append(chrom)
            infile=open(fName)
            lines=infile.readlines()
            infile.close()
            
            s=lines[0].strip('\n').split('\t')
            start_ind=2
            if len(s)<2:
                break
            end_ind=len(s)-2
            while 1.0*end_ind/(len(s)-2)>0.98:
                end_ind-=1
            for line in lines[1:]:
                s=line.strip('\n').split('\t')
                for i in range(start_ind,end_ind):
                    try:
                        count_array[mut_index[(s[0][:3],s[0][4])]]+=int(float(s[i]))
                    except KeyError:
                        next
                        #print("KeyError : "+str(s[0][:3])+" Not in Index")                    
                    
    print("CHROMS INCLUDED IN THE PLOTS : ")
    print(chromsIncluded)
    return count_array

def frequency_breakdown_continental(allele,pop,spop,start_chr, end_chr):
    count_array=np.full((row,col),0.0000001)
    chromsIncluded=[]
    for pop in pops[ThisPop]:
        for chrom in range(start_chr,end_chr):
            fName=args.inputroot+'Chr'+str(chrom)+'_'+pop+'_'+spop+'_'+allele+'.MutTypeCountPerSample.txt'
            if os.path.exists(fName):
                chromsIncluded.append(chrom)
                infile=open(fName)
                lines=infile.readlines()
                infile.close()

                s=lines[0].strip('\n').split('\t')
                start_ind=2
                if len(s)<2:
                    break
                end_ind=len(s)-2
                while 1.0*end_ind/(len(s)-2)>0.98:
                    end_ind-=1
                for line in lines[1:]:
                    s=line.strip('\n').split('\t')
                    for i in range(start_ind,end_ind):
                        try:
                            count_array[mut_index[(s[0][:3],s[0][4])]]+=int(float(s[i]))
                        except KeyError:
                            next
                            #print("KeyError : "+str(s[0][:3])+" Not in Index")                    

        return count_array
    

pop_counts=dict({})
big_pop_counts=dict({})
num_variants=dict({})
big_num_variants=dict({})
##THIS IS WHERE THE POPULATIONS ARE DEFINED
groups=['EUR','EAS','SAS','AFR','AMR']
pops=dict({})
pops['EUR']=['CEU','GBR','FIN','IBS','TSI']
pops['SAS']=['GIH','ITU','PJL','BEB','STU']
pops['EAS']=['CHB','CHS','JPT','KHV','CDX']
pops['AFR']=['YRI','MSL','LWK','ESN','GWD']
pops['AMR']=['CLM','MXL','PUR','PEL','ACB','ASW']

HomHet=['HomRef','Hetero','HomAlt']

bigpop=dict({})
for group in groups:
    for pop in pops[group]:
        bigpop[pop]=group

for h in range(len(HomHet)):
    allele=HomHet[h]
    for p in range(len(groups)):
        ThisPop = groups[p]
        for pop in pops[ThisPop]:
            pop_counts[pop]=frequency_breakdown(allele,ThisPop,pop,1, 23)
            #pop_counts[pop]=frequency_breakdown(allele,ThisPop,pop,int(args.chrom), int(args.chrom)+1)
            num_variants[pop]=pop_counts[pop].sum()

        name=dict({})
        for p in zip(pops[ThisPop],pops[ThisPop]):
            name[p[0]]=p[1]

        for pop_ind1 in range(len(pops[ThisPop])):
            refpop=pops[ThisPop][pop_ind1]
            subplot_ind=1
            for pop_ind2 in range(pop_ind1)+range(pop_ind1+1,len(pops[ThisPop])):
                ratio_list=[]
                pop=pops[ThisPop][pop_ind2]
                #print refpop, pop
                ratio_grid=np.zeros((row,col))
                sig_x,sig_y=[],[]
                for i in range(row):
                    for j in range(col):
                        chi2_results=chi2_contingency(np.array(
                                [[pop_counts[pop][i][j],num_variants[pop]],[pop_counts[refpop][i][j],num_variants[refpop]]]))
                        this_pval=chi2_results[1]
                        ratio_grid[i][j]=pop_counts[pop][i][j]*num_variants[refpop]/(num_variants[pop]*pop_counts[refpop][i][j])
                        if this_pval<0.00001:
                            sig_x.append(j+0.5)
                            sig_y.append(i+0.5)
                plt.subplot(1,len(pops[ThisPop]),subplot_ind)
                subplot_ind+=1
                plt.title(name[pop]+' v\n'+name[refpop],fontsize=10)

                if subplot_ind==2:
                    plt.xticks((0.5,1.5,2.5,3.5),('3\'-A','C','G','T'))
                    plt.yticks(tuple(ypos),tuple(ylabel))
                else:
                    plt.xticks((0.5,1.5,2.5,3.5),('A','C','G','T'))
                    plt.yticks(())
                for k in range(1,6):
                    plt.axhline(y=k*4,color='black')
                plt.pcolor(ratio_grid,vmin=0.84,vmax=1.16,cmap='seismic')
                if subplot_ind==len(pops[ThisPop]):
                    plt.colorbar()
                plt.scatter(sig_x,sig_y,marker='.')

            fig=plt.gcf()
            plt.savefig(args.outputroot+ThisPop+'_heatmap_v_'+refpop+'_'+allele+'_nosingle.pdf',format='pdf')
            plt.clf()

        ##BIG POPS
        big_pop_counts[ThisPop]=frequency_breakdown_continental(allele,ThisPop,pop,1, 23)
        big_num_variants[ThisPop]=big_pop_counts[ThisPop].sum()
        
        name=dict({})
    for p in zip(groups,['Europe','East Asia','South Asia','Africa','America']):
        name[p[0]]=p[1]
        
    for pop_ind1 in range(len(groups)):
        refpop=groups[pop_ind1]
        subplot_ind=1
        for pop_ind2 in range(pop_ind1)+range(pop_ind1+1,len(groups)):
            ratio_list=[]
            pop=groups[pop_ind2]
            print refpop, pop
            ratio_grid=np.zeros((row,col))
            sig_x,sig_y=[],[]
            for i in range(row):
                for j in range(col):
                    print big_pop_counts.keys()
                    print big_pop_counts.keys()
                    chi2_results=chi2_contingency(np.array([[big_pop_counts[pop][i][j],big_num_variants[pop]],[big_pop_counts[refpop][i][j],big_num_variants[refpop]]]))
                    this_pval=chi2_results[1]
                    ratio_grid[i][j]=big_pop_counts[pop][i][j]*big_num_variants[refpop]/(big_num_variants[pop]*big_pop_counts[refpop][i][j])
                    if this_pval<0.00001:
                        sig_x.append(j+0.5)
                        sig_y.append(i+0.5)
            plt.subplot(1,len(groups),subplot_ind)
            subplot_ind+=1
            plt.title(name[pop]+' v\n'+name[refpop],fontsize=10)

            if subplot_ind==2:
                plt.xticks((0.5,1.5,2.5,3.5),('3\'-A','C','G','T'))
                plt.yticks(tuple(ypos),tuple(ylabel))
            else:
                plt.xticks((0.5,1.5,2.5,3.5),('A','C','G','T'))
                plt.yticks(())
            for k in range(1,6):
                plt.axhline(y=k*4,color='black')
            plt.pcolor(ratio_grid,vmin=0.84,vmax=1.16,cmap='seismic')
            if subplot_ind==len(groups):
                plt.colorbar()
            plt.scatter(sig_x,sig_y,marker='.')

        fig=plt.gcf()
        plt.savefig(args.outputroot+ThisPop+'_heatmap_v_'+refpop+'_'+allele+'_contiental.pdf',format='pdf')
        plt.clf()