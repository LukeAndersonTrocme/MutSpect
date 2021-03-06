import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.mlab import PCA
import argparse

if __name__ == "__main__":
        parser = argparse.ArgumentParser(description='make heatmaps for mutation spectra')
        parser.add_argument('-i', dest='inputroot', help='should be same as outputroot from get_finescale_mut_spectra_nosingle_argparse.py')
        parser.add_argument('-out', dest='outputroot', help='output directory')
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
            ylabel.append('5\' -'+b1)
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

infile=open(args.inputroot)
lines=infile.readlines()
s=lines[0].strip('\n').split(' ')


#groups=['EUR','EAS','SAS','AFR']
groups=['EUR','EAS','SAS','AFR','AMR']
pops=dict({})
pops['EUR']=['CEU','GBR','FIN','IBS','TSI']
pops['SAS']=['GIH','ITU','PJL','BEB','STU']
pops['EAS']=['CHB','CHS','JPT','KHV','CDX','NAG']
pops['AFR']=['YRI','MSL','LWK','ESN','GWD']
pops['AMR']=['CLM','MXL','PUR','PEL','ACB','ASW']

longname=dict({})
for tup in zip(pops['EUR'],['CEU','British','Finnish','Spanish','Italian']):
    longname[tup[0]]=tup[1]
for tup in zip(pops['SAS'],['Gujarati','Telugu','Punjabi','Bengali','Tamil']):
    longname[tup[0]]=tup[1]
for tup in zip(pops['EAS'],['Han Chinese','Southern Han','Japanese','Kinh','Dai', 'New_NAG']):
    longname[tup[0]]=tup[1]
for tup in zip(pops['AFR'],['Yoruba','Mende','Luhya','Esan','Gambian']):
    longname[tup[0]]=tup[1]
for tup in zip(pops['AMR'],['Colombian','Mexican','Puerto Rican','Peruvian','Afro Caribbean','African American']):
    longname[tup[0]]=tup[1]

bigpop=dict({})
for group in groups:
    for pop in pops[group]:
        bigpop[pop]=group

indices=dict({})
for group in groups:
    indices[group]=[]
    for pop in pops[group]:
        indices[pop]=[]

pop_this_index=dict({})
for i in range(1,len(s)):
    pop_this_index[i]=s[i]
    indices[s[i]].append(i-1)
    indices[bigpop[s[i]]].append(i-1)
print pop_this_index

mut_counts=np.zeros((2*(len(s)-1),len(lines)-1))
print len(mut_counts),len(mut_counts[0])

mut_list=[]
for chrom in range(22,23):
    #infile=open('mut_type_v_allele_freq_'+pop+'_chr'+args.chrom+'_nosingle.txt')
    infile=open(args.inputroot)
    lines=infile.readlines()
    infile.close()

    for i in range(len(lines)-1):
        s=lines[i+1].strip('\n').split(' ')
        if chrom==1:
            mut_list.append(s[0])
        for j in range(len(s)-1):
#            print j,i,len(s)
            mut_counts[j][i]+=int(s[j+1])
for j in range(len(s)-1):
    der_count=mut_counts[j].sum()
    for i in range(len(mut_counts[j])):
        mut_counts[j][i]*=1.0/der_count

averaged_mut_counts=[]
for j in range((len(s)-1)/2):
    averaged_mut_counts.append([])
    for i in range(len(mut_counts[0])):
        averaged_mut_counts[-1].append(0.5*(mut_counts[2*j][i]+mut_counts[2*j+1][i]))
mut_counts=np.array(averaged_mut_counts)

for group in groups:
    #print group
    group_mut_counts=[]
    for i in indices[group]:
        #print mut_counts[i]
        group_mut_counts.append(mut_counts[i])
    group_mut_counts=np.array(group_mut_counts)
    myPCA=PCA(group_mut_counts)

    colors=['blue','green','red','purple','black','orange']
    for i in range(len(pops[group])):
        x,y=[],[]
        for ind in indices[pops[group][i]]:
            this_point=myPCA.project(mut_counts[ind])
            x.append(this_point[0])
            y.append(this_point[1])
        plt.scatter(x,y,color=colors[i],label=longname[pops[group][i]])
    plt.legend(loc='lower left',ncol=2,prop={'size':8})
    plt.xticks(())
    plt.yticks(())
    plt.xlabel('PC1 ('+str(int(100*myPCA.fracs[0]))+'% variance explained)')
    plt.ylabel('PC2 ('+str(int(100*myPCA.fracs[1]))+'% variance explained)')
    fig=plt.gcf()
    fig.set_size_inches((4.5,3.5))
    plt.savefig(args.outputroot+group+'_mut_PCA_1kg_nosingle_altlegend.pdf')
    plt.clf()

    PC_vector=myPCA.Wt

    for j in range(2):
        PC_grid=np.zeros((row,col))
        for i in range(len(mut_list)):
            mut=mut_list[i]
            PC_grid[mut_index[(mut[:3],mut[4])]]-=float(PC_vector[j][i])
        plt.xticks((0.5,1.5,2.5,3.5),('3\' -A', 'C','G','T'))
        plt.yticks(tuple(ypos),tuple(ylabel))
        for i in [4,8,12,16,20]:
            plt.axhline(y=i,color='black')
        plt.pcolor(PC_grid,cmap='seismic',vmin=-0.45,vmax=0.45)
        plt.colorbar()
        fig=plt.gcf()
        fig.set_size_inches((4,6))
        plt.title('Principal Component '+str(j+1))
        plt.tight_layout()
        plt.savefig(args.outputroot+group+'_mut_PC'+str(j)+'_nosingle.pdf',format='pdf')
        plt.clf()
