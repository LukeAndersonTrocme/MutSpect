library(ggplot2)

args = commandArgs(trailingOnly=T)
count = read.table(paste(args[1],"_derived_each_lineage_chr",args[2],"_nosingle.txt",sep=""), stringsAsFactors=F, fill=T)
#count<-read.table("/Users/luke/Documents/MutSpect/OutputFiles/ReTest/ReTest_derived_each_lineage_chr22_nosingle.txt", stringsAsFactors=F, fill=T)
repos = (args[3])
out = (args[4])
#repos<-("~/bin/smaller_mut_spectrum_pipeline/")

##Read input file with all SNPs per individual

popNames<-read.table(paste(repos,"1000genomes_phase3_sample_IDs_NAG.txt", sep=""), fill=T, stringsAsFactors=FALSE)
Names<-t(read.table(paste(repos,"1000G_NAG_masked_snp_HWE_miss_ID.Names.txt", sep="")))
NamePop<-merge(as.data.frame(t(Names)), popNames, by="V1")

populs=list()
populs[[1]]=c('CHB','JPT','CHS','CDX','KHV','NAG')
populs[[2]]=c('GIH','PJL','BEB','STU','ITU')
populs[[3]]=c('CEU','TSI','GBR','FIN','IBS')
populs[[4]]=c('CLM','MXL','PUR','PEL','ACB','ASW')
populs[[5]]=c('YRI','LWK','GWD','MSL','ESN')
big_populs=c('EAS','SAS','EUR','AMR','AFR')
names(populs)<-big_populs


for(f in 1:nrow(popNames)){
	if(NamePop$V3[f] %in% populs[['EAS']]){NamePop$V2[f]<-"EAS"}
	else if(NamePop$V3[f] %in% populs[['SAS']]){NamePop$V2[f]<-"SAS"}
	else if(NamePop$V3[f] %in% populs[['EUR']]){NamePop$V2[f]<-"EUR"}
	else if(NamePop$V3[f] %in% populs[['AMR']]){NamePop$V2[f]<-"AMR"}
	else if(NamePop$V3[f] %in% populs[['AFR']]){NamePop$V2[f]<-"AFR"}
}

count[1,1] = c("AAA_C")
colnames(count)[1]<- "Mut_type"

tot.count <- sapply(seq(2,ncol(count),2), function(i) {
  rowMeans(count[,c(i, i+1)], na.rm=T)
})

colnames(tot.count)<-NamePop$V3
total<-as.data.frame(colSums(tot.count))
total$pop<-NamePop $V2
total$names<-NamePop $V1

ggplot(total, aes(x=total$names,y=total$'colSums(tot.count)', colour=total$pop))+geom_jitter(shape=1, alpha=0.3)+ylab("Total Variants")+xlab("Individuals")+ggtitle("Total Number of Variants Per Individual and Population")+theme(plot.title = element_text(hjust = 0.5), axis.ticks.x=element_blank(), axis.text.x=element_blank())+ guides(colour=guide_legend(title="Continental Populations"))
ggsave(paste(out,"_1000Genome_NAGjapan_totVariants.jpg", sep=""), height=10, width=10)

tot.count <- prop.table(as.matrix(tot.count),2)
tot.count <- tot.count[ , apply(tot.count, 2, var) != 0]

pca.tot = prcomp(t(tot.count), scale=T)
scores.tot = as.data.frame(pca.tot$x)
ggplot(data = scores.tot, aes(x = PC1, y = PC2, label = rownames(scores.tot), colour = colnames(tot.count))) +
  geom_text(alpha = 0.8, size = 2)+
  theme_classic()+
  ggtitle("PCA using Mutation Spectrum using 1000Genome and NAG")+
  theme(plot.title = element_text(hjust= 0.5))
#ggsave("/Users/luke/Documents/MutSpect/OutputFiles/plots/1000Genome_NAGjapan_PCA.jpg", height=10, width=10)
ggsave(paste(out,"_1000Genome_NAGjapan_PCA.jpg", sep=""), height=10, width=10)
 
colnames(tot.count)<-NamePop$V3
NAG<-tot.count[, grep("NAG", colnames(tot.count))]

pca.nag = prcomp(t(NAG), scale=T)
scores.nag = as.data.frame(pca.nag$x)
ggplot(data = scores.nag, aes(x = PC1, y = PC2, label = rownames(scores.nag))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(colour = "tomato", alpha = 0.8, size = 1)+
  ggtitle("PCA using Mutation Spectrum using NAG")+
  theme_classic()+
  theme(plot.title = element_text(hjust= 0.5))
#ggsave("/Users/luke/Documents/MutSpect/OutputFiles/plots/NAGjapan_PCA.jpg", height=10, width=10)
ggsave(paste(out,"_NAG_PCA.jpg", sep=""), height=10, width=10)

colnames(tot.count)<-NamePop$V3
JPT<-tot.count[, grep("JPT", colnames(tot.count))]

pca.JPT = prcomp(t(JPT), scale=T)
scores.JPT = as.data.frame(pca.JPT$x)
ggplot(data = scores.JPT, aes(x = PC1, y = PC2, label = rownames(scores.JPT))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_point(colour = "tomato", alpha = 0.8, size = 1)+
  ggtitle("PCA using Mutation Spectrum using JPT")+
  theme_classic()+
  theme(plot.title = element_text(hjust= 0.5))
#ggsave("/Users/luke/Documents/MutSpect/OutputFiles/plots/NAGjapan_PCA.jpg", height=10, width=10)
ggsave(paste(out,"_JPT_PCA.jpg", sep=""), height=10, width=10)



colnames(tot.count)<-NamePop$V2
EAS<-tot.count[, grep("EAS", colnames(tot.count))]
colnames(EAS)<-NamePop[which(NamePop$V2=="EAS"),]$V3
EAS<-EAS[, -which(colnames(EAS) == "NAG")]

pca.eas = prcomp(t(EAS))
scores.eas = as.data.frame(pca.eas$x)
ggplot(data = scores.eas, aes(x = PC1, y = PC2, label = rownames(scores.eas), colour = colnames(EAS))) +
  geom_text(alpha = 0.8, size = 2)+
  theme_classic()+
  ggtitle("PCA using Mutation Spectrum using 1000Genome and NAG")+
  theme(plot.title = element_text(hjust= 0.5))
#ggsave("/Users/luke/Documents/MutSpect/OutputFiles/plots/NAGjapan_PCA.jpg", height=10, width=10)
ggsave(paste(out,"_EAS_PCA.jpg", sep=""), height=10, width=10)
