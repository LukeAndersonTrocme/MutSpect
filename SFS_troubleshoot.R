library(ggplot2)

##read files
no_col <- max(count.fields("~/Documents/MutSpect/data/1000Genome_filtered_JUST_JPT_freq.frq", sep = "\t"))
thou <- read.table("~/Documents/MutSpect/data/1000Genome_filtered_JUST_JPT_freq.frq",sep="\t",fill=TRUE,col.names=1:no_col, skip=1)
jap <- read.table("~/Documents/MutSpect/data/NAGJapan_filtered_freq.frq",sep="\t",fill=TRUE,col.names=1:no_col, skip=1)
toh <- read.table("~/genomes/genomes/hum0015.v1.freq.v1.txt", fill=TRUE, skip=1)
toh$V1 <- gsub("chr","",toh$V1)
thou <- read.table(args[1],sep="\t",fill=TRUE,col.names=1:no_col, skip=1)
jap <- read.table(args[2],sep="\t",fill=TRUE,col.names=1:no_col, skip=1)

##we only care about the first few columns anyways
thou<-thou[c(1:6)]
jap<-jap[c(1:6)]

join<-merge(jap, thou, all=T, by=c("X1", "X2"))
join[is.na(join)]<-0
join$X5.x=1-join$X5.x
join$X5.y=1-join$X5.y


join2<-merge(toh, thou, all=T, by.x=c("V1","V2"), by.y=c("X1", "X2"))
join2[is.na(join2)]<-0
join2$X5=1-join2$X5
join2$V7=1-join2$V7

join3 <-merge(jap, toh, all=T, by.y=c("V1","V2"), by.x=c("X1", "X2"))
join3[is.na(join3)]<-0
join3$X5=1-join3$X5
join3$V7=1-join3$V7

#make a joint frequency spectrum plot
ggplot(data=join, aes(x=join$X5.x, y=join$X5.y))+geom_bin2d(binwidth=c(0.006))+scale_fill_gradientn(limits=c(0,1000), breaks=seq(0, 1000, by=100), colours=rainbow(3))+theme_minimal()+xlab(label="NewJapanese")+ylab(label="100GenomeJapanese")+ggtitle("Site Frequency Spectrum of 1000GenomeJapanese and NewJapanese") + theme(plot.title = element_text(hjust = 0.5, size=20))

ggsave(filename=paste("/Users/luke/Documents/MutSpect/SFS_troubleshoot/SiteFrequencyJPT_NAG_0.jpg", sep=""), height=12, width=12)

ggplot(data=join2, aes(x=join2$V7, y=join2$X5))+geom_bin2d(binwidth=c(0.008))+scale_fill_gradientn(limits=c(0,1000), breaks=seq(0, 1000, by=100), colours=rainbow(3))+theme_minimal()+xlab(label="NewJapanese")+ylab(label="100GenomeJapanese")+ggtitle("SFS of 1000GenomeJapanese and Tohoku") + theme(plot.title = element_text(hjust = 0.5, size=20))

ggsave(filename=paste("/Users/luke/Documents/MutSpect/SFS_troubleshoot/SiteFrequencyJPT_TOH_0.jpg", sep=""), height=12, width=12)

ggplot(data=join3, aes(x=join3$X5, y=join3$V7))+geom_bin2d(binwidth=c(0.008))+scale_fill_gradientn(limits=c(0,1000), breaks=seq(0, 1000, by=100), colours=rainbow(3))+theme_minimal()+xlab(label="NewJapanese")+ylab(label="100GenomeJapanese")+ggtitle("SFS of Tohoku and NewJapanese") + theme(plot.title = element_text(hjust = 0.5, size=20))

ggsave(filename="/Users/luke/Documents/MutSpect/SFS_troubleshoot/SiteFrequencyTOH_NAG_0.jpg", height=12, width=12)