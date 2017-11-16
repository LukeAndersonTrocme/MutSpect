library(VennDiagram)
library(ggplot2)
##read files
no_col <- max(count.fields("~/Documents/MutSpect/OutputFiles/Aug29_/1000Genome_filtered_JUST_JPT_freq.frq", sep = "\t"))
thou <- read.table("~/Documents/MutSpect/OutputFiles/Aug29_/1000Genome_filtered_JUST_JPT_freq.frq",sep="\t",fill=TRUE,col.names=1:no_col, skip=1)

jap <- read.table("~/Documents/MutSpect/OutputFiles/Aug29_/NAGJapan_filtered_freq.frq",sep="\t",fill=TRUE,col.names=1:no_col, skip=1)

toh <- read.table("~/genomes/genomes/hum0015.v1.freq.v1.txt", fill=TRUE, skip=1)
toh$V1 <- gsub("chr","",toh$V1)

##Get columns and rows we want
tohV<-toh[which(toh$V1 == "22"), c(1,2,9)]
thouV<-thou[which(thou$X1 == "22"), c(1,2,6)]
japV<-jap[which(jap$X1 == "22"), c(1,2,6)]
colnames(tohV)=c("CHROM","POS", "FREQ")
colnames(thouV)=c("CHROM","POS", "FREQ")
colnames(japV)=c("CHROM","POS", "FREQ")


##POPULATION SIZES for bin widths
NAGpop=886
THOUpop=104
TOHOKUpop=1040

tot<-merge(tohV, thouV, all=T, by=c("CHROM", "POS"))
total<-merge(tot, japV, all=T, by=c("CHROM", "POS"))
colnames(total)=c("CHROM","POS", "Tohoku","ThouGenome","NAG")
##SINCE WE ARE USING ALT FREQUENCIES
##we assume that missing data is REF i.e. ALT freq ==0%
total[is.na(total)]<-0

##allHave are sites that are present in NAG and Tohoku
allHave <-total[which(total$NAG > 0 & total$Tohoku > 0),]
write.table(allHave[c(1,2)], file = "~/Documents/MutSpect/data/SitesInNAGandTohoku.txt", quote=F, col.names=F, row.names=F)


dat <- transform(total, min = pmin(Tohoku,ThouGenome,NAG), max= pmax(Tohoku,ThouGenome,NAG))
keep <-dat[which(dat$min==0 & dat$max>0.1),]

Tohoku = keep[which(keep$Tohoku > 0),]
ThouGenome = keep[which(keep$ThouGenome > 0),]
NAG = keep[which(keep$NAG > 0),]

Tohoku_ThouGenome = keep[which(keep$Tohoku > 0 & keep$ThouGenome > 0),]
ThouGenome_NAG = keep[which(keep$ThouGenome > 0 & keep$NAG > 0),]
NAG_Tohoku = keep[which(keep$NAG > 0 & keep$Tohoku > 0),]

draw.triple.venn(area1= nrow(Tohoku), area2= nrow(ThouGenome), area3= nrow(NAG), n12= nrow(Tohoku_ThouGenome), n23= nrow(ThouGenome_NAG), n13= nrow(NAG_Tohoku), n123=0, category = c("Tohoku","ThouGenome","NAG")  )


Tohoku.unique = keep[which(keep$Tohoku > 0.1 & keep$ThouGenome == 0 & keep$NAG == 0),]
ThouGenome.unique = keep[which(keep$ThouGenome > 0.1 & keep$NAG == 0 & keep$Tohoku == 0),]
NAG.unique = keep[which(keep$NAG > 0.1 & keep$ThouGenome == 0 & keep$Tohoku == 0),]

ggplot(Tohoku.unique, aes(x=Tohoku.unique$POS, y=0))+
	geom_point(shape=124, color='red')+
	geom_point(data= ThouGenome.unique, aes(x= ThouGenome.unique$POS, y=1), shape=124, color='blue')+
	geom_point(data= NAG.unique, aes(x= NAG.unique $POS, y=2), shape=124, color='black')+
	geom_point(data= Tohoku_ThouGenome, aes(x= Tohoku_ThouGenome $POS, y=3), shape=124, color= 'green')+
	geom_point(data= ThouGenome_NAG, aes(x= ThouGenome_NAG $POS, y=4), shape=124, color='orange')+
	geom_point(data= NAG_Tohoku, aes(x= NAG_Tohoku $POS, y=5), shape=124, color='purple')+
	theme_classic()+
	scale_y_continuous(breaks= 0:5, labels = c("Tohoku","ThouGenome","NAG", "Tohoku+ThouGenome","ThouGenome+NAG","NAG+Tohoku"))+ylab("Sites unique to Population")
	
ggsave("~/Documents/MutSpect/UniqueSitesPerPopulation_0.1.pdf", height=4, width=40)

##Write positions for each category
write.table(Tohoku.unique[c(1,2)], file = "~/Documents/MutSpect/data/SitesIn_Tohoku.unique.txt", quote=F, col.names=F, row.names=F)
write.table(ThouGenome.unique[c(1,2)], file = "~/Documents/MutSpect/data/SitesIn_ThouGenome.unique.txt", quote=F, col.names=F, row.names=F)
write.table(NAG.unique[c(1,2)], file = "~/Documents/MutSpect/data/SitesIn_NAG.unique.txt", quote=F, col.names=F, row.names=F)
write.table(Tohoku_ThouGenome[c(1,2)], file = "~/Documents/MutSpect/data/SitesIn_Tohoku_ThouGenome.txt", quote=F, col.names=F, row.names=F)
write.table(ThouGenome_NAG[c(1,2)], file = "~/Documents/MutSpect/data/SitesIn_ThouGenome_NAG.txt", quote=F, col.names=F, row.names=F)
write.table(NAG_Tohoku[c(1,2)], file = "~/Documents/MutSpect/data/SitesIn_NAG_Tohoku.txt", quote=F, col.names=F, row.names=F)


##Plot SFS of each pop 
ggplot(data=total, aes(x=total$ThouGenome, y=total$NAG))+geom_bin2d(binwidth=c(1/THOUpop, 1/NAGpop))+scale_fill_gradientn(limits=c(0,2000), breaks=seq(0, 2000, by=250), colours=rainbow(3))+theme_minimal()+ylab(label="NewJapanese")+xlab(label="100GenomeJapanese")+ggtitle("SFS of 1000GenomeJapanese and NewJapanese") + theme(plot.title = element_text(hjust = 0.5, size=20))

ggsave(filename="~/Documents/MutSpect/SFS_troubleshoot/SFS_1000GenomeJapanese_NewJapanese.jpg", height=12, width=12)

ggplot(data=total, aes(x=total$ThouGenome, y=total$Tohoku))+geom_bin2d(binwidth=c(1/THOUpop, 1/TOHOKUpop))+scale_fill_gradientn(limits=c(0,2000), breaks=seq(0, 2000, by= 250), colours=rainbow(3))+theme_minimal()+ylab(label="Tohoku")+xlab(label="1000GenomeJapanese")+ggtitle("SFS of 1000GenomeJapanese and Tohoku") + theme(plot.title = element_text(hjust = 0.5, size=20))

ggsave(filename="~/Documents/MutSpect/SFS_troubleshoot/SFS_1000GenomeJapanese_Tohoku.jpg", height=12, width=12)


ggplot(data=total, aes(x=total$Tohoku, y=total$NAG))+geom_bin2d(binwidth=c(1/TOHOKUpop, 1/NAGpop))+scale_fill_gradientn(limits=c(0,2000), breaks=seq(0, 2000, by= 250), colours=rainbow(3))+theme_minimal()+ylab(label="NewJapanese")+xlab(label="Tohoku")+ggtitle("SFS of Tohoku and NewJapanese") + theme(plot.title = element_text(hjust = 0.5, size=20))

ggsave(filename="~/Documents/MutSpect/SFS_troubleshoot/SFS_NewJapanese_Tohoku.jpg", height=12, width=12)
