args = commandArgs(trailingOnly=T)

no_col <- max(count.fields(args[1], sep = "\t"))
##read files
#thou <- read.table("~/Documents/MutSpect/OutputFiles/Oct30_Filter2.0/1000Genome_filtered_JUST_JPT_freq.frq",sep="\t",fill=TRUE,col.names=1:no_col, skip=1)
thou<-read.table(args[1], sep="\t",fill=TRUE,col.names=1:no_col, skip=1)

#jap <- read.table("~/Documents/MutSpect/OutputFiles/Oct30_Filter2.0/NAGJapan_filtered_freq.frq",sep="\t",fill=TRUE,col.names=1:no_col, skip=1)
jap<-read.table(args[2], sep="\t",fill=TRUE,col.names=1:no_col, skip=1)
out = (args[3])

##we only care about the first few columns anyways
thou<-thou[c(1:6)]
jap<-jap[c(1:6)]

join<-merge(jap, thou, all=T, by=c("X1", "X2"))
join[is.na(join)]<-0
join$diff<-(join$X5.x-join$X5.y)^2
join$NAG=1-join$X5.x
join$JPT=1-join$X5.y

MissingJPT<-join[which(join$NAG > 0.01 & join$JPT ==0),]
MissingNAG<-join[which(join$JPT > 0.01 & join$NAG ==0),]
Remove<-join[which(join$NAG > 0.1 & join$JPT == 0| join$JPT > 0.1 & join$NAG == 0),]
Remove0.01<-join[which(join$NAG > 0.01 & join$JPT == 0| join$JPT > 0.01 & join$NAG == 0),]
Remove0.05<-join[which(join$NAG > 0.05 & join$JPT == 0| join$JPT > 0.05 & join$NAG == 0),]



write.table(MissingJPT[,c(1,2)], file=paste(out,"/MissingInJPT.txt",sep=""), quote=F, row.names=F, col.names=F)
write.table(MissingNAG[,c(1,2)], paste(out,"/MissingInNAG.txt",sep=""), quote=F, row.names=F, col.names=F)
write.table(Remove[,c(1,2)], file=paste(out,"/RemoveSites_0.1.txt",sep=""), quote=F, row.names=F, col.names=F, sep = "\t")
write.table(Remove0.01[,c(1,2)], file=paste(out,"/RemoveSites_0.01.txt",sep=""), quote=F, row.names=F, col.names=F, sep = "\t")
write.table(Remove0.05[,c(1,2)], file=paste(out,"/RemoveSites_0.05.txt",sep=""), quote=F, row.names=F, col.names=F, sep = "\t")


#ggplot(data= MissingNAG, aes(x= MissingNAG $X5.x, y= MissingNAG $X5.y))+geom_bin2d(binwidth=c(0.006))+scale_fill_gradientn(limits=c(0,1000), breaks=seq(0, 1000, by=100), colours=rainbow(3))+theme_minimal()+xlab(label="NewJapanese")+ylab(label="100GenomeJapanese")+ggtitle("Site Frequency Spectrum of 1000GenomeJapanese and NewJapanese") + theme(plot.title = element_text(hjust = 0.5, size=20))
