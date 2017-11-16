##mut type per position.
args = commandArgs(trailingOnly=T)
#args[1] 	input file
#args[2]	prefix

p<-read.table(args[1])
p<-as.data.frame(sapply(p, toupper))
p<-p[nchar(as.character(p$V3)) == 1 & nchar(as.character(p$V4)) == 1,]
cat<-unique(p[order(p$V3,p$V4,p$V5),][c("V3","V4","V5")])

catNum<-list()
for(f in 1:nrow(cat)){
	this_mut<-p[which(p$V3 == cat[f,1] & p$V4 == cat[f,2] & p$V5 == cat[f,3]),]
	this_mut $end<-as.numeric(as.character(this_mut $V2)) + 1
	write.table(this_mut[c('V1','V2','end')], file = paste("~/Documents/MutSpect/MutSpectPosition/BED/",args[2],"_",cat[f,2],"_",cat[f,3],".bed",sep=""), quote=F, row.names=F, col.names=F, sep="\t")
	#print(cat[f,])
	#print(nrow(this_mut))
	catNum[[paste(cat[f,2],"_",cat[f,3],sep="")]]<-nrow(this_mut)
}

#signal<-p[which(p$V4 == "A" & p$V5 == "CCC"),]
#signal$end<-as.numeric(as.character(signal$V2)) + 1
#write.table(signal[c('V1','V2','end')], file = "/Users/luke/Documents/MutSpect/1000Genome_NAG_filtered_chr22CACtoCCC_pos.bed", quote=F, row.names=F, col.names=F, sep="\t")
