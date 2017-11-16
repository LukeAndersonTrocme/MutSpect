library(ggplot2)
library(Biostrings)
args = commandArgs(trailingOnly=T)
#read files
#p.preFilter<-read.table("~/Documents/MutSpect/MutSpectPosition/PreFilter+Context_TestVariants_position.txt")
p.preFilter<-read.table(args[1])
#make them all uppercase and only keep SNPs
p.preFilter<-as.data.frame(sapply(p.preFilter, toupper))
p.preFilter<-p.preFilter[nchar(as.character(p.preFilter$V3)) == 1 & nchar(as.character(p.preFilter$V4)) == 1,]

#Deal with the reverse compliment
r.preFilter<-p.preFilter[which(p.preFilter$V3 == "A" | p.preFilter$V3 == "C"),]

c.preFilter<-p.preFilter[which(p.preFilter$V3 == "G" | p.preFilter$V3 == "T"),]
		t<-DNAStringSet(x=c.preFilter$V5)
	c.preFilter$V5<-as.list(as.character(reverseComplement(t)))
		t<-DNAStringSet(x=c.preFilter$V4)
	c.preFilter$V4<-as.list(as.character(reverseComplement(t)))
		t<-DNAStringSet(x=c.preFilter$V3)
	c.preFilter$V3<-as.list(as.character(reverseComplement(t)))

s.preFilter<-rbind(c.preFilter,r.preFilter)

#p.postFilter<-read.table("~/Documents/MutSpect/MutSpectPosition/ONLY_MissingInNAG+Context_TestVariants_position.txt")
p.postFilter<-read.table(args[2])
p.postFilter<-as.data.frame(sapply(p.postFilter, toupper))
p.postFilter<-p.postFilter[nchar(as.character(p.postFilter$V3)) == 1 & nchar(as.character(p.postFilter$V4)) == 1,]

r.postFilter<-p.postFilter[which(p.postFilter$V3 == "A" | p.postFilter$V3 == "C"),]
c.postFilter<-p.postFilter[which(p.postFilter$V3 == "G" | p.postFilter$V3 == "T"),]
		t<-DNAStringSet(x=c.postFilter$V5)
	c.postFilter$V5<-as.list(as.character(reverseComplement(t)))
		t<-DNAStringSet(x=c.postFilter$V4)
	c.postFilter$V4<-as.list(as.character(reverseComplement(t)))
		t<-DNAStringSet(x=c.postFilter$V3)
	c.postFilter$V3<-as.list(as.character(reverseComplement(t)))

s.postFilter<-rbind(c.postFilter,r.postFilter)

#Get all categories
cat<-as.data.frame(unique(s.preFilter[order(s.preFilter$V3,s.preFilter$V4,s.preFilter$V5),][c("V3","V4","V5")]))
cat$freq<-NA
for(f in 1:nrow(cat)){
	this_mut.preFilter<-s.preFilter[which(s.preFilter$V3 == cat[f,1] & s.preFilter$V4 == cat[f,2] & s.preFilter$V5 == cat[f,3]),]
	this_mut.postFilter<-s.postFilter[which(s.postFilter$V3 == cat[f,1] & s.postFilter$V4 == cat[f,2] & s.postFilter$V5 == cat[f,3]),]
	#Get proportion of sites removed by this filter
	cat$freq[f]<-((nrow(this_mut.preFilter)-nrow(this_mut.postFilter))/nrow(this_mut.preFilter))*100
}

#Merge two columns for plotting
cat$Mut<-paste(cat$V5,cat$V4,sep="->")
#only keep the 96 sites
cat1<-cat[complete.cases(cat),]
#color the categories appropriately
for(f in 1:nrow(cat1)){
	if(substr(cat1$Mut[f],2,3)=="AC" &
		substr(cat1$Mut[f],6,6)=="C"){
		cat1$color[f]='red'
	}
	else if(substr(cat1$Mut[f],1,3)=="TAT" &
		substr(cat1$Mut[f],6,6)=="T"){
		cat1$color[f]='blue'
	}
	else{cat1$color[f]='black'}
	}
cat1$Mut<-factor(cat1$Mut, levels=cat1[order(-cat1$freq),]$Mut)

#mean<-mean(cat1$freq)
#dpois(96,mean)
ggplot(cat1, aes(x=cat1$Mut,y =cat1$freq, fill=cat1$color))+geom_bar(stat="identity")+theme_minimal()+ylab("fraction of sites removed")+xlab("Mutation Type")+ggtitle("Fraction of sites remove for each mutation type", "Filter applied : .postFilter")+theme(axis.text.x=element_text(angle = 90), plot.title =element_text(hjust=0.5), plot.subtitle=element_text(hjust=0.5))
ggsave(paste("~/Documents/MutSpect/MutSpectPosition/FractionOfSitesRemoved.jpg",sep=""))

ggsave(paste(args[3],"FractionOfSitesRemoved",".jpg",sep=""))