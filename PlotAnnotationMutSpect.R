library(ggplot2)
library(Biostrings)
library(reshape2)

p.preFilter<-read.table("/Users/luke/Documents/MutSpect/MutSpectPosition/NoFilter+Context_FINAL_position_ANNOTATED.bed")
r.preFilter<-p.preFilter[which(p.preFilter$V5 == "A" | p.preFilter$V5 == "C"),]

c.preFilter<-p.preFilter[which(p.preFilter$V5 == "G" | p.preFilter$V5 == "T"),]
		t<-DNAStringSet(x=c.preFilter$V4)
	c.preFilter$V4<-as.character(reverseComplement(t))
		t<-DNAStringSet(x=c.preFilter$V5)
	c.preFilter$V5<-as.character(reverseComplement(t))
		t<-DNAStringSet(x=c.preFilter$V6)
	c.preFilter$V6<-as.character(reverseComplement(t))

s.preFilter<-rbind(c.preFilter,r.preFilter)


p<-read.table("/Users/luke/Documents/MutSpect/MutSpectPosition/ONLY_MissingInNAG_2+Context_FINAL_position_ANNOTATED.bed")

r<-p[which(p$V5 == "A" | p$V5 == "C"),]

c<-p[which(p$V5 == "G" | p$V5 == "T"),]
		t<-DNAStringSet(x=c$V4)
	c$V4<-as.character(reverseComplement(t))
		t<-DNAStringSet(x=c$V5)
	c$V5<-as.character(reverseComplement(t))
		t<-DNAStringSet(x=c$V6)
	c$V6<-as.character(reverseComplement(t))

s<-rbind(c,r)
cat<-as.data.frame(unique(s[order(s$V4,s$V5,s$V6),][c("V4","V5","V6")]))

pl<-as.data.frame(s$V7)
pl$Mut<-paste(s$V5,s$V6,sep="-")
count<-as.data.frame(table(pl$Mut))

ggplot(pl, aes(x=pl$Mut, fill=pl$'s$V7'))+geom_bar()+geom_text(stat='count',aes(label=..count..),vjust=0)+theme_minimal()+ylab("number of sites removed")+xlab("Mutation Type")+ggtitle("Number of sites remove for each mutation type", "Filter applied : .postFilter")+theme(axis.text.x=element_text(angle = 90), plot.title =element_text(hjust=0.5), plot.subtitle=element_text(hjust=0.5))