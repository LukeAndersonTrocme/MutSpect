library(ggplot2)

args = commandArgs(trailingOnly=T)
##get number of columns
no_col <- max(count.fields(args[1], sep = "\t"))
##read files
thou <- read.table(args[1],sep="\t",fill=TRUE,col.names=1:no_col, skip=1)
jap <- read.table(args[2],sep="\t",fill=TRUE,col.names=1:no_col, skip=1)
out = (args[3])
##we only care about the first few columns anyways
thou<-thou[c(1:6)]
jap<-jap[c(1:6)]

join<-merge(jap, thou, all=T, by=c("X1", "X2"))
join[is.na(join)]<-0
join$X5.x=1-join$X5.x
join$X5.y=1-join$X5.y

##POP SIZES
NAGpop=886
THOUpop=104

#make a joint frequency spectrum plot
ggplot(data=join, aes(x=join$X5.x, y=join$X5.y))+
geom_bin2d(binwidth=c(1/NAGpop, 1/THOUpop))+
scale_fill_gradientn(limits=c(0,1000), breaks=seq(0, 1000, by=100), colours=rainbow(3))+
theme_minimal()+
xlab(label="NewJapanese")+ylab(label="100GenomeJapanese")+
ggtitle("SFS of 1000GenomeJapanese and NewJapanese") + 
theme(plot.title = element_text(hjust = 0.5, size=20))

ggsave(filename=paste(out,"_SiteFrequencyJPT_NAG_0.jpg", sep=""), height=12, width=12)
