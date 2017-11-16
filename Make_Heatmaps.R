##Make Heatmap
library(stringr)
library(data.table)
library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
setwd(args[1])
#setwd("~/Documents/MutSpect/")


##This function cleans the input data
cleanUp = function(X){
	##add up all the genotypes
	X$sum<-rowSums(X[,c(3:ncol(X))])
	total<-sum(X$sum)
	X$sum<-X$sum/total
	##we just want mutation type and total genotypes
	count<-X[,c(1,ncol(X))]
	
	##substring to separate and group by mutation
	count<-count[,c(1,ncol(count))]
	count$start<-substr(count$Mut, 1,1)
	count$snp<-substr(count$Mut, 2,2)
	count$end<-substr(count$Mut, 3,3)
	count$swap<-substr(count$Mut, 5,5)
	count$Mut<-NULL
	
	##sort the table and extract the unique snp combinations
	sort<-count[order(count$snp,count$swap),]
	uniq <- unique(sort[c("snp","swap")])
	names<-paste(uniq$snp, uniq$swap, sep="→")
	
	##prep table for plotting and grouping
	sort$mut<-paste(sort $snp, sort $swap, sep="→")
	sort$snp<-NULL
	sort$swap<-NULL
	##return it as the output of the function
	return(sort)
	}
##This function takes each pop and compares it to other pops
##It also plots and saves the plots
getDiff = function(ListOfTables, ListOfPops, BigPop){
	for(i in 1:length(listOfTables)){
		listPerPop<-list()
		for(j in 1:length(listOfTables)){
			if(i==j){next}
				
			popRef<-listOfTables[[i]]
			popComp<-listOfTables[[j]]
			
			diff<-merge(popComp,popRef, by=c("start", "end", "mut"))
			diff$ratio<-diff$sum.x/diff$sum.y
			diff$sum.x<-NULL
			diff$sum.y<-NULL
			diff$title<-paste(listOfPops[[i]],"vs",listOfPops[[j]],sep=" ")

			listOfData[[j]] <-diff[order(diff $mut),]
			}
		listPerPop[[1]]<-listOfData	
		#myplots <- lapply(listOfData, plot_data)
		#myPlot<-do.call("grid.arrange", c(myplots, nrow=1))
		#ggsave(paste("MutSpect_HeatMap_ratio",listOfPops[[i]],BigPop,
		#"chr22.jpg",sep="_"),plot=myPlot, height=20, width=20)	
		}
	return(listPerPop)	
	}	

##Function Plotting Heatmap RATIO between populations
plot_data = function (DatPlot){
    ggplot(data= DatPlot, aes(x= DatPlot$end, y= DatPlot$start, fill=DatPlot$ratio))+
		geom_tile()+
		scale_fill_gradient2(low="blue", high="red", mid="white", midpoint=1)+
		labs(y=NULL, x=NULL)+
		guides(fill=F)+
		ggtitle(DatPlot$title[1])+
		facet_grid(DatPlot$mut~.)+
		theme(axis.ticks.x=element_blank(),
			plot.margin=(unit(c(0.1,0.1,0.1,0.1),"cm")),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_blank(),
			panel.background = element_blank(),
			plot.title = element_text(hjust = 0.5),
			strip.background = element_rect(fill="white"), 
			strip.text = element_text(size=8, lineheight=0.1, angle=90))
			}	
##scale_fill_gradientn(colours = c("blue", "white", "red"), 
								#values = scales::rescale(c(0.80, 1, 1.17)) )+
##Get Pop and continental Pops

populs=list()
populs[[1]]=c('CHB','JPT','CHS','CDX','KHV','NAG')
populs[[2]]=c('GIH','PJL','BEB','STU','ITU')
populs[[3]]=c('CEU','TSI','GBR','FIN','IBS')
populs[[4]]=c('CLM','MXL','PUR','PEL','ACB','ASW')
populs[[5]]=c('YRI','LWK','GWD','MSL','ESN')
big_populs=c('EAS','SAS','EUR','AMR','AFR')
names(populs)<-big_populs

##load all file names
filenames <- as.data.frame(list.files(".",pattern="*_chr22_nosingle.txt"))

##This function get's the populations from each BigPop	
getPop = function(X){
	grep(paste(X,collapse="|"),filenames[,1], value=TRUE)}

popFiles<-lapply(populs, getPop)
names(popFiles)<-big_populs

##load all the files	in each pop					
EAS<-lapply(popFiles[['EAS']], read.table, header=T)
SAS<-lapply(popFiles[['SAS']], read.table, header=T)
EUR<-lapply(popFiles[['EUR']], read.table, header=T)
AMR<-lapply(popFiles[['AMR']], read.table, header=T)
AFR<-lapply(popFiles[['AFR']], read.table, header=T)
##Make it happen for each population
lapply(big_populs, function(X) getDiff(lapply(X, cleanUp), populs[['X']],X))

#head(getDiff(lapply(EAS, cleanUp)[[1]], populs[['EAS']], EAS)[[1]][[1]])

#test<-lapply(EAS, cleanUp)[[1]]

#plot_data(getDiff(lapply(EAS, cleanUp)[[1]], populs[['EAS']], EAS)[[1]][[1]])