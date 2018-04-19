# Author: Ana Poets
# Description: Plot Fst results
####################################################################################################################
rm(list=ls())

NAMES_results<-c("Fst_2500vs1000", "Fst_5000_vs1000", "Fst_1000vsAll","Fst_5000_vs2500")

for (i in 1:length(NAMES_results)){
	Results_all_genmap <-read.table(paste("/Users/Mia/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/", NAMES_results[i],".txt", sep=""), header=T)
	
	pdf(paste("~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Plots/", NAMES_results[i],".pdf",sep=""), width=7, height=5)
	plot(Results_all_genmap[,2],as.character(Results_all_genmap[,4]), ylab="Fst",xlab="Genetic Position", xaxt="n", cex=0.6)
	
	#Separate by chromosome
	  CHR1<-Results_all_genmap[c(which(Results_all_genmap[,1] == "1H")),]
	  CHR2<-Results_all_genmap[c(which(Results_all_genmap[,1] == "2H")),]
	  CHR3<-Results_all_genmap[c(which(Results_all_genmap[,1] == "3H")),]
	  CHR4<-Results_all_genmap[c(which(Results_all_genmap[,1] == "4H")),]
	  CHR5<-Results_all_genmap[c(which(Results_all_genmap[,1] == "5H")),]
	  CHR6<-Results_all_genmap[c(which(Results_all_genmap[,1] == "6H")),]
	  CHR7<-Results_all_genmap[c(which(Results_all_genmap[,1] == "7H")),]
	  CHR_UN<-Results_all_genmap[c(which(Results_all_genmap[,1] == "UN")),]
	
	points(CHR1[,2],as.character(CHR1[,4]), col="red" , cex=0.6)
	points(CHR2[,2],as.character(CHR2[,4]), col="BLUE" , cex=0.6)
	points(CHR3[,2],as.character(CHR3[,4]), col="GREeN" , cex=0.6)
	points(CHR4[,2],as.character(CHR4[,4]), col="ORANGE" , cex=0.6)
	points(CHR5[,2],as.character(CHR5[,4]), col="magenta" , cex=0.6)
	points(CHR6[,2],as.character(CHR6[,4]), col="cyan" , cex=0.6)
	points(CHR7[,2],as.character(CHR7[,4]), col="dark green" , cex=0.6)
	points(CHR_UN[,2],as.character(CHR_UN[,4]), col="gray" , cex=0.6)
	
	abline(h=quantile(as.numeric(as.character(Results_all_genmap[,4])), 0.975, na.rm=TRUE), col= "red", lty=2)
	
	#find mid point to set tick mark in each chromosome
	tick1<-CHR1 [round(dim(CHR1)[1]/2),2]
	tick2<-CHR2 [round(dim(CHR2)[1]/2),2]
	tick3<-CHR3 [round(dim(CHR3)[1]/2),2]
	tick4<-CHR4 [round(dim(CHR4)[1]/2),2]
	tick5<-CHR5 [round(dim(CHR5)[1]/2),2]
	tick6<-CHR6 [round(dim(CHR6)[1]/2),2]
	tick7<-CHR7 [round(dim(CHR7)[1]/2),2]
	tickUN<-CHR_UN [round(dim(CHR_UN)[1]/2),2]
	
	axis(1, at=c(tick1, tick2, tick3, tick4, tick5, tick6, tick7, tickUN), c("1H","2H","3H","4H","5H","6H","7H","UN"))
dev.off()


	# Get outlier markers
	Threshold<-quantile(as.numeric(as.character(Results_all_genmap[,4])), 0.975, na.rm=TRUE)
	Outliers<-Results_all_genmap[which(Results_all_genmap[,4] >= Threshold),]
	#sort outliers by Fst value
	Outliers_or<-Outliers[order(Outliers[,4], decreasing=T),]
	write.table(Outliers_or,paste("~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/Outliers/", NAMES_results[i],"_outliers.txt",sep=""),quote=F,row.names=F,col.names=T,sep="\t")
}
