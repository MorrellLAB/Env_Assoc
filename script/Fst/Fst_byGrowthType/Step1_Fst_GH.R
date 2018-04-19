#Author: Ana Poets
# Description: Allele frequency differenciation among the landraces at different latitudes
#################################################################################################
rm(list=ls())

library(hierfstat)

# === IMPORT INPUT FILES ===========================================================================
#Load the growth type assignations 
assig<-read.csv("/Users/Mia/Dropbox/GITHUB/Barley_Landrace_EnvAssoc/Fst/Fst_byGrowthType/784_landrace_habit.csv",header=F,sep=",")
# Import genotypes. Downloaded from Github companion of Poets et al 2015.
genotypes<-read.table("~/Dropbox/Landrace_Environmental_Assocation/Data/Poets_etal2015/Land_6152_SNPs_AB.txt",row.names=1, header=T)
# Select new geographic location, after Fumi inspected every data point. These are the phenotypes 
# Fumi used for GAPIT
Latlong<-read.csv("~/Dropbox/Landrace_Environmental_Assocation/Analyses/GWAS-GAPIT/Input/myY1.v2.csv",row.names=1)

# Import genetic map from Munñoz et al 2011
genMap<-read.table("~/Dropbox/Landrace_Environmental_Assocation/Data/Muñoz_etal2011/Consensus_iSelect.txt",header=T)

# ==== DATA TRANFORMATIONS ========================================================================
# Change genotypes to numeric values AA=2, BB=0, AB=1 ,NA=NA
CHANGE<-function(dat){
	dat[which(dat == "AA")]<-2
	dat[which(dat == "BB")]<-0
	dat[which(dat == "AB")]<-1
	return(dat)	
}
GENOTYPE_num<-as.data.frame(apply(genotypes,2, CHANGE))

### === 1. Divide samples in within the Spring and winter

#within wild range
SPRING<-assig[which(assig[,2] == 'S'),1]

WINTER<-assig[which(assig[,2] == 'W'),1]

### the next tree rows will arrange the files for Focal Fst analysis
# A. SPRING vs WINTER
as.character(SPRING)->all_pops
GENOTYPE_num[(row.names(GENOTYPE_num) %in% all_pops),]->TOTAL
GENOTYPE_num[(row.names(GENOTYPE_num) %in% WINTER),]->FOCAL


#Ready to calculate hierFst Fst
loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))
levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))

#look at the Fst at each locus
	Fst<-NULL
		for (i in 1:dim(loci)[[2]]) {
		
         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
        }
  Fst_SPRING_vs_WINTER_all<-as.data.frame(Fst) 
  row.names(Fst_SPRING_vs_WINTER_all)<-colnames(loci)
  
 #write.table(Fst_SPRING_vs_WINTER_all, "~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/GrowthHabit/Fst_SPRING_vs_WINTER_all.txt",quote=F,row.names=T,col.names=T,sep="\t")

# B. Select 80 individuals from each partition SPRING vs WINTER
# Run 10 iterations with 80 samples each
# Store the values from each iteration
Fst_SPRING_vs_WINTER_80samples_Iteration<-as.data.frame(rep(NA, dim(GENOTYPE_num)[2]))

x<-1
while (x <=10){
	x
	sample(SPRING,80,replace=F)->all_pops
	#c(SPRING)->all_pops
	GENOTYPE_num[(row.names(GENOTYPE_num) %in% all_pops),]->TOTAL
	GENOTYPE_num[(row.names(GENOTYPE_num) %in% WINTER),]->FOCAL
	
	
	#Ready to calculate hierFst Fst
	loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))
	levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))
	
	#look at the Fst at each locus
		Fst<-NULL
			for (i in 1:dim(loci)[[2]]) {
			
	         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
	        }
	  Fst_SPRING_vs_WINTER_80samples<-as.data.frame(Fst) # The outlier from this comparison is :SCRI_RS_154144
	  row.names(Fst_SPRING_vs_WINTER_80samples)<-colnames(loci)
	 Fst_SPRING_vs_WINTER_80samples_Iteration <- cbind(Fst_SPRING_vs_WINTER_80samples_Iteration, Fst_SPRING_vs_WINTER_80samples)
	 x<-x+1
}  

# Remove first column of NAs
Fst_SPRING_vs_WINTER_80samples_Iteration_resutls<-Fst_SPRING_vs_WINTER_80samples_Iteration[,-1]

#get the average Fst after10 iterations
Fst_average_80samp_10iter <-as.data.frame(apply(Fst_SPRING_vs_WINTER_80samples_Iteration_resutls, 1, mean))
row.names(Fst_average_80samp_10iter)<-row.names(Fst_SPRING_vs_WINTER_80samples_Iteration_resutls)
 # write.table(Fst_SPRING_vs_WINTER_80samples, "~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/GrowthHabit/Fst_SPRING_vs_WINTER_80samples.txt",quote=F,row.names=T,col.names=T,sep="\t")
write.table(Fst_average_80samp_10iter, "~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/GrowthHabit/Fst_SPRING_vs_WINTER_80samples_average.txt",quote=F,row.names=T,col.names=T,sep="\t")



##============= Add SNP information to the output files =============================================================

NAMES_results<-c("Fst_SPRING_vs_WINTER_all","Fst_average_80samp_10iter")
for ( i in 1:length(NAMES_results)){
	DATA<-get(NAMES_results[i])
	RESULTS<-(cbind(row.names(DATA),as.numeric(DATA[,1])))
	
	#find outliers using the 95 percent tail
	Threshold<-quantile(as.numeric(as.matrix(RESULTS[,2])), 0.972, na.rm=T)
	#abline(h=Threshold, col="red")
	
	#Find SNPs in genetic map
	genMap_results<-genMap[(genMap[,3] %in% (RESULTS[,1])),]
	#Unknown marker positions
	UNKNOWN<-(setdiff((RESULTS[,1]),genMap_results[,3]))
	
	#Separte the results in SNPs with known and Unknown positions. Then order the known position by cM order
	Results_un<-RESULTS[(RESULTS[,1] %in% UNKNOWN),]
	Results_known<-RESULTS[(RESULTS[,1] %in% genMap_results$SNP_name),]
	  Results_known_or<-Results_known[match(genMap_results$SNP_name, Results_known[,1]),]
	  if(identical(as.character(Results_known_or[,1]), as.character(genMap_results$SNP_name)) == FALSE)stop("SNP in results and genMap are in different order")
	
	# Add genetic position
	Results_known_or_genmap<-cbind(as.data.frame(genMap_results[,c(4,6)]),Results_known_or)
	
	#Add "Fake" genetic information to the unknown SNPs
	CHR_UN<-rep("UN",(dim(Results_un)[1]))
	Pos_start<-as.numeric((Results_known_or_genmap[(dim(Results_known_or_genmap)[1]),2]) +1)
	Pos_end<-Pos_start +50
	Pos_UN<- (seq(Pos_start, Pos_end, 0.02))[1:length(CHR_UN)]
	
	CHR_UN_fake<-cbind(CHR_UN,Pos_UN)
	Results_unknown_or_genmap<-cbind(as.data.frame(CHR_UN_fake),Results_un)
	colnames(Results_unknown_or_genmap)<-c("Chromosome","Accumulative_cM","1","2")
	
	# Combine all the results
	Results_all_genmap<-rbind(Results_known_or_genmap,Results_unknown_or_genmap)
	colnames(Results_all_genmap)<-c("Chromosome","Cumulative_cM","SNP","FST")
	write.table(Results_all_genmap, paste("~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/GrowthHabit/", NAMES_results[i],".txt", sep=""),quote=F,row.names=F,col.names=T,sep="\t")

}


