#Author: Ana Poets
# Description: Allele frequency differenciation among the landraces at three different elevations
#################################################################################################
rm(list=ls())

library(hierfstat)

# === IMPORT INPUT FILES ===========================================================================
#Load the genetic assignment for k=4 landraces 
assig<-read.table("/Users/Mia/Dropbox/Landrace_Environmental_Assocation/Data/Poets_etal2015/Landraces_geneticAssingment.txt",header=T)
# Import genotypes. Downloaded from Github companion of Poets et al 2015.
genotypes<-read.table("~/Dropbox/Landrace_Environmental_Assocation/Data/Poets_etal2015/Land_6152_SNPs_AB.txt",row.names=1, header=T)
# Select new geographic location, after Fumi inspected every data point. These are the phenotypes 
# Fumi used for GAPIT
Latlong<-read.csv("~/Dropbox/Landrace_Environmental_Assocation/Analyses/GWAS-GAPIT/Input/myY1.v2.csv",row.names=1)

# Import genetic map from Munñoz et al 2011
genMap<-read.table("/Users/agonzale/Dropbox/Landrace_Environmental_Assocation/Data/Muñoz_etal2011/Consensus_iSelect.txt",header=T)

# ==== DATA TRANFORMATIONS ========================================================================
# Change genotypes to numeric values AA=2, BB=0, AB=1 ,NA=NA
CHANGE<-function(dat){
	dat[which(dat == "AA")]<-2
	dat[which(dat == "BB")]<-0
	dat[which(dat == "AB")]<-1
	return(dat)	
}
GENOTYPE_num<-as.data.frame(apply(genotypes,2, CHANGE))
# Use the altitude to devide the landraces into categories at ≤1000m, 1000≥2500, >2500
m1000<-row.names(Latlong)[which(Latlong$altitude <=1000)]
m2500<-row.names(Latlong)[which(Latlong$altitude > 1000 & Latlong$altitude <=2500)]
m5000<-row.names(Latlong)[which(Latlong$altitude > 2500)]
mless2500<-row.names(Latlong)[which(Latlong$altitude <= 2500)]
mMore3000<-row.names(Latlong)[which(Latlong$altitude >3000)]
mless3000<-row.names(Latlong)[which(Latlong$altitude <= 3000)]

### the next tree rows will arrange the files for Focal Fst analysis
# A. 1000 m vs >1000
c(m2500, m5000)->all_pops
GENOTYPE_num[(row.names(GENOTYPE_num) %in% all_pops),]->TOTAL
GENOTYPE_num[(row.names(GENOTYPE_num) %in% m1000),]->FOCAL


#Ready to calculate hierFst Fst
loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))
levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))

#look at the Fst at each locus
	Fst<-NULL
		for (i in 1:dim(loci)[[2]]) {
		
         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
        }
  Fst_1000vsAll<-as.data.frame(Fst) # The outlier from this comparison is :SCRI_RS_154144
  row.names(Fst_1000vsAll)<-colnames(loci)
  
 # write.table(Fst_1000vsAll, "/Users/Mia/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/Fst_1000vsAll.txt",quote=F,row.names=T,col.names=T,sep="\t")
# ===============B. 1000m vs. 2500m

GENOTYPE_num[(row.names(GENOTYPE_num) %in% m2500),]->TOTAL
GENOTYPE_num[(row.names(GENOTYPE_num) %in% m1000),]->FOCAL


#Ready to calculate hierFst Fst
loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))
levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))

#look at the Fst at each locus
	Fst<-NULL
		for (i in 1:dim(loci)[[2]]) {
		
         Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
        }

Fst_2500vs1000<-as.data.frame(Fst) # The outlier from this comparison is :"X11_10855"      "SCRI_RS_149432" "SCRI_RS_154144"
row.names(Fst_2500vs1000)<-colnames(loci)

#write.table(Fst_2500vs1000, "/Users/Mia/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/Fst_2500vs1000.txt",quote=F,row.names=T,col.names=T,sep="\t")


#================== C. 2500m vs. 5000m ============================================

GENOTYPE_num[(row.names(GENOTYPE_num) %in% m2500),]->TOTAL
GENOTYPE_num[(row.names(GENOTYPE_num) %in% m5000),]->FOCAL


#Ready to calculate hierFst Fst
loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))
levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))

#look at the Fst at each locus
Fst<-NULL
for (i in 1:dim(loci)[[2]]) {
  
  Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
}


Fst_5000_vs2500<-as.data.frame(Fst) # The outlier from this comparison is :"SCRI_RS_146573"
row.names(Fst_5000_vs2500)<-colnames(loci)

#write.table(Fst_5000_vs2500, "/Users/Mia/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/Fst_5000_vs2500.txt",quote=F,row.names=T,col.names=T,sep="\t")


#================== D. 1000m vs. 5000m ============================================

GENOTYPE_num[(row.names(GENOTYPE_num) %in% m1000),]->TOTAL
GENOTYPE_num[(row.names(GENOTYPE_num) %in% m5000),]->FOCAL


#Ready to calculate hierFst Fst
loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))
levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))

#look at the Fst at each locus
Fst<-NULL
for (i in 1:dim(loci)[[2]]) {
  
  Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
}


Fst_5000_vs1000<-as.data.frame(Fst) # The outlier from this comparison is :"X12_11285" "X12_21479"
row.names(Fst_5000_vs1000)<-colnames(loci)

#write.table(Fst_5000_vs1000, "/Users/Mia/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/Fst_5000_vs1000.txt",quote=F,row.names=T,col.names=T,sep="\t")

#####========= E. Compare <=2500m to >2500m
c(m2500, mless2500)->all_pops

GENOTYPE_num[(row.names(GENOTYPE_num) %in% m2500),]->TOTAL
GENOTYPE_num[(row.names(GENOTYPE_num) %in% mless2500),]->FOCAL


#Ready to calculate hierFst Fst
loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))
levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))

#look at the Fst at each locus
Fst<-NULL
for (i in 1:dim(loci)[[2]]) {
  
  Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
}


Fst_less2500_vsmore2500<-as.data.frame(Fst) 
row.names(Fst_less2500_vsmore2500)<-colnames(loci)

#write.table(Fst_less2500_vsmore2500, "/Users/agonzale/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/Elevation/Fst_less2500_vsmore2500.txt",quote=F,row.names=T,col.names=T,sep="\t")

###### -===== F. Compare more than 3000m to less or equal 3000
c(mMore3000, mless3000)

GENOTYPE_num[(row.names(GENOTYPE_num) %in% mMore3000),]->TOTAL
GENOTYPE_num[(row.names(GENOTYPE_num) %in% mless3000),]->FOCAL


#Ready to calculate hierFst Fst
loci<-rbind(as.data.frame(FOCAL),as.data.frame(TOTAL))
levels<-c(rep(0,dim(FOCAL)[1]),rep(1,dim(TOTAL)[1]))

#look at the Fst at each locus
Fst<-NULL
for (i in 1:dim(loci)[[2]]) {
  
  Fst[i]<-varcomp(data.frame(levels,as.numeric(loci[,i])),diploid=FALSE)$F[1,1]   
}


Fst_less3000_vsmore3000<-as.data.frame(Fst) 
row.names(Fst_less3000_vsmore3000)<-colnames(loci)



##=============

NAMES_results<-c("Fst_2500vs1000", "Fst_5000_vs1000", "Fst_1000vsAll","Fst_5000_vs2500","Fst_less2500_vsmore2500","Fst_less3000_vsmore3000")

for ( i in 1:length(NAMES_results)){
	DATA<-get(NAMES_results[i])
	RESULTS<-(cbind(row.names(DATA),as.numeric(DATA[,1])))
	
	#find outliers using the 95 percent tail
	Threshold<-quantile(as.numeric(as.matrix(RESULTS[,2])), 0.972, na.rm=T)
	abline(h=Threshold, col="red")
	
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
	write.table(Results_all_genmap, paste("/Users/agonzale/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/Elevation/", NAMES_results[i],".txt", sep=""),quote=F,row.names=F,col.names=T,sep="\t")

}
