#Author: Ana Poets
# Description: Allele frequency differenciation among the landraces at three different elevations
#################################################################################################
rm(list=ls())

library(hierfstat)

# === IMPORT INPUT FILES ===========================================================================
#Load the growth type assignations 
GH<-read.csv("/Users/Mia/Dropbox/GITHUB/Barley_Landrace_EnvAssoc/Fst/Fst_byGrowthType/784_landrace_habit.csv",header=F,sep=",")
# Select new geographic location, after Fumi inspected every data point. These are the phenotypes 
# Fumi used for GAPIT
Latlong<-read.csv("~/Dropbox/Landrace_Environmental_Assocation/Analyses/GWAS-GAPIT/Input/myY1.v2.csv",row.names=1)


# Import genotypes. Downloaded from Github companion of Poets et al 2015.
genotypes<-read.table("~/Dropbox/Landrace_Environmental_Assocation/Data/Poets_etal2015/Land_6152_SNPs_AB.txt",row.names=1, header=T)

# Import genetic map from Munñoz et al 2011
genMap<-read.table("~/Dropbox/Landrace_Environmental_Assocation/Data/Muñoz_etal2011/Consensus_iSelect.txt",header=T)
################ Combine lat, long, altitude and GH levels #########################################
Samples_geoInfo<-intersect(GH[,1],row.names(Latlong))

GH_sh<-GH[(GH[,1] %in% Samples_geoInfo ),]
Latlong_sh<-Latlong[(row.names(Latlong) %in% Samples_geoInfo),1:3]

# sort GH_sh by latlong samples
GH_sh_or<-GH_sh[match(row.names(Latlong_sh),GH_sh[,1]),]

if (identical(row.names(Latlong_sh), as.character(GH_sh_or[,1])) == FALSE)stop('GH and LatLong have samples in different order. Cannot merge!')

# Join variables
Latlong_sh_GH<-cbind(Latlong_sh,GH_sh_or[,2])
names(Latlong_sh_GH)[4]<-"GrowthHabit"

# Remove GH other than W,S,F
Latlong_sh_GH_complete<-Latlong_sh_GH[which(Latlong_sh_GH[,4] == "S" | Latlong_sh_GH[,4] == "W"|Latlong_sh_GH[,4]=="F"),]

# Convert to levels
ConvertGH<-function(dat){
  dat[which(dat == "S")]<-"1"
  dat[which(dat == "W")]<-"2"
  dat[which(dat == "F")]<-"3"
  return(dat)
}

ConvertLatitude<-function(dat){
  dat[which(dat <30)]<-"1"
  dat[which(dat >= 30 & dat <=40)]<-"2"
  dat[which(dat >40)]<-"3"
  return(dat)
}

ConvertElevation<-function(dat){
  dat[which(dat <3000)]<-"1"
  dat[which(dat >=3000)]<-"2"
  return(dat)
}

GHLevels<-ConvertGH(as.character(Latlong_sh_GH_complete[,4]))
LatitudeLevels<-ConvertLatitude(Latlong_sh_GH_complete[,1])
ElevationLevels<-ConvertElevation(Latlong_sh_GH_complete[,3])

LEVELS<-cbind(GHLevels,LatitudeLevels,ElevationLevels)
row.names(LEVELS)<-row.names(Latlong_sh_GH_complete)
# sort levels by GH
LEVELS_or<-LEVELS[order(LEVELS[,1]),]

# Get the genotypes only for the samples for which we have all the geographic and GH info.
genotypes_lev<-genotypes[(row.names(genotypes) %in% row.names(LEVELS_or)),]
# order genotypes as in Levels
genotypes_lev_or<-genotypes_lev[match(row.names(LEVELS_or),row.names(genotypes_lev)),]

if (identical(row.names(genotypes_lev_or),row.names(LEVELS_or))==FALSE)stop("Levels and genotypes are in different order")
# ==== DATA TRANFORMATIONS ========================================================================
# Change genotypes to numeric values AA=2, BB=0, AB=1 ,NA=NA
CHANGE<-function(dat){
  dat[which(dat == "AA")]<-2
  dat[which(dat == "BB")]<-0
  dat[which(dat == "AB")]<-1
  return(dat)	
}
GENOTYPE_num<-as.data.frame(apply(genotypes_lev_or,2, CHANGE))

# Calculate Fst for Latitude WITHIN GH


#look at the Fst at each locus
Fst<-NULL
for (i in 1:dim(GENOTYPE_num)[[2]]) {
  #Fst[i]<-varcomp(data.frame(LEVELS_or[,1],as.numeric(GENOTYPE_num[,i])),diploid=FALSE)$F[1,1]
  Fst<-c(Fst,mean(test.within(data.frame(GENOTYPE_num[,i]),test=LEVELS_or[,2],within=LEVELS_or[,1],nperm=100,diploid=FALSE)$g))   
}
Fst_Latitude_withinGH<-as.data.frame(Fst) 
row.names(Fst_Latitude_withinGH)<-colnames(GENOTYPE_num)

write.table(Fst_Latitude_withinGH, "~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/GrowthHabit/Fst_Latitude_withinGH.txt",quote=F,row.names=T,col.names=T,sep="\t")

# Fst in Elevation within Growth habit
Fst2<-NULL
for (i in 1:dim(GENOTYPE_num)[[2]]) {
  #Fst[i]<-varcomp(data.frame(LEVELS_or[,1],as.numeric(GENOTYPE_num[,i])),diploid=FALSE)$F[1,1]
  Fst2<-c(Fst2,mean(test.within(data.frame(GENOTYPE_num[,i]),test=LEVELS_or[,3],within=LEVELS_or[,1],nperm=100,diploid=FALSE)$g))   
}
Fst_Elevation_withinGH<-as.data.frame(Fst2) 
row.names(Fst_Elevation_withinGH)<-colnames(GENOTYPE_num)

write.table(Fst_Elevation_withinGH, "~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/GrowthHabit/Fst_Elevation_withinGH.txt",quote=F,row.names=T,col.names=T,sep="\t")

##============= Add SNP information to the output files =============================================================

NAMES_results<-c("Fst_Elevation_withinGH")
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
  colnames(Results_all_genmap)<-c("Chromosome","Cumulative_cM","SNP","G")
  write.table(Results_all_genmap, paste("~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/GrowthHabit/", NAMES_results[i],".txt", sep=""),quote=F,row.names=F,col.names=T,sep="\t")
  
}
