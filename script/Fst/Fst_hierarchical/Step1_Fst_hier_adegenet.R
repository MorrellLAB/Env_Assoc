#Author: Ana Poets
# Description: Allele frequency differenciation in a hierarchical fashion for latitude within
# growth habit and elevation within growth habit
#################################################################################################
rm(list=ls())
#install.packages("adegenet",dep=TRUE)
library(hierfstat)
library(adegenet)
library(pegas)
# === IMPORT INPUT FILES ===========================================================================
#Load the growth type assignations 
GH<-read.csv("~/Dropbox/GITHUB/Barley_Landrace_EnvAssoc/Fst/Fst_byGrowthType/784_landrace_habit.csv",header=F,sep=",")
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
  dat[which(dat == "AA")]<-"2/2"
  dat[which(dat == "BB")]<-"1/1"
  dat[which(dat == "AB")]<-"1/2"
  return(dat)	
}
GENOTYPE_num<-as.data.frame(apply(genotypes_lev_or,2, CHANGE))



#================================= Calculate Fst for Latitude WITHIN GH ====================================================
# Combine Genotypes with levels
str<-cbind(LEVELS_or[,c(1:2)],GENOTYPE_num)
dim(str)
write.table(str,"~/Desktop/temp.txt",quote=F,row.names=T,col.names=T,sep="\t")
str<-read.table("~/Desktop/temp.txt",header=T)
strpop <- str[,c(1,2)] #pull out the strata

testStr <- df2genind(str[,-c(1,2)],sep="/",ploidy=2,pop =str$LatitudeLevels, strata = strpop)
fstat_out<-fstat(testStr)
Fst_loci<-Fst(as.loci(testStr))

write.table(fstat_out,"~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/Hierarchical_adegenet/Fstat_out_Lat_withinGH.txt",quote=F,row.names=T,col.names=T,sep="\t")
write.table(Fst_loci,"~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/Hierarchical_adegenet/FstLoci_out_Lat_withinGH.txt",quote=F,row.names=T,col.names=T,sep="\t")

#================================= Calculate Fst for Elevation WITHIN GH ====================================================
# Combine Genotypes with levels
str<-cbind(LEVELS_or[,c(1,3)],GENOTYPE_num)
dim(str)
write.table(str,"~/Desktop/temp.txt",quote=F,row.names=T,col.names=T,sep="\t")
str<-read.table("~/Desktop/temp.txt",header=T)
strpop <- str[,c(1,2)] #pull out the strata

testStr <- df2genind(str[,-c(1,2)],sep="/",ploidy=2,pop =str$ElevationLevels, strata = strpop)
fstat_out<-fstat(testStr)
Fst_loci<-Fst(as.loci(testStr))

write.table(fstat_out,"~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/Hierarchical_adegenet/Fstat_out_elev_withinGH.txt",quote=F,row.names=T,col.names=T,sep="\t")
write.table(Fst_loci,"~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/Hierarchical_adegenet/FstLoci_out_elev_withinGH.txt",quote=F,row.names=T,col.names=T,sep="\t")


##============= Add SNP information to the output files =============================================================

NAMES_results<-c("FstLoci_out_Lat_withinGH","FstLoci_out_elev_withinGH")
for ( i in 1:length(NAMES_results)){
  DATA<-read.table(paste("~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/Hierarchical_adegenet/", NAMES_results[i],".txt", sep=""), header=T)
  # combine SNP name with Fst value
  RESULTS<-(cbind(row.names(DATA),as.numeric(DATA[,2])))
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
  write.table(Results_all_genmap, paste("~/Dropbox/Landrace_Environmental_Assocation/Analyses/Fst/Results/Hierarchical_adegenet/", NAMES_results[i],"_SNPorder.txt", sep=""),quote=F,row.names=F,col.names=T,sep="\t")
  
}