#   Plot the P-values versus the allele frequency
#by Li Lei: 03/19/2018
#   We want to make beeswarm plots
library(beeswarm)

#read the files:
bio14 <- read.csv("/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/SNP_density/GAPIT_bio14_GWAS_results_physPos_MAF.csv", header=T)
head(bio14)

# Assign allele frequency 
cyc <- sapply(
  bio14$Maf,
  function(x) {
    if(x>=0 && x<=0.025) {
      return("0-0.025")
    }
    else if(x>0.025 && x<=0.05) {
      return("0.025-0.05")
    }
    else if(x>0.05 && x<=0.075) {
      return("0.05-0.075")
    }
    else if(x>0.075 && x<=0.1) {
      return("0.075-0.1")
    }
    else if(x>0.1 && x<=0.125) {
      return("0.1-0.125")
    }
    else if(x>0.125 && x<=0.15) {
      return("0.125-0.15")
    }
    else if(x>0.15 && x<=0.175) {
      return("0.15-0.175")
    }
    else if(x>0.175 && x<=0.2) {
      return("0.175-0.2")
    }
    else if(x>0.2 && x<=0.225) {
      return("0.2-0.225")
    }
    else if(x>0.225 && x<=0.25) {
      return("0.225-0.25")
    }
    else if(x>0.25 && x<=0.275) {
      return("0.25-0.275")
    }
    else if(x>0.275 && x<=0.3) {
      return("0.275-0.3")
    }
    else if(x>0.3 && x<=0.325) {
      return("0.3-0.325")
    }
    else if(x>0.325 && x<=0.35) {
      return("0.325-0.35")
    }
    else if(x>0.35 && x<=0.375) {
      return("0.35-0.375")
    }
    else if(x>0.375 && x<=0.4) {
      return("0.375-0.4")
    }
    else if(x>0.4 && x<=0.425) {
      return("0.4-0.425")
    }
    else if(x>0.425 && x<=0.45) {
      return("0.425-0.45")
    }
    else if(x>0.45 && x<=0.475) {
      return("0.45-0.475")
    }
    else {
      return("0.475-0.5")
    }
  }
)

bio14$MAF <- factor(cyc, levels=c("0-0.025", "0.025-0.05", "0.05-0.075", "0.075-0.1","0.1-0.125","0.125-0.15","0.15-0.175","0.175-0.2","0.2-0.225","0.225-0.25","0.25-0.275","0.275-0.3","0.3-0.325","0.325-0.35","0.35-0.375","0.375-0.4","0.4-0.425","0.425-0.45","0.45-0.475","0.475-0.5"))
bio14$log_p <- -log10(bio14$P.value)

str(bio14)
head(bio14)
pdf(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/maf_p_value/MAF_P_bio14.pdf",25,10)
library(beeswarm)
beeswarm(
  log_p ~ MAF,
  data=bio14,
  col=sample(colors(), 19),
  pch=19,
  cex=0.35,
  method="hex",
  ylim=c(0, 10),
  xlab="MAF",
  ylab="-log10(P-value)",
  main="BIO14",
  axes=T)
abline(h=3.30103, lty=3,col="black")

dev.off()




#break the frequency into several bins:
breaks <- seq(0, 0.5, by=0.025)
bio14.sfs <- table(cut(as.numeric(bio14$Maf), breaks=breaks, include.lowest=TRUE))/nrow(bio14)
head(bio14.sfs)

# Read the dosages
dos <- read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Results/Summaries/GP_Deleterious_TotalDosages.txt", header=TRUE)
# Random lines
rand <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Pedigrees/Random_Lines.txt", header=FALSE)$V1)
# Selected lines
selected <- as.character(read.table("/Users/tomkono/Dropbox/GitHub/Deleterious_GP/Data/Pedigrees/Selected_Lines.txt", header=FALSE)$V1)


#   Assign a color to plot the random, selected, and 'none' lines
#   They will be black, red, and grey, respectively.
type <- sapply(
  dos$LineID,
  function(x) {
    if(as.character(x) %in% selected) {
      return("Sel")
    }
    else if(as.character(x) %in% rand) {
      return("Ran")
    }
    else {
      return("Non")
    }
  }
)

dos$Type <- type

# Assign a cycle
cyc <- sapply(
  dos$LineID,
  function(x) {
    if(grepl("MS10", x)) {
      return("C1")
    }
    else if(grepl("MS11", x)) {
      return("C2")
    }
    else if(grepl("MS12", x)) {
      return("C3")
    }
    else {
      return("C0")
    }
  }
)

dos$Cycle <- factor(cyc, levels=c("C0", "C1", "C2", "C3"))

pdf(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/maf_p_value/MAF_P_bio14.pdf", height=6, width=6)

beeswarm(
  non$Dosage ~ non$Cycle,
  col="#cccccc",
  pch=19,
  cex=0.35,
  method="hex",
  ylim=c(400, 800),
  xlab="Cycle",
  ylab="Total Dosage of Deleterious Alleles",
  main="",
  axes=F)



#   Separate the data now, so we can clearly show the random and selected lines
ran <- dos[dos$Type == "Ran",]
sel <- dos[dos$Type == "Sel",]
non <- dos[dos$Type == "Non",]

pdf(file="DM_By_Cycle.pdf", height=6, width=6)
par(mar=c(4, 4, 0.1, 0.1), mgp=c(2, 1, 0))
#   Plot the not selected nor random lines
beeswarm(
  non$Dosage ~ non$Cycle,
  col="#cccccc",
  pch=19,
  cex=0.35,
  method="hex",
  ylim=c(400, 800),
  xlab="Cycle",
  ylab="Total Dosage of Deleterious Alleles",
  main="",
  axes=F)
beeswarm(
  ran$Dosage ~ ran$Cycle,
  col="#333333",
  pch=19,
  cex=0.5,
  method="hex",
  ylim=c(400, 800),
  add=TRUE,
  side=-1,
  axes=F)
beeswarm(
  sel$Dosage ~ sel$Cycle,
  col="#aa0000",
  pch=19,
  cex=0.5,
  method="hex",
  ylim=c(400, 800),
  add=TRUE,
  side=1,
  axes=F)
boxplot(
  Dosage~Cycle + Type,
  data=droplevels(rbind(ran, sel)),
  at=c(1.75, 2.75, 3.75, 2.25, 3.25, 4.25),
  boxwex=0.2,
  lwd=1,
  border=c("black", "black", "black", "red", "red", "red"),
  axes=F,
  add=TRUE)
legend("topright", pch=19, col=c("black", "red"), legend=c("Random", "Selected"))
axis(
  side=1,
  at=c(1, 2, 3, 4),
  labels=c("Parents", "C1", "C2", "C3"))
axis(
  side=2)
box()
dev.off()