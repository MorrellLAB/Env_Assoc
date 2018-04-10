#   Plot the P-values versus the allele frequency
#by Li Lei: 03/19/2018
#   We want to make beeswarm plots
library(beeswarm)

#read the files:
bio14 <- read.delim("/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/maf_p_value/GAPIT_bio14_GWAS_results_physPos_MAF_DAF.txt",sep="\t", header=T)
head(bio14)

# Assign allele frequency 
cyc <- sapply(
  bio14$DAF,
  function(x) {
    if(x>=0 && x<=0.1) {
      return("0-0.1")
    }
    else if(x>0.1 && x<=0.2) {
      return("0.1-0.2")
    }
    else if(x>0.2 && x<=0.3) {
      return("0.2-0.3")
    }
    else if(x>0.3 && x<=0.4) {
      return("0.3-0.4")
    }
    else if(x>0.4 && x<=0.5) {
      return("0.4-0.5")
    }
    else if(x>0.5 && x<=0.6) {
      return("0.5-0.6")
    }
    else if(x>0.6 && x<=0.7) {
      return("0.6-0.7")
    }
    else if(x>0.7 && x<=0.8) {
      return("0.7-0.8")
    }
    else if(x>0.8 && x<=0.9) {
      return("0.8-0.9")
    }
    else {
      return("0.9-1")
    }
  }
)

bio14$DAF <- factor(cyc, levels=c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4","0.4-0.5","0.5-0.6","0.6-0.7","0.7-0.8","0.8-0.9","0.9-1"))
bio14$log_p <- -log10(bio14$P.value)

str(bio14)
head(bio14)
pdf(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/maf_p_value/DAF_P_bio14.pdf",20,9)
library(beeswarm)
beeswarm(
  log_p ~ DAF,
  data=bio14,
  pch=21,
  cex=0.45,
  col="#626567",
  method="hex",
  ylim=c(0, 10),
  xlab="Derived allele frequency",
  ylab="-log10(P-value)",
  main="BIO14",
  axes=T)
abline(h=3.30103, lty=3,col="black")

dev.off()



