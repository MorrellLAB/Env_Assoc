# Generate a folded site frequency spectrum for 9k and exon-capture data.

# Read the data files

exonSNP <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/SFS/derived_alle_fre_60NAM_noNA.txt",header = T, sep="\t")

head(exonSNP)
arrarSNP <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/SFS/trans_myGD.V2.5749.DAF_noNA.txt",header = T, sep="\t")
head(arrarSNP)

# Bin them up in 2.5% frequency bins
breaks <- seq(0, 1, by=0.1)
exonSNP.sfs <- table(cut(as.numeric(exonSNP$DAF), breaks=breaks, include.lowest=TRUE))/nrow(exonSNP)
arrarSNP.sfs <- table(cut(as.numeric(arrarSNP$DAF), breaks=breaks, include.lowest=TRUE))/nrow(arrarSNP)
nrow(arrarSNP)
arrarSNP.sfs
# Make plot

pdf(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/SFS/Unfolded_SFS.pdf", width=10, height=16)
par(mfrow=c(2, 1))
barplot(t(arrarSNP.sfs),
        beside=TRUE,
        col=c("blue"),
        cex.lab=1.2, 
        cex.axis=1.2,
        xlab="Derived Allele Frequency",
        ylab="Proportion",
        main="Derived allele frequency spectrum in 9K SNPs",
        cex.main=1.4)
barplot(t(exonSNP.sfs),
        beside=TRUE,
        col=c("blue"),
        cex.lab=1.2, 
        cex.axis=1.2,
        xlab="Derived Allele Frequency",
        ylab="Proportion",
        ylim=c(0,0.5),
        main="Derived allele frequency spectrum in exon captured SNPs",
        cex.main=1.4)

dev.off()