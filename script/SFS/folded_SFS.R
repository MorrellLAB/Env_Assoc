# Generate a folded site frequency spectrum for 9k and exon-capture data.

# Read the data files

exonSNP <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/LD/data/OnlyLandrace_biallelic_Barley_NAM_Parents_Final_renamed.maf",header = T, sep="\t")

head(exonSNP)
arrarSNP <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Contributors/Fumi/Fumi_Env_GWAS/GWAS.landraces/trans_myGD.v2.sorted.maf",header = T, sep="\t")
head(arrarSNP)


# Bin them up in 2.5% frequency bins
breaks <- seq(0, 0.5, by=0.05)
exonSNP.sfs <- table(cut(as.numeric(exonSNP$MAF), breaks=breaks, include.lowest=TRUE))/nrow(exonSNP)
arrarSNP.sfs <- table(cut(as.numeric(arrarSNP$Maf), breaks=breaks, include.lowest=TRUE))/nrow(arrarSNP)
nrow(arrarSNP)
arrarSNP.sfs
# Make plot

pdf(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/SFS/Folded_SFS.pdf", width=10, height=16)
par(mfrow=c(2, 1))
barplot(t(arrarSNP.sfs),
        beside=TRUE,
        col=c("blue"),
        cex.lab=1.2, 
        cex.axis=1.2,
        xlab="Minor Allele Frequency",
        ylab="Proportion",
        cex.main=1.4,
        main="Folded SFS in 9K SNPs")
barplot(t(exonSNP.sfs),
        beside=TRUE,
        col=c("blue"),
        xlab="Minor Allele Frequency",
        ylab="Proportion",
        cex.lab=1.2, 
        cex.axis=1.2,
        ylim=c(0,0.5),
        cex.main=1.4,
        main="Folded SFS in exon captured SNP")

dev.off()
