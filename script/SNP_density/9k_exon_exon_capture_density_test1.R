#   Make a fancy plot of9k and exon SNPs
#by Li Lei 2018/03/15

library(ggplot2)

#   Read the exon capture density data,9k genptyping data, and exon-capture SNP data;

excap <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/SNP_density/ExomeCaptureTargets_per_Mb.txt",header = T, sep="\t")
head(excap)
nrow(excap)

arraySNP <- read.csv("/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/SNP_density/GAPIT_bio6_GWAS_results_physPos.csv", header=T)
head(arraySNP)

nrow(arraySNP)

exonSNP <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/SNP_density/sorted_11_10380_Only_landrace_biallelic_NAM_final.pos",header = T, sep="\t")
head(exonSNP)
nrow(exonSNP)

rec_rate <- read.table("/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/cM_Mb/out/9k_cM_2Mb_Smoothed.txt", header=T)
head(rec_rate)

#   Drop the unmapped chromosome
arraySNP <- arraySNP[arraySNP$Chromosome != "chrUn",]
exonSNP <- exonSNP[exonSNP$Chromosome != "chrUn",]
excap <- excap[excap$Chromosome != "chrUn",]
rec_rate <- rec_rate[rec_rate$Chromosome != "chrUn",]

get_num_excap <- function(snp) {
  chrom <- as.character(snp["Chromosome"])
  pos <- as.numeric(snp["Position"])
  sel_row <- which((excap$Chromosome == chrom) & (pos > excap$Start) & (pos <= excap$End))
  n_excap <- excap[sel_row, "NExCap"]
  if(length(n_excap) == 0) {
    return(NA)
  }
  return(n_excap)
}

exonSNP$Position <- as.numeric(exonSNP$Position)

head(exonSNP)

get_recomb_rate <- function(snp) {
  chrom <- as.character(snp["Chromosome"])
  pos <- as.numeric(snp["Position"])
  sel_row <- which((rec_rate$Chromosome == chrom) & (pos > rec_rate$LeftBP) & (pos < rec_rate$RightBP))
  rrate <- rec_rate[sel_row, "Smoothed_cMMb"]
  #   Set recombination rates that are too high to NA. Same with those that
  #   do not match any given interval
  if(length(rrate) == 0) {
    return(NA)
  }
  if(is.na(rrate)) {
    return(NA)
  }
  if(rrate > 10) {
    return(NA)
  }
  return(rrate)
}


## creat a full picture fo all chromosomes! Adjust a little! 

pdf(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/SNP_density/exon_9k_SNPs_2M.pdf", 10, 6)

ggplot(exonSNP) +
geom_vline(aes(xintercept=Position/1000000), size=0.02, alpha=0.1, color="#a6cee3") +
  geom_line(aes(x=(Start+End)/2000000, y=NExCap/10), data=excap, size=0.75, color="#1f78b4", alpha=0.7) +
  geom_line(aes(x=(LeftBP+RightBP)/2000000, y=Smoothed_cMMb), data=rec_rate, color="#9900ff", size=1, alpha=1) +
  geom_point(aes(x=PhysPos_2016/1000000, y=15), data=arraySNP, shape=17, color="#ff4000", size=1, alpha=0.45) +
  facet_grid(Chromosome~.) +
  scale_y_continuous(limits=c(0, 15), breaks=seq(0,15,by=5)) +
  scale_x_continuous(limits=c(0, 725), breaks=seq(0, 725, by=50)) +
  theme_bw() +
  theme(
    strip.background=element_blank(),
    strip.text.y=element_text(size=10, colour="black", angle=0),
    axis.ticks.y=element_blank()) +
  labs(y="", x="Physical Position (Mb)")

dev.off()
       


#geom_point(aes(x=PhysPos_2016/1000000,y=-log10(P.value)), data=arraySNP, color="#ff4000", size=0.50, alpha=0.45) +
###We need to split chr3H out and plot it as the main text figure and then put the rest of them as supplemental figures:

#   extract 3H only:
arraySNP_3H <- arraySNP[arraySNP$Chromosome == "chr3H",]
exonSNP_3H <- exonSNP[exonSNP$Chromosome == "chr3H",]
excap_3H <- excap[excap$Chromosome == "chr3H",]
rec_rate_3H <- rec_rate[rec_rate$Chromosome == "chr3H",]
  
pdf(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/SNP_density/exon_9k_SNPs_3H#.pdf", 10, 3)

ggplot(exonSNP_3H) +
  geom_vline(aes(xintercept=Position/1000000), size=0.02, alpha=0.1, color="#a6cee3") +
  geom_line(aes(x=(Start+End)/2000000, y=NExCap/10), data=excap_3H, size=0.75, color="#1f78b4", alpha=0.7) +
  geom_line(aes(x=(LeftBP+RightBP)/2000000, y=Smoothed_cMMb), data=rec_rate_3H, color="#9900ff", size=1, alpha=1) +
  geom_point(aes(x=PhysPos_2016/1000000,y=-log10(P.value)), data=arraySNP_3H, color="#ff4000", size=0.50, alpha=0.45) +
  geom_hline(yintercept = 3.30103,size=0.25,color="#000000",linetype="dotted") +
  facet_grid(Chromosome~.) +
  scale_y_continuous(limits=c(0, 10), breaks=c(0, 5, 10)) +
  scale_x_continuous(limits=c(0, 725), breaks=seq(0, 725, by=50)) +
  theme_bw() +
  theme(
    strip.background=element_blank(),
    strip.text.y=element_text(size=10, colour="black", angle=0),
    axis.ticks.y=element_blank()) +
  labs(y="", x="Physical Position (Mb)")

dev.off()

#   extract the rest of chrs only:
arraySNP_rest <- arraySNP[arraySNP$Chromosome != "chr3H",]
exonSNP_rest <- exonSNP[exonSNP$Chromosome != "chr3H",]
excap_rest <- excap[excap$Chromosome != "chr3H",]
rec_rate_rest <- rec_rate[rec_rate$Chromosome != "chr3H",]

pdf(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/SNP_density/exon_9k_SNPs_rest.pdf#", 10, 6)

ggplot(exonSNP_rest) +
  geom_vline(aes(xintercept=Position/1000000), size=0.02, alpha=0.1, color="#a6cee3") +
  geom_line(aes(x=(Start+End)/2000000, y=NExCap/10), data=excap_rest, size=0.75, color="#1f78b4", alpha=0.7) +
  geom_line(aes(x=(LeftBP+RightBP)/2000000, y=Smoothed_cMMb), data=rec_rate_rest, color="#9900ff", size=1, alpha=1) +
  geom_point(aes(x=PhysPos_2016/1000000,y=-log10(P.value)), data=arraySNP_rest, color="#ff4000", size=0.50, alpha=0.45) +
  geom_hline(yintercept = 3.30103,size=0.25,color="#000000",linetype="dotted") +
  facet_grid(Chromosome~.) +
  scale_y_continuous(limits=c(0, 10), breaks=c(0, 5, 10)) +
  scale_x_continuous(limits=c(0, 725), breaks=seq(0, 725, by=50)) +
  theme_bw() +
  theme(
    strip.background=element_blank(),
    strip.text.y=element_text(size=10, colour="black", angle=0),
    axis.ticks.y=element_blank()) +
  labs(y="", x="Physical Position (Mb)")

dev.off()
