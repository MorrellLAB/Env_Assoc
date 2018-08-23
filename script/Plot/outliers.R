#   Make a fancy plot of outliers for Fst and GWAS results:
#by Li Lei 2018/05/07

library(ggplot2)

#   Read the exon capture density data,9k genptyping data, and exon-capture SNP data;

Elv <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/Elevation/Outliers/Fst_less3000_vsmore3000_outliers_99th_sorted_physPos.txt",header = T, sep="\t")
Elv $Chromosome <- as.character(Elv$Chr_2016)
Elv <- na.omit(Elv)
head(Elv)
nrow(Elv)

Lat_30_40_above <- read.delim("/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/Latitude/Outliers/Fst_wild_range30_40_vs_higherLat40_outliers_99th_sorted_physPos.txt", header=T)
Lat_30_40_above$Chromosome <- as.character(Lat_30_40_above$Chr_2016)
Lat_30_40_above <- na.omit(Lat_30_40_above)
head(Lat_30_40_above)
str(Lat_30_40_above)
nrow(Lat_30_40_above)

Lat_30_40_less <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/Latitude/Outliers/Fst_moreoreq30to40vsless30_no_NA_outliers_99th_sorted_physPos.txt",header = T, sep="\t")
Lat_30_40_less$Chromosome <- as.character(Lat_30_40_less$Chr_2016)
Lat_30_40_less <- na.omit(Lat_30_40_less)
head(Lat_30_40_less)
nrow(Lat_30_40_less)

habit <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/GrowthHabit/Outliers/outlier_Fst_SPRING_vs_WINTER_80samples_average_noNA_PhyPos_99th.txt",header = T, sep="\t")
habit$Chromosome <- as.character(habit$Chr_2016)
habit <- na.omit(habit)
head(habit)
nrow(habit)
#habit <- na.omit(habit)
#threshold <- quantile(habit$Fst, probs = seq(0, 1, by= 0.005))[196] #take 97.5% as threhold
#0.1059391
#outliers <- habit[(habit$Fst >0.1059391),]
#write.table(outliers, file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/Fst/Results/AdegenetIndividual/FstLoci_Habit_outliers.txt",quote = F,sep = "\t", na = "NA",row.names = F, col.names = TRUE)


#habit <- habit[(habit$Fst >0.1059391),]



GWAS <- read.delim(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/GWAS-GAPIT/factor_categories/matrix_physPos.txt",header = T, sep="\t")
GWAS <- na.omit(GWAS)
GWAS$Chromosome <- as.character(GWAS$Chr_2016)
GWAS <- GWAS[GWAS$Chromosome != "chrUn",]
head(GWAS)

preci <- GWAS[(GWAS$precipitation==1),]
temp <- GWAS[(GWAS$temperature==1),]
geo <- GWAS[(GWAS$geographic==1),]

head(preci)
str(habit)
#habit$PhysPos_2016 <- as.numeric(habit$Chr_2016)
#   Drop the unmapped chromosome
Elv  <- Elv[Elv$Chromosome != "chrUn",]
Lat_30_40_above <- Lat_30_40_above[Lat_30_40_above$Chromosome != "chrUn",]
Lat_30_40_less <- Lat_30_40_less[Lat_30_40_less$Chromosome != "chrUn",]
habit <- habit[habit$Chromosome != "chrUn",]
GWAS <- GWAS[GWAS$Chromosome != "chrUn",]

head(GWAS)
#exonSNP$Position <- as.numeric(exonSNP$Position)

###read the genes

HEout <- read.delim("/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/outliers/High_elevation.txt",header=T,sep="\t")
head(HEout)
str(HEout)

HLout <- read.delim("/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/outliers/High_latitude.txt",header=T,sep="\t")
head(HLout)
str(HLout)

LLout <- read.delim("/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/outliers/Low_latitude.txt",header=T,sep="\t")
head(LLout)
str(LLout)

GHout <- read.delim("/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/outliers/Growth_habit.txt",header=T,sep="\t")
head(GHout)

str(GHout)

geosig <- read.delim("/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/outliers/geo_ass.txt",header=T,sep="\t")
head(geosig)
str(geosig)

tempsig <- read.delim("/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/outliers/temp_ass.txt",header=T,sep="\t")
head(tempsig)
str(tempsig)

humsig <- read.delim("/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/outliers/hum_ass.txt",header=T,sep="\t")
head(humsig)
str(humsig)

##Adjust a little! 

pdf(file="/Users/lilei/Dropbox (Morrell Lab)/Landrace_Environmental_Assocation/Analyses/outliers/outlier_99th.pdf", 10, 6)
ggplot(Elv) +
  geom_point(aes(x=PhysPos_2016/1000000, y=5), data=Elv, shape=21,fill="#0000FF", color="#0000FF", size=1.5, alpha=0.5)+
  geom_point(aes(x=PhysPos_2016/1000000, y=10), data=Lat_30_40_above, shape=21,fill="#0000FF", color="#0000FF", size=1.5, alpha=0.5) +
  geom_point(aes(x=PhysPos_2016/1000000, y=15), data=Lat_30_40_less, shape=21,fill="#0000FF", color="#0000FF", size=1.5, alpha=0.5) +
  geom_point(aes(x=PhysPos_2016/1000000, y=20), data=habit, shape=21,fill="#0000FF", color="#0000FF", size=1.5, alpha=0.5) +
  geom_point(aes(x=PhysPos_2016/1000000, y=25), data=preci, shape=21, fill="#0000FF", color="#0000FF", size=1.5, alpha=0.5) +
  geom_point(aes(x=PhysPos_2016/1000000, y=30), data=temp, shape=21, fill="#0000FF", color="#0000FF", size=1.5, alpha=0.5) +
  geom_point(aes(x=PhysPos_2016/1000000, y=35), data=geo, shape=21, fill="#0000FF", color="#0000FF", size=1.5, alpha=0.5) +
  geom_point(aes(x=PhysPos_2016/1000000, y=5), data=HEout, shape=21, fill="#17202A", color="#17202A", size=1.5, alpha=1) +
  geom_point(aes(x=PhysPos_2016/1000000, y=10), data=HLout, shape=21,fill="#17202A", color="#17202A", size=1.5, alpha=1) +
  geom_point(aes(x=PhysPos_2016/1000000, y=15), data=LLout, shape=21,fill="#17202A", color="#17202A", size=1.5, alpha=1) +
  geom_point(aes(x=PhysPos_2016/1000000, y=20), data=GHout, shape=21,fill="#17202A", color="#17202A", size=1.5, alpha=1) +
  geom_point(aes(x=PhysPos_2016/1000000, y=25), data=humsig, shape=21, fill="#17202A", color="#17202A", size=1.5, alpha=1) +
  geom_point(aes(x=PhysPos_2016/1000000, y=30), data=tempsig, shape=21, fill="#17202A", color="#17202A", size=1.5, alpha=1) +
  geom_point(aes(x=PhysPos_2016/1000000, y=35), data=geosig, shape=21, fill="#17202A", color="#17202A", size=1.5, alpha=1) +
  facet_grid(Chromosome~.) +
  scale_y_continuous(limits=c(0, 40), breaks=seq(0,40,by=5)) +
  scale_x_continuous(limits=c(0, 725), breaks=seq(0, 725, by=50)) +
  theme_bw() +
  theme(
    strip.background=element_blank(),
    axis.text.y = element_blank(),
    strip.text.y=element_text(size=10, colour="black", angle=0)
  ) +
  labs(y="", x="Physical Position (Mb)")

dev.off()




