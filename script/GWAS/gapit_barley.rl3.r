#### 151211. myY's 774 rows were not the same as MyGD.
#### redo it.

#setwd("C:\\run_space\\myGAPIT\\barley.landraces")

#### prepare the data for GAPIT
## phenotype data
coor.dat = read.delim('coor.dat.txt',header=T,row.names=1,as.is=c(2,6))

## dealing with missing pheno data
coor.x.dat = coor.dat[!complete.cases(coor.dat) & !is.na(coor.dat$Latitude),]
# check on Google Earth and change the coordinates
#write.table(coor.x.dat[,1:5],file='missing.pheno.data.point.csv',sep=',',quote=F)
coor.ch = read.csv('coor.fix.dat.csv',header=T,row.names=1)
coor.dat[!complete.cases(coor.dat) & !is.na(coor.dat$Latitude),] = coor.ch

## ICA
library(ica)
coor1.dat = coor.dat[complete.cases(coor.dat),-c(1:5)]

## standardize and look at the distribution of each bioclim
env.ds = apply(coor1.dat, 2, function(x){
  mean.v = mean(x)
  sd.v = sd(x)
  (x - mean.v)/sd.v
})
ic.env = icaimax(env.ds,3)
i.mat1 = ic.env$R
env.icx = cbind(env.ds,ic.env$S)
round(cor(env.icx)^2,digits=3)

coor.f.dat = data.frame(Taxa=rownames(coor.dat[complete.cases(coor.dat),]),
                        coor.dat[complete.cases(coor.dat),], round(ic.env$S,digits=3))
sum(complete.cases(coor.f.dat))
myY1 = coor.f.dat
myY1 = myY1[,-c(2,6)]
colnames(myY1)[24:26]=c('IC1','IC2','IC3')
write.table(myY1,'myY1.csv',sep=',',quote=F,row.names=F)


## genotype data
## SNPs with unknown positions are assigned to chr 8 with arbitrary postions
geno.d = read.delim('Land_6152_SNPs_AB.txt',header=T,row.names=1)
gmap.d = read.delim('GeneticMap_iSelect_9k.txt',header=T)
geno.n = colnames(geno.d)
gmap.n = as.character(gmap.d[,'SNP'])
gn.com = intersect(geno.n,gmap.n)
gn.uncom = setdiff(geno.n, gn.com) 

rownames(gmap.d)=gmap.d[,'SNP']
gmap.d1 = gmap.d[gn.com,c('Chro','cm','index')]
gmap.d1 = gmap.d1[order(gmap.d1[,'index']),]

gmap.du = data.frame(Chro=paste(rep(8,length(gn.uncom)),'H',sep=''),
                     cm=1:length(gn.uncom)/50,
                     index=(nrow(gmap.d1)+1):(nrow(gmap.d1)+length(gn.uncom)))
rownames(gmap.du)=gn.uncom
gmap.d1 = rbind(gmap.d1,gmap.du)

myGM = data.frame(Name = rownames(gmap.d1), 
                     Chromosome=as.numeric(gmap.d1[,'Chro']),
                     Position = gmap.d1[,'cm']*1000)
gn.com = rownames(gmap.d1)

geno.dx = as.matrix(geno.d[,gn.com])
table(as.character(geno.dx), useNA='ifany')
#      AA      AB      BB    <NA> 
# 2192030    5246 2722471   20309 

geno.dx1 = matrix(NA,ncol=ncol(geno.dx),nrow=nrow(geno.dx))
geno.dx1[geno.dx=='AA']=0
geno.dx1[geno.dx=='BB']=2
colnames(geno.dx1)=colnames(geno.dx)
rownames(geno.dx1)=rownames(geno.dx)
geno.dx1 = geno.dx1[as.character(myY1[,1]),]
table(as.character(geno.dx1), useNA='ifany')
#        0       2    <NA> 
#  2139846 2658623   24699 

### Remove the SNPs with > 20% NA
good.cols = apply(geno.dx1,2,function(x) sum(is.na(x)) < 0.2* length(x))
geno.dx1 = geno.dx1[,good.cols]

### Remove the samples with > 20% NA
good.rows = apply(geno.dx1,1,function(x) sum(is.na(x)) < 0.2* length(x))
geno.dx1 = geno.dx1[good.rows,]

### Remove the SNPs with the same genotype. Priorities:
### (1) mapped; (2) fewer NAs
skip.list = c()
keep.list = c()
agg.rep = list()
snp.names = colnames(geno.dx1)
for (i in snp.names){
  if (i %in% skip.list) next
  skip.list = c(skip.list, i)
  if (length(skip.list)==length(snp.names)){
    keep.list = c(keep.list, i)
    next
  }
  list.q = setdiff(snp.names, skip.list)
  geno.q = geno.dx1[,i]
  if (length(list.q)==1){
    g.comp = geno.q == list.q
    if (sum(g.comp==F, na.rm=T)){
      int.list = F
    } else {
      int.list = T
    }
  } else {
    int.list = apply(geno.dx1[,list.q],2,function(x){
      g.comp = geno.q == x
      if (sum(g.comp==F, na.rm=T)) {
        return(F)
      } else {
        return(T)
      }
    })
  }
  if (sum(int.list)==0) {
    keep.list = c(keep.list, i)
    next
  }
  overlap.list = cbind(geno.dx1[,i], (geno.dx1[,list.q])[,int.list])
  if (sum(int.list)==1){
    colnames(overlap.list)=c(i,colnames(geno.dx1[,list.q])[int.list])
  } else {
    colnames(overlap.list)[1]=i
  }
  oll.na = is.na(overlap.list)
  if (sum(oll.na)==0){
    keep.list = c(keep.list, i)
    skip.list = c(skip.list, colnames(overlap.list)[-1])
    agg.rep[[i]] = list(inp=colnames(overlap.list), outp=i)
    next
  }
  oll.na.co = apply(oll.na, 1, sum)
  oll.na = oll.na[as.logical(oll.na.co),]
  if (class(oll.na)=='logical'){
    oll.na.n = names(oll.na)
    oll.na = as.numeric(oll.na)
    names(oll.na) = oll.na.n
    oll.na = sort(oll.na)
    keep.list = c(keep.list, names(oll.na)[1])
    skip.list = c(skip.list, colnames(overlap.list)[-1])
    agg.rep[[i]] = list(inp=colnames(overlap.list), names(oll.na)[1])
    next
  } 
  oll.na.ord = apply(oll.na, 2, sum)    
  oll.na = oll.na[,order(oll.na.ord)]
  k1.list = colnames(oll.na)[1]
  for (j in 2:ncol(oll.na)){
    colq = oll.na[,j]
    colq.k = sapply(1:length(colq), function(x){
      if (colq[x]==T){
        return(F)
      } else {
        if (sum(!oll.na[x,1:(j-1)])){
          return(F)
        } else {
          return(T)
        }
      }
    })
    if (sum(colq.k)) {
      k1.list = c(k1.list, colnames(oll.na)[j])
    }
  }
  keep.list = c(keep.list, k1.list)
  skip.list = c(skip.list, colnames(overlap.list)[-1])
  agg.rep[[i]] = list(inp=colnames(overlap.list), outp=k1.list)
}
date() #~14 min

### sanity check
sum(!(keep.list %in% snp.names))
length(keep.list)
#[1] 5818
length(snp.names)
#[1] 6152
## 334 snps were removed

sum(keep.list %in% gn.uncom)
#[1] 1548
length(gn.uncom)
#[1] 1625
## among removed 77 were unmapped

sum(is.na(geno.dx1[,keep.list]))/prod(dim(geno.dx1[,keep.list]))
#[1] 0.005173292
sum(is.na(geno.dx1))/prod(dim(geno.dx1))
#[1] 0.005120908
## almost no change in NA percentage

snp.order = 1:length(snp.names)
names(snp.order) = snp.names
# up to 4527 are mapped
to.check = list()
for (i in names(agg.rep)){
  q.list = agg.rep[[i]]$inp
  q.numb = snp.order[q.list]
  q.numb = q.numb[q.numb <=4527]
  if (length(q.numb) < 2) next
  q.chr = myGM[q.numb,'Chromosome']
  q.pos = myGM[q.numb, 'Position']
  if (length(unique(q.chr)) > 1 | max(q.pos) - min(q.pos) > 2000){
    to.check[[i]] = myGM[q.numb,]
  }
}

to.check
#$X12_21131
#Name Chromosome Position
#249       X12_21131          1    60140
#252  SCRI_RS_145336          1    61220
#445       X12_20429          1   142740
#1761      X12_10680          3   105980
#3617      X12_11181          6    59990
#3665      X12_11321          6    65380
#3680      X12_10392          6    70890
#4141 SCRI_RS_169639          7    75740

#$X12_21396
#Name Chromosome Position
#1023 X12_21396          2   126230
#1521 X12_20863          3    62970

#$X12_30232
#Name Chromosome Position
#2314      X12_30232          4    90060
#2366      X12_11194          4   107770
#3518 SCRI_RS_137706          6    58810
#4143      X12_21167          7    76060

#$SCRI_RS_150768
#Name Chromosome Position
#4279 SCRI_RS_150768          7    80470
#4308 SCRI_RS_207246          7    83790

## "X12_21131","X12_21396","X12_30232" are problematic
# "X12_21131", 249, 445, 1761, 3617, 3665, 3680, 4141
comp.set.n = c(249, 445, 1761, 3617, 3665, 3680, 4141)
comp.set = geno.dx1[,comp.set.n]
cor(comp.set)
table(comp.set[,1]==comp.set[,2],useNA='ifany')
apply(comp.set,2,function(x) sum(x==0,na.rm=T))
apply(comp.set,2,function(x) sum(x==2,na.rm=T))
## these are non-polymorphic, should be removed

# "X12_21396", 1023, 1521
comp.set.n = c(1023, 1521)
comp.set = geno.dx1[,comp.set.n]
cor(comp.set)
table(comp.set[,1]==comp.set[,2],useNA='ifany')
apply(comp.set,2,function(x) sum(x==0,na.rm=T))
apply(comp.set,2,function(x) sum(x==2,na.rm=T))
## one minor allele conincides

# "X12_30232" 2314, 3518, 4143
comp.set.n = c(2314, 3518, 4143)
comp.set = geno.dx1[,comp.set.n]
cor(comp.set)
table(comp.set[,1]==comp.set[,2],useNA='ifany')
apply(comp.set,2,function(x) sum(x==0,na.rm=T))
apply(comp.set,2,function(x) sum(x==2,na.rm=T))
## non-polymorphic, should be removed

### filter for at least 2 minor alleles
geno.dx2 = geno.dx1[,keep.list]
k2.list = apply(geno.dx2, 2, function(x){
  al0 = sum(x==0,na.rm=T)
  al2 = sum(x==2,na.rm=T)
  if (al0+al2 - abs(al0-al2) < 3) {
    return(F)
  } else {
    return(T)
  }
})
sum(!k2.list)
#[1] 18 would be removed
geno.dx2 = geno.dx2[,k2.list]
dim(geno.dx2)
#[1]  784 5800

myGD = data.frame(taxa=rownames(geno.dx2),geno.dx2)
myGM1 = myGM[myGM[,'Name'] %in% colnames(myGD),]

## check
setequal(as.character(myGD[,1]),as.character(myY1[,1]))
setequal(colnames(myGD)[-1],as.character(myGM1[,1]))

#save(myGD, myGM1, myY1, file='cleaned.myGD.fixed.RData')


rm(list=ls())
load('cleaned.myGD.fixed.RData')

# GAPIT - Genomic Association and Prediction Integrated Tool
# Designed by Zhiwu Zhang
# Written by Zhiwu Zhang, Alex Lipka, Feng Tian and You Tang
# Last update: September 15, 2015

#Install packages (Do this section only for new installation of R)
#-------------------------------------------------------------------------------
#source("http://www.bioconductor.org/biocLite.R") 
#biocLite("multtest")
#install.packages("gplots")
#install.packages("scatterplot3d")#The downloaded link at: http://cran.r-project.org/package=scatterplot3d

#Step 0: Import library and GAPIT functions run this section each time to start R)
#######################################################################################
library('MASS') # required for ginv
library(multtest)
library(gplots)
library(compiler) #required for cmpfun
library("scatterplot3d")

source("http://www.zzlab.net/GAPIT/emma.txt")
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")

#source("/Users/Zhiwu/Dropbox/Current/revolutionr/gapit/gapit_functions.txt")
#############################################################################################
## since PCA cannot take any NA, convert NA to 1 in myGD
myGDc = as.matrix(myGD[,-1])
myGDc[is.na(myGDc)]=1
myGDc = data.frame(myGDc)
myGDc = data.frame(taxa=rownames(myGDc),myGDc)

#save.image('myGD.cleaning.fixed.RData')

#Step 2: Run GAPIT
myGAPIT <- GAPIT(
Y=myY1,
GD=myGDc,
GM=myGM1,
PCA.total=3,
)

# results are stored in /PCA3.fixed.all

## save myGD, myGM1, myY1
write.table(myY1,'myY1.v2.csv',sep=',',quote=F,row.names=F)
write.table(myGD,'myGD.v2.csv',sep=',',quote=F,row.names=F)
write.table(myGM1,'myGM1.v2.csv',sep=',',quote=F,row.names=F)

### compile results from /PCA3.fixed.all
setwd('./PCA3.fixed.all')
all.files = dir()
gwas.files = all.files[grep('GWAS.Results',all.files)]
#compile
compiled.file = c()
for (i in gwas.files){
  g.f = read.csv(i, header=T, as.is=1)
  pheno.n = sub('GAPIT\\.\\.(.+)\\.GWAS.+$','\\1',i)
  g.f = cbind(rep(pheno.n,nrow(g.f)),g.f)
  compiled.file = rbind(compiled.file, g.f)
}

sel.snp = compiled.file[compiled.file$P.value < 5e-4 & compiled.file$maf > 0.01,]
colnames(sel.snp)[1]='Phenotype'
write.table(sel.snp,'compiled.5e_4.0.01.v2.csv',sep=',',quote=F,row.names=F)

##### above, 12/14/2015



#Tutorial 12: Compare to Power against FDR for GLM,MLM,CMLM,ECMLM,SUPER(PLINK)
#Hint:Program runing time is more than 24 hours for repetition 100 times.
#Run description:Please refer to page 34,35 of the User manual on GAPIT website.
#----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myGD <-read.table("mdp_numeric.txt", head = TRUE)
myGM <-read.table("mdp_SNP_information.txt", head = TRUE)
myKI <- read.table("KSN.txt", head = FALSE)
myG <- read.table("mdp_genotype_test.hmp.txt", head = FALSE)

#Step 2: Run GAPIT
#GAPIT.Power.compare.plink
GAPIT.Power.compare(
myG=myG,
myGD=myGD,
myGM=myGM,
myKI=myKI,
rel=100,
h2=0.9,
NQTN=5
)


#Tutorial 13: Genetic Prediction one time by cross validation
#-----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY<-read.table("mdp_traits.txt", head = TRUE)
myK<-read.table("KSN.txt", head = FALSE)

#Step 2: Run GAPIT
OnePre<-GAPIT.Prediction(
myK=myK,
y<-myY[,c(1,3)],
##y=y[,1:2],
num=5
)


#Tutorial 14: Compare accuracy to different folds by replicate times
#-----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY<- read.table("mdp_traits.txt", head = TRUE)
myGD <-read.table("mdp_numeric.txt", head = TRUE)

#Step 2: Run GAPIT
myCross<-GAPIT.cross_validation.compare(
myGD=myGD,
y=myY,
#y<-y[,c(1,3)],
rel=100,
tc<-c(2,5,10,20)  ##input compare to folds num
)


#Tutorial 15: Marker density and decade of linkage disequilibrium over distance
#-----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files

myG <- read.delim("mdp_genotype_test.hmp.txt", head = FALSE)

#Step 2: Run GAPIT

myGenotype<-GAPIT.Genotype.View(
myG=myG,
#chr=1,
#w1_start=30,
#w1_end=230,
#mav1=10
)

#Tutorial 16: Statistical distributions of phenotype
#-----------------------------------------------------------------------------------------
#Step 1: Set data directory and import files
myY  <- read.table("mdp_traits.txt", head = TRUE)

myPhenotype<-GAPIT.Phenotype.View(
myY=myY
)



