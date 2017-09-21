setwd('C:\\run_space\\srtm')

library(raster)
library(sp)
library(rgdal)

### load sample location coordinates and bio values
load('worldclim.plot.RData')

## load genotype data
geno.d = read.delim('Land_6152_SNPs_AB.txt',header=T,row.names=1)

#### visualization of top Fst markers for each bioclim data
gmap.d = read.delim('GeneticMap_iSelect_9k.txt',header=T)
rownames(gmap.d)=gmap.d[,'SNP']
bound.v = tapply(gmap.d$Cumulative,factor(gmap.d$Chromosome),max)

#### import the SNP marker and bioc list
snp.bio = read.csv('compiled.5e_4.0.01.v2.csv', header=T)


pdf('compiled.5e_4.0.01.v2.pdf',width=7.5,height=10)
opar=par(las=1, mfrow=c(2,1))
for (to.map in 1:nrow(snp.bio)){
  bioc = as.character(snp.bio[to.map,'Phenotype'])
  if (!grepl('bio',bioc)) next
  bio.n = toupper(bioc)
  mk.n = as.character(snp.bio[to.map,'SNP'])
  if (mk.n %in% mapped.mk){
    sub.t = paste(mk.n, ', Chro: ', gmap.d[mk.n,'Chro'], 
                  ', ',gmap.d[mk.n,'cm'], 'cM, P-value: ',
                  sprintf('%.3e',snp.bio[to.map,'P.value']),
                  sep='')
  } else {
    sub.t = paste(mk.n, ', Chro location unavailable, P-value: ',
                  sprintf('%.3e',snp.bio[to.map,'P.value']),
                  sep='')
  }
  plot(wc.dat[[bioc]], xlim=c(-20,150),ylim=c(-10,65),
       main=paste(bio.n,get(bio.n)),
       sub=sub.t)
  geno.r = geno.d[,mk.n]
  is.A=grepl('A',geno.r)
  is.B=grepl('B',geno.r)*16
  AB.pch = is.A + is.B
  AB.pch[AB.pch==17]=13
  AB.pch[AB.pch==0]=NA
  names(AB.pch)=rownames(geno.d)
  AB.col = AB.pch
  AB.col[AB.col==1]='red'
  AB.col[AB.col==13|AB.col==16]='blue'
  fake=sapply(rownames(geno.d),function(x){
    points(coor.dat[coor.dat$Accession.ID==x,'Longitude'],
           coor.dat[coor.dat$Accession.ID==x,'Latitude'],
           pch=AB.pch[x], col=AB.col[x], cex=0.7)
  })
}
dev.off()

