library(raster)
library(sp)
library(rgdal)

### download worldclim, bioclim data
wc.dat = getData('worldclim',var='bio',res=2.5)
# http://www.worldclim.org/bioclim
BIO1 = 'Annual Mean Temperature'
BIO2 = 'Mean Diurnal Range (Mean of monthly (max temp - min temp))'
BIO3 = 'Isothermality (BIO2/BIO7) (* 100)'
BIO4 = 'Temperature Seasonality (standard deviation *100)'
BIO5 = 'Max Temperature of Warmest Month'
BIO6 = 'Min Temperature of Coldest Month'
BIO7 = 'Temperature Annual Range (BIO5-BIO6)'
BIO8 = 'Mean Temperature of Wettest Quarter'
BIO9 = 'Mean Temperature of Driest Quarter'
BIO10 = 'Mean Temperature of Warmest Quarter'
BIO11 = 'Mean Temperature of Coldest Quarter'
BIO12 = 'Annual Precipitation'
BIO13 = 'Precipitation of Wettest Month'
BIO14 = 'Precipitation of Driest Month'
BIO15 = 'Precipitation Seasonality (Coefficient of Variation)'
BIO16 = 'Precipitation of Wettest Quarter'
BIO17 = 'Precipitation of Driest Quarter'
BIO18 = 'Precipitation of Warmest Quarter'
BIO19 = 'Precipitation of Coldest Quarter'

## plotting test
plot(wc.dat[[1]], xlim=c(-20,150),ylim=c(-10,65))

## load 'cleaned' coordinate data
coor.dat = read.delim('location.altitude_60fixed_cleaned.txt',header=T,row.names=1)

for (bioc in 1:19)){
  dat = wc.dat[[bioc]]
  dat = setMinMax(dat)
  tile.loc = dat@extent
  lonmin = tile.loc@xmin
  lonmax = tile.loc@xmax
  latmin = tile.loc@ymin
  latmax = tile.loc@ymax
  
  bioc.vals = c()
  dat1 = as.matrix(dat)
  for (loc.d in 1:nrow(coor.dat)){
    sel.t = coor.dat[loc.d,]
    if (is.na(sel.t$Longitude)) {
      bioc.vals = c(bioc.vals, NA)
      next
    }
    longi = sel.t[1,'Longitude']
    lati = sel.t[1,'Latitude']
    x.coor = round( (longi - lonmin) / (lonmax - lonmin) * dat@ncols )
    y.coor = round( (latmax - lati) / (latmax - latmin) * dat@nrows )
    bioc.vals = c(bioc.vals, dat1[y.coor, x.coor])
  }
  coor.dat[,names(wc.dat)[bioc]]=bioc.vals
}

save(coor.dat,file='coor.dat.RData')
write.table(coor.dat,file='coor.dat.txt',sep='\t',quote=F,row.names=F)

#### R code to compute Hudson Fst estimator
#### modified from https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CB0QFjAAahUKEwi6w6iF7LPIAhUCGD4KHQ47B5c&url=http%3A%2F%2Ffiles.figshare.com%2F1428100%2FText_S1.doc&usg=AFQjCNHgFbiyGqD-TLTe4iAejKj8U0YBkQ&sig2=GIRE73T1dZ4377LSMXs-lQ&cad=rja
# input data frame pop1 is a N x 4 matrix
# where N is the number of SNPs
# row names correspond to the SNP name
# MAF represent the minor allele frequency
# NCHROBS represent the number of chromosome observed (2 x sample size)
# A1 common allele
# A2 variant allele
# example
#> head(pop1,5)
#A1  A2  MAF	NCHROBS
#rs3094315	G	A	0.18590	156
#rs3131972	A	G	0.18350	158
#rs3115860	C	A	0.13160	152
#rs12562034	A	G	0.09615	156
#rs12124819	G	A	0.20950	148
#rs2980300	A	G	0.13290	158
# similarly for pop2

Hudson.Fst <- function(pop1,pop2,call.rate = 0.95,top.number = 10){
  # remove the SNPs that are not in common between the 2 populations
  snp.to.keep <- intersect(row.names(pop1),row.names(pop2))
  if (length(snp.to.keep) == 0){print("Error: no SNP in common");return(NULL)}
  pop1.k <- pop1[snp.to.keep,]
  pop2.k <- pop2[snp.to.keep,]
  # change the reference allele if is not concordant between the 2 populations
  if (sum(pop1.k$A1 == pop2.k$A1) != length(snp.to.keep)){
    idx <- which(pop1.k$A1 != pop2.k$A1)
    idx.rev <- which(pop1.k$A1 != pop2.k$A1 & pop1.k$A1 == pop2.k$A2)
    idx.rm  <- which(pop1.k$A1 != pop2.k$A1 & pop1.k$A1 != pop2.k$A2)
    if(length(idx.rev) > 0){
      provv <- pop1.k$A1[idx.rev]
      pop1.k$A1[idx.rev] <- pop1.k$A2[idx.rev]
      pop1.k$A2[idx.rev] <- provv
      pop1.k$MAF[idx.rev] <- 1 - pop1.k$MAF[idx.rev]
    }
    if(length(idx.rm) > 0){      
      pop1.k <- pop1.k[-idx.rm,]
      pop2.k <- pop2.k[-idx.rm,]
    }
  }
  mk.names = rownames(pop1.k)
  # remove SNPs with low call rate in one or both populations
  N1 <- pop1.k$NCHROBS
  N2 <- pop1.k$NCHROBS
  idx.rm.pop1 <- which(N1 < max(N1)*call.rate)
  idx.rm.pop2 <- which(N2 < max(N2)*call.rate)
  idx.rm.all <- union(idx.rm.pop1,idx.rm.pop2)
  if (length(idx.rm.all)) {
    pop1.k <- pop1.k[-idx.rm.all,]
    pop2.k <- pop2.k[-idx.rm.all,]
    mk.names = mk.names[-idx.rm.all]
  }
  # compute Hudson SNP_Fst and global Fst estimators
  p1 <- pop1.k$MAF
  p2 <- pop2.k$MAF
  n1 <- pop1.k$NCHROBS
  n2 <- pop2.k$NCHROBS
  fst.N <- (p1 - p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
  fst.D <- p1*(1-p2) + p2*(1-p1)
  Fst.v <- fst.N/fst.D
  names(Fst.v) <- mk.names
  Fst.o <- Fst.v[order(Fst.v,decreasing=TRUE)]
  q.Fst = quantile(Fst.o,probs=c(0.5,0.75,0.95,0.975),na.rm=T)
  mu1 <- mean(fst.N)
  mu2 <- mean(fst.D)
  se1 <- sd(fst.N)/sqrt(length(fst.N))
  se2 <- sd(fst.D)/sqrt(length(fst.D))
  F.global <- mu1/mu2
  se.F <- sqrt(se1^2+se2^2)
  output <- list()
  output[[1]] <- c(F.global,se.F,sum(!is.na(Fst.v)),
                   q.Fst[1],q.Fst[2],q.Fst[3],q.Fst[4])
  names(output[[1]]) <- c("Hudson.Fst.mean","Hudson.Fst.sem","Total.marker.number",
                          "Hudson.Fst.50","Hudson.Fst.75","Hudson.Fst.95","Hudson.Fst.97.5")
  output[[2]] <- data.frame(Fst.o[1:top.number])
  names(output[[2]]) <- c("Hudson.Fst")
  return(output)
}


### calculate Hudson Fst
## load genotype data
geno.d = read.delim('Land_6152_SNPs_AB.txt',header=T,row.names=1)

coor.d1 = coor.dat[!is.na(coor.dat$Latitude),]
fst.rall = list()
for (bioc in 7:ncol(coor.d1)){
  bioc.d0 = coor.d1[,bioc]
  bioc.d = bioc.d0[!is.na(bioc.d0)]
  names(bioc.d) = coor.d1[!is.na(bioc.d0),'Accession.ID']
  cutoffs = quantile(bioc.d, probs = seq(0.05, 0.95, 0.05))
  cutoffs[cutoffs==min(bioc.d)]=NA
  cutoffs[cutoffs==max(bioc.d)]=NA
  fst.r = list()
  for (co in names(cutoffs)){
    cutoff = cutoffs[co]
    if (is.na(cutoff)) {
      fst.r[[co]]=NA
      next
    }
    g1 = bioc.d > cutoff
    g1.loc.n = names(bioc.d)[g1]
    g2.loc.n = names(bioc.d)[!g1]
    g1.mks = as.matrix(geno.d[g1.loc.n,])
    g2.mks = as.matrix(geno.d[g2.loc.n,])
    g1.pop = t(apply(g1.mks, 2, function(x){
      mk.str = paste(x, collapse='')
      mk.vec = unlist(strsplit(mk.str,split=''))
      A.numb = length(grep('A',mk.vec))
      to.al = length(mk.vec)
      A.T.rat = A.numb/to.al
      return(c(A.T.rat,to.al))
    }))
    colnames(g1.pop) = c('MAF','NCHROBS')
    g1.pop.t = data.frame(A1='B',A2='A',g1.pop)
    g2.pop = t(apply(g2.mks, 2, function(x){
      mk.str = paste(x, collapse='')
      mk.vec = unlist(strsplit(mk.str,split=''))
      A.numb = length(grep('A',mk.vec))
      to.al = length(mk.vec)
      A.T.rat = A.numb/to.al
      return(c(A.T.rat,to.al))
    }))
    colnames(g2.pop) = c('MAF','NCHROBS')
    g2.pop.t = data.frame(A1='B',A2='A',g2.pop)
    
    fst.r[[co]]= Hudson.Fst(g1.pop.t,g2.pop.t,call.rate = 0.95,top.number = 10)
  }
  fst.rall[[colnames(coor.d1)[bioc]]] = fst.r
}

save(fst.rall,file='fst.rall.RData')

#### visualization of top Fst markers for each bioclim data
gmap.d = read.delim('GeneticMap_iSelect_9k.txt',header=T)
rownames(gmap.d)=gmap.d[,'SNP']
bound.v = tapply(gmap.d$Cumulative,factor(gmap.d$Chromosome),max)

pdf('bioclim.19.HFst.top10.gmap.pdf',width=7.5,height=10)
opar=par(las=1, mfrow=c(2,1))
for (bioc in names(fst.rall)) {
  fst.r = fst.rall[[bioc]]
  bio.n = toupper(bioc)
  
  fst.d = data.frame()
  fst.mean = c()
  fst.sd = c()
  for (cutoff in names(fst.r)){
    if (is.na(fst.r[[cutoff]][1])) {
      f.d = data.frame(marker.name=NA, ch.position=NA,
                       Hudson.Fst=NA,cutoff=cutoff)
      fst.mean = c(fst.mean, NA)
      fst.sd = c(fst.sd, NA)
    } else {
      fst.mean = c(fst.mean, fst.r[[cutoff]][[1]][1])
      fst.sd = c(fst.sd, fst.r[[cutoff]][[1]][2])
      #mk.numb = fst.r[[cutoff]][[1]][3]
      top10.fst = fst.r[[cutoff]][[2]]
      mk.names = rownames(top10.fst)
      ch.pos = gmap.d[mk.names,'Cumulative']
      f.d = data.frame(marker.name=mk.names,
                       ch.position=ch.pos,
                       top10.fst,cutoff=cutoff)
    }
    fst.d = rbind(fst.d, f.d)
  }
  col1 = cut(fst.d$Hudson.Fst,breaks=c(0,0.3,0.4,0.5,1),
             labels=c('blue','green','orange','red'))
  
  stripchart(ch.position ~ cutoff,data=fst.d[col1=='blue',],
             col='blue',
             xlim=c(0,1430),ylim=c(0,21.5),
             xlab='chromosomal position', ylab='Cutoff percentile', 
             main=get(bio.n))
  fake = sapply(1:19, function(x) {
    lines(c(0,1113),c(x,x),col='gray')
  })
  stripchart(ch.position ~ cutoff,data=fst.d[col1=='green',],
             col='green', add=T)
  stripchart(ch.position ~ cutoff,data=fst.d[col1=='orange',],
             col='orange', add=T)
  stripchart(ch.position ~ cutoff,data=fst.d[col1=='red',],
             col='red', add=T)
  abline(v=c(0,bound.v),lty=3,col='skyblue')
  centr.pos = ( c(0, bound.v[-7]) + bound.v )/2
  text(centr.pos,rep(20.5,7),c('I','II','III','IV','V','VI','VII'),
       col='skyblue')
  
  bioc.d0 = coor.d1[,bioc]
  bioc.d = bioc.d0[!is.na(bioc.d0)]
  cutoffs = quantile(bioc.d, probs = seq(0, 1, 0.05))
  
  fake = sapply(1:19, function(x){
    text(c(1125,1230,1335),c(x,x),
         labels=as.character(round(c(cutoffs[x+1],fst.mean[x],fst.sd[x]),digits=4)),
         cex=0.6, pos=4)
  })
  text(c(1125,1125),c(0,20),paste(cutoffs[c(1,21)],c(' min',' max'),sep=''),
       cex=0.6, pos=4)
  text(c(1125,1230,1335),c(21.2,21.2,21.2),
       labels=c('Cutoff', 'Fst mean', 'Fst SD'),cex=0.6, pos=4)
  
  
}
par(opar)
dev.off()

#### do the same with the standardized Fst
pdf('bioclim.19.HFst.top10.norm.gmap.pdf',width=7.5,height=10)
opar=par(las=1, mfrow=c(2,1))
for (bioc in names(fst.rall)) {
  fst.r = fst.rall[[bioc]]
  bio.n = toupper(bioc)
  
  fst.d = data.frame()
  fst.mean = c()
  fst.sd = c()
  fst.75 = c()
  for (cutoff in names(fst.r)){
    if (is.na(fst.r[[cutoff]][1])) {
      f.d = data.frame(marker.name=NA, ch.position=NA,
                       Hudson.Fst=NA,cutoff=cutoff)
      fst.mean = c(fst.mean, NA)
      fst.sd = c(fst.sd, NA)
      fst.75 = c(fst.75, NA)
    } else {
      norm.f = fst.r[[cutoff]][[1]][5]
      # 75 percentile normalization
      fst.mean = c(fst.mean, fst.r[[cutoff]][[1]][1]/norm.f)
      fst.sd = c(fst.sd, fst.r[[cutoff]][[1]][2]/norm.f)
      fst.75 = c(fst.75, norm.f)
      top10.fst = fst.r[[cutoff]][[2]]/norm.f
      mk.names = rownames(top10.fst)
      ch.pos = gmap.d[mk.names,'Cumulative']
      f.d = data.frame(marker.name=mk.names,
                       ch.position=ch.pos,
                       top10.fst,cutoff=cutoff)
    }
    fst.d = rbind(fst.d, f.d)
  }
  col1 = cut(fst.d$Hudson.Fst,breaks=c(0,4.5,5.25,6,10),
             labels=c('blue','green','orange','red'))
  
  stripchart(ch.position ~ cutoff,data=fst.d[col1=='blue',],
             col='blue',
             xlim=c(0,1535),ylim=c(-0.2,21.4),
             xlab='Cumulative chromosomal position', ylab='Cutoff percentile', 
             main=get(bio.n))
  fake = sapply(1:19, function(x) {
    lines(c(0,1113),c(x,x),col='gray')
  })
  stripchart(ch.position ~ cutoff,data=fst.d[col1=='green',],
             col='green', add=T)
  stripchart(ch.position ~ cutoff,data=fst.d[col1=='orange',],
             col='orange', add=T)
  stripchart(ch.position ~ cutoff,data=fst.d[col1=='red',],
             col='red', add=T)
  abline(v=c(0,bound.v),lty=3,col='skyblue')
  centr.pos = ( c(0, bound.v[-7]) + bound.v )/2
  text(centr.pos,rep(20.5,7),paste(1:7,'H',sep=''),
       col='skyblue')
  
  bioc.d0 = coor.d1[,bioc]
  bioc.d = bioc.d0[!is.na(bioc.d0)]
  cutoffs = quantile(bioc.d, probs = seq(0, 1, 0.05))
  
  fake = sapply(1:19, function(x){
    text(c(1125,1230,1335,1440),x,
         labels=as.character(round(c(cutoffs[x+1],fst.mean[x],fst.sd[x],fst.75[x]),digits=4)),
         cex=0.6, pos=4)
  })
  text(c(1125,1125),c(0,20),paste(cutoffs[c(1,21)],c(' min',' max'),sep=''),
       cex=0.6, pos=4)
  text(c(1125,1230,1335,1440),21.2,
       labels=c('Cutoff', 'Fst mean', 'Fst SD','Norm F'),cex=0.6, pos=4)
  legend(0,0.45,c('<4.5','4.5-5.25','5.25-6','>6'),pch=0,cex=0.6,pt.cex=0.8,
         col=c('blue','green','orange','red'),ncol=4)
  
  
}
par(opar)
dev.off()

#### mapping consistent top markers on the bioclim map
### marker identification for each bioclim
### criteia
### in top10 in atleast 3 percentile cutoffs and at least one red or 3 orange
### OR in top10 in at least 5 percentile cutoffs and at least one orange or red
### OR in top10 in at least 7 percentile cutoffs and at least two green

sel.mk.l = list()
for (bioc in names(fst.rall)) {
  fst.r = fst.rall[[bioc]]
  
  fst.d = data.frame()
  for (cutoff in names(fst.r)){
    if (is.na(fst.r[[cutoff]][1])) {
      f.d = data.frame(marker.name=NA, ch.position=NA,
                       Hudson.Fst=NA,cutoff=cutoff)
    } else {
      norm.f = fst.r[[cutoff]][[1]][5]
      top10.fst = fst.r[[cutoff]][[2]]/norm.f
      mk.names = rownames(top10.fst)
      ch.pos = gmap.d[mk.names,'Cumulative']
      f.d = data.frame(marker.name=mk.names,
                       ch.position=ch.pos,
                       top10.fst,cutoff=cutoff)
    }
    fst.d = rbind(fst.d, f.d)
  }
  ta.m = table(fst.d$marker.name)
  ap3 = names(ta.m)[ta.m>=3]
  ap5 = names(ta.m)[ta.m>=5]
  ap7 = names(ta.m)[ta.m>=7]
  sel.mk.n = c()
  for (mk.n in ap3){
    if (max(fst.d[fst.d$marker.name == mk.n,'Hudson.Fst'],na.rm=T) > 6|
          sort(fst.d[fst.d$marker.name == mk.n,'Hudson.Fst'],decreasing=T,na.last=T)[3]>5.25) {
      sel.mk.n = c(sel.mk.n, mk.n)
    }
  }
  for (mk.n in ap5) {
    if (max(fst.d[fst.d$marker.name == mk.n,'Hudson.Fst'],na.rm=T) > 5.25){
      sel.mk.n = c(sel.mk.n, mk.n)
    }
  }
  for (mk.n in ap7) {
    if (sort(fst.d[fst.d$marker.name == mk.n,'Hudson.Fst'],decreasing=T,na.last=T)[2]>4.5) {
      sel.mk.n = c(sel.mk.n, mk.n)
    }
  }
  sel.mk.n = as.character(unique(sel.mk.n))
  sel.mk.l[[bioc]]=sel.mk.n
}
sel.mk.ll = unlist(sapply(1:length(sel.mk.l),function(x) {
  rep(names(sel.mk.l)[x],length(sel.mk.l[[x]]))
}))
sel.mk.lt = data.frame(marker.name=unlist(sel.mk.l),
                       bioclim=sel.mk.ll,stringsAsFactors=F)

mapped.mk = rownames(gmap.d)

pdf('bioclim.19.HFst.sel.marker.alleles.gmap.pdf',width=10,height=7.5)

for (bioc in names(sel.mk.l)) {
  fst.r = fst.rall[[bioc]]
  bio.n = toupper(bioc)
  sel.mk.n = sel.mk.l[[bioc]]  
  for (mk.n in sel.mk.n){
    bioc.l = sel.mk.lt[sel.mk.lt$marker.name==mk.n,'bioclim']
    bioc.l = sub('bio(\\d+)','\\1',bioc.l)
    if (mk.n %in% mapped.mk){
      sub.t = paste(mk.n, ', Chro: ', gmap.d[mk.n,'Chro'], 
                    ', ',gmap.d[mk.n,'cm'], 'cM, Bioclim: ',
                    paste(bioc.l,collapse=','),sep='')
    } else {
      sub.t = paste(mk.n, ', Chro location unavailable, Bioclim: ',
                    paste(bioc.l,collapse=','),sep='')
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
}

dev.off()

save.image('worldclim.miner.RData')

########### ver 3
#### calculate the allele frequency (fit to a sigmoid) for eveyr marker
#### visualize the distribution of the profile along the envrionmental data value
#### identify those significantly different in the sigmoid parameters (99 percentile?)

### to allow easy comparisons
### the more frequent allele as 0
### scale the env values to 0-1
### flip the env value, if needed, to make the rarer allele associate with higher env value

### calculate the distribution of allele frequencies along the env variable
### for evey marker (call.rate=0.95)
call.rate = 0.95
sig.cutoff = 0.001 # significanc cutoff for the slope in linear regression

all.mk.names = colnames(geno.d)

geno.byn = apply(geno.d,2,function(x){
  is.A=grepl('A',x)
  is.B=grepl('B',x)
  if (sum(is.B) >= sum(is.A)){
    AB.ind = is.A
    AB.flip = 0
  } else {
    AB.ind = is.B
    AB.flip = 1
  }
  AB.ind[is.A+is.B ==0]=NA
  AB.ind[is.A+is.B ==2]=0.5
  return(list(AB.ind = AB.ind, AB.flip=AB.flip))
})
geno.byn.df = data.frame(sapply(geno.byn,'[',1))
rownames(geno.byn.df)=rownames(geno.d)
colnames(geno.byn.df)=colnames(geno.d)
geno.byn.fl = unlist(sapply(geno.byn,'[',2))
names(geno.byn.fl)=colnames(geno.byn.df)
coor.dat.allmk = data.frame(coor.d1,geno.byn.df[coor.d1$Accession.ID,])
coor.dat.allmk[,names(wc.dat)] = apply(coor.dat.allmk[,names(wc.dat)],2,function(x){
  max.bioc = max(x,na.rm=T)
  min.bioc = min(x,na.rm=T)
  return((x - min.bioc)/(max.bioc - min.bioc))
})

### fist a linear regression - remove those insignificant and determine the slope sign
slope.rll = list()
int.rll = list()
for (bioc in names(wc.dat)) {
  slope.rl = c()
  int.rl = c()
  for (mk.n in colnames(geno.byn.df)) {
    sub.t = coor.dat.allmk[,c(bioc,mk.n)]
    sub.t = sub.t[complete.cases(sub.t),]
    if (nrow(sub.t) < call.rate * nrow(coor.dat.allmk)){
      slope.rl = c(slope.rl, NA)
      int.rl = c(int.rl, NA)
      next
    }
    if (sum(sub.t[mk.n]) == 0){
      slope.rl = c(slope.rl, 0)
      int.rl = c(int.rl, 0)
      next
    }
    lm.out = lm(sub.t[,mk.n] ~ sub.t[,bioc])
    if (summary(lm.out)$coef[2,'Pr(>|t|)'] > sig.cutoff) {
      slope.rl = c(slope.rl,0)
      int.rl = c(int.rl, summary(lm.out)$coef[1,'Estimate'])
      next
    } 
    slope.rl = c(slope.rl,summary(lm.out)$coef[2,'Estimate']) 
    int.rl = c(int.rl, summary(lm.out)$coef[1,'Estimate'])
  }
  names(slope.rl)=colnames(geno.byn.df)
  names(int.rl)=colnames(geno.byn.df)
  slope.rll[[bioc]]=slope.rl
  int.rll[[bioc]] = int.rl
}

### second logistic regression
logis.paras.rll = list()
for (bioc in names(wc.dat)) {
  logis.para.rl = list()
  mk.names = c()
  for (mk.n in colnames(geno.byn.df)) {
    sub.t = coor.dat.allmk[,c(bioc,mk.n)]
    sub.t = sub.t[complete.cases(sub.t),]
    sl = slope.rll[[bioc]][mk.n]
    if (is.na(sl) | sl==0) next
    if (sl < 0) {
      sub.t[,bioc]= -sub.t[,bioc]
    }
    ## sliding window mean to approximate the allele frequency
    ## as nls is finicky
    ww = 50 # window size by the number of samples
    start.p = 1:(nrow(sub.t)-ww+1)
    sub.t.ord = sub.t[order(sub.t[,bioc]),]
    win.means = sapply(start.p, function(x){
      mean(sub.t.ord[x:(x+ww-1),mk.n])
    })
    win.x =  (start.p + (ww-1)/2)/nrow(sub.t)
    win.xq = quantile(sub.t[,bioc],prob=win.x)
    
    plot(sub.t[,bioc], sub.t[,mk.n] + rnorm(nrow(sub.t),0,0.02), 
         xlab=paste(bio.n,get(bio.n)),ylab='Allele frequency: AA=0,BB=1',
         cex = 0.4,
         cex.lab=0.7,cex.axis=0.5)
    lines(win.xq,win.means, col='red')
    glm.out = glm(win.means ~ win.xq,
                  family=binomial(logit))
    lines(sort(win.xq),glm.out$fitted[order(win.xq)],col='blue')
    slope = coef(glm.out)[2]
    mid.p = -coef(glm.out)[1]/coef(glm.out)[2]
    s.glm.out = summary(glm.out)
    sig.l = 1-pchisq(s.glm.out$deviance,s.glm.out$df.residual)
    
    nls.out = nls(win.means ~ maxv / ( 1 + exp(-sl *(win.xq-midv) )) +minv,
                  start=list(sl=abs(sl),midv=0.7,
                             minv=0.2, maxv=1))
                                        
    if (summary(lm.out)$coef[2,'Pr(>|t|)'] > sig.cutoff) {
      slope.rl = c(slope.rl,0)
      int.rl = c(int.rl, summary(lm.out)$coef[1,'Estimate'])
      next
    } 
    slope.rl = c(slope.rl,summary(lm.out)$coef[2,'Estimate']) 
    int.rl = c(int.rl, summary(lm.out)$coef[1,'Estimate'])
  }
  names(slope.rl)=colnames(geno.byn.df)
  names(int.rl)=colnames(geno.byn.df)
  slope.rll[[bioc]]=slope.rl
  int.rll[[bioc]] = int.rl
}

























sig.lv = list()
slope.v = list()
max.vv = list()
min.vv = list()
for (bioc in names(sel.mk.l)) {
  sig.l1=c()
  slope1 = c()
  min.v = c()
  max.v = c()
  mk.n1 = c()
  for (mk.n in samp.mk.names) {
    sub.t = coor.dat.allmk[,c(bioc,mk.n)]
    sub.t = sub.t[complete.cases(sub.t),]
    glm.out = glm(sub.t[,mk.n] ~ sub.t[,bioc],
                  family=binomial(logit))
    s.glm.out = summary(glm.out)
    sig.l1 = c(sig.l1, 1-pchisq(s.glm.out$deviance,s.glm.out$df.residual))
    slope1 = c(slope1, coef(glm.out)[2])
    range.v = range(fitted(glm.out))
    min.v = c(min.v, range.v[1])
    max.v = c(max.v, range.v[2])
    mk.n1 = c(mk.n1, mk.n)
  }
  names(sig.l1)=names(slope1)=names(min.v)=names(max.v)=mk.n1
  sig.lv[[bioc]]=sig.l1
  slope.v[[bioc]]=slope1
  max.vv[[bioc]]=max.v
  min.vv[[bioc]]=min.v
}

pdf('bioclim.19.sel.marker.allele.freq2.pdf',width=7.5,height=10)
opar=par(mfrow=c(5,4),las=1,mar=c(2.6,2.2,0.2,0.4), mgp=c(1.2,0,0))

for (bioc in names(sel.mk.l)) {
  bio.n = toupper(bioc)
  sel.mk.n = sel.mk.l[[bioc]]  
  
  for (mk.n in sel.mk.n){
    sub.t = coor.dat.allmk[,c(bioc,mk.n)]
    sub.t = sub.t[complete.cases(sub.t),]
    bioc.l = sel.mk.lt[sel.mk.lt$marker.name==mk.n,'bioclim']
    bioc.l = sub('bio(\\d+)','\\1',bioc.l)
    if (mk.n %in% mapped.mk){
      mk.d = paste(mk.n, ', Chro: ', gmap.d[mk.n,'Chro'], 
                    ', ',gmap.d[mk.n,'cm'], 'cM',sep='')
      other.b = paste('Bioclim: ',
                    paste(bioc.l,collapse=','),sep='')
    } else {
      mk.d = paste(mk.n, ', Chro loc unavailable',sep='')
      other.b = paste('Bioclim: ',
                    paste(bioc.l,collapse=','),sep='')
    }
    
    ## sliding window mean to approximate the allele frequency
    ww = 30 # window size by the number of samples
    start.p = 1:(nrow(sub.t)-ww+1)
    sub.t.ord = sub.t[order(sub.t[,bioc]),]
    win.means = sapply(start.p, function(x){
      mean(sub.t.ord[x:(x+ww-1),mk.n])
    })
    win.x =  (start.p + (ww-1)/2)/nrow(sub.t)
    win.xq = quantile(sub.t[,bioc],prob=win.x)
    
    plot(sub.t[,bioc], sub.t[,mk.n] + rnorm(nrow(sub.t),0,0.02), 
         xlab=paste(bio.n,get(bio.n)),ylab='Allele frequency: AA=0,BB=1',
         cex = 0.4,
         cex.lab=0.7,cex.axis=0.5)
    lines(win.xq,win.means, col='red')
    glm.out = glm(sub.t[,mk.n] ~ sub.t[,bioc],
                  family=binomial(logit))
    lines(sort(sub.t[,bioc]),glm.out$fitted[order(sub.t[,bioc])],col='red')
    slope = coef(glm.out)[2]
    mid.p = -coef(glm.out)[1]/coef(glm.out)[2]
    s.glm.out = summary(glm.out)
    sig.l = 1-pchisq(s.glm.out$deviance,s.glm.out$df.residual)

    xmax = max(sub.t[,bioc])
    xmin = min(sub.t[,bioc])
    yl = 0.9
    if (slope <0) {
      xl = xmax
      pos.v =2
    } else {
      xl = xmin
      pos.v =4
    }
    text(xl,c(yl-0:4*0.07),c(
      mk.d, other.b,
      paste('slope = ',round(slope,digits = 3)),
      paste('mid.point = ', round(mid.p, digits = 1)),
      paste('significance = ', round(sig.l, digits = 4))),
      col='blue', pos=pos.v, cex=0.6)
  }
}
par(opar)
dev.off()

first.s = max.vv[[1]]>0.75 & min.vv[[1]]<0.25
sum(first.s)
#[1] 61
second.s = first.s & sig.lv[[1]] >0.1
sum(second.s)
#[1] 8

