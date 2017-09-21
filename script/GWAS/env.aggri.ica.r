### aggregate the bioclim data
### using ICA

## load packages
library(ica)
library(RColorBrewer)

## load data
coor.dat = read.delim('coor.dat.txt',header=T,row.names=1)
coor.dat = coor.dat[complete.cases(coor.dat),] #29 uncertain locations removed
env.d = as.matrix( coor.dat[,6:24] )

## standardize and look at the distribution of each bioclim
env.ds = apply(env.d, 2, function(x){
  mean.v = mean(x)
  sd.v = sd(x)
  (x - mean.v)/sd.v
})

plot(0,0, xlim=c(min(env.ds) ,max(env.ds)),
     ylim = c(0,1.3), type='n')
col1 = brewer.pal(8, 'Set2')
for (i in 1:ncol(env.ds)){
  e.d = env.ds[,i]
  points(density(e.d), type='l', col=col1[i %% 8 +1],
         lty=i)
}
points(density(env.ds), type='l', lwd=2)
# no common distribution patterns
# some are bimodal or more
# just use the standardized

## PCA first
pc.env = prcomp(env.ds)
plot(pc.env$sdev^2,type='o')
sum(pc.env$sdev[1:4]^2)/sum(pc.env$sdev^2)
# 88% by PCs1-4
sum(pc.env$sdev[1:2]^2)/sum(pc.env$sdev^2)
# 66% by PCs1-2
# let's take PCs1-6

env.pcx = cbind(env.ds,pc.env$x[,1:6])
round(cor(env.pcx)^2,digits=3)
# none of the PCs1-6 explain well for bio19,18, particularly
# while PC1~bio5,10, PC2~bio4,7
# it doesn't make sense to use PCs nor those closest to PCs

## ICA
ic.nv = c()
for (nc in 2:8){
  ic.env = icaimax(env.ds,nc)
  ic.nv = c(ic.nv, sum(ic.env$vafs))
}
names(ic.nv) = 2:8
ic.nv
#2         3         4         5         6         7         8 
#0.6601323 0.7877981 0.8761619 0.9224002 0.9534895 0.9757805 0.9901049 

pc.nv = c()
pc.nv[1] = pc.env$sdev[1]^2
for (pc in 2:8){
  pc.nv = c(pc.nv, pc.nv[length(pc.nv)]+pc.env$sdev[pc]^2)
}
names(pc.nv)=1:8
pc.nv = pc.nv/sum(pc.env$sdev^2)
pc.nv
#1         2         3         4         5         6         7         8 
#0.3807936 0.6601323 0.7877981 0.8761619 0.9224002 0.9534895 0.9757805 0.9901049 
# no difference in explained variation between PCA and ICA??
# this means that ICA opperates in the PC space!!

ic.env = icaimax(env.ds,3)
i.mat1 = ic.env$R
env.icx = cbind(env.ds,ic.env$S)
round(cor(env.icx)^2,digits=3)

ic.env = icaimax(env.ds,4)
i.mat1 = ic.env$R
env.icx = cbind(env.ds,ic.env$S)
round(cor(env.icx)^2,digits=3)

## ICA with nc=3 makes sense: 
## IC1, temperature; IC2, precipitation dry/cold; IC3, precipitation wet/warm
## then add bio2,3,8,15.
## however, I think that this would lose a lot of information.
## association map separately then compare is the best way.
ic.env = icaimax(env.ds,3)
i.mat1 = ic.env$R
env.icx = cbind(env.ds,ic.env$S)
round(cor(env.icx)^2,digits=3)

hist(ic.env$S[,1])
hist(ic.env$S[,2])
hist(ic.env$S[,3])

coor.ic.dat = data.frame(coor.dat,round(ic.env$S,digits=3))
colnames(coor.ic.dat)[25:27]=c('IC1','IC2','IC3')

write.table(coor.ic.dat,file='coor.ic.dat.txt',sep='\t',quote=F)

