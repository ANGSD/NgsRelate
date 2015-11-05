## --------------------------
## Make test data set
## --------------------------

library(Relate)
source("simrelatedatafuns.R")
source("ngsrelateHMM_270812.R")

set.seed(25)

nind=50
nsnp=50
nchr=2
dat <-simData(nind=nind,nsnp=nsnp,nchr=nchr,k0=0.75,k1=0.25)
geno<-dat$geno
pos<-dat$pos
chr<-dat$chr
freq<-dat$truefreq
af<-apply(geno-1,2,mean,na.rm=TRUE)/2   
genolikes<-rbind(t(getEpsilonLikes_old(geno[1,]-1)),t(getEpsilonLikes_old(geno[2,]-1)))


## ------------------------------------------------------------------------------------------
## Test 2: does Anders relate and my ngs relate agree??
## (First one run one real genos (epsilon=0.01) and the other one run with epsiolon 0.01 genolikes)
## ------------------------------------------------------------------------------------------

# Do we get same values when par is fixed?
ra1<-runHmmld(geno,pair=c(1,2),pos=pos/1e6,chr=chr,back=0,ld_adj=F,epsilon=0.01,par=c(0.02,0,0.25,0.75),calc.a=T)
ra2<-runHmmld(geno,pair=c(1,2),pos=pos/1e6,chr=chr,back=0,ld_adj=F,epsilon=0.01,par=c(0.02,0,0.25,0.75),calc.a=F)

ri1<-relateHMM(geno-1,pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,0.15),epsilon=0.01,par=c(0.02,0,0.25,0.75),calc.a=T,af=af)
ri2<-relateHMM(geno-1,pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,0.15),epsilon=0.01,par=c(0.02,0,0.25,0.75),calc.a=F,af=af)

ri1x<-ngsrelateHMM(genolikes,pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,0.15),epsilon=0.01,par=c(0.02,0,0.25,0.75),calc.a=T,af=af)
ri2x<-ngsrelateHMM(genolikes,pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,0.15),epsilon=0.01,par=c(0.02,0,0.25,0.75),calc.a=F,af=af)

# - same llh ? 
print(c(ra1$kResult,ra1$kLike))
print(c(ra2$kResult,ra2$kLike))
print(c(ri1$k,ri1$k.loglike))
print(c(ri2$k,ri2$k.loglike))
print(c(ri1x$k,ri1x$k.loglike))
print(c(ri2x$k,ri2x$k.loglike))


# same viterbi path?
all(ra2$path==ri1x$decode$path)
all(ra2$path==ri2x$decode$path)
# all true :) 

# same posterior?
max(ra2$post[2,]-ri2x$decode$post[,2])
#[1] 5.684342e-14
max(ra2$post[3,]-ri2x$decode$post[,3])
#[1] 1.920686e-14

max(ra2$post[2,]-ri1x$decode$post[,2])
#[1] 5.684342e-14
max(ra2$post[3,]-ri1x$decode$post[,3])
#[1] 1.920686e-14







