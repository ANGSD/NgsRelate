## --------------------------
## Make test data set
## --------------------------

library(Relate)
source("simrelatedatafuns.R")
source("ngsrelateHMM_270812.R")

set.seed(25)

nind=50
nsnp=500
nchr=2
dat <-simData(nind=nind,nsnp=nsnp,nchr=nchr,k0=0.75,k1=0.25)
geno<-dat$geno
pos<-dat$pos
chr<-dat$chr

## ------------------------------------------------------
## Test 1: does Anders' relate and my relate agree??
## (Both run with out LD correction and on true genos)
## ------------------------------------------------------

# Do we reach more or less same ML estimates
ra1<-runHmmld(geno,pair=c(1,2),pos=pos/1e6,chr=chr,back=0,ld_adj=F,alim=c(0.001,0.05),epsilon=0.01,timesToRun=50,timesToConv=20)
ra2<-runHmmld(geno,pair=c(1,2),pos=pos/1e6,chr=chr,back=0,ld_adj=F,alim=c(0.001,0.05),epsilon=0.01,timesToRun=50,timesToConv=20,calc.a=F)

ri1<-relateHMM(geno-1,pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,0.05),epsilon=0.01,opti=50,calc.a=T)
ri2<-relateHMM(geno-1,pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,0.05),epsilon=0.01,opti=50,calc.a=F)


ro<-function(x){sapply(x,function(y){round(y,2)})}
print(ro(c(ra1$a,ra1$kResult,ra1$kLike)))
print(ro(c(ra2$a,ra2$kResult,ra2$kLike)))
print(ro(c(ri1$a,ri1$k,ri1$k.loglike)))
print(ro(c(ri2$a,ri2$k,ri2$k.loglike)))
# more or less...

# Do we get same values when par is fixed?
ra1<-runHmmld(geno,pair=c(1,2),pos=pos/1e6,chr=chr,back=0,ld_adj=F,epsilon=0.01,par=c(ra2$a,ra2$kResult))
ra2<-runHmmld(geno,pair=c(1,2),pos=pos/1e6,chr=chr,back=0,ld_adj=F,epsilon=0.01,par=c(ra2$a,ra2$kResult),calc.a=F)

ri2<-relateHMM(geno-1,pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,0.15),epsilon=0.01,par=c(ra2$a,ra2$kResult),calc.a=F)
ri1<-relateHMM(geno-1,pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,0.15),epsilon=0.01,par=c(ra2$a,ra2$kResult),calc.a=T)

print(c(ra1$kResult,ra1$kLike))
print(c(ra2$kResult,ra2$kLike))
print(c(ri1$k,ri1$k.loglike))
print(c(ri2$k,ri2$k.loglike))

# Yep :)

ra1<-runHmmld(geno,pair=c(1,2),pos=pos/1e6,chr=chr,back=0,ld_adj=F,epsilon=0.01,par=c(0.02,0,0.25,0.75))
ra2<-runHmmld(geno,pair=c(1,2),pos=pos/1e6,chr=chr,back=0,ld_adj=F,epsilon=0.01,par=c(0.02,0,0.25,0.75),calc.a=F)

ri2<-relateHMM(geno-1,pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,0.15),epsilon=0.01,par=c(0.02,0,0.25,0.75),calc.a=F)
ri1<-relateHMM(geno-1,pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,0.15),epsilon=0.01,par=c(0.02,0,0.25,0.75),calc.a=T)

# - same llh ? 
print(c(ra1$kResult,ra1$kLike))
print(c(ra2$kResult,ra2$kLike))
print(c(ri1$k,ri1$k.loglike))
print(c(ri2$k,ri2$k.loglike))
# yep :)

# same viterbi path?
all(ra2$path==ri1$decode$path)
all(ra2$path==ri2$decode$path)
all(ra2$path==ra1$path)
# all true :) 

# same posterior?
max(ra2$post[2,]-ri2$decode$post[,2])
#[1] 4.547474e-13
max(ra2$post[3,]-ri2$decode$post[,3])
#[1] 4.547474e-13

max(ra2$post[3,]-ri1$decode$post[,3])
#[1] 4.547474e-13
max(ra2$post[2,]-ri1$decode$post[,2])
#[1] 4.547474e-13






