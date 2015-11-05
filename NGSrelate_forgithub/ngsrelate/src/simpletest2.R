library(Relate)
# test 2: Rettet sim
set.seed(1)
n<-100 #number of individuals
snp=1000 #number of SNPs
freq<-runif(snp,min=0.3,max=0.7) #the population allele frequency of the SNPs (assumed to be uniform for SNP chip data)
data<-matrix(rbinom(n*snp*10,2,rep(freq,10*n)),ncol=snp*10,byrow=T)+1
r<-replicate(10,sim_chr(snp,freq=freq, min=0.5, max=0.95, k=c(0.95,0.05,0), a=0.06, number_per_cm=5 ),simplify = F)
pos<-unlist(sapply(r,function(x) x$pos))
geno<-unlist(sapply(r,function(x) x$geno))
data<-rbind(rbind(as.vector(geno[1:snp,]),as.vector(geno[1:snp+snp,]))+1,data)

t3<-runHmmld(data,pair=c(1,2),pos=pos/1e6,chr=rep(1:10,each=1000))#estimate tracts of relatedness
plot(t3)
t3$kResult
t3$a


