library(Relate)
set.seed(1)

# test 1: Anders basis kode (har kun aendret et par parameter vaerdier)
n<-100 #number of individuals
snp=1000 #number of SNPs

nchr <- 1

freq<-runif(snp,min=0.05,max=0.95) #the population allele frequency of the SNPs (assumed to be uniform for SNP chip data)
data<-matrix(rbinom(n*snp*nchr,2,rep(freq,nchr)),ncol=snp*nchr)+1

r<-replicate(nchr,sim_chr(snp,freq=freq, min=0.5, max=0.95, k=c(0.95,0.05,0), a=0.06, number_per_cm=5 ),simplify = F)
pos<-unlist(sapply(r,function(x) x$pos))
geno<-unlist(sapply(r,function(x) x$geno))
data<-rbind(rbind(as.vector(geno[1:snp,]),as.vector(geno[1:snp+snp,]))+1,data)


t<-runHmmld(data,pair=c(1,2),pos=pos/1e6,chr=rep(1:nchr,each=snp))#estimate tracts of relatedness
plot(t)
t$kResult

# Kommentar: Resultatet ser fint ud men der er noget der undrer mig nemlig af for det simulrede dataset:
#> dim(data)
#[1]   108 10000
af = apply(data-1,2,mean)/2
range(af)
#[1] 0.3333333 0.6372549
par(mfrow=2:1)
hist(af)
hist(freq)
af[1:10]
# [1] 0.5441176 0.5245098 0.4558824 0.5147059 0.5098039 0.5441176 0.4362745
# [8] 0.4705882 0.5000000 0.4117647
freq[1:10]
# [1] 0.2889578 0.3849115 0.5655680 0.8673870 0.2315137 0.8585507 0.9002077
# [8] 0.6447180 0.6162026 0.1056076
plot(af,rep(freq,nchr))

# NB: de er langfra de samme som freq og ligger ml 0.33 og 0.63 saa jeg mistaenker der er noget galt med simkoden

# test 2: Rettet sim
set.seed(1)
n<-100 #number of individuals
snp=1000 #number of SNPs
freq<-runif(snp,min=0.05,max=0.95) #the population allele frequency of the SNPs (assumed to be uniform for SNP chip data)
data<-matrix(rbinom(n*snp*nchr,2,rep(freq,nchr*n)),ncol=snp*nchr,byrow=T)+1
r<-replicate(nchr,sim_chr(snp,freq=freq, min=0.5, max=0.95, k=c(0.95,0.05,0), a=0.06, number_per_cm=5 ),simplify = F)
pos<-unlist(sapply(r,function(x) x$pos))
geno<-unlist(sapply(r,function(x) x$geno))
data<-rbind(rbind(as.vector(geno[1:snp,]),as.vector(geno[1:snp+snp,]))+1,data)

# Nu ser af mere rigtig ud:
#> dim(data)
#[1]   108 10000
af = apply(data-1,2,mean)/2
range(af)
#[1] 0.02450980 0.98039216                                       
par(mfrow=2:1)
hist(af)
hist(freq)
af[1:10]
# [1] 0.3333333 0.4313725 0.5735294 0.8725490 0.2107843 0.8235294 0.9068627
# [8] 0.6176471 0.6225490 0.1176471
freq[1:10]
# [1] 0.2889578 0.3849115 0.5655680 0.8673870 0.2315137 0.8585507 0.9002077
# [8] 0.6447180 0.6162026 0.1056076
plot(af,rep(freq,nchr))

# men inferens resltaterne ser mindre fine ud...
t2<-runHmmld(data,pair=c(1,2),pos=pos/1e6,chr=rep(1:nchr,each=1000))#estimate tracts of relatedness
plot(t2)
t2$kResult
#    IBD=2     IBD=1     IBD=0
#0.0000000 0.3869699 0.6130301
t2$a
#[1] 0.15

