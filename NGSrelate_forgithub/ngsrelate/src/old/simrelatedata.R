library(Relate)

# find par for grandgrandson-grandgrandparent 
k2=0
k1=1/2^2
k0=1-k2-k1

calculate.a<-function(k0,k1,k2,phi=0.013){
    if(k2==0){
       ma<-1-log(k1)/log(2)
       mb=0
    }else{
       xa<-(k1+2*k2+sqrt((k1+2*k2)^2-4*k2))/2
       xb<-k2/xa
       ma<-1-log(xa)/log(2)
       mb<-1-log(xb)/log(2)
    }
    m<-ma+mb
    a<--m*log(1-phi)
    #cat("m=",m," ma=",ma," mb=",mb," a=",a," k0=",k0," k1=",k1," k2=",k2,"\n")
     return(a)
}

a <- calculate.a(k0,k1,k2)


n<-100 #number of individuals
snp=1000 #number of SNPs
freq<-runif(snp,min=0.025,max=0.975) #the population allele frequency of the SNPs (assumed to be uniform for SNP chip data)
data<-matrix(rbinom(n*snp*10,2,rep(freq,10)),ncol=snp*10)+1 # sim 10 times more snps than snp to get 10 chr for each ind witn snp snps on each
r<-replicate(10,sim_chr(snp,freq=freq, min=0.5, max=0.95, k=c(k0,k1,k2), a=a, number_per_cm=5 ),simplify = F)
pos<-unlist(sapply(r,function(x) x$pos))
geno<-unlist(sapply(r,function(x) x$geno))

data<-rbind(rbind(as.vector(geno[1:snp,]),as.vector(geno[1:snp+snp,]))+1,data)

#t<-runHmmld(data,pair=c(1,2),pos=pos/1e6,chr=rep(1:10,each=1000))#estimate tracts of relatedness using real genotypes
t<-runHmmld(data,pair=c(1,2),pos=pos/1e6,chr=rep(1:10,each=1000),back=0,ld_adj=F)
plot(t)


getLikes<-function(geno,dep,e=0.01,norm=FALSE){
  d<-dep
  n<-length(geno)
  dep<-rpois(n,d)
  nA<-rbinom(n,dep,c(e,0.5,1-e)[geno+1])
  res<-rbind(dbinom(nA,dep,e),dbinom(nA,dep,0.5),dbinom(nA,dep,1-e))
  if(norm)
    res<-t(t(res)/colSums(res))
  res
}

getPerfectLikes<-function(geno){
  res <- matrix(0,length(geno),3)
  res[,geno+1]=1
  res
}


perfectlike1<-getPerfectLikes(as.vector(geno[1:snp,]))
perfectlike2<-getPerfectLikes(as.vector(geno[1:snp+snp,]))

like1<-getLikes(as.vector(geno[1:snp,]),dep=2,norm=T)
like2<-getLikes(as.vector(geno[1:snp+snp,]),dep=10,norm=T)

ML1<-apply(like1,2,which.max)
ML2<-apply(like2,2,which.max)
data_likebased<-rbind(ML1,ML2,data)
t2<-runHmmld(data_likebased,pair=c(1,2),pos=pos/1e6,chr=rep(1:10,each=1000),back=0,ld_adj=F)#estimate tracts of relatedness based on called genos
par(mfrow=2:1)
plot(t)
plot(t2)


## measures of comparison
s = snp*10

## Using post:
## - fraction of correctly inferred sites (using post>0.95) -- 1 is best
for(postthres in c(0.75,0.9,0.95,0.99)){
cor = sum(unlist(lapply(r,function(x){x$state}))==apply(t$post ,2,function(x){if(max(x)>postthres){3-which.max(x)}else{NA}}),na.rm=T)/s

cor2 = sum(unlist(lapply(r,function(x){x$state}))==apply(t2$post,2,function(x){if(max(x)>postthres){3-which.max(x)}else{NA}}),na.rm=T)/s


## - fraction of incorrectly inferred sites (using post>0.95) -- 0 is best
incor = sum(unlist(lapply(r,function(x){x$state}))!=apply(t$post ,2,function(x){if(max(x)>postthres){3-which.max(x)}else{NA}}),na.rm=T)/s

incor2 = sum(unlist(lapply(r,function(x){x$state}))!=apply(t2$post,2,function(x){if(max(x)>postthres){3-which.max(x)}else{NA}}),na.rm=T)/s

## - fraction of NAs (using post>0.95) -- 0 is best
nas = sum(is.na(apply(t$post ,2,function(x){if(max(x)>postthres){3-which.max(x)}else{NA}})))/s

nas2 = sum(is.na(apply(t2$post,2,function(x){if(max(x)>postthres){3-which.max(x)}else{NA}})))/s

print("")
print(paste("Postthres",postthres))
print(c(cor,incor,nas))
print(c(cor2,incor2,nas2))
}

## Using viterbi
cor = sum(unlist(lapply(r,function(x){x$state}))==t$path)/s
cor2 = sum(unlist(lapply(r,function(x){x$state}))==t2$path)/s


## - fraction of incorrectly inferred sites (using post>0.95) -- 0 is best
incor = sum(unlist(lapply(r,function(x){x$state}))!=t$path)/s
incor2 = sum(unlist(lapply(r,function(x){x$state}))!=t2$path)/s

print("")
print("Viterbi")
print(c(cor,incor,0))
print(c(cor2,incor2,0))

