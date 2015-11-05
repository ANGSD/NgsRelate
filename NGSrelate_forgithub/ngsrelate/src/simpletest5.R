library(Relate)
set.seed(1)

# test 1: Anders basis kode (har kun aendret et par parameter vaerdier)

ffun<-function(k0=0.95,k1=0.05,k2=0,a=0.06){

  n <- 100
  snp <- 1000 #number of SNPs
  freq<-runif(snp,min=0.01,max=0.99) #the population allele frequency of the SNPs (assumed to be uniform for SNP chip data)
  data<-matrix(rbinom(n*snp,2,freq),ncol=snp,byrow=T)+1
  s<- sim_chr(snp,freq=freq, min=0.5, max=0.95, k=c(k0,k1,k2), a=a, number_per_cm=5 )#simulate chromosomes for a sib pair
  data<-rbind(t(s$geno)+1,data)
  return(list(data=data,pos=s$pos,freq=freq))
 }


pos<-c()
data<-c()
freq=c()
for(tal in 1:10){
   f<-ffun()
   pos<-c(pos,f$pos)
   freq <- c(freq,f$freq)
   data<-cbind(data,f$data)

 }
chr<-rep(1:10,each=1000)

# af blive nu fine
af = apply(data-1,2,mean)/2
plot(af,freq)

# men det goer resultaterne ikke...
t<-runHmmld(data,pair=c(1,2),pos=pos/1e6,chr=chr)#estimate tracts of relatedness
plot(t)
t$kResult


