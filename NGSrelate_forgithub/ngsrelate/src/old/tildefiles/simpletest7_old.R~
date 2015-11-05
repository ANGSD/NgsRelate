set.seed(1)
library(Relate)
source("sim_chr.R")

# test 1: Anders basis kode (har kun aendret et par parameter vaerdier)

ffun<-function(){
  n<-100 #number of individuals

  snp=1000 #number of SNPs

  freq<-runif(snp,min=0.05,max=0.95) #the population allele frequency of the SNPs (assumed to be uniform for SNP chip data)

  data<-matrix(rbinom(n*snp,2,freq),ncol=snp,by=TRUE)

  s<- sim_chr_new(snp,freq=freq, min=0.5, max=0.95, k=c(0.95,0.05,0), a=0.06, number_per_cm=5 )#simulate chromosomes for a sib pair

  data<-rbind(t(s$geno)+1,data+1)

  return(list(data=data,pos=s$pos,freq=freq,states=s$state))

}



pos<-c()
data<-c()
freq <- c()
state <- c()
statelist <- list()

for(tal in 1:10){
  f <- ffun()
  pos<-c(pos,f$pos)
  posx<-c(pos,f$pos+rev(pos)[1])
  data<-cbind(data,f$data)
  freq <- c(freq,f$freq)
  state <- c(state,f$state)
  statelist[[tal]] <- f$states
}
chr<-rep(1:10,each=1000)

# af blive nu fine
af = apply(data-1,2,mean)/2
#plot(af,freq)

# men det goer resultaterne ikke...
t<-runHmmld(data,pair=c(1,2),pos=pos/1e6,chr=chr)#estimate tracts of relatedness
t$kResult
pdf("test7.pdf")
plot(t)
nsnp=1000
for(chr in 1:10){
  plot(t,chr=chr)
  points(pos[((chr-1)*nsnp+1):(chr*nsnp)]/1e6,rep(.5,nsnp),col=-statelist[[chr]]+3,pch=20)
}
graphics.off()
