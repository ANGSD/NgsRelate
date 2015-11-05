source("ngsrelate/src/sim_chr.R")
source("ngsrelate/src/simrelatedatafuns.R")
source("ida2.R")



getSampleks <- function(x){
  dat <- x
  pos <- dat$pos
  
  tp<-dat$truepath
  
  splits<-which(tp[-1]!=tp[-length(tp)])+1
  diffs <- diff(c(pos[c(1,splits)],tail(pos,1)+1))
  
  r <- cbind(c(tp[splits-1],tp[length(tp)]),diffs)
  
  span <- diff(range(pos))+1
  
  ks <- c(sum(r[r[,1]==0,2]),sum(r[r[,1]==1,2]),sum(r[r[,1]==2,2]))/span
  names(ks)<-paste("k",0:2)
  ks
}



readStuff <- function(file="test3.beagle.bres.gz"){
  ff <- gzfile(file,"rb")
  strLen <- readBin(ff,integer(),1)
  chrNam <- readBin(ff,character(),n=strLen,size=1)
  nSites <- readBin(ff,integer(),1)
  pars <- readBin(ff,numeric(),4)
  liks <- readBin(ff,numeric(),3)
  pos <- readBin(ff,integer(),nSites)
  post <- t(matrix(readBin(ff,numeric(),nSites*3),nrow=3,byrow=T))
  forward <- matrix(readBin(ff,numeric(),nSites*3),nrow=3,byrow=T)
  backward <- matrix(readBin(ff,numeric(),nSites*3),nrow=3,byrow=T)
  emis <- matrix(readBin(ff,numeric(),nSites*3),nrow=3,byrow=T)
  trans <- matrix(readBin(ff,numeric(),(nSites+1)*9),nrow=9,byrow=T)
  viterbi <- as.integer(readBin(ff,what="raw",nSites))


  dec=list(path=viterbi,post=post[,3:1],forward=forward,backward=backward)
  class(dec) <- "decode"
  close(ff)
  
  relatedness<-pars[4]+pars[3]/2
  a <- pars[1]
  k <- rev(pars[-1])
  obj <- list(k=k,k.loglike=liks[1],relatedness=relatedness,a=pars[1],u.loglike=liks[2],po.loglike=liks[3],snp=nSites,position=pos/1e6,logE=emis,Tp=trans,viterbi=viterbi,chr=rep(1,nSites))
#conv=conv,LD=LD,t=t,alim=alim,af=af,chr=chr,conv=conv,geno=ind1
  class(obj) <- "relateHMM"
  obj$decode = dec
  return(obj)
}

makePerfect <- function(x){
  tmp <- getPerfectLikes(x$dat$geno-1)
  t(matrix(tmp,nrow=nrow(x$dat$geno)*3))
}



doStuff <- function(nind=20,nsnp=5e3,nchr=1,dep=8,err=0.005,file=NULL,k0=0.75,k1=0.25,k2=0.0,nsnpprcm=15,...){
  ##simulate data
  dat <-simData(nind=nind,nsnp=nsnp,nchr=nchr,k0=k0,k1=k1,nsnpprcm=nsnpprcm,...)
  ##simulate GLS
  gl <- getLikes(dat$geno-1,dep=dep,e=err,norm=T,...)
  gls <- t(matrix(gl,nrow=nrow(dat$geno)*3))

  ##function to generate a beagle file
  get.beagle <- function(dat,gls){
    nams <- paste(dat$chr,dat$pos,sep="_")
    major <- rep(0,nsnp)
    minor <- rep(1,nsnp)
    d <- data.frame(nams,major,minor,gls)
    colnames(d) <- c(colnames(d)[1:3],1:(3*nrow(dat$geno)))
    return(d)
  }

  ##dump a beagle file if file is specified
  if(!is.null(file)){
    cat("-> Writing file: ",file,"\n")
    write.table(get.beagle(dat,gls),file=file,col.names=TRUE,quote=FALSE,row.names=F)
  }
  

  freq <- apply(dat$geno-1,2,mean)/2
  ks <- getSampleks(dat)
  return(list(dat=dat,gls=gls,freq=freq,sampleTrue=ks))
 
}

set.seed(0)
d <- doStuff(nind=50,nsnp=50e3,file="ida.beagle",k0=0.75,k1=0.25,k2=0.0,nsnpprcm=1000)

#run like

system("./NGSrelate -beagle ida.beagle")
system("./NGSrelate -beagle ida.beagle -freqfile ida.beagle.freq")
cstuff <- readStuff(file="test3.beagle.bin.gz.bres.gz")


##run old versions

library(Relate)
##genolikesPerfect <- t(makePerfect(d)[,1:6])
genolikesNonPerfect <- t(d$gls[,1:6])


#ra1<-runHmmld(d$dat$geno,pair=c(1,2),pos=d$dat$pos/1e6,chr=d$dat$chr,back=0,ld_adj=F,epsilon=0.01,par=c(0.02,0,0.25,0.75),calc.a=F,giveCrap=T )
af<-scan("test3.beagle.freq")
source("ida3.R")


ri1x<-ngsrelateHMM(genolikesNonPerfect,pair=c(1,2),pos=d$dat$pos/1e6,chr=d$dat$chr,alim=c(0.001,0.15),epsilon=0.01,par=c(0.046,0,0.21,0.78),calc.a=F,af=af)


cstuff
ri1x


d <- doStuff(nind=2,nsnp=150e3,file=NULL,k0=0.75,k1=0.25,k2=0.0)

getSampleks(d)

