source("ngsrelate/src/sim_chr.R")
source("ngsrelate/src/simrelatedatafuns.R")
source("ida2.R")
library(multicore)
getSampleks <- function(x){
  dat <- x
  posAll <- dat$pos
  chr <- dat$chr
  tpAll<-dat$truepath
  res <- c()
  for(i in unique(chr)){
    tp <- tpAll[chr==i]
    splits<-which(tp[-1]!=tp[-length(tp)])+1
    pos <- posAll[chr==i]
    diffs <- diff(c(pos[c(1,splits)],tail(pos,1)+1))
    r <- cbind(c(tp[splits-1],tp[length(tp)]),diffs)
    span <- diff(range(pos))+1
  
    ks <- c(sum(r[r[,1]==0,2]),sum(r[r[,1]==1,2]),sum(r[r[,1]==2,2]))/span
    res <- rbind(res,ks)
  }
  colnames(res)<-paste("k",0:2)
  rownames(res) <- unique(chr)
  res
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



doStuff <- function(nind=20,nsnp=5e3,nchr=1,dep=8,err=0.005,file=NULL,ks=NULL,k2=0.0,nsnpprcm=15,...){
  if(is.null(ks)||length(ks)!=3||sum(ks)!=1)
    stop("ks should be a vector of length3 summing to one")

  k0=ks[1]
  k1=ks[2]
  k2=ks[3]

  
  
  ##simulate data
  dat <-simData(nind=nind,nsnp=nsnp,nchr=nchr,k0=k0,k1=k1,nsnpprcm=nsnpprcm,...)
  dat$truepaths <- as.numeric(t(dat$truepaths))
  freq <- apply(dat$geno-1,2,mean)/2
  
  ##remove invar
  if(TRUE){
    keep <- freq>0&freq<1
    dat$geno <- dat$geno[,keep]
    dat$pos <- dat$pos[keep]
    dat$chr <- dat$chr[keep]
    dat$truepaths <- dat$truepaths[keep]
    freq <- freq[keep]
  }
  ##simulate GLS
  gl <- getLikes(dat$geno-1,dep=dep,e=err,norm=T,...)
  gls <- t(matrix(gl,nrow=nrow(dat$geno)*3))

  
  ks <- getSampleks(dat)
  res <- list(dat=dat,gls=gls,freq=freq,sampleTrue=ks)



  
  ##function to generate a beagle file
  get.beagle <- function(dat,gls){
    nams <- paste(dat$chr,dat$pos,sep="_")
    major <- rep(0,length(dat$chr))
    minor <- rep(1,length(dat$chr))
    d <- data.frame(nams,major,minor,gls)
    colnames(d) <- c(colnames(d)[1:3],1:(3*nrow(dat$geno)))
    return(d)
  }

  ##dump a beagle file if file is specified
  if(!is.null(file)){
    aName <- paste(file,".beagle",sep="")
    cat("-> Writing file: ",aName,"\n")
    write.table(get.beagle(dat,gls),file=aName,col.names=TRUE,quote=FALSE,row.names=F)
  }
  ##dump the RData object
  if(!is.null(file)){
    aName <- paste(file,".RData",sep="")
    cat("-> Writing file: ",aName,"\n")
    save(res,file=aName)
  }
  ##dump the freqs
  if(!is.null(file)){
    aName <- paste(file,".freq",sep="")
    cat("-> Writing file: ",aName,"\n")
    write.table(freq,file=aName,col.names=FALSE,quote=FALSE,row.names=F,eol=' ')
  }  
  
  return(res)
 
}

if(FALSE){
  set.seed(0)
  #source("batch2.R")
  nChr <- 22
  nRep <- 20
  dname <- "batch2/"

  po <- c(0,1,0)
  fs <- c(0.25,0.5,0.25)
  fc <- c(0.75,0.25,0)
  sc <- c(0.9375,0.0625,0)
  simArg <- rbind(po,fs,fc,sc)
  
  for(bname in rownames(simArg)){
    d <- mclapply(1:nRep,function(x) doStuff(nind=50,nsnp=5e3,nchr=nChr,file=paste(dname,bname,x,sep=''),ks=simArg[bname,],nsnpprcm=30)$sampleTrue)
    write.table(t(sapply(d,function(x) apply(x,2,mean))),file=paste(dname,bname,".true",sep=''))
  }

  if(FALSE){
    system("./NGSrelate -beagle batch2/po.beagle -freqfile batch2/po.freq -k2 0")
    system("./NGSrelate -beagle batch2/fs.beagle -freqfile batch2/fs.freq -calcA 0")
    system("./NGSrelate -beagle batch2/fc.beagle -freqfile batch2/fc.freq -k2 0")
    system("./NGSrelate -beagle batch2/sc.beagle -freqfile batch2/sc.freq -k2 0")
    system("./NGSrelate -beagle batch2/sc.beagle -freqfile batch2/sc.freq -k2 0 -pair1 3 -pair2 4")
  }

}

if(FALSE){

  calcA <- grep("calc",grep("all",list.files("batch2","*.res",full.names=T),val=T),val=T)
  all<- grep("calc",grep("all",list.files("batch2","*.res",full.names=T),val=T),val=T,invert=T)
  TT <- list.files("batch2","*.true",include.dirs=T,full.names=T)

  calcA <- lapply(calcA,function(x) as.matrix(read.table(x)[,2:4]))
  all <- lapply(all,function(x) as.matrix(read.table(x)[,2:4]))
  TT <- lapply(TT,function(x) as.matrix(read.table(x)))
  TT[[5]]<-matrix(rep(c(1,0,0),20),ncol=3,byrow=T)


  nams <- c("fc","fs","po","sc","u")
  tmp1<- lapply(1:5,function(x) {tmp<-TT[[x]]-calcA[[x]];colnames(tmp)<-paste("k",0:2);tmp})
  tmp2<- lapply(1:5,function(x) {tmp<-TT[[x]]-all[[x]];colnames(tmp)<-paste("k",0:2);tmp})

  brewer()
  par(mfrow=c(2,5))
  lapply(1:5,function(x) boxplot(tmp1[[x]],main=nams[x],col=1:3))
  lapply(1:5,function(x) boxplot(tmp2[[x]],main=nams[x],col=1:3))
  

  
}
