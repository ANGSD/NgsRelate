source("ngsrelate/src/sim_chr.R")
source("ngsrelate/src/simrelatedatafuns.R")
source("ida2.R")

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



doStuff <- function(nind=20,nsnp=5e3,nchr=1,dep=8,err=0.005,file=NULL,k0=0.75,k1=0.25,k2=0.0,nsnpprcm=15,...){
  if(k0+k1+k2!=1)
    stop("sum of ks should be one")
  
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
  #source("batch1.R")

  
  d <- doStuff(nind=50,nsnp=5e3,nchr=22,file="batch1/po",k0=0,k1=1,k2=0.0,nsnpprcm=30)
  system("./NGSrelate -beagle batch1/po.beagle -freqfile batch1/po.freq")

  d <- doStuff(nind=50,nsnp=5e3,nchr=22,file="batch1/fs",k0=0.25,k1=0.50,k2=0.25,nsnpprcm=30)
  system("./NGSrelate -beagle batch1/fs.beagle -freqfile batch1/fs.freq")

  d <- doStuff(nind=50,nsnp=5e3,nchr=22,file="batch1/fc",k0=0.75,k1=0.25,k2=0.0,nsnpprcm=30)
  system("./NGSrelate -beagle batch1/fc.beagle -freqfile batch1/fc.freq")

  d <- doStuff(nind=50,nsnp=5e3,nchr=22,file="batch1/sc",k0=0.9375,k1=0.0625,k2=0,nsnpprcm=30)
  system("./NGSrelate -beagle batch1/sc.beagle -freqfile batch1/sc.freq")

  d <- doStuff(nind=50,nsnp=5e3,nchr=22,file="batch1/u",k0=1,k1=0,k2=0,nsnpprcm=30)
  system("./NGSrelate -beagle batch1/u.beagle -freqfile batch1/u.freq")

}

if(FALSE){
  nams <- c("po","fs","fc","sc","u")
  nams <- nams[-5]
  fun <- function(x){
    load(paste("batch1/",x,".RData"))
    sampl <- res$sampleTrue
    est <- read.table(paste("batch1/",x,".beagle.res"))[,-c(1,5)]
    res <- sampl-est
    colnames(res) <- c("k0","k1","k2")
    invisible(boxplot(res,ylim=c(-0.5,0.5),main=x,col=1:3))
    0
  }
  pdf("batch1.pdf")
  par(mfrow=c(2,2))
  lapply(nams,fun)
  dev.off()
}
