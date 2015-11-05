`pruning` <-  function(ld,snp,LD="rmisc",back=100,prune=0.2,abs=ifelse(LD=="D",TRUE,FALSE)){
  mea<-c()
  mea<-cbind(mea,t(ld[[LD]])[(snp-1):1,])
  if(abs)
    mea<-abs(mea)
  tal<-1
  keep<-TRUE
  if(back>1){
    for(i in 2:snp){
      if(sum(mea[tal,]>prune)>0){
        mea<-mea[-tal,]
        keep[i]<-FALSE
        if(tal>=dim(mea)[1])
          break
        if(back>2)
          for(i in 2:min(back-1,dim(mea)[1]-tal+2))
            mea[tal+i-2,i:back]<-c(mea[tal+i-2,(i+1):back],0)
        mea[min(tal+back-2,dim(mea)[1]),back]<-0
      }
      else{
        tal=tal+1
        keep[i]<-TRUE
      }
      if(tal>dim(mea)[1])
        break
    }
  }
  else{
    keep<-c(TRUE,mea<=prune)
  }
  return(keep)
}
