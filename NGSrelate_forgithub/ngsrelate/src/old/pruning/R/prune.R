 
prune<-function(data,prune=0.1,LD="D",back=50){
  
  #error messages
  if(class(data)!="matrix")
    stop("data must be a matrix")
  snp<-dim(data)[2]
  d<-as.integer(data)
  d[is.na(d)]<-0
  d<-as.integer(d)
  data<-matrix(d,ncol=snp)


  
  LD_var <- ifelse(LD=="rsq2",TRUE,FALSE)
 
  res<-.Call(interface,data,prune,as.integer(LD_var),as.numeric(back),PACKAGE="pruning")
  snp<-dim(res$data)[2]
  d<-as.integer(res$data)
  d[d==0]<-NA
  d<-as.integer(d)
  res$data<-matrix(d,ncol=snp)
  
  class(res)<-"pruning"
  return(res)
}
