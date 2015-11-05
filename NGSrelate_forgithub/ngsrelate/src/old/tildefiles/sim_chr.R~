 
sim_chr <- function(snp,freq=0.2, min=0.5, max=0.95, k=k, a=3.854, number_per_cm=1){
  CM=1e6
  unrelated <- function(snp, freq, min, max){
    if (!is.null(freq))
      p <- sample(freq, snp, replace=TRUE)
    else 
      p <- runif(snp, min, max)
    prob <- t(matrix( c(p,1-p), nrow=length(p)))
    apply(prob, 2, sample, x=c(1,0), size=1, replace=TRUE)
  }

  recom<-function(pos,r=3.854,k,a=NULL){
    if(is.null(a))
      a<-0.1
  
    k0<-k[1]
    k1<-k[2]
    k2<-k[3]
    em<-function(state,t,k0,k1,k2,a){
      IBD01<-(1-exp(-a*t))*k1
      IBD02<-exp(-a*k1*t)*k2/(k1-1)+exp(-a*t)*k1+exp(-a*t)*k0*k1/(k1-1)+k2
      IBD02[IBD02<1e-15]<-0 #underflow problem loesning
      IBD00<-1-IBD01-IBD02
      if(k1==0){
        IBD10<-IBD01
        IBD12<-IBD01
      }
      else{
        IBD10<-IBD01*k0/k1
        IBD12<-IBD01*k2/k1
      }
      IBD11<-1-IBD10-IBD12
      IBD21<-IBD01
      IBD20<-exp(-a*k1*t)*k0/(k1-1)+exp(-a*t)*k1+exp(-a*t)/(k1-1)*k2*k1-(1-k1)*k0/(k1-1)
      IBD20[IBD20<1e-15]<-0 #underflow problem loesning
      IBD22<-1-IBD21-IBD20
      if(state==0)
        prob<-c(IBD00,IBD01,IBD02)
      else if(state==1)
        prob<-c(IBD10,IBD11,IBD12)
      else
        prob<-c(IBD20,IBD21,IBD22)
      return(prob)
    }
  t<-diff(pos)
  state<-sample(0:2,1,prob=c(k0,k1,k2)+0.0000001)
  for(tal in 2:length(pos))  
    state[tal]<-sample(0:2,1,prob=em(state[tal-1],t[tal-1],k0,k1,k2,a)+0.00000001)
  return(state)
  }








  hap11 <- unrelated(snp,freq,min,max)
  hap12 <- unrelated(snp,freq,min,max)
  hap21 <- unrelated(snp,freq,min,max)
  hap22 <- unrelated(snp,freq,min,max)
  pos <- as.integer(rnorm(snp, mean=CM, sd=CM/10)) ## interval between positions
  pos<-round(pos/number_per_cm,0)
  pos<-cumsum(pos)
  state=recom(pos/CM,k=k,a=a)

  hap11[state==2]<-hap21[state==2]
  hap12[state==2]<-hap22[state==2]
  hap11[state==1]<-hap21[state==1]
  geno <- cbind(hap11+hap12,hap21+hap22)
  obs <- list(geno=geno, state=state, k=k, a=a, freq=freq,
              number_per_cm=number_per_cm, min=min, max=max,
              snp=snp, pos=pos)
  class(obs) = "sim_chr"
  return(obs)
}
