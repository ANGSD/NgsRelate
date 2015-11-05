#####################################
# Output functions
######################################

print.relateHMM<-function(x){
	if (!inherits(x, "relateHMM")) 
        stop("Not an object of class relateHMM!")
	print(structure(list(k=round(x$k,2),k.loglike=round(x$k.loglike,3),
						 relatedness=signif(x$relatedness,3),a=signif(x$a,3),
						 u.loglike=round(x$u.loglike,1),po.loglike=round(x$po.loglike,1)),
						 class="power.htest"))
}

plot.relateHMM<-function(x,col=1:3,lwd=2,chr=NULL,...){
post<-x$decode$post
path<-x$decode$path
pos<-x$position
if(is.null(x$chr)){
  plot(pos[1:(x$snp-1)],post[2:(x$snp),1],type="l",ylim=c(0,1),col=col[1],lwd=lwd,...)
  lines(pos[1:(x$snp-1)],post[2:(x$snp),2],col=col[2],lwd=lwd)
  lines(pos[1:(x$snp-1)],post[2:(x$snp),3],col=col[3],lwd=lwd)
  points(pos[1:(x$snp-1)],rep(1,x$snp-1),col=3-path[-(x$snp-1)],pch="|")
}
else{
  if(!is.null(chr)){
    pos<-pos[x$chr==chr]
    post<-post[x$chr==chr,]
    plot(pos[-length(pos)],post[-1,1],type="l",ylim=c(0,1),col=col[1],lwd=lwd,...)
    lines(pos[-length(pos)],post[-1,2],col=col[2],lwd=lwd)
    lines(pos[-length(pos)],post[-1,3],col=col[3],lwd=lwd)
    
  }
  else{
    C<-names(table(x$chr))
    m<-c(0,cumsum(tapply(pos,x$chr,max)))
    pos2<-rep(NA,length(pos))
    for(tal in 1:length(C))
      pos2[x$chr==C[tal]]<-pos[x$chr==C[tal]]+m[tal]
    plot(pos2[1:(x$snp-1)],post[2:(x$snp),1],type="l",ylim=c(0,1),col=col[1],lwd=lwd,...)
    lines(pos2[1:(x$snp-1)],post[2:(x$snp),2],col=col[2],lwd=lwd)
    lines(pos2[1:(x$snp-1)],post[2:(x$snp),3],col=col[3],lwd=lwd)
    abline(v=m)
    for(tal in 1:length(C))
      text(m[tal]+diff(m)[tal]/2,0.5,C[tal],col="gray")
  }
}
invisible(cbind(post,pos))
}

######################################
# Help functions
######################################

# Function for calculating a (from Albrechtsen et al. 2009)
# ------------------------------------------------------------

calculate.a<-function(k0,k1,k2,phi=0.013){
  if(k2==0){
    ma<-1-log(k1)/log(2)
    mb=0
  }
  else{
    xa<-(k1+2*k2+sqrt((k1+2*k2)^2-4*k2))/2
    xb<-k2/xa
    ma<-1-log(xa)/log(2)
    mb<-1-log(xb)/log(2)
  }
  m<-ma+mb
  a<--m*log(1-phi)
  return(a)
}

# Function for making an error table (from Albrechtsen et al. 2009)
# --------------------------------------------------------------------

error_old<-function(obsgeno,epsilon,snp){
	ep<-matrix(NA,ncol=3,nrow=snp)
	
	ep[obsgeno==0,1]<-(1-epsilon)^2
	ep[obsgeno==1,1]<-(1-epsilon)*epsilon
	ep[obsgeno==2,1]<-(epsilon)^2
	
	ep[obsgeno==0,2]<-2*(1-epsilon)*epsilon
	ep[obsgeno==1,2]<-(1-epsilon)^2+epsilon^2
	ep[obsgeno==2,2]<-2*(1-epsilon)*epsilon
	
	ep[obsgeno==0,3]<-(epsilon)^2
	ep[obsgeno==1,3]<-(1-epsilon)*epsilon
	ep[obsgeno==2,3]<-(1-epsilon)^2
	
	return(ep)
}

error_new<-function(obsgeno,epsilon,snp){
	ep<-matrix(NA,ncol=3,nrow=snp)
	
	ep[obsgeno==0,1]<-(1-epsilon)^2
	ep[obsgeno==1,1]<-(1-epsilon)*epsilon
	ep[obsgeno==2,1]<-(epsilon)^2
	
	ep[obsgeno==0,2]<-2*(1-epsilon)*epsilon
	ep[obsgeno==1,2]<-(1-epsilon)^2+epsilon^2
	ep[obsgeno==2,2]<-2*(1-epsilon)*epsilon
	
	ep[obsgeno==0,3]<-(epsilon)^2
	ep[obsgeno==1,3]<-(1-epsilon)*epsilon
	ep[obsgeno==2,3]<-(1-epsilon)^2
	
	return(ep)
}


# Functions for calculating emissions
# --------------------------------------

emission_relateHMM_old<-function(freq,ep1,ep2){
# Input:   freq is freqs for all loci in the pop, 
#	   ep1 and ep2 are matrices with P(G_obs|G_real,epsilon) (a col per posible G_real, a row per observed locus) 
# Output:  a matrix E with a col per posible hidden state X (# alleles IBD), a row perlocus 
#          E[i,j] = P(G_obs in i|X=(j-1)) = sum over G_real P(G_real|X=(j-1))P(G_obs|G_real,epsilon)

	nsnps = length(freq)
	freqA = freq
	freqa = 1-freq	

    	# Summing over over all 9 possible combinations of (G^ind1_real,G^ind2_real) 
	E<-cbind(freqA^4,freqA^3,freqA^2)*ep1[,1]*ep2[,1]				# G_real=(AA,AA) 
	E<-E+cbind(4*freqA^3*freqa,2*freqA^2*freqa,0)*ep1[,1]*ep2[,2]			# G_real=(AA,Aa)
	E<-E+cbind(2*freqA^2*freqa^2,0,0)*ep1[,1]*ep2[,3]				# G_real=(AA,aa)
	E<-E+cbind(4*freqA^3*freqa,2*freqA^2*freqa,0)*ep1[,2]*ep2[,1]			# G_real=(Aa,AA)
	E<-E+cbind(4*freqA^2*freqa^2,freqA*freqa,2*freqA*freqa)*ep1[,2]*ep2[,2]		# G_real=(Aa,Aa)
	E<-E+cbind(4*freqA*freqa^3,2*freqA*freqa^2,0)*ep1[,2]*ep2[,3]			# G_real=(Aa,aa)
	E<-E+cbind(2*freqA^2*freqa^2,0,0)*ep1[,3]*ep2[,1]				# G_real=(aa,AA)
	E<-E+cbind(4*freqA*freqa^3,2*freqA*freqa^2,0)*ep1[,3]*ep2[,2]			# G_real=(aa,Aa)
	E<-E+cbind(freqa^4,freqa^3,freqa^2)*ep1[,3]*ep2[,3]				# G_real=(aa,aa)

        # Fixing bug...
        Etmp<-E
        E[,1]<-Etmp[,3]
        E[,3]<-Etmp[,1]
    	return(t(E))
}



emission_relateHMM_new<-function(freq,ep1,ep2){
# Input:   freq is freqs for all loci in the pop, 
#	   ep1 and ep2 are matrices with P(G_obs|G_real,epsilon) (a col per posible G_real, a row per observed locus) 
# Output:  a matrix E with a col per posible hidden state X (# alleles IBD), a row perlocus 
#          E[i,j] = P(G_obs in i|X=(j-1)) = sum over G_real P(G_real|X=(j-1))P(G_obs|G_real,epsilon)
  #print(length(freq))
##  stop("yoyo")
	nsnps = length(freq)
	freqA = freq
	freqa = 1-freq	
	
	# Summing over over all 9 possible combinations of (G^ind1_real,G^ind2_real) 
	E<-cbind(freqA^4,freqA^3,freqA^2)*ep1[,1]*ep2[,1]				# G_real=(AA,AA) 
	E<-E+cbind(2*freqA^3*freqa,1*freqA^2*freqa,0)*ep1[,1]*ep2[,2]			# G_real=(AA,Aa)
	E<-E+cbind(1*freqA^2*freqa^2,0,0)*ep1[,1]*ep2[,3]				# G_real=(AA,aa)
	E<-E+cbind(2*freqA^3*freqa,1*freqA^2*freqa,0)*ep1[,2]*ep2[,1]			# G_real=(Aa,AA)
	E<-E+cbind(4*freqA^2*freqa^2,freqA*freqa,2*freqA*freqa)*ep1[,2]*ep2[,2]		# G_real=(Aa,Aa)
	E<-E+cbind(2*freqA*freqa^3,1*freqA*freqa^2,0)*ep1[,2]*ep2[,3]			# G_real=(Aa,aa)
	E<-E+cbind(1*freqA^2*freqa^2,0,0)*ep1[,3]*ep2[,1]				# G_real=(aa,AA)
	E<-E+cbind(2*freqA*freqa^3,1*freqA*freqa^2,0)*ep1[,3]*ep2[,2]			# G_real=(aa,Aa)
	E<-E+cbind(freqa^4,freqa^3,freqa^2)*ep1[,3]*ep2[,3]				# G_real=(aa,aa)

	# Fixing bug...
        Etmp<-E
        E[,1]<-Etmp[,3]
        E[,3]<-Etmp[,1]
    	return(t(E))
}

emission_ngsrelateHMM_new<-function(freq,l1,l2){
# Input:   freq is freqs for all loci in the pop, 
#	   l1 and l2 are matrices with P(G_obs|G_real) (a col per posible G_real, a row per observed locus) (these are genotype likes!) 
# Output:  a matrix E with a col per posible hidden state X (# alleles IBD), a row perlocus 
#          E[i,j] = P(G_obs in i|X=(j-1)) = sum over G_real P(G_real|X=(j-1))P(G_obs|G_real)
	
	nsnps = length(freq)
	freqA = freq
	freqa = 1-freq	
	
	# Summing over over all 9 possible combinations of (G^ind1_real,G^ind2_real) 
	E<-cbind(freqA^4,freqA^3,freqA^2)*l1[,1]*l2[,1]						# G_real=(AA,AA) 
	E<-E+cbind(2*freqA^3*freqa,1*freqA^2*freqa,0)*l1[,1]*l2[,2]				# G_real=(AA,Aa)
	E<-E+cbind(1*freqA^2*freqa^2,0,0)*l1[,1]*l2[,3]						# G_real=(AA,aa)
	E<-E+cbind(2*freqA^3*freqa,1*freqA^2*freqa,0)*l1[,2]*l2[,1]				# G_real=(Aa,AA)
	E<-E+cbind(4*freqA^2*freqa^2,freqA*freqa,2*freqA*freqa)*l1[,2]*l2[,2]			# G_real=(Aa,Aa)
	E<-E+cbind(2*freqA*freqa^3,1*freqA*freqa^2,0)*l1[,2]*l2[,3]				# G_real=(Aa,aa)
	E<-E+cbind(1*freqA^2*freqa^2,0,0)*l1[,3]*l2[,1]						# G_real=(aa,AA)
	E<-E+cbind(2*freqA*freqa^3,1*freqA*freqa^2,0)*l1[,3]*l2[,2]				# G_real=(aa,Aa)
	E<-E+cbind(freqa^4,freqa^3,freqa^2)*l1[,3]*l2[,3]					# G_real=(aa,aa)
	
    	return(t(E))
}



# Function for calculating transitions (from Albrechtsen et al 2009)
# --------------------------------------------------------------------

transition_relateHMM<-function(a,t,k0,k1,k2){

	# -- trans x=0 -> 
	T01<-(1-exp(-a*t))*k1  
	T02<-exp(-a*k1*t)*k2/(k1-1)+exp(-a*t)*k1+exp(-a*t)*k0*k1/(k1-1)+k2 # follows T 
	if(k1==1) # unless we know k2=0, in which case it is 0  
		T02<-rep(0,length(T02))
	T02[T02<1e-15]<-0 # solving underflow problem
	T00<-1-T01-T02  

	# -- trans x=1 -> 
	if(k1==0){
		T10<-T01 # not relevant since x!=1 so set to 0
		T12<-T01 # not relevant since x!=1 so set to 0
	}else{
		T10<-T01*k0/k1 
		T12<-T01*k2/k1 
	}   
	T11<-1-T10-T12 

	# -- trans x=2 -> 
	T21<-T01 
	T20<-exp(-a*k1*t)*k0/(k1-1)+exp(-a*t)*k1+exp(-a*t)*k2*k1/(k1-1)-(1-k1)*k0/(k1-1) # why this trick?
	T20[T20<1e-15]<-0 #underflow problem loesning
	if(k1==1)
		T20<-rep(0,length(T20))
	T22<-1-T21-T20

        # To avoid underflow errors
        llim=1e-15
        T00[abs(T00)<llim]=0
        T01[abs(T01)<llim]=0
        T02[abs(T02)<llim]=0
        T10[abs(T10)<llim]=0
        T11[abs(T11)<llim]=0
        T12[abs(T12)<llim]=0
        T20[abs(T20)<llim]=0
        T21[abs(T21)<llim]=0
        T22[abs(T22)<llim]=0
        
	Tp = list(T00=T00,T01=T01,T02=T02,T10=T10,T11=T11,T12=T12,T20=T20,T21=T21,T22=T22)	
	return(Tp)
}	


# Calculating likelihood (using the forward algorithm)
# --------------------------------------------------------------------

forward_relateHMM<-function(k0,k1,k2,logE,Tp){
	
	# - getting parameter values ready
	snp = dim(logE)[2]
	T00 = Tp$T00
	T01 = Tp$T01
	T02 = Tp$T02
	T10 = Tp$T10	
	T11 = Tp$T11
	T12 = Tp$T12
	T20 = Tp$T20
	T21 = Tp$T21
	T22 = Tp$T22	

	# - init 
	logf0<-c()
	logf1<-c()
	logf2<-c()
	logf0[1]<- log(k0)+logE[3,1]
	logf1[1]<- log(k1)+logE[2,1]
	logf2[1]<- log(k2)+logE[1,1]
     #   cat("f0=",logf0[1],"f1=",logf1[1],"f2=",logf2[1],"\n")
	m<-0
        
	# - recursion
	for (i in (2:snp)){
		f0m<-exp(logf0[i-1]-m)
		f1m<-exp(logf1[i-1]-m)
		f2m<-exp(logf2[i-1]-m)
                
		logf0[i]<-m+log(T00[i]*f0m+T10[i]*f1m+T20[i]*f2m)+logE[3,i]
		logf1[i]<-m+log(T01[i]*f0m+T11[i]*f1m+T21[i]*f2m)+logE[2,i]
		logf2[i]<-m+log(T02[i]*f0m+T12[i]*f1m+T22[i]*f2m)+logE[1,i]
		m<-max(logf0[i],logf1[i],logf2[i])
                ##print(m)

	}

	# - termination
	loglike = m+log(exp(logf0[snp]-m) + exp(logf1[snp]-m)+ exp(logf2[snp]-m))
        cat("loglike in forward=",loglike,"\n")
	# NB: the -m+m is a trick (add protect) to avoid underflow 
	
	# - return both the final like and the values of the forward variables
	return(list(loglike=loglike,logf0=logf0,logf1=logf1,logf2=logf2))
}

fastlike_relateHMM<-function(k0,k1,k2,logE,Tp){
#  cat("pars givine k0: ",k0,"k1: ",k1,"k2: ",k2,"\n");
	# - getting parameter values ready
	snp = dim(logE)[2]
	T00 = Tp$T00
	T01 = Tp$T01
	T02 = Tp$T02
	T10 = Tp$T10	
	T11 = Tp$T11
	T12 = Tp$T12
	T20 = Tp$T20
	T21 = Tp$T21
	T22 = Tp$T22	

	# - init
	log_IBD0<- log(k0)
	log_IBD1<- log(k1)
	log_IBD2<- log(k2)
  #cat("pars init: ",log_IBD0,"k1: ",log_IBD1,"k2: ",log_IBD2,"\n");
	m<-0
	
	# - recursion
  for (i in (2:snp)){
##  for (i in (2:14)){
          logp0<-exp(log_IBD0+logE[3,i-1]-m)
          logp1<-exp(log_IBD1+logE[2,i-1]-m)
          logp2<-exp(log_IBD2+logE[1,i-1]-m)
 ##         tmp <- c(logp0,logp1,logp2)
  ##        print(tmp)

          log_IBD0<-m+log(T00[i]*logp0+T10[i]*logp1+T20[i]*logp2)
          log_IBD1<-m+log(T01[i]*logp0+T11[i]*logp1+T21[i]*logp2)
          log_IBD2<-m+log(T02[i]*logp0+T12[i]*logp1+T22[i]*logp2)
          m<-max(log(logp1)+m,log(logp2)+m,log(logp0)+m)
    ##      cat("m i: ",i-1,"val: ",m,"\n")
          ##             stop()
          if(is.na(m)|is.na(log_IBD0)){
            cat("error in likelihood\n")
            cat("pars givine k0: ",k0,"k1: ",k1,"k2: ",k2,"\n");
            stop("asdfasdfasdf")
            print(delta)
            break
          }
	}
	logp0 = logE[3,snp] + log_IBD0
	logp1 = logE[2,snp] + log_IBD1
	logp2 = logE[1,snp] + log_IBD2
	
	m = max(logp0, logp1, logp2)
        print(m)
	loglike = m+log(exp(logp0-m) + exp(logp1-m)+ exp(logp2-m))
	
	return(list(loglike=loglike))
}

# Calculating likelihood (using the backward algorithm)
# --------------------------------------------------------------------

backward_relateHMM<-function(k0,k1,k2,logE,Tp){

	# - getting parameter values ready
	snp = dim(logE)[2]
	T00 = Tp$T00
	T01 = Tp$T01
	T02 = Tp$T02
	T10 = Tp$T10	
	T11 = Tp$T11
	T12 = Tp$T12
	T20 = Tp$T20
	T21 = Tp$T21
	T22 = Tp$T22	
	
	# - init
	logb0<-c()
	logb1<-c()
	logb2<-c()
	logb0[snp]<-0
	logb1[snp]<-0
	logb2[snp]<-0	
	m<-0

	# - recursion 
	for (i in ((snp-1):1)){
          p0m<-exp(logE[3,i+1]+logb0[i+1]-m)
          p1m<-exp(logE[2,i+1]+logb1[i+1]-m)
          p2m<-exp(logE[1,i+1]+logb2[i+1]-m)
##          print(c(i, p0m,p1m,p2m))

          logb0[i]<-m+log(p2m*T02[i+1]+p1m*T01[i+1]+p0m*T00[i+1])
          logb1[i]<-m+log(p2m*T12[i+1]+p1m*T11[i+1]+p0m*T10[i+1])
          logb2[i]<-m+log(p2m*T22[i+1]+p1m*T21[i+1]+p0m*T20[i+1])

  ##        print(c(i,logb0[i],logb1[i],logb2[i])
          
          m<-max(logb0[i],logb1[i],logb2[i])

         # print(c(i,m))
    #      if(i<4990)
    #        stop("")
	}

	# - termination 
	logp0 =  logb0[1]+logE[3,1]+log(k0)
	logp1 =  logb1[1]+logE[2,1]+log(k1)
	logp2 =  logb2[1]+logE[1,1]+log(k2)
	m = max(logp0, logp1, logp2)
	loglike = m+log(exp(logp0-m)+exp(logp1-m)+exp(logp2-m))
	# NB: the -m+m is a trick (add protect) to avoid underflow 
	cat("loglike in backward: ",loglike,"\n")
	# - return both the final like and the values of the backward variables	
	return(list(loglike=loglike,logb0=logb0,logb1=logb1,logb2=logb2))
}



######################################
# Decode functions
######################################

# ------------------------------------
# Function for viterbi decoding
# ------------------------------------

viterbidecode_relateHMM <- function(x){

	# prepare for decoding
	# ---------------------
    
	# - get distances and number of snps
  t<-x$t
  snp<-x$snp

	# - get parameter values
  a <- x$a
  k0 <- x$k[3]
  k1 <- x$k[2]
  k2 <- x$k[1]

	# - get emission probabilities
  logE<-x$logE
	
	# - get transition probabilities
  Tp = x$Tp
  logT00<-log(Tp$T00)
  logT01<-log(Tp$T01)
  logT02<-log(Tp$T02)
  logT10<-log(Tp$T10)
  logT11<-log(Tp$T11)
  logT12<-log(Tp$T12)
  logT20<-log(Tp$T20)
  logT21<-log(Tp$T21)
  logT22<-log(Tp$T22)
  
	# run viterbi
	# ---------------------
    
	# - make container for result
  path<-c()
	
	# - init: P(X0=x)P(g0|X0=x) for the 3 diff possible x vals
  v0<-log(k0)+logE[3,1] 
  v1<-log(k1)+logE[2,1]
  v2<-log(k2)+logE[1,1]
  
	# - recursion
  ptr0<--1
  ptr1<--1
  ptr2<--1
  for (i in (2:snp)){
    v0[i]<-logE[3,i]+max(c(logT00[i]+v0[i-1],logT10[i]+v1[i-1],logT20[i]+v2[i-1]))
    v1[i]<-logE[2,i]+max(c(logT01[i]+v0[i-1],logT11[i]+v1[i-1],logT21[i]+v2[i-1]))
    v2[i]<-logE[1,i]+max(c(logT02[i]+v0[i-1],logT12[i]+v1[i-1],logT22[i]+v2[i-1]))
    
    ptr0[i]<-which.max(c(v0[i-1]+logT00[i],v1[i-1]+logT10[i],v2[i-1]+logT20[i]))-1
    ptr1[i]<-which.max(c(v0[i-1]+logT01[i],v1[i-1]+logT11[i],v2[i-1]+logT21[i]))-1
    ptr2[i]<-which.max(c(v0[i-1]+logT02[i],v1[i-1]+logT12[i],v2[i-1]+logT22[i]))-1
  }
  

    # - termination
    l = max(c(v0[snp],v1[snp],v2[snp]))
  cat("like in viterbi= ",l,"\n")
    pi<-rep(0,snp)
    pi[snp]<-which.max(c(v0[snp],v1[snp],v2[snp]))-1

    # - traceback
    for(i in snp:2){
      if(pi[i]==2){
        pi[i-1]<-ptr2[i]
      }else{
        if(pi[i]==1){
          pi[i-1]<-ptr1[i]
        }else{
          pi[i-1]<-ptr0[i]
        }
      }
    }
  
    # - return most probable path
    return(pi)
}


viterbidecode_relateHMM_tsk<- function(x){

	# prepare for decoding
	# ---------------------
    
	# - get distances and number of snps
  t<-x$t
  snp<-x$snp

	# - get parameter values
  a <- x$a
  k0 <- x$k[3]
  k1 <- x$k[2]
  k2 <- x$k[1]

	# - get emission probabilities
  logE<-x$logE
	
	# - get transition probabilities
  Tp = x$Tp
  logT00<-log(Tp$T00)
  logT01<-log(Tp$T01)
  logT02<-log(Tp$T02)
  logT10<-log(Tp$T10)
  logT11<-log(Tp$T11)
  logT12<-log(Tp$T12)
  logT20<-log(Tp$T20)
  logT21<-log(Tp$T21)
  logT22<-log(Tp$T22)
  
	# run viterbi
	# ---------------------
    
	# - make container for result
  path<-c()
	
	# - init: P(X0=x)P(g0|X0=x) for the 3 diff possible x vals
  v0<-log(k0)+logE[3,1] 
  v1<-log(k1)+logE[2,1]
  v2<-log(k2)+logE[1,1]
  print(c(v0,v1,v2))
	# - recursion
  ptr0<--1
  ptr1<--1
  ptr2<--1

  m <- 0
  for (i in (2:snp)){
    v0[i]<-logE[3,i]+max(c(logT00[i]+v0[i-1],logT10[i]+v1[i-1],logT20[i]+v2[i-1]))
    v1[i]<-logE[2,i]+max(c(logT01[i]+v0[i-1],logT11[i]+v1[i-1],logT21[i]+v2[i-1]))
    v2[i]<-logE[1,i]+max(c(logT02[i]+v0[i-1],logT12[i]+v1[i-1],logT22[i]+v2[i-1]))
    ptr0[i]<-which.max(c(v0[i-1]+logT00[i],v1[i-1]+logT10[i],v2[i-1]+logT20[i]))-1
    ptr1[i]<-which.max(c(v0[i-1]+logT01[i],v1[i-1]+logT11[i],v2[i-1]+logT21[i]))-1
    ptr2[i]<-which.max(c(v0[i-1]+logT02[i],v1[i-1]+logT12[i],v2[i-1]+logT22[i]))-1
  }
  

    # - termination
    l = max(c(v0[snp],v1[snp],v2[snp]))
  cat("like in viterbi_tsk= ",l,"\n")
    pi<-rep(0,snp)
    pi[snp]<-which.max(c(v0[snp],v1[snp],v2[snp]))-1

    # - traceback
    for(i in snp:2){
      if(pi[i]==2){
        pi[i-1]<-ptr2[i]
      }else{
        if(pi[i]==1){
          pi[i-1]<-ptr1[i]
        }else{
          pi[i-1]<-ptr0[i]
        }
      }
    }
  
    # - return most probable path
    return(pi)
}

# ------------------------------------
# Function for posterior decoding
# ------------------------------------
backs <- c()
posteriordecode_relateHMM<-function(x){

	# prepare for decoding
	# ---------------------

    # - get distances
  t<-x$t
	
	# - get parameter values
  a <- x$a
  k0 <- x$k[3]
  k1 <- x$k[2]
  k2 <- x$k[1]

	# - get emissions 
  logE<-x$logE
	
	# - calc transition probabilities
  Tp = x$Tp
	
	# - run forward algorithm
  resf = forward_relateHMM(k0,k1,k2,logE,Tp)
  forwards<<-resf
	# - run backward algorithm
  resb = backward_relateHMM(k0,k1,k2,logE,Tp)
  backs<<-resb
	
	# - Calculation of posterior
  l = resf$loglike
  post<-cbind(exp(resf$logf2+resb$logb2-l),
              exp(resf$logf1+resb$logb1-l),
              exp(resf$logf0+resb$logb0-l))

	# - sanity check: does post add to 1?
  apply(post,1,sum)

    # return post
  return(post)	
}

# ----------------------------------------------------------
# Function for doing both viterbi and posterior decoding
# ----------------------------------------------------------

decode_relateHMM<-function(x){

  path<-viterbidecode_relateHMM_tsk(x)
  post<-posteriordecode_relateHMM(x)	
  obj=list(path=path,post=post)
  class(obj)="decode"

  return(obj)
}


###################################################################
# Main function (ML estimates parameters adnt hen does decoding)
###################################################################

`relateHMM`<-function(geno,pair=c(1,2),pos,af=NULL,par=NULL,
	minmaf=0,alim=c(0.01,0.5),opti=1,chr=NULL,method=NULL,
	calc.a=F,fix.a=NULL,fix.k2=NULL,back=0,LD="rsq2",prune=NULL,
	abs=ifelse(LD=="D",TRUE,FALSE),epsilon=0.01,start=NULL,phi=0.013){

	# Simple checks/reformatting of input data
	# - make sure there is exactly one pos for each locus with a genotype
	if(dim(geno)[2]!=length(pos))
		stop("number of SNP does not match the number of positions")

	# - calc af if necessary
	if(is.null(af))
		af<-apply(geno,2,mean,na.rm=TRUE)/2

	# - remove loci with maf<minmaf of missing in any of the two inds
	keep<-!is.na(geno[pair[1],])&!is.na(geno[pair[2],])&af>=minmaf&af<=1-minmaf # NB changed (added check of ind2)
	keep<-!is.na(geno[pair[1],])&af>=minmaf&af<=1-minmaf
	af<-af[keep]
	geno<-geno[,keep]
    	pos<-pos[keep]
	if(!is.null(chr))
    	chr<-chr[keep]

    	# - prune for ld if needed
    	if(!is.null(prune)){
		snp<-dim(geno)[2]
		ld<-ld.snp2(geno[,snp:1],depth=back)
		keep<-pruning(ld,snp,LD,back,prune,abs)
		af<-af[keep]
		geno<-geno[,keep]
		pos<-pos[keep]
		if(!is.null(chr))
			chr<-chr[keep]
		cat("pruning done\n")
	}

	# - extract extra info
	t<-diff(pos)
	t[t<0]<-1e20
    	t<-c(0,t,0)
    	nind<-dim(geno)[1]
    	snp<-dim(geno)[2]
	
    	ind1<-geno[pair[1],]
    	ind2<-geno[pair[2],]
    	af<-1-af

    	# - precalc emissions
	ep1<-error(ind1,epsilon,snp)
        ep2<-error(ind2,epsilon,snp)
        logE=log(emission_relateHMM_old(af,ep1,ep2)) # NB changed (to relatedness emissions)
	
	# - define internal function for calculating minusloglike			 
	minusloglike_relateHMM<-function(delta,logE,t,alim=c(0.01,100),fix.a=NULL,fix.k2=NULL,calc.a=FALSE,phi=0.013){

        # Getting parameters values:				 
		# - delta is adjusted according to special cases
          if(calc.a)
            fix.a<-alim[1]
          if(!is.null(fix.a))
            delta<-c(fix.a,delta)
          if(!is.null(fix.k2))
            delta<-c(delta[1],fix.k2,delta[-1])
				 
		# - parameter values are extracted from delta
          a<-delta[1]
          k2 <- delta[2] # NB changed (switched between k2 ond k0!!) 
          k1 <- delta[3]
          k0 <- 1-(k2+k1)
                
		# - parameters are checked 
                # -- are k vals valid for calculate.a? (if so like=1e20 returned)
                if(calc.a&&((k1+2*k2)^2<4*k2)){
                  cat("None valid k vector\n")
                  return(1e20)
		}
		 
		# -- is  any of the ks are<0 or>1, if a<mina or a>maxa? (if so like=1e20 returned)
		if(sum(c(k0,k2,k1)<0,a<alim[1])>0|sum(c(k0,k1,k2,a-alim[2]+1)>1)>0)
			return(1e20)
				 
		# - if calc.a a is calculated
		if(calc.a){
			a<-calculate.a(k0,k1,k2,phi=phi)
                        if(is.na(a)|a<alim[1]|a>alim[2])
                          return(1e20)
                }
                        
		# - transition probabilities are calculated
		Tp = transition_relateHMM(a,t,k0,k1,k2)
		
		# - likelihood is calculated
		loglike = fastlike_relateHMM(k0,k1,k2,logE,Tp)$loglike
		return(-loglike)
	} 
	
	# - get loglike for extremes 
#	loglike_parentoffspring<--minusloglike_relateHMM(delta=c(calculate.a(0,1,0),0,1),logE=logE,t=t,alim=alim) # NB changed (from 1,0 to 0,1) 
#	loglike_unrelated<--minusloglike_relateHMM(delta=c(calculate.a(0,0,1),0,0),logE=logE,t=t,alim=alim)
	loglike_parentoffspring<--minusloglike_relateHMM(delta=c(alim[1],0,1),logE=logE,t=t,alim=alim)
	loglike_unrelated<--minusloglike_relateHMM(delta=c(alim[1],0,0),logE=logE,t=t,alim=alim) 
        conv<-c()
   
    # - if no method is chosen BFGS is used
    if(is.null(method))
      method<-"BFGS"

    # - if parameter values are not specified as input they are estimated using ML 
    if(is.null(par)){
		# - make up a "good" starting value for the optim function
		cat("Start opti to find ML parameter values...\n")
		conv$value<-Inf
		for(tal in 1:opti){
			print(paste("Doing optim run",tal,"out of",opti))
			# - construct starting vals (triple consisting of (a,k2,k1)) 
			# -- first sampling 8 values "close to" 0.1
			temp.start<-diff(c(0,sort(runif(8)),1))[1:8]
			# -- then transform them into two numbers 1-four of them,sum of the last two
			temp.start<-c(1-sum(temp.start[c(1,5)])-sum(temp.start[c(3,6)]),
						  sum(temp.start[c(3,6)])) 
			# -- add a 3rd random number between min alim and max alim
			temp.start <- c(runif(1,min=alim[1],max=alim[2]),temp.start)
			if(tal==1&&!is.null(start)){
                          temp.start=start[-4]
                        }
			# - adjust format of delta to special cases
			# -- if fix.k2 then the k2 start val is rm
			if(!is.null(fix.k2))
				temp.start<-temp.start[-2]  
			# -- if fix.a or calc.a then the "a" start val is rm
			if(!is.null(fix.a)|calc.a)
				temp.start<-temp.start[-1]

			# - find ML with temp.start as init value
			conv_temp<-optim(temp.start,minusloglike_relateHMM,logE=logE,t=t,
							 alim=alim,fix.a=fix.a,fix.k2=fix.k2,
							 method=method,calc.a=calc.a,phi=phi)
			print(paste("Current opti is",conv$value,". Suggested opti is",conv_temp$value),sep="")
			print(paste("The corresponding parameters are",conv_temp$par[1],conv_temp$par[2],
						conv_temp$par[3],"which gives a likelihood of",conv_temp$value))
			if(conv_temp$value<conv$value)
				conv<-conv_temp
                      }

		par<-conv$par
                print("ML pars")
                print(sapply(par,function(x){round(x,2)}))
                # once the ML values are found extra values for special cases are inserted 
		if(calc.a)
			fix.a<-1  
		if(!is.null(fix.a))
			par<-c(fix.a,par)
		if(!is.null(fix.k2))
			par<-c(par[1],fix.k2,par[-1])
		if(calc.a) # NB changed (I added this function! from relateAnette on brutus)
			par<-c(calculate.a(1-sum(par[2:3]),par[3],par[2]),par[-1])  
      
		a <- par[1]
		k<-c(par[2:3],1-sum(par[2:3]))
		k.loglike<--conv$value
              }else{
                print("User given pars:")
                print(par)
		a<-par[1]
		k<-par[2:4]
		k.loglike<--minusloglike_relateHMM(c(a,par[2:3]),logE,t=t,alim=alim)
		conv<-NULL
	}

	Tp = transition_relateHMM(a,t,k[3],k[2],k[1])
	relatedness<-k[1]+k[2]/2 

	# - Perform decode
    	obj<-list(k=k,k.loglike=k.loglike,relatedness=relatedness,a=a,u.loglike=loglike_unrelated,
			  po.loglike=loglike_parentoffspring,conv=conv,LD=LD,t=t,snp=snp,position=pos,
			  logE=logE,Tp=Tp,alim=alim,af=af,chr=chr,conv=conv,geno=ind1)
	class(obj)="relateHMM"
	obj$decode<-decode_relateHMM(obj)
        names(obj$k)<-c("IBD=2","IBD=1","IBD=0")
	return(obj)
}





`ngsrelateHMM`<-function(geno,pair=c(1,2),pos,af=NULL,par=NULL,
                         minmaf=0,alim=c(0.01,0.5),opti=1,chr=NULL,method=NULL,
                         calc.a=F,fix.a=NULL,fix.k2=NULL,back=0,LD="rsq2",prune=NULL,
                         abs=ifelse(LD=="D",TRUE,FALSE),epsilon=0.01,start=NULL,phi=0.013){

	# Simple checks/reformatting of input data
	# - make sure there is exactly one pos for each locus with a genotype
	if(dim(geno)[2]!=length(pos))
		stop("number of SNP does not match the number of positions")

	# - calc af if necessary
      if(FALSE){
        if(is.null(af))
		af<-apply(geno,2,mean,na.rm=TRUE)/2
      }else{
        if(is.null(af)){
          print("Can only be run with af not equal to NULL!")
          return(NA)
        }
      }
        
	# - remove loci with maf<minmaf of missing in any of the two inds
      if(FALSE){
        keep<-!is.na(geno[pair[1],])&!is.na(geno[pair[2],])&af>=minmaf&af<=1-minmaf # NB changed (added check of ind2)
	keep<-!is.na(geno[pair[1],])&af>=minmaf&af<=1-minmaf
      }else{
        keep<-af>=minmaf&af<=1-minmaf
      }
        cat("we remove: ",sum(keep==0),"\n");
        af<-af[keep]
	geno<-geno[,keep]
    	pos<-pos[keep]
	if(!is.null(chr))
    	chr<-chr[keep]

        if(FALSE){
    	# - prune for ld if needed
    	if(!is.null(prune)){
		snp<-dim(geno)[2]
		ld<-ld.snp2(geno[,snp:1],depth=back)
		keep<-pruning(ld,snp,LD,back,prune,abs)
		af<-af[keep]
		geno<-geno[,keep]
		pos<-pos[keep]
		if(!is.null(chr))
			chr<-chr[keep]
		cat("pruning done\n")
	}
      }
	# - extract extra info
	t<-diff(pos)
	t[t<0]<-1e20
    	t<-c(0,t,0)
    	snp<-dim(geno)[2]

      if(FALSE){
    	# - precalc emissions
    	ind1<-geno[pair[1],]
    	ind2<-geno[pair[2],]
	ep1<-error(ind1,epsilon,snp)
        ep2<-error(ind2,epsilon,snp)
      }else{
        # - extract genolikes
        ep1<-t(geno[((pair[1]-1)*3+1):((pair[1]-1)*3+3),])
        ep2<-t(geno[((pair[2]-1)*3+1):((pair[2]-1)*3+3),])
      }
      af<-1-af
      logE=log(emission_relateHMM_new(af,ep1,ep2)) # NB changed (to relatedness emissions)
        eem<<-(emission_relateHMM_new(af,ep1,ep2))
	# - define internal function for calculating minusloglike			 
               minusloglike_relateHMM<-function(delta,logE,t,alim=c(0.01,100),fix.a=NULL,fix.k2=NULL,calc.a=FALSE,phi=0.013){
	
        # Getting parameters values:				 
		# - delta is adjusted according to special cases
		if(calc.a)
			fix.a<-alim[1]
		if(!is.null(fix.a))
			delta<-c(fix.a,delta)
		if(!is.null(fix.k2))
			delta<-c(delta[1],fix.k2,delta[-1])
				 
		# - parameter values are extracted from delta
		a<-delta[1]
		k2 <- delta[2] # NB changed (switched between k2 ond k0!!) 
		k1 <- delta[3]
		k0 <- 1-(k2+k1)
                
		# - parameters are checked 
                # -- are k vals valid for calculate.a? (if so like=1e20 returned)
                if(calc.a&&(k1+2*k2)^2<4*k2){
                  cat("None valid k vector\n")
                  return(1e20)
		}
		
		# -- is  any of the ks are<0 or>1, if a<mina or a>maxa? (if so like=1e20 returned)
		if(sum(c(k0,k2,k1)<0,a<alim[1])>0|sum(c(k0,k1,k2,a-alim[2]+1)>1)>0)
			return(1e20)
				 
		# - if calc.a a is calculated
		if(calc.a){
			a<-calculate.a(k0,k1,k2,phi=phi)
                        if(is.na(a)|a<alim[1]|a>alim[2])
                          return(1e20)
                }
                        
		# - transition probabilities are calculated
		Tp = transition_relateHMM(a,t,k0,k1,k2)
		
		# - likelihood is calculated
		loglike = fastlike_relateHMM(k0,k1,k2,logE,Tp)$loglike
		return(-loglike)
	} 
	
	# - get loglike for extremes 
#	loglike_parentoffspring<--minusloglike_relateHMM(delta=c(calculate.a(0,1,0),0,1),logE=logE,t=t,alim=alim) # NB changed (from 1,0 to 0,1) 
#	loglike_unrelated<--minusloglike_relateHMM(delta=c(calculate.a(0,0,1),0,0),logE=logE,t=t,alim=alim)
	loglike_parentoffspring<--minusloglike_relateHMM(delta=c(alim[1],0,1),logE=logE,t=t,alim=alim) # NB changed (from 1,0 to 0,1) 
	loglike_unrelated<--minusloglike_relateHMM(delta=c(alim[1],0,0),logE=logE,t=t,alim=alim) # same as Relate if c(xx,1,0) -- Error in Relate??
        conv<-c()
   
    # - if no method is chosen BFGS is used
    if(is.null(method))
      method<-"BFGS"

    # - if parameter values are not specified as input they are estimated using ML 
    if(is.null(par)){
		# - make up a "good" starting value for the optim function
      cat("Start opti to find ML parameter values...\n")
      conv$value<-Inf
      for(tal in 1:opti){
        print(paste("Doing optim run",tal,"out of",opti))
        ## - construct starting vals (triple consisting of (a,k2,k1)) 
        ## -- first sampling 8 values "close to" 0.1
        temp.start<-diff(c(0,sort(runif(8)),1))[1:8]
        ## -- then transform them into two numbers 1-four of them,sum of the last two
        temp.start<-c(1-sum(temp.start[c(1,5)])-sum(temp.start[c(3,6)]), sum(temp.start[c(3,6)])) 
        ## -- add a 3rd random number between min alim and max alim
        temp.start <- c(runif(1,min=alim[1],max=alim[2]),temp.start)
        if(tal==1&&!is.null(start))
          temp.start=start[-4]
        ##	temp.start=start

                                        ## - adjust format of delta to special cases
        ## -- if fix.k2 then the k2 start val is rm
        if(!is.null(fix.k2))
          temp.start<-temp.start[-2]  
			## -- if fix.a or calc.a then the "a" start val is rm
        if(!is.null(fix.a)|calc.a)
          temp.start<-temp.start[-1]
        
			## - find ML with temp.start as init value
        conv_temp<-optim(temp.start,minusloglike_relateHMM,logE=logE,t=t,
							 alim=alim,fix.a=fix.a,fix.k2=fix.k2,
							 method=method,calc.a=calc.a,phi=phi)
			print(paste("Current opti is",conv$value,". Suggested opti is",conv_temp$value),sep="")
			print(paste("The corresponding parameters are",conv_temp$par[1],conv_temp$par[2],conv_temp$par[3]))
			if(conv_temp$value<conv$value)
				conv<-conv_temp
                      }

		par<-conv$par
                print("ML pars")
                print(sapply(par,function(x){round(x,2)}))
                # once the ML values are found extra values for special cases are inserted 
		if(calc.a)
			fix.a<-1  
		if(!is.null(fix.a))
			par<-c(fix.a,par)
		if(!is.null(fix.k2))
			par<-c(par[1],fix.k2,par[-1])
		if(calc.a) # NB changed (I added this function! from relateAnette on brutus)
			par<-c(calculate.a(1-sum(par[2:3]),par[3],par[2]),par[-1])  
      
		a <- par[1]
		k<-c(par[2:3],1-sum(par[2:3]))
		k.loglike<--conv$value
              }else{
                print("User given pars:")
                print(sapply(par,function(x){round(x,2)}))
		a<-par[1]
		k<-par[2:4]
		k.loglike<--minusloglike_relateHMM(c(a,par[2:3]),logE,t=t,alim=alim)
		conv<-NULL
              }

	Tp = transition_relateHMM(a,t,k[3],k[2],k[1])
	relatedness<-k[1]+k[2]/2 

	# - Perform decode
    	obj<-list(k=k,k.loglike=k.loglike,relatedness=relatedness,a=a,u.loglike=loglike_unrelated,
			  po.loglike=loglike_parentoffspring,conv=conv,LD=LD,t=t,snp=snp,position=pos,
			  logE=logE,Tp=Tp,alim=alim,af=af,chr=chr,conv=conv)
	class(obj)="relateHMM"
	obj$decode<-decode_relateHMM(obj)
        names(obj$k)<-c("IBD=2","IBD=1","IBD=0")
	return(obj)
}



