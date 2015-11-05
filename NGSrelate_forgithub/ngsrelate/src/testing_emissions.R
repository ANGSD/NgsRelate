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

error<-function(obsgeno,epsilon,snp){
	ep<-matrix(NA,ncol=3,nrow=snp)
	
	ep[obsgeno==0,1]<-(1-epsilon)^2
	ep[obsgeno==1,1]<-2*(1-epsilon)*epsilon
	ep[obsgeno==2,1]<-(epsilon)^2
	
	ep[obsgeno==0,2]<-(1-epsilon)*epsilon
	ep[obsgeno==1,2]<-(1-epsilon)^2+epsilon^2
	ep[obsgeno==2,2]<-(1-epsilon)*epsilon
	
	ep[obsgeno==0,3]<-(epsilon)^2
	ep[obsgeno==1,3]<-2*(1-epsilon)*epsilon
	ep[obsgeno==2,3]<-(1-epsilon)^2
	
	return(ep)
}





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


epsilon= 0
geno1=rep(c(0,0,0,1,1,1,2,2,2),3)
geno2=rep(c(0,1,2),9)
ep1 = error(geno1,epsilon,81)
ep2 = error(geno2,epsilon,81)
freq = rep(0.2,81)
e_old = emission_relateHMM_old(freq,ep1,ep2)
apply(e_old[,1:9],1,sum)
apply(e_old[,10:18],1,sum)
apply(e_old[,19:27],1,sum)
apply(e_old[,28:36],1,sum)
apply(e_old[,37:45],1,sum)
apply(e_old[,46:54],1,sum)
apply(e_old[,55:63],1,sum)
apply(e_old[,64:72],1,sum)
apply(e_old[,73:81],1,sum)

e_new = emission_relateHMM_new(freq,ep1,ep2)
apply(e_new[,1:9],1,sum)
apply(e_new[,10:18],1,sum)
apply(e_new[,19:27],1,sum)
apply(e_new[,28:36],1,sum)
apply(e_new[,37:45],1,sum)
apply(e_new[,46:54],1,sum)
apply(e_new[,55:63],1,sum)
apply(e_new[,64:72],1,sum)
apply(e_new[,73:81],1,sum)
