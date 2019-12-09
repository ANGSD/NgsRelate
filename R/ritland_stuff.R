library('data.table')

toGeno <- function(x,data){
    colSums(matrix(as.numeric(data[x,]),2))
}

d<-fread("plink.tped",data.table=F)

d<-as.matrix(d[,-c(1:4)])
nind <- ncol(d)/2

g<-t(sapply(1:nrow(d),toGeno,data=d))
##d,g is now sites x inds

freq1 <- rowSums(d)/nind/2
if(FALSE){
    ##this blocks checks that the freqs is the same for us as in plink
   
    freq2 <- rowMeans(g)/2
    freq3<-read.table("plink.frq",he=T)$MAF
    plot(apply(cbind(freq1,1-freq1),1,min),freq3)
}
##freq is frequency of ones (1) not zeroes(1)
i1<-g[,1]
i2<-g[,2]

freq0=1-freq1

myfun <- function(i1,i2,freq0,freq1){
    if(length(table(sapply(list(i1,i2,freq0,freq1),length)))!=1)
        stop("Has to be same length")

    s1<-rep(NA,length(i1))
    s2<-s1
    ##wang 2010 doi:10.1534/genetics.110.124438
    ##check for no sharing
    s1[i1==0 & i2==2] <- 0
    s1[i1==2 & i2==0] <- 0
    s2[i1==0 & i2==2] <- 0
    s2[i1==2 & i2==0] <- 0
    ##check for double sharing
    s1[i1==0 & i2==0] <- 1
    s2[i1==0 & i2==0] <- 0
    s1[i1==2 & i2==2] <- 0
    s2[i1==2 & i2==2] <- 1
    ##check for sharing exactly 1 allele
    s1[i1==1 & i2==1] <- 0.25
    s2[i1==1 & i2==1] <- 0.25
    ##check for sharing one (het x hom)
    s1[i1==1 & i2==2] <- 0     ##s1[i1==1 & i2==2] <- 0.5
    s2[i1==1 & i2==2] <- 0.5
    s1[i1==1 & i2==0] <- 0.5
    s2[i1==1 & i2==0] <- 0     ##s2[i1==1 & i2==0] <- 0.5
    ##check for sharing one (hom x het)
    s1[i1==0 & i2==1] <- 0.5
    s2[i1==0 & i2==1] <- 0      ##s2[i1==0 & i2==1] <- 0.5 
    s1[i1==2 & i2==1] <- 0	##s1[i1==2 & i2==1] <- 0.5
    s2[i1==2 & i2==1] <- 0.5

    rhats <- (s1/freq0+s2/freq1-1) * 2

    ##rhats are site specific, we obtain the multilocus estimator b weighting the site-specific by the inverse of their sampling variance
    ##wang 2017 heredity 2017 119,302-313
    ## sum(rhat*(nalleles_i-1))/sum(nalleles_i-1)
    ## for diallelic, this is the arithmic mean
    return(mean(rhats))
}

res<-c()
for(i in 1:(ncol(g)-1))
    for(j in (i+1):ncol(g))
        res<-rbind(res,c(i,j,myfun(g[,i],g[,j],freq0,freq1)))

write.table(res, file="res.txt", quote=F, row.names=F, col.names=F)


