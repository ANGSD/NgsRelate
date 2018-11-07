library(parallel)##for mclapply,

##function to normalize e.g sum to one.
norm <- function(x) x/sum(x)

##collapse the 'allele' vectors for one indivdual.
rs2 <- function(x)
    x[1,1]$gam+x[1,2]$gam


#get the true jacquard from haplotypes.
getjac <- function(a,b){
    i1h1 <- a[2,1]$h
    i1h2 <- a[2,2]$h
    i2h1 <- b[2,1]$h
    i2h2 <- b[2,2]$h
    j <- rep(NA,9)
    ##           -                 _            |,         ,|           \            /
    j[1] <-      sum(i1h1==i1h2 & i2h1==i2h2 & i1h1==i2h1 & i1h2==i2h2 & i1h1==i2h2 & i1h2==i2h1 )
    j[2] <-      sum(i1h1==i1h2 & i2h1==i2h2 & i1h1!=i2h1 & i1h2!=i2h2 & i1h1!=i2h2 & i1h2!=i2h1 )
    j[3] <-      sum(i1h1==i1h2 & i2h1!=i2h2 & i1h1==i2h1 & i1h2!=i2h2 & i1h1!=i2h2 & i1h2==i2h1 )
    j[3] <- j[3]+sum(i1h1==i1h2 & i2h1!=i2h2 & i1h1!=i2h1 & i1h2==i2h2 & i1h1==i2h2 & i1h2!=i2h1 )
    j[4] <-      sum(i1h1==i1h2 & i2h1!=i2h2 & i1h1!=i2h1 & i1h2!=i2h2 & i1h1!=i2h2 & i1h2!=i2h1 )
    j[5] <-      sum(i1h1!=i1h2 & i2h1==i2h2 & i1h1==i2h1 & i1h2!=i2h2 & i1h1==i2h2 & i1h2!=i2h1 )
    j[5] <- j[5]+sum(i1h1!=i1h2 & i2h1==i2h2 & i1h1!=i2h1 & i1h2==i2h2 & i1h1!=i2h2 & i1h2==i2h1 )
    j[6] <-      sum(i1h1!=i1h2 & i2h1==i2h2 & i1h1!=i2h1 & i1h2!=i2h2 & i1h1!=i2h2 & i1h2!=i2h1 )
    j[7] <-      sum(i1h1!=i1h2 & i2h1!=i2h2 & i1h1==i2h1 & i1h2==i2h2 & i1h1!=i2h2 & i1h2!=i2h1 )
    j[7] <- j[7]+sum(i1h1!=i1h2 & i2h1!=i2h2 & i1h1!=i2h1 & i1h2!=i2h2 & i1h1==i2h2 & i1h2==i2h1 )
    j[8] <-      sum(i1h1!=i1h2 & i2h1!=i2h2 & i1h1==i2h1 & i1h2!=i2h2 & i1h1!=i2h2 & i1h2!=i2h1 )
    j[8] <- j[8]+sum(i1h1!=i1h2 & i2h1!=i2h2 & i1h1!=i2h1 & i1h2==i2h2 & i1h1!=i2h2 & i1h2!=i2h1 )
    j[8] <- j[8]+sum(i1h1!=i1h2 & i2h1!=i2h2 & i1h1!=i2h1 & i1h2!=i2h2 & i1h1==i2h2 & i1h2!=i2h1 )
    j[8] <- j[8]+sum(i1h1!=i1h2 & i2h1!=i2h2 & i1h1!=i2h1 & i1h2!=i2h2 & i1h1!=i2h2 & i1h2==i2h1 )
    j[9] <-      sum(i1h1!=i1h2 & i2h1!=i2h2 & i1h1!=i2h1 & i1h2!=i2h2 & i1h1!=i2h2 & i1h2!=i2h1 )
    return(j)
}

#get the inbreeding for both indivudals based on the condensed jacquard coefficients
whatisf <- function(x)
    return(c(sum(x[1]+x[2],x[3],x[4]),sum(x[1],x[2],x[5],x[6])))

##get the relatedness coefficients,based on the condencsed jadquards coefficients
whatisrxy <- function(x)
    return(x[1]+x[7]+3/4*(x[3]+x[5])+x[8]*0.5)

##generate a gamet by sampling either of the two alleles.
makegamet2 <- function(a,b){
    if(length(a)!=2||length(b)!=2)
        stop("sampling of two haplotypes require two haplotypes")
    if(any(c(names(a)!=names(b),names(a)!=c("gam","h"))))
        stop("Dance straight")

    gam <- rep(NA,length(a$gam))
    h <- gam
    ss <- sample.int(2,length(a$gam),rep=T)
    for(s in 1:length(ss)){
        gam[s] = ifelse(ss[s]==1,a$gam[s],b$gam[s])
        h[s] = ifelse(ss[s]==1,a$h[s],b$h[s])
    }
    return(list(gam=gam,h=h))
}


## calculate emissions assumption here is that i1, and i2 each are matrices with 2column each and all entries are either 0/1 such that the genotypes for either individual is 0,1,2. NOT USING GLS

emis9 <-  function(i1,i2,freq){
    g1 <- rowSums(i1)
    g2 <- rowSums(i2)
    freq1 <- freq
    freq0 <- 1-freq1

    emis <- matrix(NA,ncol=length(freq),nrow=9)

    ##00&00 S1
    keep <- g1+g2==0
    tmp0 <- freq0[keep]
    emis[,keep] <- rbind(tmp0,tmp0^2,tmp0^2, tmp0^3,tmp0^2,tmp0^3, tmp0^2,tmp0^3,tmp0^4)
    ##11&11 S1
    keep <- g1+g2==4
    tmp1 <- freq1[keep]
    emis[,keep] <- rbind(tmp1,tmp1^2,tmp1^2, tmp1^3,tmp1^2,tmp1^3, tmp1^2,tmp1^3,tmp1^4)

    ##11&00 S2 ##difference between table1 in BG milligan Maximum-likelihood estimation of realtedness
    ##genetics 2002, and anderson,weier a maximums likelihood method fr the estimation of pairwise relatedness in structured popuations, genetics 207
    keep <- g1==2&g2==0
    tmp0 <- freq0[keep]
    tmp1 <- freq1[keep]
    emis[,keep] <- rbind(0,tmp0*tmp1,0,    tmp0^2*tmp1,0,tmp1^2*tmp0,  0,0,tmp1^2*tmp0^2)
    ##00&11 S2
    keep <- g1==0&g2==2
    tmp0 <- freq0[keep]
    tmp1 <- freq1[keep]
    emis[,keep] <- rbind(0,tmp0*tmp1,0,    tmp0*tmp1^2,0,tmp0^2*tmp1  ,0,0,tmp1^2*tmp0^2)

    ##11&01 S3
    keep <- g1==2&g2==1
    tmp0 <- freq0[keep]
    tmp1 <- freq1[keep]
    emis[,keep] <- rbind(0,0,tmp0*tmp1, 2*tmp1^2*tmp0,0,0, 0,tmp1^2*tmp0,2*tmp1^3*tmp0)

    ##00&01 S3
    keep <- g1==0&g2==1
    tmp0 <- freq0[keep]
    tmp1 <- freq1[keep]
    emis[,keep] <- rbind(0,0,tmp0*tmp1, 2*tmp0^2*tmp1,0,0, 0,tmp0^2*tmp1,2*tmp0^3*tmp1)


    ##01&11 S5
    keep <- g1==1&g2==2
    tmp0 <- freq0[keep]
    tmp1 <- freq1[keep]
    emis[,keep] <- rbind(0,0,0, 0,tmp0*tmp1,2*tmp1^2*tmp0, 0,tmp1^2*tmp0,2*tmp1^3*tmp0)

    ##01&00 S5
    keep <- g1==1&g2==0
    tmp0 <- freq0[keep]
    tmp1 <- freq1[keep]
    emis[,keep] <- rbind(0,0,0, 0,tmp0*tmp1,2*tmp0^2*tmp1, 0,tmp0^2*tmp1,2*tmp0^3*tmp1)

    ##01&01 S7
    keep <- g1==1&g2==1
    tmp0 <- freq0[keep]
    tmp1 <- freq1[keep]
    emis[,keep] <- rbind(0,0,0, 0,0,0, 2*tmp0*tmp1,tmp1*tmp0,4*tmp1^2*tmp0^2)


    emis
}

##do emstep
emStep <- function(x,emis){
    tmp <- x*emis
    tmp <- t(tmp)/colSums(tmp)
    nextx <- colMeans(tmp)
    if(any(nextx<0))
        stop("emStep gives negative")
    nextx
}
##calculate llh
llh <- function(x,emis){
    if(any(x<0|x>1))
        stop("pars outside of parsspace")
    -sum(log(colSums(x*emis)))
}

##normal em algorithm
em <- function(x,emis,niter=1500,tol=1e-6){
    lastp <- x
    lastllh <- llh(lastp,emis)
    iter <- 1
    while(iter <niter){
        newpar <- emStep(lastp,emis) ##this is the maximation of the EM
        newllh <- llh(newpar,emis)
        barplot(newpar,main=paste0("iter: ",iter," llh: ",round(newllh,4)," diff: ",round(newllh-lastllh,4) ),ylim=c(0,1))
        if(newllh>lastllh)
            stop("Problem new llh is worse than last llh")
        if(abs(lastllh-newllh)<tol){
            cat("breaking at ",iter,"\n")
            lastp <- newpar
            lastllh <- newllh
            break
        }
        lastp <- newpar
        lastllh <- newllh
        iter <- iter+1

    }
    return(list(par = lastp, value.objfn = lastllh, iter = iter, convergence = ifelse(iter>=niter,FALSE,TRUE)))
}

##acceleratedemalgorithm,
em2 <- function (x,emis,niter=1500,tol=1e-6)
{
    x <- norm(x)
    step.min <- 1
    step.max0 <- 1
    step.max <- 1
    mstep <- 4
    iter <- 1
    p <-x
    lold <- llh(p, emis)
    leval <- 1
    feval <- 0
    while (feval < niter) {
        p1 <- emStep(p, emis)
        feval <- feval + 1
        q1 <- p1 - p
        sr2 <- crossprod(q1)
        if (sqrt(sr2) < tol){
            p <- p1
            break
        }
        p2 <- emStep(p1, emis)
        feval <- feval + 1
        q2 <- p2 - p1
        sq2 <- sqrt(crossprod(q2))
        if (sq2 < tol){
            p <- p2
            break
        }
        sv2 <- crossprod(q2 - q1)
        srv <- crossprod(q1, q2 - q1)
        alpha <- sqrt(sr2/sv2)
        alpha <- max(step.min, min(step.max, alpha))
        p.new <- p + 2 * alpha * q1 + alpha^2 * (q2 - q1)
        if(any(p.new<0|p.new>1)){
            ##            cat("Problem p.new is out of parameterspace will reseat the interpolation step with the older emstep\n")
            ##           print(p2)
            p.new <- p2
        }
        if (abs(alpha - 1) > 0.01) {
            p.new <- emStep(p.new, emis)
            feval <- feval + 1
        }

        lnew <- llh(p.new, emis)
        leval <- leval + 1

        if (alpha == step.max)
            step.max <- mstep * step.max
        if (step.min < 0 & alpha == step.min)
            step.min <- mstep * step.min
        p <- p.new
        if (!is.nan(lnew))
            lold <- lnew
        iter <- iter + 1
    }

    lold <- llh(p, emis)
    leval <- leval + 1
    names(p) <- NULL
    return(list(par = p, value.objfn = lold, iter = iter, fpevals = feval,
        objfevals = leval, convergence =  ifelse(feval >= niter,FALSE,TRUE) ))
}


##7sep
##some gl simulation stuff
##input vector with entries from 0,1,2. Output are the gls inmatrix. each column are the gls for different sites for one indivitual
getLikes<-function(geno,dep=8,e=0.01,norm=FALSE){
  d<-dep
  n<-length(geno)
  dep<-rpois(n,d)
  nA<-rbinom(n,dep,c(e,0.5,1-e)[geno+1])
  res<-rbind(dbinom(nA,dep,e),dbinom(nA,dep,0.5),dbinom(nA,dep,1-e))
  if(norm)
    res<-t(t(res)/colSums(res))
  res
}

## same as above but likelihoods are now 0,0,1;0,1,0 or 1,0,0
getPerfectLikes<-function(genovec){
  res <- matrix(0,ncol=length(genovec),nrow=3)
  res[1,which(genovec==0)]=1
  res[2,which(genovec==1)]=1
  res[3,which(genovec==2)]=1
  res
}

##wrapper for above functions.
toGL <- function(x,e=0.01,...){
    if(e==0)
        getPerfectLikes(x)
    else
        getLikes(x,...)
}


#gl <- getLikes(dat$geno-1,dep=dep,e=err,norm=T,...)
#gls <- t(matrix(gl,nrow=nrow(dat$geno)*3))

##freq is the frequency of 1=one not 0=zero
##calculate emission matrix assuming input are gls.
##this calculates the 3 'non inbreed' ibds.
## The orginal ngsRelate, ASSUMNG GLS
emis3.gl <-  function(gl1,gl2,freq){
    ##gl[1,] is 00
    ##gl[2,] is 01
    ##gl[3,] is 11
    freq1 <- freq
    freq0 <- 1-freq1


    ##00&00 S1
    emis <- cbind(freq0^2,freq0^3,freq0^4)*gl1[1,]*gl2[1,]
    ##11&11 S1
    emis <- emis + cbind(freq1^2,freq1^3,freq1^4)*gl1[3,]*gl2[3,]
    ##11&00 S2
    emis <- emis + cbind(0,0,freq1^2*freq0^2)*gl1[3,]*gl2[1,]
    ##00&11 S2
    emis <- emis + cbind(0,0,freq1^2*freq0^2)*gl1[1,]*gl2[3,]
    ##11&01 S3
    emis <- emis + cbind(0,freq1^2*freq0,2*freq1^3*freq0)*gl1[3,]*gl2[2,]
    ##00&01 S3
    emis <- emis + cbind(0,freq0^2*freq1,2*freq0^3*freq1)*gl1[1,]*gl2[2,]
    ##01&11 S3
    emis <- emis + cbind(0,freq1^2*freq0,2*freq1^3*freq0)*gl1[2,]*gl2[3,]
    ##01&00
    emis <- emis + cbind(0,freq0^2*freq1,2*freq0^3*freq1)*gl1[2,]*gl2[1,]
    ##01&01
    emis <- emis + cbind(2*freq0*freq1,freq1*freq0,4*freq1^2*freq0^2)*gl1[2,]*gl2[2,]


    t(emis)
}

##calculate emissions matrix for single individual inbreeding
##ASSMUING GLS
emis2.gl <-  function(i1,freq){
    tmp0 <- 1-freq
    tmp1 <- 1-tmp0
    emis <- cbind(tmp0^2,tmp0)*i1[1,]
    emis <- emis + cbind(2*tmp0*tmp1,0)*i1[2,]
    emis <- emis + cbind(tmp1^2,tmp1)*i1[3,]
    t(emis)
}

##calculate emissions matrix for all 9 ibd pattern for the9 condensed jacquard coefficient assuming gls for our two samples
##ASSUMING GLS
emis9.gl <-  function(gl1,gl2,freq){
    freq1 <- freq
    freq0 <- 1-freq1



    ##00&00 S1
    emis<- cbind(freq0,freq0^2,freq0^2, freq0^3,freq0^2,freq0^3, freq0^2,freq0^3,freq0^4)*gl1[1,]*gl2[1,]
    ##11&11 S1
    emis <- emis + cbind(freq1,freq1^2,freq1^2, freq1^3,freq1^2,freq1^3, freq1^2,freq1^3,freq1^4)*gl1[3,]*gl2[3,]
    ##11&00
    emis <- emis + cbind(0,freq0*freq1,0,    freq0^2*freq1,0,freq1^2*freq0,  0,0,freq1^2*freq0^2)*gl1[3,]*gl2[1,]
    ##00&11 S2
    emis <- emis + cbind(0,freq0*freq1,0,    freq0*freq1^2,0,freq0^2*freq1  ,0,0,freq1^2*freq0^2)*gl1[1,]*gl2[3,]
    ##11&01 S3
    emis <- emis + cbind(0,0,freq0*freq1, 2*freq1^2*freq0,0,0, 0,freq1^2*freq0,2*freq1^3*freq0)*gl1[3,]*gl2[2,]
    ##00&01 S3
    emis <- emis + cbind(0,0,freq0*freq1, 2*freq0^2*freq1,0,0, 0,freq0^2*freq1,2*freq0^3*freq1)*gl1[1,]*gl2[2,]
    ##01&11 S5
    emis <- emis + cbind(0,0,0, 0,freq0*freq1,2*freq1^2*freq0, 0,freq1^2*freq0,2*freq1^3*freq0)*gl1[2,]*gl2[3,]
    ##01&00 S5
    emis <- emis + cbind(0,0,0, 0,freq0*freq1,2*freq0^2*freq1, 0,freq0^2*freq1,2*freq0^3*freq1)*gl1[2,]*gl2[1,]
    ##01&01 S7
    emis <- emis + cbind(0,0,0, 0,0,0, 2*freq0*freq1,freq1*freq0,4*freq1^2*freq0^2)*gl1[2,]*gl2[2,]
    t(emis)
}


get_likelihood <- function(emis, par){
    return(llh(par, emis))
}

tol <- 1e-10  ## tolerance
e <- 1e-3     ## error rate
iter <- 10000 ## iterations
ndip <- 10 ## number of diploid samples
nsites <- 1e4 ## 10k sites
maf=0.05

set.seed(0)
freq <- runif(nsites,min=maf,max=1-maf)
haps <- mclapply(freq,function(x) sample(c(0,1),2*ndip,rep=T,prob=c(1-x,x)))
haps <- matrix(unlist(haps),nsites,byrow=T);
haps <- mclapply(1:ncol(haps),function(x) list(gam=haps[,x],h=rep(x,nrow(haps))))
mm <- cbind(makegamet2(haps[[1]],haps[[2]]),makegamet2(haps[[3]],haps[[4]]))
mf <- cbind(makegamet2(haps[[5]],haps[[6]]),makegamet2(haps[[7]],haps[[8]]))
fm <- cbind(makegamet2(haps[[9]],haps[[10]]),makegamet2(haps[[11]],haps[[12]]))
ff <- cbind(makegamet2(haps[[13]],haps[[14]]),makegamet2(haps[[15]],haps[[16]]))
m <-  cbind(makegamet2(mm[,1],mm[,2]),makegamet2(mf[,1],mf[,2]))
m2 <- cbind(makegamet2(mm[,1],mm[,2]),makegamet2(mf[,1],mf[,2]))
f <-  cbind(makegamet2(ff[,1],ff[,2]),makegamet2(fm[,1],fm[,2]))
c1_inbreed <- cbind(makegamet2(m2[,1],m2[,2]),makegamet2(mf[,1], mf[,2]))
b1 <- cbind(makegamet2(f[,1],f[,2]),makegamet2(m[,1],m[,2]))
ind1 <- b1
ind2 <- c1_inbreed
jac <- norm(getjac(ind1,ind2))
ind1.gl <- toGL(rs2(ind1),e=e)
ind2.gl <- toGL(rs2(ind2),e=e)

glsPair <- rbind(ind1.gl, ind2.gl)
outname=as.character("test")
con <- file(paste0(outname,".glf"), 'wb')
writeBin(as.vector(log(glsPair)),con)
close(con)
system(paste0("gzip -f ", paste0(outname,".glf")))
write.table(freq,paste0(outname,".freq"),row.names=F,col.names=F,quote=F)
system(paste0("gzip -f ", paste0(outname,".freq")))

## print("INBREEDING: ind1")
## emis.inbred <- emis2.gl(ind1.gl,freq)
## print(em2(runif(2,min=0.01, max=0.99),emis.inbred,iter,tol=tol))
## print("INBREEDING: ind2")
## emis.inbred <- emis2.gl(ind2.gl,freq)
## print(em2(runif(2,min=0.01, max=0.99),emis.inbred,iter,tol=tol))


emis <- emis9.gl(ind1.gl,ind2.gl,freq)
d <- em2(runif(9,min=0.01, max=0.99),emis,iter,tol=tol)
est <- d$par

print(round(rbind(jac,est),3))
print(c("SSE: ", sum((jac-est)**2)))
print("Inbreeding (true, estimate): ")
print(round(rbind(whatisf(jac), whatisf(est)),3))
print("Relatedness (true, estimate): ")
print(round(rbind(whatisrxy(jac), whatisrxy(est)),3))
