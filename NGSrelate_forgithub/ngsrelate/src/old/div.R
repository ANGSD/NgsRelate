# testing emission prob functions
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




  perfectlike1<-getPerfectLikes(as.vector(geno[1:snp,]))
perfectlike2<-getPerfectLikes(as.vector(geno[1:snp+snp,]))

like1<-getLikes(as.vector(geno[1:snp,]),dep=2,norm=T)
like2<-getLikes(as.vector(geno[1:snp+snp,]),dep=10,norm=T)

ML1<-doMLgenoCalling(like1)
ML2<-doMLgenoCalling(like1)
data_likebased<-rbind(ML1,ML2,data)
t2<-runHmmld(data_likebased,pair=c(1,2),pos=pos/1e6,chr=rep(1:10,each=1000),back=0,ld_adj=F)#estimate tracts of relatedness based on called genos
par(mfrow=2:1)
plot(t)
plot(t2)


## measures of comparison
s = snp*10

## Using post:
## - fraction of correctly inferred sites (using post>0.95) -- 1 is best
for(postthres in c(0.75,0.9,0.95,0.99)){
cor = sum(unlist(lapply(r,function(x){x$state}))==apply(t$post ,2,function(x){if(max(x)>postthres){3-which.max(x)}else{NA}}),na.rm=T)/s

cor2 = sum(unlist(lapply(r,function(x){x$state}))==apply(t2$post,2,function(x){if(max(x)>postthres){3-which.max(x)}else{NA}}),na.rm=T)/s


## - fraction of incorrectly inferred sites (using post>0.95) -- 0 is best
incor = sum(unlist(lapply(r,function(x){x$state}))!=apply(t$post ,2,function(x){if(max(x)>postthres){3-which.max(x)}else{NA}}),na.rm=T)/s

incor2 = sum(unlist(lapply(r,function(x){x$state}))!=apply(t2$post,2,function(x){if(max(x)>postthres){3-which.max(x)}else{NA}}),na.rm=T)/s

## - fraction of NAs (using post>0.95) -- 0 is best
nas = sum(is.na(apply(t$post ,2,function(x){if(max(x)>postthres){3-which.max(x)}else{NA}})))/s

nas2 = sum(is.na(apply(t2$post,2,function(x){if(max(x)>postthres){3-which.max(x)}else{NA}})))/s

print("")
print(paste("Postthres",postthres))
print(c(cor,incor,nas))
print(c(cor2,incor2,nas2))
}

## Using viterbi
cor = sum(unlist(lapply(r,function(x){x$state}))==t$path)/s
cor2 = sum(unlist(lapply(r,function(x){x$state}))==t2$path)/s


## - fraction of incorrectly inferred sites (using post>0.95) -- 0 is best
incor = sum(unlist(lapply(r,function(x){x$state}))!=t$path)/s
incor2 = sum(unlist(lapply(r,function(x){x$state}))!=t2$path)/s

print("")
print("Viterbi")
print(c(cor,incor,0))
print(c(cor2,incor2,0))

