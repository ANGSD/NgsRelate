
#not for ML and low epsilon but otherwise correct pars highr higher or equal llh
k0real=(table(dat$truepaths)/(dim(dat$truepaths)[2]*nchr))[1]
k1real=(table(dat$truepaths)/(dim(dat$truepaths)[2]*nchr))[2]
k2real=1-k0real-k1real

relateHMM(geno_likebased-1,pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,alim2),opti=opti,epsilon=0.01,calc.a=F,af=af,par=c(a,k2real,k1real,k0real)) 
ri2
ri20.15k2is0eps0.01

relateHMM(geno_likebased-1,pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,alim2),opti=opti,epsilon=0.05,calc.a=F,af=af,par=c(a,0,1,0)) 
ri20.15k2is0eps0.05

relateHMM(geno_likebased-1,pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,alim2),opti=opti,epsilon=0.1,calc.a=F,af=af,par=c(a,0,1,0)) 
ri20.15k2is0eps0.10

relateHMM(geno_likebased-1,pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,alim2),opti=opti,epsilon=0.15,calc.a=F,af=af,par=c(a,0,1,0)) 
ri20.15k2is0eps0.15

ngsrelateHMM(rbind(like1,like2),pair=c(1,2),pos=pos/1e6,chr=chr,alim=c(0.001,alim2),opti=opti,epsilon=0.01,calc.a=F,af=af,par=c(a,0,1,0))
ri2x
ri2x0.15k2is0

par(mfrow=3:1)
plot(ri2)
plot(ri2x)





for(t in c(ri2)){
for(postthres in c(0.5,0.75,0.9,0.95,0.99)){
  print(paste("Proportion of genome where post of sharing 1 allele is >", postthres," is ",mean(t$decode$post[,2]>postthres,na.rm=T), sep=""))
}
print(paste("Proportion of genome where viterbi path has state 1 is ",mean(t$decode$path==1,na.rm=T),sep=""))
print(paste("ML k1 is",t$k[2]))
print(paste("Proportion of genome where 1 allele is shared ibd is",k1real))  


s = nsnp*nchr
