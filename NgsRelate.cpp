/*
  ida@binf.ku.dk thorfinn@binf.ku.dk
  20 june, 2015

  thorfinn@binf.ku.dk
  ida@binf.ku.dk

  http://www.popgen.dk/software/
  https://github.com/ANGSD/NgsRelate/


  mod for inbreed anders
  g++ NgsRelateV3.cpp -O3 -lz -o ngsrelate3

*/

#include <vector>
#include <cstring>
#include <zlib.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <map>

#define LENS  4096
int refToInt[256] = {
  0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//15
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//31
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//47
  0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//63
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//79
  4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//95
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//111
  4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//127
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//143
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//159
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//175
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//191
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//207
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//223
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//239
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4//255
};
//this is as close to the bound we will allow
float emTole=1e-12;
double TINY=1e-12;
void normalize(double *tmp,int len){
  double s=0;
  for(int i=0;i<len;i++)
    s += tmp[i];
  for(int i=0;i<len;i++)
    tmp[i] /=s;
}


void print(FILE *fp,int x,int y,double **m){
  for(int i=0;i<x;i++){
    for(int j=0;j<y;j++)
      fprintf(fp,"%f ",m[i][j]);
    fprintf(fp,"\n");
  }

}

void print(FILE *fp,int x,double *m){
  for(int i=0;i<x;i++)
    fprintf(fp,"%f ",m[i]);
  fprintf(fp,"\n");
}

double loglike(double *p,double **emis,int len){
  //  fprintf(stderr,"%f %f %f\n",p[0],p[1],p[2]);
  double ret =0;
  
  for(int i=0;i<len;i++){
    double tmp = 0;
    for(int j=0;j<3;j++)
      tmp += p[j]*emis[i][j];
    ret +=log(tmp);
  }
  return ret;
}

double loglikeInbreed(double *p,double **emis,int len){
  //  fprintf(stderr,"%f %f %f\n",p[0],p[1],p[2]);
  double ret =0;
  
  for(int i=0;i<len;i++){
    double tmp = 0;
    for(int j=0;j<2;j++)
      tmp += p[j]*emis[i][j];
    ret +=log(tmp);
  }
  return ret;
}



void emStep1_inbreed(double *pre,double **emis,double *post,int len){
  // fprintf(stderr,"%f %f\n",pre[0],pre[1]);
  double inner[2];
  for(int x=0;x<2;x++)
    post[x] =0.0;
  
  for(int i=0;i<len;i++){
    for(int x=0;x<2;x++){
      inner[x] = pre[x]*emis[i][x];
      //fprintf(stderr,"%f %f\n",emis[i][0],emis[i][1]);
    }
  
    normalize(inner,2);
    for(int x=0;x<2;x++)
      post[x] += inner[x];
    //   fprintf(stderr,"%f %f %f\n",post[0],post[1],post[2]);
  }
  //set bounds
  normalize(post,2);
  //  fprintf(stderr,"%f %f %f\n",post[0],post[1],post[2]);
}


void emStep1(double *pre,double **emis,double *post,int len){
  //  fprintf(stderr,"%f %f %f\n",pre[0],pre[1],pre[2]);
  double inner[3];
  for(int x=0;x<3;x++)
    post[x] =0.0;
  
  for(int i=0;i<len;i++){
    for(int x=0;x<3;x++)
      inner[x] = pre[x]*emis[i][x];
  
   normalize(inner,3);
   for(int x=0;x<3;x++)
     post[x] += inner[x];
   //   fprintf(stderr,"%f %f %f\n",post[0],post[1],post[2]);
  }
  //set bounds
#if 0
  for(int i=0;i<3;i++){
    if(post[i]<emTole)
      post[i] = emTole;
    if(post[i]>(1-emTole))
      post[i] = emTole;
  }
#endif
  normalize(post,3);
  //  fprintf(stderr,"%f %f %f\n",post[0],post[1],post[2]);
}

void minus(double fst[3],double sec[3],double res[3]){
  for(int i=0;i<3;i++)
    res[i] = fst[i]-sec[i];
}
void minus2(double fst[2],double sec[2],double res[2]){
  for(int i=0;i<2;i++)
    res[i] = fst[i]-sec[i];
}


double sumSquare(double mat[3]){
  double tmp=0;
  for(size_t i=0;i<3;i++){
    //    fprintf(stderr,"%f \n",mat[i]);
    tmp += mat[i]*mat[i];
  }

  return tmp;
}

double sumSquare2(double mat[2]){
  double tmp=0;
  for(size_t i=0;i<2;i++){
    //    fprintf(stderr,"%f \n",mat[i]);
    tmp += mat[i]*mat[i];
  }

  return tmp;
}



int emAccel(double *F,double **emis,double *F_new,int len){
  //  fprintf(stderr,"calling emaccel \n");
   double ttol=0.0000001;

  //  fprintf(stderr,"tol:%f\n",tol);
  //maybe these should be usersettable?
  double stepMin =1;
  double stepMax0 = 1;
  static double stepMax=stepMax0;
  double mstep=4;
  double objfnInc=1;


  double F_em1[3];
  double F_diff1[3];
  double F_em2[3];
  double F_diff2[3];
  double F_diff3[3];
  double F_tmp[3];

  emStep1(F,emis,F_em1,len);
  minus(F_em1,F,F_diff1);
  double sr2 = sumSquare(F_diff1);
  
  if(sqrt(sr2)<ttol){
    //    fprintf(stderr,"sr2 break:%f\n",sr2);
    return 0;
    //break;
  }
  emStep1(F_em1,emis,F_em2,len);
  minus(F_em2,F_em1, F_diff2);

  double sq2 = sumSquare(F_diff2);
  if(sqrt(sq2)<ttol){
    //fprintf(stderr,"sq2\n");
    return 0;
    //    break;
  }


  minus(F_diff2,F_diff1, F_diff3);
  
  double sv2 = sumSquare(F_diff3);
  
  double alpha = sqrt(sr2/sv2);
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  for(size_t i=0;i<3;i++)
      F_new[i] = F[i]+2*alpha*F_diff1[i]+alpha*alpha*F_diff3[i];

  if (fabs(alpha - 1) > 0.01){
    emStep1(F_new,emis,F_tmp,len);
    for(int i=0;i<3;i++)
      std::swap(F_new[i],F_tmp[i]);
  }

  double lnew = 1;
  if ((alpha - stepMax) > -0.001) {
    stepMax = mstep*stepMax;
  }
  //  print(stderr,3,F_new);
  //  fprintf(stderr,"alpha %f stepMax %f\n",alpha,stepMax);

  // fprintf(stderr,"calling emaccel \n");
  return 1;

  
}


int emAccel_inbreed(double *F,double **emis,double *F_new,int len){
  //  fprintf(stderr,"calling emaccel \n");
  double ttol=0.0000001;

  //  fprintf(stderr,"tol:%f\n",tol);
  //maybe these should be usersettable?
  double stepMin =1;
  double stepMax0 = 1;
  static double stepMax=stepMax0;
  double mstep=4;
  double objfnInc=1;


  double F_em1[2];
  double F_diff1[2];
  double F_em2[2];
  double F_diff2[2];
  double F_diff3[2];
  double F_tmp[2];

  emStep1_inbreed(F,emis,F_em1,len);
  minus2(F_em1,F,F_diff1);
  double sr2 = sumSquare2(F_diff1);
  
  if(sqrt(sr2)<ttol){
    //    fprintf(stderr,"sr2 break:%f\n",sr2);
    return 0;
    //break;
  }
  emStep1_inbreed(F_em1,emis,F_em2,len);
  minus2(F_em2,F_em1, F_diff2);

  double sq2 = sumSquare2(F_diff2);
  if(sqrt(sq2)<ttol){
    //fprintf(stderr,"sq2\n");
    return 0;
    //    break;
  }


  minus2(F_diff2,F_diff1, F_diff3);
  
  double sv2 = sumSquare2(F_diff3);
  
  double alpha = sqrt(sr2/sv2);
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  for(size_t i=0;i<2;i++)
      F_new[i] = F[i]+2*alpha*F_diff1[i]+alpha*alpha*F_diff3[i];

  if (fabs(alpha - 1) > 0.01){
    emStep1_inbreed(F_new,emis,F_tmp,len);
    for(int i=0;i<2;i++)
      std::swap(F_new[i],F_tmp[i]);
  }

  double lnew = 1;
  if ((alpha - stepMax) > -0.001) {
    stepMax = mstep*stepMax;
  }
  //  print(stderr,3,F_new);
  //  fprintf(stderr,"alpha %f stepMax %f\n",alpha,stepMax);

  // fprintf(stderr,"calling emaccel \n");
  return 1;

  
}



int em1(double *sfs,double  **emis,double tole,int maxIter,int len,int verbose){
  double oldLik,lik;
  oldLik = loglike(sfs,emis,len);
  if(verbose)
    fprintf(stderr,"startlik=%f %f %f %f\n",oldLik,sfs[0],sfs[1],sfs[2]);
  fflush(stderr);

  double tmp[3];
  int it;
  for(it=0;it<maxIter;it++) {
    emStep1(sfs,emis,tmp,len);
    for(int i=0;i<3;i++)
      sfs[i]= tmp[i];
    lik = loglike(sfs,emis,len);

    if(verbose)
      fprintf(stderr,"[%d] lik=%f diff=%g\n",it,lik,fabs(lik-oldLik));

    if(fabs(lik-oldLik)<tole){
     
      oldLik=lik;
      break;
    }
    oldLik=lik;
  }
  return it;
}


int em_inbreed(double *pars,double  **emis,double tole,int maxIter,int len,int model,int verbose){
  double oldLik,lik;
  oldLik = loglikeInbreed(pars,emis,len);
  if(verbose){
    fprintf(stderr,"startlik=%f %f %f\n",oldLik,pars[0],pars[1]);
    fflush(stderr);
  }

  double tmp[2];
  int it;
  for(it=0;it<maxIter;it++) {
    if(model==0)
      emStep1_inbreed(pars,emis,tmp,len);
    else
      emAccel_inbreed(pars,emis,tmp,len);
    for(int i=0;i<2;i++)
      pars[i]= tmp[i];
    lik = loglikeInbreed(pars,emis,len);
    if(verbose)
      fprintf(stderr,"[%d] lik=%f diff=%g\n",it,lik,fabs(lik-oldLik));

    if(fabs(lik-oldLik)<tole){
     
      oldLik=lik;
      break;
    }
    oldLik=lik;
  }
  return it;
}



int em2(double *sfs,double  **emis,double tole,int maxIter,int len,int verbose){
  double oldLik,lik;
  oldLik = loglike(sfs,emis,len);
  if(verbose){
    fprintf(stderr,"startlik=%f\n",oldLik);
    fflush(stderr);
  }

  double tmp[3];
  int it;
  for(it=0;it<maxIter;it++) {
    emAccel(sfs,emis,tmp,len);

    for(int i=0;i<3;i++)
      sfs[i]= tmp[i];
    lik = loglike(sfs,emis,len);
    if(verbose)
      fprintf(stderr,"[%d] lik=%f diff=%e\n",it,lik,fabs(lik-oldLik));

    if(fabs(lik-oldLik)<tole){
      // fprintf(stderr,"breaking\n");
      oldLik=lik;
      break;
    }
    oldLik=lik;
  }
  
  return it;
}

int em3(double *sfs,double  **emis,double tole,int maxIter,int len,int verbose){
  //  exit(0);
  double oldLik,lik;
  oldLik = loglike(sfs,emis,len);
  if(verbose){
    fprintf(stderr,"em3startlik=%f\n",oldLik);
    fflush(stderr);
  }

  double tmp[3];
  int it;
  int speedy=1;
  for(it=0;it<maxIter;it++) {
    if(speedy)
      emAccel(sfs,emis,tmp,len);
    else
      emStep1(sfs,emis,tmp,len);
    lik = loglike(tmp,emis,len);

    fprintf(stderr,"[%d] lik=%f diff=%e\n",it,lik,fabs(lik-oldLik));
    if(std::isnan(lik)||lik<oldLik){
      fprintf(stderr,"Problem llh is now bigger or nan, will go back and use regular em\n");
      fprintf(stderr,"This is offending pars: %f %f %f\n",sfs[0],sfs[1],sfs[2]);
      speedy=0;
      continue;
    }
    for(int i=0;i<3;i++)
      sfs[i]= tmp[i];
  
    if(fabs(lik-oldLik)<tole){
      // fprintf(stderr,"breaking\n");
      oldLik=lik;
      break;
    }
    oldLik=lik;
    
  }
  
  return it;
}





void emission_ngsrelate(double *freq,double **l1,double **l2,double **emis,int len){

  for(int i=0;i<len;i++){
    double freqA=freq[i];
    double freqa=1-freqA;
      //1 E<-cbind(freqA^4,freqA^3,freqA^2)*l1[,1]*l2[,1];			      // G_real=(AA,AA)
      emis[i][0] = pow(freqA,4)*l1[i][0]*l2[i][0];
      emis[i][1] = pow(freqA,3)*l1[i][0]*l2[i][0];
      emis[i][2] = pow(freqA,2)*l1[i][0]*l2[i][0];
      //2   E<-E+cbind(2*freqA^3*freqa,1*freqA^2*freqa,0)*l1[,1]*l2[,2] ;		# G_real=(AA,Aa)
      emis[i][0] += 2*pow(freqA,3)*freqa*l1[i][0]*l2[i][1];
      emis[i][1] += pow(freqA,2)*freqa*l1[i][0]*l2[i][1];
      emis[i][2] += 0;
      //3   E<-E+cbind(1*freqA^2*freqa^2,0,0)*l1[,1]*l2[,3]				# G_real=(AA,aa)
      emis[i][0] += pow(freqA,2)*pow(freqa,2)*l1[i][0]*l2[i][2];
      emis[i][1] += 0;
      emis[i][2] += 0;
      //4   E<-E+cbind(2*freqA^3*freqa,1*freqA^2*freqa,0)*l1[,2]*l2[,1]			# G_real=(Aa,AA)
      emis[i][0] += 2*pow(freqA,3)*freqa*l1[i][1]*l2[i][0];
      emis[i][1] += pow(freqA,2)*freqa*l1[i][1]*l2[i][0];
      emis[i][2] += 0;
      //5     E<-E+cbind(4*freqA^2*freqa^2,freqA*freqa,2*freqA*freqa)*l1[,2]*l2[,2]	# G_real=(Aa,Aa)
      emis[i][0] += 4*pow(freqA,2)*pow(freqa,2)*l1[i][1]*l2[i][1];
      emis[i][1] += freqA*freqa*l1[i][1]*l2[i][1];
      emis[i][2] += 2*freqA*freqa *l1[i][1]*l2[i][1];
      //6  E<-E+cbind(2*freqA*freqa^3,1*freqA*freqa^2,0)*l1[,2]*l2[,3]			# G_real=(Aa,aa)      
      emis[i][0] += 2*freqA*pow(freqa,3)*l1[i][1]*l2[i][2];
      emis[i][1] += freqA*pow(freqa,2)*l1[i][1]*l2[i][2];
      emis[i][2] += 0;
      //7  E<-E+cbind(1*freqA^2*freqa^2,0,0)*l1[,3]*l2[,1]				# G_real=(aa,AA)
      emis[i][0] += pow(freqA,2)*pow(freqa,2)*l1[i][2]*l2[i][0];
      emis[i][1] += 0;
      emis[i][2] += 0;
      //8  E<-E+cbind(2*freqA*freqa^3,1*freqA*freqa^2,0)*l1[,3]*l2[,2]		        # G_real=(aa,Aa)
      emis[i][0] += 2*freqA*pow(freqa,3)*l1[i][2]*l2[i][1];
      emis[i][1] += freqA*pow(freqa,2)*l1[i][2]*l2[i][1];
      emis[i][2] += 0;
      //9  E<-E+cbind(freqa^4,freqa^3,freqa^2)*l1[,3]*l2[,3]				# G_real=(aa,aa)
      emis[i][0] += pow(freqa,4)*l1[i][2]*l2[i][2];
      emis[i][1] += pow(freqa,3)*l1[i][2]*l2[i][2];
      emis[i][2] += pow(freqa,2)*l1[i][2]*l2[i][2];

      
      if(i==295&&0)
	fprintf(stderr,"emis[%d]:freq:%f %f %f %f\n",i,freqa,emis[i][0],emis[i][1],emis[i][2]);
      //      exit(0);
    }

}











void emission_ngsInbreed(double *freq,double **l1,double **emis,int len){

    for(int i=0;i<len;i++){
    double freqA=freq[i];
    double freqa=1-freqA;
    // G_real=(AA)
    emis[i][0] = pow(freqA,2)*l1[i][0];
    emis[i][1] = freqA*l1[i][0];

    // G_real=(Aa)
    emis[i][0] += 2*freqa*freqA*l1[i][1];
    emis[i][1] += 0;

    // G_real=(aa)
    emis[i][0] += pow(freqa,2)*l1[i][2];
    emis[i][1] += freqa*l1[i][0];

    //if(i==295&&0)
    //    fprintf(stderr,"emis[%d]:freq:%f %f %f\n",i,freqa,emis[i][0],emis[i][1]);
      //      exit(0);
    }

}







 double **getGL(const char *fname,int sites, int nInd){
   gzFile gz=Z_NULL;
   if(((gz=gzopen(fname,"rb")))==Z_NULL){
     fprintf(stderr,"Problem opening file:%s\n",fname);
     exit(0);
   }
   

   double **ret=new double*[sites+10];
   
   int i=0;
   while(1){
     //   for(int i=0;i<sites;i++){
     ret[i] = new double[3*nInd];
     int nbit = gzread(gz,ret[i],sizeof(double)*nInd*3);
     if(nbit==0)
       break;
     if(sizeof(double)*nInd*3!=nbit){
       fprintf(stderr,"\t-> Problem reading full chunk\n");
       exit(0);
     }
     for(int g =0;g<3*nInd;g++)
       ret[i][g] = exp(ret[i][g]);
#if 0
     for(int g=0;g<nInd;g++){
       double ts = 0;
       for(int gg=0;gg<3;gg++)
	 ts += ret[i][g*3+gg];
       for(int gg=0;gg<3;gg++)
	 ret[i][g*3+gg] /= ts;

     }
#endif
     i++;
     if(i>sites){
       fprintf(stderr,"too many sites in glf file. Looks outof sync\n");
       exit(0);
     }
   }
   if(i!=sites){
     fprintf(stderr,"nsites: %d assumed but %d read\n",sites,i);
     exit(0);
   }
   gzclose(gz);
   return ret;
 }


 std::vector<double> getDouble(const char *fname) {
   gzFile gz=Z_NULL;
   if(((gz=gzopen(fname,"rb")))==Z_NULL){
     fprintf(stderr,"Problem opening file:%s\n",fname);
     exit(0);
   }
   char *buf = new char[10000];
   std::vector<double> ret;
   while(gzgets(gz,buf,10000)){
     ret.push_back(atof(strtok(buf,"\n\t\r ")));
     char *tok=NULL;
     while(((tok=strtok(NULL,"\n\t\r ")))){
       ret.push_back(1-atof(tok));

     }

   }
   fprintf(stderr,"frequency file:%s contain %lu number of sites\n",fname,ret.size());
   gzclose(gz);
   delete [] buf;
   return ret;
 }

void print_info(FILE *fp){
  fprintf(fp, "\n");
  fprintf(fp, "Usage main analyses: ./ngsrelate  [options] \n");
  fprintf(fp, "Options:\n");
  fprintf(fp, "   -f <filename>       Name of file with frequencies\n");
  fprintf(fp, "   -m <INTEGER>        model 0=normalEM 1=acceleratedEM\n");
  fprintf(fp, "   -i <UINTEGER>       Maximum number of EM iterations\n");
  fprintf(fp, "   -t <FLOAT>          Tolerance for breaking EM\n");
  fprintf(fp, "   -r <FLOAT>          Seed for rand\n");
  fprintf(fp, "   -g gfile            Name of genotypellh file\n");
  fprintf(fp, "   -c <INT>            Should call genotypes instead?\n");
  fprintf(fp, "   -s <INT>            Should you swich the freq with 1-freq?\n");
  fprintf(fp, "   -v <INT>            Verbose. print like per iteration\n");
  fprintf(fp, "   -F <INT>            Estimate inbreeding instead of relatedness\n");
  fprintf(fp, "   -e <INT>            Errorrates when calling genotypes?\n");
  fprintf(fp, "   -a <INT>            First individual used for analysis? (zero offset)\n");
  fprintf(fp, "   -b <INT>            Second individual used for analysis? (zero offset)\n");
  fprintf(fp, "   -n <INT>            Number of samples in glf.gz\n");
  fprintf(fp, "   -l <INT>            minMaf or 1-Maf filter\n");
  fprintf(fp, "\n");
  fprintf(fp,"Or\n ./ngsrelate extract_freq pos.glf.gz plink.bim plink.freq\n");
  exit(0);
}

void callgenotypesEps(double **gls,int len,double eps){
  fprintf(stderr,"\t-> Call genotypes: %f\n",eps);
  double g00 = (1-eps)*(1-eps);
  double g01 = 2*(1-eps)*eps;
  double g02 = eps*eps;
  double g10 = (1-eps)*eps;
  double g11 = (1-eps)*(1-eps)+eps*eps;

  for(int s=0;s<len;s++){
    int whmax=0;
    for(int i=1;i<3;i++){
      //      fprintf(stderr,"%d %d | %f %f\n",s,i,gls[s][i],gls[s][whmax]);
      if(gls[s][i]>gls[s][whmax])
	whmax=i;
    }
    //    fprintf(stdout,"i\t%d\n",whmax);
    if(whmax==0){
      gls[s][0] = g00;gls[s][1]=g01;gls[s][2]=g02;
    }else if(whmax==1){
      gls[s][0]=gls[s][2]=g10;
      gls[s][1]=g11;
    }else if(whmax==2){
      gls[s][0] = g02;gls[s][1]=g01;gls[s][2]=g00;
    }else{
      fprintf(stderr,"never happens\n");
    }
      
  }

}

void callgenotypesHwe(double **gls,int len,double eps,double *freq){
  fprintf(stderr,"\t-> Call genotypesHwe: %f\n",eps);

  for(int s=0;s<len;s++){
    gls[s][0] = gls[s][0]*freq[s]*freq[s];
    gls[s][1] = 2*gls[0][1]*(1-freq[s])*freq[s];
    gls[s][2] = gls[s][2]*(1-freq[s])*(1-freq[s]);

    int whmax=0;
    for(int i=1;i<3;i++){
      //      fprintf(stderr,"%d %d | %f %f\n",s,i,gls[s][i],gls[s][whmax]);
      if(gls[s][i]>gls[s][whmax])
	whmax=i;
    }
    //    fprintf(stdout,"i\t%d\n",whmax);
    for(int i=1;i<3;i++)
      gls[s][i] = 0;
    gls[s][whmax]=1;
  }

}

struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};

typedef struct{
  char *chr;
  int pos;
}gpos;


struct ltstr2
{
  bool operator()(const gpos &s1, const gpos &s2) const
  {
    if(strcmp(s1.chr,s2.chr)==0)
      return s1.pos< s2.pos;
    else
      return strcmp(s1.chr, s2.chr) < 0;
  }
};



typedef struct{
  int major;
  int minor;
  double freq;
}datum;



typedef std::map<const char *,datum,ltstr> rsMap;


typedef std::map<const gpos,datum,ltstr2> posMap;



posMap getBim(char *bname,char *fname){
  //first generate rsnumber to freq
  char *buf = new char[LENS];
  gzFile gz = Z_NULL;
  gz = gzopen(fname,"rb");
  assert(gz!=Z_NULL);
  gzgets(gz,buf,LENS);//header
  
  rsMap rMap;
  while(  gzgets(gz,buf,LENS)){
    strtok(buf,"\t\n ");
    char *rs = strdup(strtok(NULL,"\t\n "));
    int A1 = refToInt[strtok(NULL,"\t\n ")[0]];
    int A2 = refToInt[strtok(NULL,"\t\n ")[0]];
    char *tok = strtok(NULL,"\t\n ");
    if(strcmp(tok,"NA")==0)
      continue;
    double freq = atof(tok);
    datum d;
    d.minor = A1;//http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#freq
    d.major = A2;//always forget which one if major and minor
    d.freq = freq;
    assert(rMap.find(rs)==rMap.end());
    rMap[rs]= d;
  }
  fprintf(stderr,"nfreqs:%lu read from plink.frq:\'%s\' file\n",rMap.size(),fname);
  gzclose(gz);gz=Z_NULL;
  gz = gzopen(bname,"rb");
  assert(gz!=Z_NULL);

  posMap pm;
  int linenr=0;
  while(  gzgets(gz,buf,LENS)){
    linenr++;
    gpos gp;
    gp.chr = strdup(strtok(buf,"\t\n "));
    char *rs = strtok(NULL,"\t\n ");
    strtok(NULL,"\t\n ");
    gp.pos = atoi(strtok(NULL,"\t\n "));
    //check that position is new
    assert(pm.find(gp)==pm.end());

    rsMap::iterator rit = rMap.find(rs);
    if(rit == rMap.end()){
      fprintf(stderr,"\t-> rsnumber:%s from bimfile[%d]:\'%s\' doesn't exists in freqfile: \'%s\' will continue filereading but discarding the site \n",rs,linenr,bname,fname);
      //exit(0);
    }
    pm[gp] = rit->second;
    
  }
  fprintf(stderr,"nsites:%lu read from plink.bim file\n",pm.size());


  return pm;
}


int extract_freq(int argc,char **argv){
  ++argv;
  char *pfile,*bfile,*ffile;
  pfile=bfile=ffile=NULL;
  pfile =*argv++;
  bfile =*argv++;
  ffile =*argv++;
  fprintf(stderr,"posfile:%s bimfile:%s ffile:%s\n",pfile,bfile,ffile);
  assert(pfile &&bfile &&ffile);
  posMap pm = getBim(bfile,ffile);
  //
  char *buf = new char[LENS];
  gzFile gz = Z_NULL;
  gz = gzopen(pfile,"rb");
  assert(gz!=Z_NULL);
  while(gzgets(gz,buf,LENS)){
    gpos gp;
    gp.chr = strtok(buf,"\t\n ");
    gp.pos = atoi(strtok(NULL,"\t\n "));
    //    fprintf(stderr,"chr:%s pos:%d\n",gp.chr,gp.pos);
    char major = refToInt[strtok(NULL,"\t\n ")[0]];
    char minor = refToInt[strtok(NULL,"\t\n ")[0]];
    posMap::iterator it = pm.find(gp);
    if(it==pm.end()){
      fprintf(stderr,"\t-> Problem finding chr:%s pos:%d from pos.glf.gz in plinkfiles\n",gp.chr,gp.pos);
      exit(0);
    }
    if(major!=it->second.major&&major!=it->second.minor){
      fprintf(stderr,"\t-> major from glf.pos.gz is not defined properly in majorminor from plink at site: (%s,%d) \n",gp.chr,gp.pos);
      exit(0);
    }
    if(minor!=it->second.major&&major!=it->second.minor){
      fprintf(stderr,"\t-> minor from glf.pos.gz is not defined properly in majorminor from plink at site: (%s,%d)\n",gp.chr,gp.pos);
      exit(0);
    }
    if(major!=it->second.major)
      fprintf(stdout,"%f\n",1-it->second.freq);
    else
      fprintf(stdout,"%f\n",it->second.freq);
  }
  
  
  return 0;
}


int main(int argc, char **argv){
  if(argc==1)
    print_info(stderr);

  if(strcasecmp(argv[1],"extract_freq")==0)
    return extract_freq(--argc,++argv);
  
  char *freqname=NULL;
  char *gname=NULL;
  int maxIter =10;
  double tole =1e-6;
  int n=-1;
  int seed=100;
  int model =1;
  int gc =0;
  double errate = 0.005;
  int pair1 =-1;
  int pair2 =-1;
  int nind =2;
  int doInbreed=0;
  int swichMaf = 0;
  int verbose = 0;
  double minMaf =0.05;
  while ((n = getopt(argc, argv, "f:i:t:r:g:m:v:s:F:c:e:a:b:n:l:")) >= 0) {
    switch (n) {
    case 'f': freqname = strdup(optarg); break;
    case 'i': maxIter = atoi(optarg); break;
    case 't': tole = atof(optarg); break;
    case 'r': seed = atoi(optarg); break;
    case 'g': gname = strdup(optarg); break;
    case 'm': model = atoi(optarg); break;      
    case 'v': verbose = atoi(optarg); break;      
    case 's': swichMaf = atoi(optarg); break;      
    case 'F': doInbreed = atoi(optarg); break;      
    case 'c': gc = atoi(optarg); break;      
    case 'a': pair1 = atoi(optarg); break;      
    case 'b': pair2 = atoi(optarg); break;      
    case 'n': nind = atoi(optarg); break;      
    case 'e': errate = atof(optarg); break;      
    case 'l': minMaf = atof(optarg); break;      
    default: {fprintf(stderr,"unknown arg:\n");return 0;}
      print_info(stderr);
    }
  }
  srand48(seed);
  if(nind==-1||freqname==NULL||gname==NULL){
    fprintf(stderr,"\t-> Must supply -n -f -g parameters (%d,%s,%s)\n",nind,freqname,gname);
    return 0;
  }
  std::vector<double> freq = getDouble(freqname);
  if(swichMaf){
    fprintf(stderr,"swiching frequencies\n");
    for(int i=0;i<freq.size();i++)
      freq[i] = 1 - freq[i];
  }
 
  double **gls = getGL(gname,freq.size(),nind);
#if 0
  print(stdout,freq.size(),3*nind,gls);
  exit(0);
#endif
   double **l1,**l2;l1=l2=NULL;
   double *newfreq = NULL;
  l1=new double*[freq.size()];
  l2=new double*[freq.size()];
  newfreq =new double[freq.size()];
  double **emis=new double*[freq.size()];

  for(int i=0;i<freq.size();i++){
    l1[i] = new double[3];
    l2[i] = new double[3];
    emis[i] = new double[3];
  }
  if(doInbreed)
    fprintf(stdout,"Pair\tZ=0\tZ=1\tloglh\tnIter\tcoverage\n");
  else
    fprintf(stdout,"Pair\tk0\tk1\tk2\tloglh\tnIter\tcoverage\n");

  if(doInbreed){
    for(int a=0;a<nind;a++){
      int nkeep=0;

      for(int i=0;i<freq.size();i++){
	for(int j=0;j<3;j++)
	  l1[nkeep][j] = gls[i][a*3+j];
	if(l1[nkeep][0]==l1[nkeep][1]&&l1[nkeep][0]==l1[nkeep][2])
	  continue;
	if(freq[i]<minMaf||(1-freq[i])<minMaf)
	  continue;
	newfreq[nkeep]=freq[i];
	nkeep++;
      }
      fprintf(stdout,"(%d,%d)\t",a,a);
      //call genotypes
      if(gc){
	if(gc>1)
	  callgenotypesHwe(l1,nkeep,errate,newfreq);
	if(gc>0)
	  callgenotypesEps(l1,nkeep,errate);
      }

      emission_ngsInbreed(newfreq,l1,emis,nkeep);
      double pars[2];
      pars[0]=drand48();
      pars[1]=1-pars[0];

      int niter;
      
      niter=em_inbreed(pars,emis,tole,maxIter,nkeep,model,verbose);
      double lopt= loglikeInbreed(pars,emis,nkeep);
      fprintf(stdout,"%f\t%f\t%f\t%d\t%f\n",pars[0],pars[1],lopt,niter,(1.0*nkeep)/(1.0*freq.size()));
    }  
  }
  else{
  for(int a=0;a<nind;a++){
    for(int b=a+1;b<nind;b++){
      //      fprintf(stderr,"a:%d b:%d pair1:%d pair2:%d\n",a,b,pair1,pair2);
      int nkeep=0;
      if(pair1!=-1)
	a=pair1;
      if(pair2!=-1)
	b=pair2;

      for(int i=0;i<freq.size();i++){
	//copoy data into l1, this might be overwritten at next iteration if, if either is missing or freq<minmaf
	for(int j=0;j<3;j++){
	  l1[nkeep][j] = gls[i][a*3+j];
	  l2[nkeep][j] = gls[i][b*3+j];
	}
	
	if(l1[nkeep][0]==l1[nkeep][1]&&l1[nkeep][0]==l1[nkeep][2])
	  continue;
	if(l2[nkeep][0]==l2[nkeep][1]&&l2[nkeep][0]==l2[nkeep][2])
	  continue;
	if(freq[i]<minMaf||(1-freq[i])<minMaf)
	  continue;
	
	newfreq[nkeep]=freq[i];
	
	nkeep++;
	
      }
      //print(stdout,freq.size(),6,gls);
      //return 0;
      //      fprintf(stdout,"(%d,%d)\t%f\t",a,b,nkeep/(1.0*freq.size()));
      fprintf(stdout,"(%d,%d)\t",a,b);
      if(gc){
	if(gc>1){
	  callgenotypesHwe(l1,nkeep,errate,newfreq);
	  callgenotypesHwe(l2,nkeep,errate,newfreq);
	}

	if(gc>0){
	  callgenotypesEps(l1,nkeep,errate);
	  callgenotypesEps(l2,nkeep,errate);
	}
      }
      emission_ngsrelate(newfreq,l1,l2,emis,nkeep);
      //      print(stdout,freq.size(),3,emis);return 0;
      double pars[3];
      pars[0] = drand48();
      pars[1] = drand48()*(1-pars[0]);
      pars[2] = 1-pars[0]-pars[1];
      //      fprintf(stdout,"loglike:%f\t",loglike(pars,emis,nkeep));return 0;
      int niter;
      if(model==0)
	niter=em1(pars,emis,tole,maxIter,nkeep,verbose);
      else if(model==1)
	niter=em2(pars,emis,tole,maxIter,nkeep,verbose);
      else//below might not work
	niter=em3(pars,emis,tole,maxIter,nkeep,verbose);
      
      double p100[3]={1-TINY,TINY/2.0,TINY/2.0};
      double p010[3]={TINY/2.0,1-TINY,TINY/2.0};
      double p001[3]={TINY/2.0,TINY/2.0,1-TINY};
      double l100=loglike(p100,emis,nkeep);
      double l010=loglike(p010,emis,nkeep);
      double l001=loglike(p001,emis,nkeep);
      double lopt= loglike(pars,emis,nkeep);
      //  fprintf(stderr,"%f %f %f\n",l100,l010,l001);
      
      double likes[4] = {l100,l010,l001,lopt};
      int best = 0;
      for(int i=0;i<4;i++){
	if(likes[i]>likes[best])
	  best=i;
      }
      // fprintf(stderr,"100:%f 010:%f 001:%f lopt:%f\n",l100,l010,l001,lopt);
      if(best==3)
	fprintf(stdout,"%f\t%f\t%f\t%f\t%d\t%f\n",pars[0],pars[1],pars[2],lopt,niter,(1.0*nkeep)/(1.0*freq.size()));
      if(best==0)
	fprintf(stdout,"%f\t%f\t%f\t%f\t%d\t%f\n",p100[0],p100[1],p100[2],l100,-1,(1.0*nkeep)/(1.0*freq.size()));
      if(best==1)
	fprintf(stdout,"%f\t%f\t%f\t%f\t%d\t%f\n",p010[0],p010[1],p010[2],l010,-1,(1.0*nkeep)/(1.0*freq.size()));
      if(best==2)
	fprintf(stdout,"%f\t%f\t%f\t%f\t%d\t%f\n",p001[0],p001[1],p001[2],l001,-1,(1.0*nkeep)/(1.0*freq.size()));
      
      if(pair1!=-1||pair2!=-1){
	//	fprintf(stderr,"BREAKING\n");
	break;
      }
    }
    if(pair1!=-1||pair2!=-1){
      //	fprintf(stderr,"BREAKING\n");
	break;
    }
  }
  }
for(int i=0;i<freq.size();i++){
  delete [] gls[i];
  delete [] l1[i];
    delete [] l2[i];
    delete [] emis[i];
  }
  delete [] gls;
  delete [] l1;
  delete [] l2;
  delete [] emis;
  free(freqname);
  free(gname);
  return 0;
}
