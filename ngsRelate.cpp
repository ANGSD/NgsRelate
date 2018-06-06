/*
  NgsRelateV2
  http://www.popgen.dk/software/
  https://github.com/ANGSD/NgsRelate/

  g++ NgsRelate.cpp -O3 -o ngsrelate -lz -lpthread

*/

#include <vector>
#include <cstring>
#include <zlib.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <map>
#include <libgen.h>  // basename
#include <pthread.h> // threading
#ifdef __WITH_BCF__
#include "vcf.h"
#endif

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
double TINY=1e-6;
double p100000000[9] = {1 - TINY,   TINY / 8.0, TINY / 8.0,
                        TINY / 8.0, TINY / 8.0, TINY / 8.0,
                        TINY / 8.0, TINY / 8.0, TINY / 8.0};

double p010000000[9] = {TINY / 8.0, 1 - TINY,   TINY / 8.0,
                        TINY / 8.0, TINY / 8.0, TINY / 8.0,
                        TINY / 8.0, TINY / 8.0, TINY / 8.0};

double p001000000[9] = {TINY / 8.0, TINY / 8.0, 1 - TINY,
                        TINY / 8.0, TINY / 8.0, TINY / 8.0,
                        TINY / 8.0, TINY / 8.0, TINY / 8.0};

double p000100000[9] = {TINY / 8.0, TINY / 8.0, TINY / 8.0,
                        1 - TINY,   TINY / 8.0, TINY / 8.0,
                        TINY / 8.0, TINY / 8.0, TINY / 8.0};

double p000010000[9] = {TINY / 8.0, TINY / 8.0, TINY / 8.0,
                        TINY / 8.0, 1 - TINY,   TINY / 8.0,
                        TINY / 8.0, TINY / 8.0, TINY / 8.0};

double p000001000[9] = {TINY / 8.0, TINY / 8.0, TINY / 8.0,
                        TINY / 8.0, TINY / 8.0, 1 - TINY,
                        TINY / 8.0, TINY / 8.0, TINY / 8.0};

double p000000100[9] = {TINY / 8.0, TINY / 8.0, TINY / 8.0,
                        TINY / 8.0, TINY / 8.0, TINY / 8.0,
                        1 - TINY,   TINY / 8.0, TINY / 8.0};

double p000000010[9] = {TINY / 8.0, TINY / 8.0, TINY / 8.0,
                        TINY / 8.0, TINY / 8.0, TINY / 8.0,
                        TINY / 8.0, 1 - TINY,   TINY / 8.0};

double p000000001[9] = {TINY / 8.0, TINY / 8.0, TINY / 8.0,
                        TINY / 8.0, TINY / 8.0, TINY / 8.0,
                        TINY / 8.0, TINY / 8.0, 1 - TINY};

double p10[2]={1-TINY,TINY};
double p01[2]={TINY,1-TINY};

int num_threads = 4;
char *freqname=NULL;
char *gname=NULL;
int maxIter =5000;
double tole =1e-6;
int n=-1;
int seed=100;
int model =1;
int gc =0;
double errate = 0.005;
int pair1 =-1;
int pair2 =-1;
int nind =2;
int do_inbred=0;
int switchMaf = 0;
int verbose = 0;
double minMaf =0.05;
int hasDef = 0;
double ttol=1e-6;
std::vector<char *> ids;

float emTole=1e-12;

void stayin(double *post){
  for(int i=0;i<9;i++){
    if(post[i]<emTole)
      post[i] = emTole;
    if(post[i]>(1-emTole))
      post[i] =1- emTole;
  }
}

double access_genotype(double **gls, const int & site, const int & indi, const int & geno){
  return gls[site][indi * 3 + geno];
}

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
  double ret =0;
  for(int i=0;i<len;i++){
    double tmp = 0;
    for(int j=0;j<9;j++)
      tmp += p[j]*emis[i][j];
    ret +=log(tmp);
  }
  return ret;
}

// megatest

double loglike_inbred(double *p,double **emis,int len){
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



void emStep1_inbred(double *pre,double **emis,double *post,int len){
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
  for(int i=0;i<9;i++){
    if(pre[i]<0||pre[i]>1){
      fprintf(
          stderr,
          "Problem with gues in emStep1: pres:(%f,%f,%f,%f,%f,%f,%f,%f,%f)\n",
          pre[0], pre[1], pre[2], pre[3], pre[4], pre[5], pre[6], pre[7],
          pre[8]);
      exit(0);
    }
  }
  double inner[9];
  for(int x=0;x<9;x++){
    post[x] =0.0;
  }
  for(int i=0;i<len;i++){
    for(int x=0;x<9;x++){
      inner[x] = pre[x]*emis[i][x];
    }
    normalize(inner,9);
    for(int x=0;x<9;x++){
      post[x] += inner[x];
    }

  }
  //set bounds
  // fprintf(stderr,"halli haloo\n");
#if 0
  stayin(post);
#endif
  normalize(post,9);
  for(int i=0;i<9;i++){
    if(post[i]<0||post[i]>1){
      fprintf(stderr,
              "Probable after normalizing: pres:(%f,%f,%f,%f,%f,%f,%f,%f,%f) "
              "post:(%f,%f,%f,%f,%f,%f,%f,%f,%f)\n",
              pre[0], pre[1], pre[2], pre[3], pre[4], pre[5], pre[6], pre[7],
              pre[8], post[0], post[1], post[2], post[3], post[4], post[5],
              post[6], post[7], post[8]);
    }
  }
}

void minus(double * fst,double * sec,double * res){
  for(int i=0;i<9;i++)
    res[i] = fst[i]-sec[i];
}
void minus2(double fst[2],double sec[2],double res[2]){
  for(int i=0;i<2;i++)
    res[i] = fst[i]-sec[i];
}

double sumSquare(double * mat){
  double tmp=0;
  for(size_t i=0;i<9;i++){
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

int emAccel(double *F,double **emis,double *F_new,int len, int & niter){
  //  maybe these should be usersettable?

  double stepMin =1;
  double stepMax0 = 1;
  static double stepMax=stepMax0;
  double mstep=4;
  //  double objfnInc=1;


  double F_em1[9];
  double F_diff1[9];
  double F_em2[9];
  double F_diff2[9];
  double F_diff3[9];
  double F_tmp[9];
  niter++;
  emStep1(F, emis, F_em1, len);
  // stayin(F_em1);
  
  minus(F_em1, F, F_diff1);

  // for (int ix=0; ix<9;ix++)
  //   fprintf(stderr,"F_diff1: %d: %lf, %lf, %lf\n",ix,F_diff1[ix], F_em1[ix], F[ix]);

  double sr2 = sumSquare(F_diff1);
  
  if(sqrt(sr2)<ttol){
    // fprintf(stderr,"sr2 break: %e\n", sqrt(sr2));
    return 0;
  }
  niter++;
  emStep1(F_em1, emis, F_em2, len);
  minus(F_em2, F_em1, F_diff2);

  double sq2 = sumSquare(F_diff2);

  if(sqrt(sq2)<ttol){
    // fprintf(stderr,"sq2 break: %e\n", sqrt(sq2));
    return 0;
  }

  minus(F_diff2,F_diff1, F_diff3);

  double sv2 = sumSquare(F_diff3);

  double alpha = sqrt(sr2/sv2);
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  for(size_t i=0;i<9;i++){
    F_new[i] = F[i]+2*alpha*F_diff1[i]+alpha*alpha*F_diff3[i];
  }
  
  // fprintf(stderr,"%d: F_new (this the linear jump: (%f,%f,%f,%f,%f,%f,%f,%f,%f)\n", niter,F_new[0],F_new[1],F_new[2],F_new[3],F_new[4],F_new[5],F_new[6],F_new[7],F_new[8]);
  
  int outofparspace =0;
  for(int i=0;i<9;i++){
    if(F_new[i]<0||F_new[i]>1){
      outofparspace++;
      // break;
    }
  }
  if(outofparspace){
    fprintf(stderr,"outofparspace will use second emstep as jump\n");
    for(int i=0;i<9;i++)
      F_new[i] = F_em2[i];
  }

  if (fabs(alpha - 1) > 0.01){
    niter++;
    emStep1(F_new,emis,F_tmp,len);
    for(int i=0;i<9;i++)
      std::swap(F_new[i],F_tmp[i]);
  }

  if ((alpha - stepMax) > -0.001) {
    stepMax = mstep*stepMax;
  }

  return 1;
}


int emAccel_inbred(double *F,double **emis,double *F_new,int len){
  //  fprintf(stderr,"calling emaccel \n");
  double ttol=0.0000001;

  //  fprintf(stderr,"tol:%f\n",tol);
  //maybe these should be usersettable?
  double stepMin =1;
  double stepMax0 = 1;
  static double stepMax=stepMax0;
  double mstep=4;
  //  double objfnInc=1;


  double F_em1[2];
  double F_diff1[2];
  double F_em2[2];
  double F_diff2[2];
  double F_diff3[2];
  double F_tmp[2];

  emStep1_inbred(F,emis,F_em1,len);
  minus2(F_em1,F,F_diff1);
  double sr2 = sumSquare2(F_diff1);
  
  if(sqrt(sr2)<ttol){
    //    fprintf(stderr,"sr2 break:%f\n",sr2);
    return 0;
    //break;
  }
  emStep1_inbred(F_em1,emis,F_em2,len);
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
    emStep1_inbred(F_new,emis,F_tmp,len);
    for(int i=0;i<2;i++)
      std::swap(F_new[i],F_tmp[i]);
  }

  //  double lnew = 1;
  if ((alpha - stepMax) > -0.001) {
    stepMax = mstep*stepMax;
  }
  //  print(stderr,3,F_new);
  //  fprintf(stderr,"alpha %f stepMax %f\n",alpha,stepMax);

  // fprintf(stderr,"calling emaccel \n");
  return 1;
 
}


int em1(double *sfs,double  **emis, int len){
  int niter = 0;
  double oldLik,lik;
  oldLik = loglike(sfs,emis,len);
  if(verbose)
    fprintf(stderr,"startlik=%f est: %f %f %f %f %f %f %f %f %f\n",oldLik,sfs[0],sfs[1],sfs[2],sfs[3],sfs[4],sfs[5],sfs[6],sfs[7],sfs[8]);
  fflush(stderr);

  double tmp[9];
  int it;
  for(it=0;niter<maxIter;it++) {
    niter++;
    emStep1(sfs,emis,tmp,len);
    for(int i=0;i<9;i++)
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
  return niter;
}

//   niter=em_inbred(pars,emis,tole,maxIter,nkeep,model,verbose);
int em_inbred(double *pars,double  **emis,double tole,int maxIter,int len,int model,int verbose){
  double oldLik,lik;
  oldLik = loglike_inbred(pars,emis,len);
  if(verbose){
    fprintf(stderr,"startlik=%f %f %f\n",oldLik,pars[0],pars[1]);
    fflush(stderr);
  }

  double tmp[2];
  int it;
  for(it=0;it<maxIter;it++) {
    if(model==0)
      emStep1_inbred(pars,emis,tmp,len);
    else
      emAccel_inbred(pars,emis,tmp,len);
    for(int i=0;i<2;i++)
      pars[i]= tmp[i];
    lik = loglike_inbred(pars,emis,len);
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


int em2(double *sfs,double  **emis, int len){
  int niter=0;
  double oldLik,lik;
  oldLik = loglike(sfs,emis,len);
  if(verbose){
    fprintf(stderr,"startlik=%f\n",oldLik);
    fflush(stderr);
  }

  double tmp[9];
  int it;
  for(it=0;niter<maxIter;it++) {
    emAccel(sfs,emis,tmp,len, niter);

    for(int i=0;i<9;i++)
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
  return niter;
}

int em3(double *sfs,double  **emis, int len){
  //  exit(0);
  int niter = 0;
  double oldLik,lik;
  oldLik = loglike(sfs,emis,len);
  if(verbose){
    fprintf(stderr,"em3startlik=%f (%f,%f,%f)\n",oldLik,sfs[0],sfs[1],sfs[2]);
    fflush(stderr);
  }

  double tmp[9];
  int it;
  int speedy=1;
  for(it=0;niter<maxIter;it++) {
    if(speedy)
      emAccel(sfs,emis,tmp,len, niter);
    else
      niter++;
      emStep1(sfs,emis,tmp,len);
    lik = loglike(tmp,emis,len);

    fprintf(stderr,"[%d]:%d (%f,%f,%F) lik=%f diff=%e\n",it,speedy,sfs[0],sfs[1],sfs[2],lik,fabs(lik-oldLik));
    if(std::isnan(lik)||lik<oldLik){
      fprintf(stderr,"Problem llh is now bigger or nan, will go back and use regular em\n");
      fprintf(stderr,"This is offending pars: %f %f %f\n",sfs[0],sfs[1],sfs[2]);
      speedy=0;
      continue;
    }
    for(int i=0;i<9;i++)
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

void emission_ngsrelate9(std::vector<double> * freq, double **gls, double **emis, int *keeplist, int & nkeep, int & ind1, int & ind2 ){
  // access_genotype(td->gls, i, td->a, 0);
  int i;
  for(int x=0;x<nkeep;x++){
    i = keeplist[x];
    double freqa=freq->at(i);  // alternative allele frequency
    double freqA=1-freqa;
    
    // i == freqA == freq0
    // j == freqa == freq1
    // emis<- cbind(freq0,freq0^2,freq0^2, freq0^3,freq0^2,freq0^3, freq0^2,freq0^3,freq0^4)*gl1[1,]*gl2[1,]
    // ##00&00
    // G_real=(AA,AA)
    double AAAA = access_genotype(gls, i, ind1, 0) * access_genotype(gls, i, ind2, 0);
    emis[x][0] = pow(freqA, 4) * AAAA;
    emis[x][1] = pow(freqA, 3) * AAAA;
    emis[x][2] = pow(freqA, 2) * AAAA;
    emis[x][3] = pow(freqA, 3) * AAAA;
    emis[x][4] = pow(freqA, 2) * AAAA;
    emis[x][5] = pow(freqA, 3) * AAAA;
    emis[x][6] = pow(freqA, 2) * AAAA;
    emis[x][7] = pow(freqA, 2) * AAAA;
    emis[x][8] = freqA * AAAA;

    // emis <- emis + cbind(freq1,freq1^2,freq1^2, freq1^3,freq1^2,freq1^3, freq1^2,freq1^3,freq1^4)*gl1[3,]*gl2[3,]
    // ##11&11
    // G_real=(aa,aa)
    double aaaa = access_genotype(gls, i, ind1, 2) * access_genotype(gls, i, ind2, 2);
    emis[x][0] += pow(freqa, 4) * aaaa;
    emis[x][1] += pow(freqa, 3) * aaaa;
    emis[x][2] += pow(freqa, 2) * aaaa;
    emis[x][3] += pow(freqa, 3) * aaaa;
    emis[x][4] += pow(freqa, 2) * aaaa;
    emis[x][5] += pow(freqa, 3) * aaaa;
    emis[x][6] += pow(freqa, 2) * aaaa;
    emis[x][7] += pow(freqa, 2) * aaaa;
    emis[x][8] += freqa * aaaa;

    // emis <- emis + cbind(0,freq0*freq1,0,    freq0*freq1^2,0,freq0^2*freq1  ,0,0,freq1^2*freq0^2)*gl1[1,]*gl2[3,]
    // ##00&11
    // G_real=(AA,aa)
    double AAaa = access_genotype(gls, i, ind1, 0) * access_genotype(gls, i, ind2, 2);
    emis[x][0] += pow(freqA, 2) * pow(freqa, 2) * AAaa;
    emis[x][1] += 0;
    emis[x][2] += 0;
    emis[x][3] += pow(freqA, 2) * freqa * AAaa;
    emis[x][4] += 0;
    emis[x][5] += freqA * pow(freqa,2) * AAaa;
    emis[x][6] += 0;
    emis[x][7] += freqA * freqa * AAaa;
    emis[x][8] += 0;

    // emis <- emis + cbind(0,freq0*freq1,0,    freq0^2*freq1,0,freq1^2*freq0,  0,0,freq1^2*freq0^2)*gl1[3,]*gl2[1,]
    // ##11&00
    // G_real=(aa,AA)
    double aaAA = access_genotype(gls, i, ind1, 2) * access_genotype(gls, i, ind2, 0);
    emis[x][0] += pow(freqA, 2) * pow(freqa, 2) * aaAA;
    emis[x][1] += 0;
    emis[x][2] += 0;
    emis[x][3] += freqA * pow(freqa, 2) * aaAA;
    emis[x][4] += 0;
    emis[x][5] += pow(freqA, 2) * freqa * aaAA;
    emis[x][6] += 0;
    emis[x][7] += freqA * freqa * aaAA;
    emis[x][8] += 0;

    // emis <- emis + cbind(0,0,freq0*freq1, 2*freq0^2*freq1,0,0, 0,freq0^2*freq1,2*freq0^3*freq1)*gl1[1,]*gl2[2,]
    // ##00&01
    // G_real=(AA,Aa)
    double AAAa = access_genotype(gls, i, ind1, 0) * access_genotype(gls, i, ind2, 1);
    emis[x][0] += 2 * pow(freqA, 3) * freqa * AAAa;
    emis[x][1] += pow(freqA, 2) * freqa * AAAa;
    emis[x][2] += 0;
    emis[x][3] += 0;
    emis[x][4] += 0;
    emis[x][5] += 2 * pow(freqA, 2) * freqa * AAAa;
    emis[x][6] += freqA * freqa * AAAa;
    emis[x][7] += 0;
    emis[x][8] += 0;

    // emis <- emis + cbind(0,0,freq0*freq1, 2*freq1^2*freq0,0,0, 0,freq1^2*freq0,2*freq1^3*freq0)*gl1[3,]*gl2[2,]
    // ##11&01
    // # G_real=(aa,Aa)
    double aaAa = access_genotype(gls, i, ind1, 2) * access_genotype(gls, i, ind2, 1);
    emis[x][0] += 2 * freqA * pow(freqa, 3) * aaAa;
    emis[x][1] += freqA * pow(freqa, 2) * aaAa;
    emis[x][2] += 0;
    emis[x][3] += 0;
    emis[x][4] += 0;
    emis[x][5] += 2 * freqA * pow(freqa, 2) * aaAa;
    emis[x][6] += freqA * freqa * aaAa;
    emis[x][7] += 0;
    emis[x][8] += 0;

    // emis <- emis + cbind(0,0,0, 0,freq0*freq1,2*freq0^2*freq1, 0,freq0^2*freq1,2*freq0^3*freq1)*gl1[2,]*gl2[1,]
    // ##01&00
    // G_real=(Aa,AA)
    double AaAA = access_genotype(gls, i, ind1, 1) * access_genotype(gls, i, ind2, 0);
    emis[x][0] += 2 * pow(freqA, 3) * freqa * AaAA;
    emis[x][1] += pow(freqA, 2) * freqa * AaAA;
    emis[x][2] += 0;
    emis[x][3] += 2 * pow(freqA, 2) * freqa * AaAA;
    emis[x][4] += freqA * freqa * AaAA;
    emis[x][5] += 0;
    emis[x][6] += 0;
    emis[x][7] += 0;
    emis[x][8] += 0;

    // emis <- emis + cbind(0,0,0, 0,freq0*freq1,2*freq1^2*freq0, 0,freq1^2*freq0,2*freq1^3*freq0)*gl1[2,]*gl2[3,]
    // ##01&11
    // G_real=(Aa,aa)
    double Aaaa = access_genotype(gls, i, ind1, 1) * access_genotype(gls, i, ind2, 2);
    emis[x][0] += 2 * freqA * pow(freqa, 3) * Aaaa;
    emis[x][1] += freqA * pow(freqa, 2) * Aaaa;
    emis[x][2] += 0;
    emis[x][3] += 2 * freqA * pow(freqa, 2) * Aaaa;
    emis[x][4] += freqA * freqa * Aaaa;
    emis[x][5] += 0;
    emis[x][6] += 0;
    emis[x][7] += 0;
    emis[x][8] += 0;

    // emis <- emis + cbind(0,0,0, 0,0,0, 2*freq0*freq1,freq1*freq0,4*freq1^2*freq0^2)*gl1[2,]*gl2[2,]
    // ##01&01 S7
    // G_real=(Aa,Aa)
    double AaAa = access_genotype(gls, i, ind1, 1) * access_genotype(gls, i, ind2, 1);
    emis[x][0] += 4 * pow(freqA, 2) * pow(freqa, 2) * AaAa;
    emis[x][1] += freqA * freqa * AaAa;
    emis[x][2] += 2 * freqA * freqa * AaAa;
    emis[x][3] += 0;
    emis[x][4] += 0;
    emis[x][5] += 0;
    emis[x][6] += 0;
    emis[x][7] += 0;
    emis[x][8] += 0;

  }

  
#if 0
  // this is just to test compare the emission matrix of the first x snps.
  // https://overiq.com/c-programming/101/fwrite-function-in-c/
  FILE *fp;
  fp = fopen("new.txt", "wb");
  if(fp == NULL){
    printf("Error opening file\n");
    exit(1);
  }
  for(int i=0;i<nkeep;i++){
    for(int x=0;x<9;x++){
      fwrite(&emis[i][x], sizeof(double), 1, fp);
    }
  }
  fclose(fp);
  
  // for(int i=0;i<nkeep;i++){
  //   fprintf(stdout, "%f", emis[i][0]);
  //   for (int x=1;x<9;x++){
  //     fprintf(stdout, "\t%f", emis[i][x]);
  //   }
  //   fprintf(stdout, "\n");
  // }
  exit(0);
#endif
}

void emission_ngs_inbred(std::vector<double> * freq, double **gls, double **emis, int *keeplist, int & nkeep, int & ind1){
  int i;
  double AA, Aa, aa, freqa, freqA;
  for(int x=0;x<nkeep;x++){
    i = keeplist[x];
    freqa=freq->at(i);  // alternative allele frequency
    freqA=1-freqa;
    // i might have to flip it: this was the old approach in ngsRelate:
    // double freqA=freq[i];
    // double freqa=1-freqA;

    AA = access_genotype(gls, i, ind1, 0);
    Aa = access_genotype(gls, i, ind1, 1);
    aa = access_genotype(gls, i, ind1, 2);

    // G_real=(AA)
    emis[i][0] = pow(freqA,2)*AA;
    emis[i][1] = freqA*AA;

    // G_real=(Aa)
    emis[i][0] += 2*freqa*freqA*Aa;
    emis[i][1] += 0;

    // G_real=(aa)
    emis[i][0] += pow(freqa,2)*aa;
    emis[i][1] += freqa*aa;

    //if(i==295&&0)
    //    fprintf(stderr,"emis[%d]:freq:%f %f %f\n",i,freqa,emis[i][0],emis[i][1]);
      //      exit(0);
    }

}


double **getGL(const char *fname, int sites, int nInd) {
  gzFile gz = Z_NULL;
  if (((gz = gzopen(fname, "rb"))) == Z_NULL) {
    fprintf(stderr, "\t-> Problem opening file: \'%s\'\n", fname);
    exit(0);
  }

  double **ret = new double *[sites + 10];

  int i = 0;
  while (1) {
    //   for(int i=0;i<sites;i++){
    ret[i] = new double[3 * nInd];
    unsigned nbit = gzread(gz, ret[i], sizeof(double) * nInd * 3);
    if (nbit == 0)
      break;
    if (sizeof(double) * nInd * 3 != nbit) {
      fprintf(stderr, "\t-> Problem reading full chunk\n");
      exit(0);
    }
    for (int g = 0; g < 3 * nInd; g++)
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
    if (i > sites) {
      fprintf(stderr, "\t-> Too many sites in glf file. Looks outof sync, or "
                      "make sure you supplied correct number of individuals "
                      "(-n)\n");
      exit(0);
    }
  }
  if (i != sites) {
    fprintf(stderr, "nsites: %d assumed but %d read\n", sites, i);
    exit(0);
  }
  gzclose(gz);
  return ret;
}

size_t getDouble(const char *fname,std::vector<double> &ret) {
  assert(ret.size()==0);
  gzFile gz = Z_NULL;
  if (((gz = gzopen(fname, "r"))) == Z_NULL) {
    fprintf(stderr, "[%s]\t-> Problem opening file:%s\n",__FUNCTION__, fname);
    exit(0);
  }
  int nbytes=10000;
  char *buf = new char[nbytes];
  while (gzgets(gz, buf, nbytes)) {
    // reads line by line
    ret.push_back(atof(buf));
  }
  fprintf(stderr, "\t-> Frequency file: \'%s\' contain %lu number of sites\n",
          fname, ret.size());
  gzclose(gz);
  delete[] buf;
  return ret.size();
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
  fprintf(fp, "   -F <INT>            Estimate inbreeding instead of the nine jacquard coefficients\n");
  fprintf(fp, "   -v <INT>            Verbose. print like per iteration\n");
  fprintf(fp, "   -e <INT>            Errorrates when calling genotypes?\n");
  fprintf(fp, "   -a <INT>            First individual used for analysis? (zero offset)\n");
  fprintf(fp, "   -b <INT>            Second individual used for analysis? (zero offset)\n");
  fprintf(fp, "   -n <INT>            Number of samples in glf.gz\n");
  fprintf(fp, "   -l <INT>            minMaf or 1-Maf filter\n");
  fprintf(fp, "   -z <INT>            Name of file with IDs (optional)\n");
  fprintf(fp, "\n");
  fprintf(fp,"Or\n ./ngsrelate extract_freq_bim pos.glf.gz plink.bim plink.freq\n");
  fprintf(fp,"Or\n ./ngsrelate extract_freq .mafs.gz .pos.glf.gz [-rmTrans]\n");
  fprintf(fp,"Or\n ./ngsrelate -h my.bcf [DEVELOPMENT]\n");
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
    gls[s][1] = 2*gls[s][1]*(1-freq[s])*freq[s];
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
  if(argc==1){
    fprintf(stderr,"\t-> supply mafs.gz glf.gz files\n");
    return 0;
  }
  int rmTrans = 0;
  ++argv;
  char *mfile,*gfile;
  mfile=gfile=NULL;
  mfile =*argv++;
  if(strcasecmp(mfile,"-rmTrans")==0){
    fprintf(stderr,"\t-> Will remove transitions by setting allele frequency to zero\n");
    rmTrans=1;
    mfile=*argv++;
  }
  gfile =*argv++;
  if(strcasecmp(gfile,"-rmTrans")==0){
    fprintf(stderr,"\t-> Will remove transitions by setting allele frequency to zero\n");
    rmTrans=1;
    gfile=*argv++;
  }
  if(*argv&& strcasecmp(*argv,"-rmTrans")==0){
    rmTrans=1;
    fprintf(stderr,"\t-> Will remove transitions by setting allele frequency to zero\n");
  }
  if(*argv&& strcasecmp(*argv,"-rmTrans")!=0)
    fprintf(stderr,"\t-> Unrecognized parameter supplied: \'%s\', only -rmTrans is implemented\n",*argv);

  fprintf(stderr,"\t-> .mafs.gz file:\'%s\' .glf.pos.gz file:\'%s\' \n",mfile,gfile);
  //exit(0);
  assert(mfile &&gfile);

  char *buf = new char[LENS];
  gzFile gz = Z_NULL;
  gz = gzopen(mfile,"rb");
  assert(gz!=Z_NULL);
  gzgets(gz,buf,LENS);//header

  posMap pm;
  while(  gzgets(gz,buf,LENS)){
    gpos gp;
    gp.chr = strdup(strtok(buf,"\t\n "));
    gp.pos = atoi(strtok(NULL,"\t\n "));
    //check that position is new
    posMap::iterator mit = pm.find(gp);
    if(mit!=pm.end()){
      fprintf(stderr,"\t-> Duplicate position detected in freqfile: %s:%d will set freq to zero and thereby discarding site from analysis\n",gp.chr,gp.pos);
      mit->second.freq = 0;
      continue;
    }

    int A1 = refToInt[strtok(NULL,"\t\n ")[0]];
    int A2 = refToInt[strtok(NULL,"\t\n ")[0]];
    double freq = atof(strtok(NULL,"\t\n "));
    if(rmTrans==1){
      if((A1==0&&A2==2)||(A1==2&&A2==0))
        freq=0;
      if((A1==1&&A2==3)||(A1==3&&A2==1))
        freq=0;
    }

    datum d;
    d.minor = A1;//http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#freq
    d.major = A2;//always forget which one if major and minor
    d.freq = freq;
    pm[gp]= d;
  }
  fprintf(stderr,"\t-> nfreqs:%lu read from mafs.gz:\'%s\' file\n",pm.size(),mfile);
  gzclose(gz);gz=Z_NULL;
  gz = gzopen(gfile,"rb");
  assert(gz!=Z_NULL);

  int linenr=0;
  while(  gzgets(gz,buf,LENS)){
    linenr++;
    gpos gp;
    gp.chr = strdup(strtok(buf,"\t\n "));
    gp.pos = atoi(strtok(NULL,"\t\n "));
    //check that position exists
    assert(pm.find(gp)!=pm.end());

    posMap::iterator pit = pm.find(gp);

    int A1 = refToInt[strtok(NULL,"\t\n ")[0]];
    int A2 = refToInt[strtok(NULL,"\t\n ")[0]];
    if((A2==pit->second.major)&&(A1==pit->second.minor))
      fprintf(stdout,"%f\n",pit->second.freq);
    else{
      fprintf(stderr,"\t-> Problem with majorminor types for position %s:%d will set freq to zero and thereby discard the site\n",gp.chr,gp.pos);
      fprintf(stdout,"%f\n",0.0);
    }

  }
  fprintf(stderr,"\t-> number of lines in output: %d\n",linenr);
  return 0;
}

int extract_freq_bim(int argc,char **argv){
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

void readids(std::vector<char *> &ids,char *fname){
  FILE *fp = NULL;
  fp=fopen(fname,"rb");
  if(fp==NULL){
    fprintf(stderr,"\t-> Problem reading file: \'%s\'\n",fname);
    exit(0);
  }
  char *buf= new char[LENS];
  while(fgets(buf,LENS,fp)){
    char *tmp = strdup(basename(buf));
    tmp[strlen(tmp)-1]='\0';
    ids.push_back(tmp);

  }
  fprintf(stderr,"\t-> Number of file names read:%lu\n",ids.size());
  fclose(fp);
}

struct worker_args {
  int thread_id, nkeep, a,b, niter, best;
  double **gls;
  std::vector<double> * freq;
  double ll, bestll;
  double pars[9];
  worker_args(int & id_a, int & id_b, std::vector<double> * f, double ** gls_arg ){
    a=id_a;
    b=id_b;
    nkeep=0;
    best=0;
    bestll=0.0;
    bestll=0.0;    
    freq = f;
    gls=gls_arg;
  }
};


void * do_work(void *threadarg){
  // https://www.tutorialspoint.com/cplusplus/cpp_multithreading.htm
  struct worker_args * td;
  td = ( worker_args * ) threadarg;
#if 0
  fprintf(stderr,"ID:%d THREAD:%d\n",td->id, td->thread_id);
#endif

  // init all in each thread
  int *keeplist = new int [td->freq->size()];
  for (size_t i = 0; i < td->freq->size(); i++) {
    // skip missing data
    if (access_genotype(td->gls, i, td->a, 0) == access_genotype(td->gls, i, td->a, 1) && 
        access_genotype(td->gls, i, td->a, 0) == access_genotype(td->gls, i, td->a, 2)){
      continue;
    }
    if (access_genotype(td->gls, i, td->b, 0) == access_genotype(td->gls, i, td->b, 1) && 
        access_genotype(td->gls, i, td->b, 0) == access_genotype(td->gls, i, td->b, 2)){
      continue;
    }

    // removing minor allele frequencies
    if (td->freq->at(i) < minMaf || (1 - td->freq->at(i)) < minMaf)
      continue;

    keeplist[td->nkeep] = i;
    td->nkeep++;
  }
  // fprintf(stderr, "\t-> keeping %d sites for downstream analyses", td->nkeep++);
  double **emis = new double *[td->nkeep];
  
  for (int i = 0; i < td->nkeep; i++) {
    emis[i] = new double[9];
  }

#if 0   // l1 + l2 -> gls not fixed yet
  if (gc) {
    if (gc > 1) {
      callgenotypesHwe(l1, td->nkeep, errate, newfreq);
      callgenotypesHwe(l2, td->nkeep, errate, newfreq);
    }

    if (gc > 0) {
      callgenotypesEps(l1, td->nkeep, errate);
      callgenotypesEps(l2, td->nkeep, errate);
    }
  }
#endif

  emission_ngsrelate9(td->freq, td->gls, emis, keeplist, td->nkeep, td->a, td->b);
  if (model == 0)
    td->niter = em1(td->pars, emis, td->nkeep);
  else if (model == 1)
    td->niter = em2(td->pars, emis, td->nkeep);
  else // below might not work
    td->niter = em3(td->pars, emis, td->nkeep);

  double l100000000 = loglike(p100000000, emis, td->nkeep);
  double l010000000 = loglike(p010000000, emis, td->nkeep);
  double l001000000 = loglike(p001000000, emis, td->nkeep);
  double l000100000 = loglike(p000100000, emis, td->nkeep);
  double l000010000 = loglike(p000010000, emis, td->nkeep);
  double l000001000 = loglike(p000001000, emis, td->nkeep);
  double l000000100 = loglike(p000000100, emis, td->nkeep);
  double l000000010 = loglike(p000000010, emis, td->nkeep);
  double l000000001 = loglike(p000000001, emis, td->nkeep);
  double lopt = loglike(td->pars, emis, td->nkeep);

  td->ll = lopt;

  double likes[10] = {l100000000, l010000000, l001000000, l000100000,
                      l000010000, l000001000, l000000100, l000000010,
                      l000000001, lopt};
  td->best = 0;
  td->bestll = likes[0];
  for (int i = 1; i < 10; i++) {
    if (likes[i] > likes[td->best]){
      td->best = i;
      td->bestll = likes[i];
    }
  }

  // double newpars[9] = {0.0625, 0.3125, 0.21875,
  //                      0.03125, 0.125, 0.03125,
  //                      0.125, 0.03125, 0.0625};

  // double newlopt = loglike(newpars, emis, td->nkeep);
  // fprintf(stdout, "%f is the likelihood of the correct estimates\n", newlopt);

  // end of work. Teardown
  for (int i = 0; i < td->nkeep; i++) {
    delete[] emis[i];
  }
  delete[] emis;
  pthread_exit(NULL);
}

void * do_work_inbred(void *threadarg){
  // https://www.tutorialspoint.com/cplusplus/cpp_multithreading.htm
  struct worker_args * td;
  td = ( worker_args * ) threadarg;
#if 0
  fprintf(stderr,"ID:%d THREAD:%d\n",td->id, td->thread_id);
#endif

  // init all in each thread
  int *keeplist = new int [td->freq->size()];
  for (size_t i = 0; i < td->freq->size(); i++) {
    // skip missing data
    if (access_genotype(td->gls, i, td->a, 0) == access_genotype(td->gls, i, td->a, 1) && 
        access_genotype(td->gls, i, td->a, 0) == access_genotype(td->gls, i, td->a, 2)){
      continue;
    }
    // removing minor allele frequencies
    if (td->freq->at(i) < minMaf)
      continue;

    keeplist[td->nkeep] = i;
    td->nkeep++;
  }
  // fprintf(stderr, "\t-> keeping %d sites for downstream analyses", td->nkeep++);
  double **emis = new double *[td->nkeep];
  for (int i = 0; i < td->nkeep; i++) {
    emis[i] = new double[2];
  }
  
#if 0   // l1 + l2 -> gls not fixed yet
  if (gc) {
    if (gc > 1) {
      callgenotypesHwe(l1, td->nkeep, errate, newfreq);
    }
    
    if (gc > 0) {
      callgenotypesEps(l1, td->nkeep, errate);
    }
  }
#endif
  emission_ngs_inbred(td->freq,td->gls, emis, keeplist, td->nkeep, td->a);
  td->niter=em_inbred(td->pars,emis,tole,maxIter,td->nkeep,model,verbose);
  td->ll = loglike_inbred(td->pars,emis,td->nkeep);
  double l01= loglike_inbred(p01,emis,td->nkeep);
  double l10= loglike_inbred(p10,emis,td->nkeep);
  double likes[3] ={l10,l01,td->ll};
  td->best = 0;
  td->bestll = likes[0];
  for(int i=1;i<3;i++){
    if(likes[i]>likes[td->best]){
      td->best=i;
      td->bestll = likes[i];
    }
  }
  for (int i = 0; i < td->nkeep; i++) {
    delete[] emis[i];
  }
  delete[] emis;
  pthread_exit(NULL);
}


int main(int argc, char **argv){
  if(argc==1)
    print_info(stderr);

  if(strcasecmp(argv[1],"extract_freq_bim")==0)
    return extract_freq_bim(--argc,++argv);
  if(strcasecmp(argv[1],"extract_freq")==0)
    return extract_freq(--argc,++argv);

  char *htsfile=NULL;
  while ((n = getopt(argc, argv, "f:i:t:r:g:m:v:s:F:c:e:a:b:n:l:z:p:h:")) >= 0) {
    switch (n) {
    case 'f': freqname = strdup(optarg); break;
    case 'i': maxIter = atoi(optarg); break;
    case 't': tole = atof(optarg); break;
    case 'r': seed = atoi(optarg); break;
    case 'g': gname = strdup(optarg); break;
    case 'm': model = atoi(optarg); break;
    case 'v': verbose = atoi(optarg); break;
    case 's': switchMaf = atoi(optarg); break;
    case 'F': do_inbred = atoi(optarg); break;
    case 'c': gc = atoi(optarg); break;
    case 'a': pair1 = atoi(optarg); break;
    case 'b': pair2 = atoi(optarg); break;
    case 'n': {nind = atoi(optarg); hasDef=1; break;}
    case 'p': num_threads = atoi(optarg);break;
    case 'e': errate = atof(optarg); break;
    case 'l': minMaf = atof(optarg); break;
    case 'h': htsfile = strdup(optarg); break;
    case 'z': readids(ids,optarg); break;
    default: {fprintf(stderr,"unknown arg:\n");return 0;}
      print_info(stderr);
    }
  }
#ifndef __WITH_BCF__
  if(htsfile){
    fprintf(stderr,"\t-> Must compile with -D__WITH_BCF__ for using htsfiles\n");
    return 0;
  }
#endif
    

  if (hasDef == 0&&htsfile==NULL) {
  fprintf(stderr, "\t-> -n parameter has not been supplied. Will assume that "
    "file contains 2 samples...\n");
  }
  if (ids.size() != 0 && ids.size() != nind) {
    fprintf(stderr, "\t-> Number of names doesnt match the -n parameter\n");
  }
  srand48(seed);
  if ((nind == -1 || freqname == NULL || gname == NULL)&&htsfile==NULL) {
    fprintf(stderr, "\t-> Must supply -n -f -g parameters (%d,%s,%s) OR -h file.bcf\n", nind,
            freqname, gname);
    return 0;
  }

  pthread_t threads[num_threads];

  std::vector<double> freq;
  double **gls=NULL;


  if(htsfile==NULL){
    getDouble(freqname,freq);
    gls = getGL(gname, freq.size(), nind);
  }
#ifdef __WITH_BCF__
  if(htsfile){
    std::vector<double *> tmpgl;
    nind=getgls(htsfile,tmpgl,freq,2,0.05);
    gls=new double *[tmpgl.size()];
    for(int i=0;i<tmpgl.size();i++)
      gls[i] = tmpgl[i];
  }
  fprintf(stderr,"\t\t-> NIND:%d\n",nind);
  for(int i=0;0&&i<freq.size();i++){
    fprintf(stdout,"%f",freq[i]);
    for(int ii=0;ii<3*nind;ii++)
      fprintf(stdout,"\t%f",gls[i][ii]);
    fprintf(stdout, "\n");
  }
  // exit(0);
#endif

  double total_sites = (1.0 * freq.size());
  if (switchMaf) {
    fprintf(stderr, "\t-> switching frequencies\n");
    for (size_t i = 0; i < freq.size(); i++)
      freq[i] = 1 - freq[i];
  }


#if 0
  print(stdout,freq.size(),3*nind,gls);
  exit(0);
#endif
  if(do_inbred){
    fprintf(stdout,"Pair\tZ=0\tZ=1\tloglh\tnIter\tcoverage\n");
    int comparison_ids_inbred = 0;
    std::vector<worker_args> all_args_inbred;
    int fake_person = -1;
    for(int a=0;a<nind;a++){
      if (pair1 != -1)
        a = pair1;
      worker_args td_args_inbred(a, fake_person, &freq, gls );        
      td_args_inbred.pars[0]=drand48();
      td_args_inbred.pars[1]=1-td_args_inbred.pars[0];
      all_args_inbred.push_back(td_args_inbred);
      comparison_ids_inbred++;
      if (pair1 != -1){
        //	fprintf(stderr,"BREAKING\n");
        break;
      }
    } // end ind b
    int cnt_inbred=0;
    while(cnt_inbred<comparison_ids_inbred){
      int nTimes_inbred;
      if(comparison_ids_inbred-cnt_inbred-num_threads>=0)
        nTimes_inbred = num_threads;
      else
        nTimes_inbred = comparison_ids_inbred-cnt_inbred;
#if 0
      for(int i=0;1 && i<nTimes_inbred;i++){
        fprintf(stderr,"cnt:%d i:%d\n",cnt_inbred+i,i);
      }
#endif
      for(int i=0;i<nTimes_inbred;i++){
        all_args_inbred[cnt_inbred + i].thread_id = i;
        pthread_create(&threads[i],NULL,do_work_inbred,&all_args_inbred[cnt_inbred+i]);
      }
      for(int i=0;i<nTimes_inbred;i++){
        pthread_join(threads[i], NULL);
        worker_args * td_out = &all_args_inbred[cnt_inbred + i];
        if(td_out->best==0)
          fprintf(stdout,"%f\t%f\t%f\t%d\t%f\n",p10[0],p10[1],td_out->bestll,-1,((double)td_out->nkeep)/((double)freq.size()));
        if(td_out->best==1)
          fprintf(stdout,"%f\t%f\t%f\t%d\t%f\n",p01[0],p01[1],td_out->bestll,-1,((double)td_out->nkeep)/((double)freq.size()));
        if(td_out->best==2)
          fprintf(stdout,"%f\t%f\t%f\t%d\t%f\n",td_out->pars[0],td_out->pars[1],td_out->bestll,td_out->niter,((double)td_out->nkeep)/((double)freq.size()));
      fflush(stdout);
      }
    }
  } else {
    if (ids.size()){
      fprintf(stdout,
              "a\tb\tida\tidb\tnSites\ts9\ts8\ts7\ts6\ts5\ts4\ts3\ts2\ts1\trab\tFa\tFb\ttheta\tloglh\tnIter\tcoverage\n");
    } else {
      fprintf(stdout, "a\tb\tnSites\ts9\ts8\ts7\ts6\ts5\ts4\ts3\ts2\ts1\trab\tFa\tFb\ttheta\tloglh\tnIter\tcoverage\n");
    }
    int comparison_ids = 0;
    std::vector<worker_args> all_args;
    for (int a = 0; a < nind; a++) {
      for (int b = a + 1; b < nind; b++) {
        if (pair1 != -1)
          a = pair1;
        if (pair2 != -1)
          b = pair2;

        worker_args td_args(a, b, &freq, gls );        
        
        double parval = 0.0, parsum = 0.0;
        for (int i = 0; i < 9; i++) {
          parval = drand48();
          td_args.pars[i] = parval;
          parsum += parval;
        }
        for (int i = 0; i < 9; i++) {
          td_args.pars[i] /= parsum;
          // fprintf(stderr, "par: %d, %f\n", i, td_args.pars[i]);
        }
        
        all_args.push_back(td_args);
        comparison_ids++;
        if (pair1 != -1 || pair2 != -1) {
          //	fprintf(stderr,"BREAKING\n");
          break;
        }
      } // end ind b
      if (pair1 != -1 || pair2 != -1) {
        //	fprintf(stderr,"BREAKING\n");
        break;
      }
    } // end ind b
    
    int cnt=0;
    while(cnt<comparison_ids){
      int nTimes;
      if(comparison_ids-cnt-num_threads>=0)
        nTimes = num_threads;
      else
        nTimes = comparison_ids-cnt;
#if 0
      for(int i=0;1 && i<nTimes;i++){
        fprintf(stderr,"cnt:%d i:%d\n",cnt+i,i);
      }
#endif
      for(int i=0;i<nTimes;i++){
        all_args[cnt + i].thread_id = i;
        pthread_create(&threads[i],NULL,do_work,&all_args[cnt+i]);
      }
      
      for(int i=0;i<nTimes;i++){
        pthread_join(threads[i], NULL);
        worker_args * td_out = &all_args[cnt + i];
        // printing results for threads[i]
        if (ids.size()) {
          fprintf(stdout, "%d\t%d\t%s\t%s\t%d", td_out->a,
                  td_out->b, ids[td_out->a],
                  ids[td_out->b], td_out->nkeep);
        } else {
          fprintf(stdout, "%d\t%d\t%d", td_out->a, td_out->b,
                  td_out->nkeep);
        }
        
        if (td_out->best == 9) {
          for (int j = 0; j < 9; j++) {
            fprintf(stdout, "\t%f", td_out->pars[j]);
          }
          double rab = (td_out->pars[8]+td_out->pars[2]) +
            (0.75 * td_out->pars[6]+td_out->pars[4]) +
            td_out->pars[1]*0.5;
          fprintf(stdout, "\t%f", rab);
          // Fa
          // sum(x[1]+x[2],x[3],x[4])
          double Fa = td_out->pars[8]+td_out->pars[7]+td_out->pars[6]+td_out->pars[5];
          fprintf(stdout, "\t%f", Fa);
          // Fb
          // sum(x[1],x[2],x[5],x[6])
          double Fb = td_out->pars[8]+td_out->pars[7]+td_out->pars[4]+td_out->pars[3];
          fprintf(stdout, "\t%f", Fb);
          // theta, coancestry / kinship coefficient
          double theta = td_out->pars[8]+0.5*(td_out->pars[6]+td_out->pars[4]+td_out->pars[2]) + 0.25 * td_out->pars[1];
          fprintf(stdout, "\t%f", theta);
          
          fprintf(stdout, "\t%f\t%d\t%f\n", td_out->ll,
                  td_out->niter,
                  (1.0 * td_out->nkeep) / total_sites);
        } else {
          for (int j = 0; j < 9; j++) {
            fprintf(stdout, "\t%f", td_out->pars[j]);
          }
          // 1 = 8;
          // 2 = 7;
          // 3 = 6;
          // 4 = 5;        
          // 5 = 4;
          // 6 = 3;
          // 7 = 2;
          // 8 = 1;
          // 9 = 0;
          // rab
          // whatisrxy <- function(x)
          // return(x[1]+x[7]+3/4*(x[3]+x[5])+x[8]*0.5)
          double rab = (td_out->pars[8]+td_out->pars[2]) +
            (0.75 * td_out->pars[6]+td_out->pars[4]) +
            td_out->pars[1]*0.5;
          fprintf(stdout, "\t%f", rab);
          // Fa
          // sum(x[1]+x[2],x[3],x[4])
          double Fa = td_out->pars[8]+td_out->pars[7]+td_out->pars[6]+td_out->pars[5];
          fprintf(stdout, "\t%f", Fa);
          // Fb
          // sum(x[1],x[2],x[5],x[6])
          double Fb = td_out->pars[8]+td_out->pars[7]+td_out->pars[4]+td_out->pars[3];
          fprintf(stdout, "\t%f", Fb);
          // theta, coancestry / kinship coefficient
          
          double theta = td_out->pars[8]+ 0.5*(td_out->pars[6]+td_out->pars[4]+td_out->pars[2]) +
            0.25 * td_out->pars[1];
          fprintf(stdout, "\t%f", theta);
          
          fprintf(stdout, "\t%f;s%d_%f\t%d\t%f\n", td_out->ll,
                  9 - td_out->best, td_out->bestll, -1,
                  (1.0 * td_out->nkeep) / total_sites);
          
        }
        
      }
      cnt += nTimes;
      fprintf(stderr, "\t-> Processed %d out of %d\r", cnt, comparison_ids);
    }
  }
  fflush(stdout);
  for (size_t i = 0; i < freq.size(); i++) {
    delete[] gls[i];
  }
  delete[] gls;
  free(freqname);
  free(gname);
  free(htsfile);
  return 0;
}
