/*
  NgsRelateV2
  http://www.popgen.dk/software/
  https://github.com/ANGSD/NgsRelate/
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
#include <time.h>    // time
#include <limits>
#include <string>
#include <algorithm> //shuffle
#include <unistd.h>
#include "filereaders.h"

#ifdef __WITH_BCF__
#include "vcf.h"
#endif


double total_sites;
int FINISHED=0;
std::vector<FILE *> spillfiles;
std::vector<char *> spillfilesnames;

typedef struct{
  int a;
  int b;//<-if inbreding a=b
  char *res;//<-results as char*
}mypair;

std::vector<mypair> mp;


bool myfunction (mypair i,mypair j) { return i.a==j.a?i.b<j.b:i.a<j.a;}

//this is as close to the bound we will allow
double TINY=1e-8;

double p10[2]={1-TINY,TINY};
double p01[2]={TINY,1-TINY};

int num_threads = 4;
char *freqname=NULL;
char *gname=NULL;
char *beaglefile=NULL;

int maxIter =10000;
double tole =1e-8;
int n=-1;
int nBootstrap = 0;//counter for howmany bootstraps
int seed=std::numeric_limits<int>::max();
int ntimes =1;
int model =1;
int gc =0;
double errate = 0.005;
int pair1 =-1;
int pair2 =-1;
int nind =-1;
int nsites_nofreqfile = 0;
size_t overall_number_of_sites = 0;
int do_2dsfs_only = 0;
int do_inbred=0;
int do_simple=0;
int switchMaf = 0;

double minMaf =0.05;
int hasDef = 0;
double ttol=1e-6;
char *vcf_format_field = strdup("PL"); // can take PL or GT
char *vcf_allele_field = strdup("AFngsrelate"); // can take any tag value e.g. AF AF1 etc

std::vector<char *> ids;

float emTole=1e-12;

// https://en.cppreference.com/w/c/numeric/math/isnan
bool is_nan(double x) { return x != x; }

double emFrequency_from_data(double *like,int numInds, int iter,double start){
     
  float W0;
  float W1;
  float W2;
  // fprintf(stderr,"start=%f\n",start);
  float p=(float)start;
  float temp_p=(float)start;
  double accu=0.00001;
  double accu2=0;
  float sum;

  
  int it=0;
  
  for(it=0;it<iter;it++){
    sum=0;
    for(int i=0;i<numInds;i++){
      W0=like[i*3+0]*pow(1-p,2);
      W1=like[i*3+1]*2*p*(1-p);
      W2=like[i*3+2]*(pow(p,2));
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      if(std::isnan(sum))
	fprintf(stderr,"PRE[%d]:gls:(%f,%f,%f) W(%f,%f,%f) sum=%f\n",i,like[i*3],like[i*3+1],like[i*3+2],W0,W1,W2,sum);
    }
    
    p=sum/numInds;
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
    temp_p=p;
  }
  
  if(std::isnan(p)){
    fprintf(stderr,"[%s] caught nan will not exit\n",__FUNCTION__);
    fprintf(stderr,"Like (3*nInd). nInd=%d\n",numInds);
    //print_array(stderr,like,3*numInds);
    fprintf(stderr,"keepList (nInd)\n");
    //print_array(stderr,keep,numInds);
    fprintf(stderr,"used Like (3*length(keep))=%d\n",numInds);
    
    for(int ii=0;1&&ii<numInds;ii++){
      for(int gg=0;gg<3;gg++)
        fprintf(stderr,"%f\t",like[ii*3+gg]);
      fprintf(stderr,"\n");
    }
    sum=0;
    for(int i=0;i<numInds;i++){
      W0=like[i*3+0]*pow(1-p,2);
      W1=like[i*3+1]*2*p*(1-p);
      W2=like[i*3+2]*(pow(p,2));
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      fprintf(stderr,"[%s.%s():%d] p=%f W %f\t%f\t%f sum=%f like: %f\n",__FILE__,__FUNCTION__,__LINE__,p,W0,W1,W2,sum,like[i*3+2]*pow(1-p,2));
      break;
    }
    p=-999;
    assert(p!=999);
    return p;
  }

  return(p);
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

void sample48(double *ary,int dim){
  for(int i=0;i<dim;i++)
    ary[i] =TINY+drand48()*(1-2*TINY);
  
  normalize(ary,dim);
}
double loglike(double *p,double **emis,int len,int dim){
  double ret =0;

  for(int i=0;i<len;i++){
    double tmp = 0;
    for(int j=0;j<dim;j++)
      tmp += p[j]*emis[i][j];
    ret += log(tmp);
  }
  return ret;
}

void emStep(double *pre,double **emis,double *post,int len,int len2){
  for(int i=0;i<len2;i++){
    if(pre[i]<0||pre[i]>1||is_nan(pre[i])){
      fprintf(stderr,"Problem with guess in emStep: ");
      for(int j=0;j<len2;j++)
	fprintf(stderr," %f ",pre[j]);
      fprintf(stderr,"\n");
      exit(0);
    }
  }
  double inner[len2];
  for(int x=0;x<len2;x++){
    post[x] =0.0;
  }
  for(int i=0;i<len;i++){
    for(int x=0;x<len2;x++){
      inner[x] = pre[x]*emis[i][x];
    }

    normalize(inner,len2);

#ifdef DB_MP
    for(int x=0;x<len2;x++){
      if(is_nan(inner[x])){
        fprintf(stderr, "pre: %d ", i);
        for(int x=0;x<len2;x++){
          fprintf(stderr, "%f ", pre[x]);
        }
        fprintf(stderr, "\n");
        
        fprintf(stderr, "emis: %d ", i);
        for(int x=0;x<len2;x++){
          fprintf(stderr, "%f ", emis[i][x]);
        }
        fprintf(stderr, "\n");
        
        fprintf(stderr, "site: %d ", i);
        for(int x=0;x<len2;x++){
          fprintf(stderr, "%f ", inner[x]);
        }
        fprintf(stderr, "\n");
        exit(0);
      }
    }
#endif
    
    
    for(int x=0;x<len2;x++){
      post[x] += inner[x];
#ifdef DB_MP
      if(is_nan(post[x])){
        fprintf(stderr, "%d %f\n", i, inner[x]);
        exit(0);
      }
#endif

    }
  }

#ifdef DB_MP
  fprintf(stderr, "post->sites: %d, params: %d ", len, len2);
  for(int i=0; i<len2; i++){
    fprintf(stderr, " %f", post[i]);
  }
  fprintf(stderr, "\n");
#endif
  
  normalize(post,len2);

}

void minus(double * fst,double * sec,double * res,int dim){
  for(int i=0;i<dim;i++)
    res[i] = fst[i]-sec[i];
}

double sumSquare(double * mat,int dim){
  double tmp=0;
  for(size_t i=0;i<dim;i++){
    tmp += mat[i]*mat[i];
  }
  return tmp;
}

int emAccel(double *F,double **emis,double *F_new,int len, int & niter,int dim){
  //  maybe these should be usersettable?
 
  double stepMin =1;
  double stepMax0 = 1;
  static double stepMax=stepMax0;
  double mstep=4;
  //  double objfnInc=1;
#ifdef DB_MP
  fprintf(stderr, "iter: %d input:  ", niter);
  for(int i=0; i<dim; i++){
    fprintf(stderr, " %f", F[i]);
  }
  fprintf(stderr, "\n");
#endif

  double F_em1[dim];
  double F_diff1[dim];
  double F_em2[dim];
  double F_diff2[dim];
  double F_diff3[dim];
  double F_tmp[dim];
  niter++;
  emStep(F, emis, F_em1, len,dim);
  // stayin(F_em1);
#ifdef DB_MP
  fprintf(stderr, "iter: %d emstep1:  ", niter);
  for(int i=0; i<dim; i++){
    fprintf(stderr, " %f", F_em1[i]);
  }
  fprintf(stderr, "\n");
#endif

  
  
  minus(F_em1, F, F_diff1,dim);

  double sr2 = sumSquare(F_diff1,dim);
  
  if(sqrt(sr2)<ttol){
#ifdef DB_MP
    fprintf(stderr,"sr2 break1: %e\n", sqrt(sr2));
#endif
    if (niter==2){
      // This is required as breaking in first round will result in empty F_new
      // copying F_em1 to F_new every time literally disables accelerated.
      for(int i=0;i<dim;i++)
        F_new[i]  = F_em1[i];
    }


    for(int i=0;0&&i<dim;i++)
      F_new[i]  = F_em1[i];
    return 1;
  }
  niter++;
  emStep(F_em1, emis, F_em2, len,dim);
  minus(F_em2, F_em1, F_diff2,dim);
#ifdef DB_MP
  fprintf(stderr, "iter: %d emstep2:  ", niter);
  for(int i=0; i<dim; i++){
    fprintf(stderr, " %f", F_em2[i]);
  }
  fprintf(stderr, "\n");
#endif

  double sq2 = sumSquare(F_diff2,dim);

  if(sqrt(sq2)<ttol){
#ifdef DB_MP
    fprintf(stderr,"sr2 break2: %e\n", sqrt(sr2));
#endif

    if (niter==3){
      // This is required as breaking in first round will result in empty F_new
      // copying F_em1 to F_new every time literally disables accelerated.      
      for(int i=0;i<dim;i++)
        F_new[i]  = F_em2[i];
    }
    

    for(int i=0;0&&i<dim;i++)
      F_new[i]  = F_em2[i];
    return 1;
  }

  minus(F_diff2,F_diff1, F_diff3,dim);

  double sv2 = sumSquare(F_diff3,dim);

  double alpha = sqrt(sr2/sv2);
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  for(size_t i=0;i<dim;i++){
    //  F_new[i] = std::min(1-ttol, std::max(ttol,F[i]+2*alpha*F_diff1[i]+alpha*alpha*F_diff3[i]));
    F_new[i] = F[i]+2*alpha*F_diff1[i]+alpha*alpha*F_diff3[i];    
  }

  int outofparspace =0;
  for(int i=0;i<dim;i++){
    if(F_new[i]<0||F_new[i]>1){
      outofparspace++;
      // break;
    }
  }
  if(outofparspace){
#ifdef DB_MP
  fprintf(stderr, "iter: %d outofspace:  ", niter);
  for(int i=0; i<dim; i++){
    fprintf(stderr, " %f", F[i]);
  }
  fprintf(stderr, "\n");
#endif
    for(int i=0;i<dim;i++)
      F_new[i] = F_em2[i];
  }

  if (fabs(alpha - 1) > 0.01){
    niter++;
    emStep(F_new,emis,F_tmp,len,dim);
    for(int i=0;i<dim;i++)
      std::swap(F_new[i],F_tmp[i]);
  }

  if ((alpha - stepMax) > -0.001) {
    stepMax = mstep*stepMax;
  }

  return 1;
}

//std em along with control logic for breaking
int em(double *sfs,double  **emis, int len, int dim){
  int niter = 0;
  double oldLik,lik;
  //  sample48(sfs,dim);//<- moved to outside 

  oldLik = loglike(sfs,emis,len,dim);

  double tmp[dim];

  int it;
  
  for(it=0;niter<maxIter;it++) {
    niter++;
    if(model==0)
      emStep(sfs,emis,tmp,len,dim);
    else
      emAccel(sfs,emis,tmp,len, niter,dim);      
    for(int i=0;i<dim;i++)
      sfs[i]= tmp[i];
    
    lik = loglike(sfs,emis,len,dim);

    if(fabs(lik-oldLik)<tole){

      oldLik=lik;
      break;
    }
    oldLik=lik;
  }
  return niter;
}

void emislike_2dsfs_gen(double **gls, double **emislike_2dsfs, int *keeplist, int & nkeep, int & ind1, int & ind2 ){
  int i;
  for(int x=0;x<nkeep;x++){
    i = keeplist[x];
    emislike_2dsfs[x][0] = access_genotype(gls, i, ind1, 0) * access_genotype(gls, i, ind2, 0);
    emislike_2dsfs[x][1] = access_genotype(gls, i, ind1, 0) * access_genotype(gls, i, ind2, 1);
    emislike_2dsfs[x][2] = access_genotype(gls, i, ind1, 0) * access_genotype(gls, i, ind2, 2);
    emislike_2dsfs[x][3] = access_genotype(gls, i, ind1, 1) * access_genotype(gls, i, ind2, 0);
    emislike_2dsfs[x][4] = access_genotype(gls, i, ind1, 1) * access_genotype(gls, i, ind2, 1);
    emislike_2dsfs[x][5] = access_genotype(gls, i, ind1, 1) * access_genotype(gls, i, ind2, 2);
    emislike_2dsfs[x][6] = access_genotype(gls, i, ind1, 2) * access_genotype(gls, i, ind2, 0);
    emislike_2dsfs[x][7] = access_genotype(gls, i, ind1, 2) * access_genotype(gls, i, ind2, 1);
    emislike_2dsfs[x][8] = access_genotype(gls, i, ind1, 2) * access_genotype(gls, i, ind2, 2);
#ifdef DB_EMIS
    fprintf(stderr,"emis2dsfs[%d]:\t",x);
    for(int j=0;j<9;j++)
      fprintf(stderr,"%f ",emislike_2dsfs[x][j]);
    fprintf(stderr,"\n");
#endif
  }
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

#ifdef myDEBUGemis9
    double test =0;
    for (int xx=0;xx<9;xx++){
      test+=emis[x][xx];
    }
    if(test<1e-30){
      fprintf(stderr, "%d %d g1:", i, x);
      for (int geno=0; geno<3; geno++)
        fprintf(stderr, " %f",  access_genotype(gls, i, ind1, geno));
      fprintf(stderr, " g2:");
      for (int geno=0; geno<3; geno++)
        fprintf(stderr, " %f",  access_genotype(gls, i, ind2, geno));
      fprintf(stderr,"\n");
    }
#endif
    
  }

  
}

void emission_ngs_inbred(std::vector<double> * freq, double **gls, double **emis, int *keeplist, int & nkeep, int & ind1){
  int i;
  double AA, Aa, aa, freqa, freqA;
  for(int x=0;x<nkeep;x++){
    i = keeplist[x];
    freqa=freq->at(i);  // alternative allele frequency
    freqA=1-freqa;
    AA = access_genotype(gls, i, ind1, 0);
    Aa = access_genotype(gls, i, ind1, 1);
    aa = access_genotype(gls, i, ind1, 2);
    // G_real=(AA)
    emis[x][0] = pow(freqA,2)*AA;
    emis[x][1] = freqA*AA;

    // G_real=(Aa)
    emis[x][0] += 2*freqa*freqA*Aa;
    emis[x][1] += 0;

    // G_real=(aa)
    emis[x][0] += pow(freqa,2)*aa;
    emis[x][1] += freqa*aa;
  }
}



void print_info(FILE *fp){
  fprintf(fp, "\n");
  fprintf(fp, "Usage main analyses: ./ngsrelate  [options] \n");
  fprintf(fp, "Options:\n");
  fprintf(fp, "   -f <filename>       Name of file with frequencies\n");
  fprintf(fp, "   -O <filename>       Output filename\n");
  fprintf(fp, "   -L <INT>            Number of genomic sites. Must be provided if -f (allele frequency file) is NOT provided \n");
  fprintf(fp, "   -m <INTEGER>        model 0=normalEM 1=acceleratedEM\n");
  fprintf(fp, "   -i <UINTEGER>       Maximum number of EM iterations\n");
  fprintf(fp, "   -t <FLOAT>          Tolerance for breaking EM\n");
  fprintf(fp, "   -r <FLOAT>          Seed for rand\n");
  fprintf(fp, "   -R <chr:from-to>    Region for analysis (only for bcf)\n");
  fprintf(fp, "   -g gfile            Name of glf (compressed binary) file\n");
  fprintf(fp, "   -G gfile            Name of beagle (compressed) file\n");  
  fprintf(fp, "   -p <INT>            threads (default 4)\n");
  fprintf(fp, "   -c <INT>            Should call genotypes instead?\n");
  fprintf(fp, "   -s <INT>            Should you swich the freq with 1-freq?\n");
  fprintf(fp, "   -F <INT>            Estimate inbreeding instead of estimating the nine jacquard coefficients\n");
  fprintf(fp, "   -o <INT>            estimating the 3 jacquard coefficient, assumming no inbreeding\n");
  fprintf(fp, "   -v <INT>            Verbose. print like per iteration\n");
  fprintf(fp, "   -e <FLOAT>          Errorrates when calling genotypes?\n");
  fprintf(fp, "   -a <INT>            First individual used for analysis? (zero offset)\n");
  fprintf(fp, "   -b <INT>            Second individual used for analysis? (zero offset)\n");
  fprintf(fp, "   -B <INT>            Number of bootstrap replicates for (only for single pairs)\n");
  fprintf(fp, "   -N <INT>            How many times to start each pair with random seed?\n");
  fprintf(fp, "   -n <INT>            Number of samples in glf.gz\n");
  fprintf(fp, "   -l <FLOAT>          minMaf or 1-Maf filter (default: 0.05)\n");  
  fprintf(fp, "   -z <INT>            Name of file with IDs (optional)\n");
  fprintf(fp, "   -T <STRING>         For -h vcf use PL (default) or GT tag\n");
  fprintf(fp, "   -A <STRING>         For -h vcf use allele frequency TAG e.g. AFngsrelate (default)\n");  
  fprintf(fp, "   -P <filename>       plink name of the binary plink file (excluding the .bed)\n");
  fprintf(fp, "\n");
  fprintf(fp,"Or\n ./ngsrelate extract_freq_bim pos.glf.gz plink.bim plink.freq\n");
  fprintf(fp,"Or\n ./ngsrelate extract_freq .mafs.gz .pos.glf.gz [-rmTrans]\n");
#ifdef __WITH_BCF__
  fprintf(fp,"Or\n ./ngsrelate -h my.bcf [DEVELOPMENT ONLY]\n");
#endif
  exit(0);
}

int is_missing(double *ary){
  if(fabs(ary[0] - ary[1])<1e-6 && fabs(ary[0] - ary[2])<1e-6 && fabs(ary[1] - ary[2])<1e-6)
    return 1;
  else if(is_nan(ary[0]) || is_nan(ary[1]) || is_nan(ary[2]))
    return 1;
  else
    return 0;
}

int is_missing_old(double *ary){
  if(ary[0]==ary[1]&&ary[0]==ary[2]&&ary[1]==ary[2])
    return 1;
  else
    return 0;
}


void callgenotypesEps(double **gls,int nsites,int nind,double eps){

  double g00 = (1-eps)*(1-eps);
  double g01 = 2*(1-eps)*eps;
  double g02 = eps*eps;
  double g10 = (1-eps)*eps;
  double g11 = (1-eps)*(1-eps)+eps*eps;

  for(int off=0;off<nind;off++){
    for(int s=0;s<nsites;s++){
      if(!is_missing(&gls[s][off*3])){
	int whmax=0;
	for(int i=1;i<3;i++)
	  if(gls[s][off*3+i]>gls[s][off*3+whmax])
	    whmax=i;
	
	if(whmax==0){
	  gls[s][off*3+0] = g00;
	  gls[s][off*3+1] = g01;
	  gls[s][off*3+2] = g02;
	}else if(whmax==1){
	  gls[s][off*3+0]=gls[s][off*3+2]=g10;
	  gls[s][off*3+1]=g11;
	}else if(whmax==2){
	  gls[s][off*3+0] = g02;
	  gls[s][off*3+1] = g01;
	  gls[s][off*3+2] = g00;
	}else{
	  assert(0!=1);
	}
      }else
	gls[s][off*3+0]=gls[s][off*3+1]=gls[s][off*3+2]=1;//how should we treat missing in the case of called genotypes?
    }
  }
}

void callgenotypesHwe(double **gls,int nsites,int nind,std::vector<double> freq){

  for(int off=0;off<nind;off++){
    for(int s=0;s<nsites;s++){
      gls[s][3*off+0] *= freq[s]*freq[s];
      gls[s][3*off+1] *= (1-freq[s])*freq[s];
      gls[s][3*off+2] *= (1-freq[s])*(1-freq[s]);
      
      int whmax=0;
      for(int i=1;i<3;i++)
	if(gls[s][3*off+i]>gls[s][3*off+whmax])
	  whmax=i;

      for(int i=0;i<3;i++)
	gls[s][3*off+i] = 0;
      gls[s][3*off+whmax]=1;
    }
  }
}

struct worker_args_t {
  int thread_id, nkeep, a,b, niter, best, niter_2dsfs;
  double **gls;
  std::vector<double> * freq;
  size_t nsites;
  double ll, bestll, ll_2dsfs;
  double pars[9], pars_2dsfs[9];
  int *keeplist;
  double **emis;
  int *bootindex;
  worker_args_t(int & id_a, int & id_b, std::vector<double> * f, double ** gls_arg, size_t & s ){
    a=id_a;
    b=id_b;
    nkeep=0;
    best=0;
    bestll=0.0;
    freq = f;
    gls=gls_arg;
    nsites = s;
    bootindex =NULL;
    keeplist = new int[nsites];
    emis = new double*[s];
    if(nBootstrap>0)
      bootindex = new int[nsites];
    for(int i=0;i<s;i++)
      emis[i] = new double[9];//<- will hold, 2dsfs,F and 9jacq emissions depending on context
  }
  ~worker_args_t(){
    for(int i=0;i<nsites;i++)
      delete [] emis[i];
    delete [] emis;
    delete [] keeplist;
    delete [] bootindex;
  }
};

typedef struct worker_args_t worker_args;

// function will update pk_keeplist and return the number of sites that should be retained for analysis
int populate_keeplist(int pk_a,int pk_b,int pk_nsites,double **pk_gls,int pk_minmaf,std::vector<double> *pk_freq,int *pk_keeplist,int *tmp){//last is workspace to avoid allocation and deallocation;

  int *keeplist = pk_keeplist;
  int nkeep=0;
  for (size_t i = 0; i < pk_nsites; i++) {

    if(is_missing(&pk_gls[i][3*pk_a]))
      continue;
    if(is_missing(&pk_gls[i][3*pk_b]))
      continue;

    

    // removing minor allele frequencies
    // if ( (!do_2dsfs_only) &&
    //      ((*pk_freq)[i] < minMaf || (1 - (*pk_freq)[i] < minMaf)))
    //   continue;
    if ((*pk_freq)[i] < minMaf || (1 - (*pk_freq)[i] < minMaf))
      continue;

#ifdef DB_GL
    for (int x=0; x<3;x++)
      fprintf(stderr, "%lu %f %f\n", i, pk_gls[i][3*pk_a+x], pk_gls[i][3*pk_b+x]);
#endif
    keeplist[nkeep] = i;//dont forget
    nkeep++;
  }
  
  if(tmp){//indicater for if we should bootstrap
    for(int i=0;i<nkeep;i++){
      tmp[i] = keeplist[lrand48() % nkeep];
      //      fprintf(stderr,"tmp:%d\n",tmp[i]);
    }
    std::sort(tmp,tmp+nkeep);
    for(int i=0;i<nkeep;i++)
      keeplist[i] = tmp[i];
  }
  
  
  return nkeep;
}

//this one does both inbreeding and the j9, the j3 is obtained by setting the j1-6 to zero
int analyse_jaq(double *pk_pars,std::vector<double> *pk_freq,double **pk_gls,int *pk_keeplist,double **pk_emis,int pk_nkeep,int pk_a,int pk_b,double &pk_ll,int &pk_best,double &pk_bestll,int &pk_niter, int ntimes){
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

  if(do_inbred==0)
    emission_ngsrelate9(pk_freq, pk_gls, pk_emis, pk_keeplist, pk_nkeep, pk_a, pk_b);
  else
    emission_ngs_inbred(pk_freq,  pk_gls, pk_emis, pk_keeplist, pk_nkeep, pk_a);

#ifdef DB_EMIS
  FILE * emis_fp = fopen("db_emis.txt", "w");
    for(int i=0;i<pk_nkeep;i++){
      fprintf(emis_fp, "%f", pk_emis[i][0]);
      for(int x=1;x<9;x++){
        fprintf(emis_fp, " %f", pk_emis[i][x]);
      }
      fprintf(emis_fp, "\n");      
    }
    fclose(emis_fp);
#endif
  //will not be used before values has been plugged in;
  double tmp_pk_pars[9];
  double tmp_pk_ll;
  int tmp_pk_niter;

  pk_ll=log(0);

#ifdef DB_MP
  fprintf(stderr, "jacquard\n");
#endif
  
  for(int n=0;n<ntimes;n++) {
    sample48(tmp_pk_pars,do_inbred?2:9);
    
    for(int i=3;do_simple&&do_inbred==0 && i<9;i++)//setting it to the old
      tmp_pk_pars[i] = 0;
    tmp_pk_niter = em(tmp_pk_pars, pk_emis, pk_nkeep,do_inbred?2:9);

    tmp_pk_ll = loglike(tmp_pk_pars, pk_emis, pk_nkeep,do_inbred?2:9);
    
    if(n==0||tmp_pk_ll>pk_ll){
      if(0&&n>0)
	fprintf(stderr,"n:%d better estimate: %f %f\n",n,tmp_pk_ll,pk_ll);
      for(int i=0;i<(do_inbred?2:9);i++){
	pk_pars[i] = tmp_pk_pars[i];
      }
      pk_ll=tmp_pk_ll;
      pk_niter=tmp_pk_niter;
    }
    
  }
  
  if(do_inbred==0){
    double l100000000 = loglike(p100000000, pk_emis, pk_nkeep,9);//0
    double l010000000 = loglike(p010000000, pk_emis, pk_nkeep,9);
    double l001000000 = loglike(p001000000, pk_emis, pk_nkeep,9);
    double l000100000 = loglike(p000100000, pk_emis, pk_nkeep,9);
    double l000010000 = loglike(p000010000, pk_emis, pk_nkeep,9);
    double l000001000 = loglike(p000001000, pk_emis, pk_nkeep,9);//5
    double l000000100 = loglike(p000000100, pk_emis, pk_nkeep,9);//6
    double l000000010 = loglike(p000000010, pk_emis, pk_nkeep,9);//7
    double l000000001 = loglike(p000000001, pk_emis, pk_nkeep,9);//8
    double likes[10] = {l100000000, l010000000, l001000000, l000100000,
                        l000010000, l000001000, l000000100, l000000010,
                        l000000001, pk_ll};
    pk_best = 0;
    pk_bestll = likes[0];
    int stop= do_simple==0?10:3;
    for (int i = 1; i < stop; i++) {
      if (likes[i] > likes[pk_best]){
        pk_best = i;
        pk_bestll = likes[i];
      }
    }
    if(stop==3&&likes[9]>likes[pk_best]){
      pk_best=9;
      pk_bestll = likes[pk_best];
    }

  }else{
    double l01= loglike(p01,pk_emis,pk_nkeep,2);
    double l10= loglike(p10,pk_emis,pk_nkeep,2);
    double likes[3] ={l10,l01,pk_ll};
    pk_best = 0;
    pk_bestll = likes[0];
    for(int i=1;i<3;i++){
      if(likes[i]>likes[pk_best]){
	pk_best=i;
	pk_bestll = likes[i];
      }
    }
  }
  return 0;
}


void anal1(int a,int b,worker_args * td,double minMaf, bool & nosites){
  
  assert(td->nsites>0);
  td->nkeep = populate_keeplist(a,b,td->nsites,td->gls,minMaf,td->freq,td->keeplist,td->bootindex);
  
  if (td->nkeep==0){
    fprintf(stderr, "\t-> sample index %d and %d have no overlapping sites with data. Pair will not be analyzed\n", a, b);
    nosites=true;
    return ;
  }
  
  //if(!do_2dsfs_only)
  analyse_jaq(td->pars,td->freq,td->gls,td->keeplist,td->emis,td->nkeep,a,b,td->ll,td->best,td->bestll,td->niter,ntimes);    


  if(do_inbred==0){

    emislike_2dsfs_gen(td->gls, td->emis,td->keeplist, td->nkeep, a, b);


    double tmp_pars_2dsfs[9];
    int tmp_niter_2dsfs;
    double tmp_ll_2dsfs;
#ifdef DB_MP
  fprintf(stderr, "sfs\n");
#endif
   
    td->ll_2dsfs=log(0);
    for(int n=0;n<ntimes;n++){
      sample48(tmp_pars_2dsfs,9);
      tmp_niter_2dsfs = em(tmp_pars_2dsfs, td->emis, td->nkeep,9);
      tmp_ll_2dsfs = loglike(tmp_pars_2dsfs, td->emis, td->nkeep,9);
      if(n==0||tmp_ll_2dsfs>td->ll_2dsfs){
	if(0&&n>0)
	  fprintf(stderr,"[sfs] n:%d better estimate: %f %f diff:%.3e\n",n,tmp_ll_2dsfs,td->ll_2dsfs,tmp_ll_2dsfs-td->ll_2dsfs);
	for(int i=0;i<9;i++){
	  td->pars_2dsfs[i] = tmp_pars_2dsfs[i];
	}
	td->ll_2dsfs=tmp_ll_2dsfs;
	td->niter_2dsfs = tmp_niter_2dsfs;
      }
    }
  }
}

char *formatoutputnosites(int a, int b){
  char retbuf[4096];
  if (ids.size()) {
    snprintf(retbuf,4096, "%d\t%d\t%s\t%s\t%d", a,
	     b, ids[a], ids[b], -1);
  } else {
    snprintf(retbuf,4096, "%d\t%d\t%d", a, b, -1);
  }

  for(int val=0; val<30; val++){
    snprintf(retbuf+strlen(retbuf), 4096, "\t%d", -1);
  }
  
  snprintf(retbuf+strlen(retbuf), 4096, "\n");
  
  return strdup(retbuf);  
}

char *formatoutput(int a, int b,worker_args *td_out,double total_sites){
  char retbuf[4096];
  if (ids.size()) {
    snprintf(retbuf,4096, "%d\t%d\t%s\t%s\t%d", a,
	     b, ids[a],
	     ids[b], td_out->nkeep);
  } else {
    snprintf(retbuf,4096, "%d\t%d\t%d", a, b,
	     td_out->nkeep);
  }

 /////////////////////////
 // Jacquard = ArrayPos //
 //        1 = 8;       //
 //        2 = 7;       //
 //        3 = 6;       //
 //        4 = 5;       //
 //        5 = 4;       //
 //        6 = 3;       //
 //        7 = 2;       //
 //        8 = 1;       //
 //        9 = 0;       //
 /////////////////////////
 
 if (td_out->best == 9) {
   for (int j = 0; j < 9; j++) {
     snprintf(retbuf+strlen(retbuf),4096, "\t%f", td_out->pars[j]);
   }
 } else {
   for (int j = 0; j < 9; j++){
     if (j==td_out->best){
       snprintf(retbuf+strlen(retbuf),4096, "\t%f", 1 - TINY);
     } else {
       snprintf(retbuf+strlen(retbuf),4096, "\t%f", TINY / 8.0);
     }
   }
 }

 //////////////////////////////////////////////////////////////////////
 // Measuring Relatedness between Inbred Individuals by Hedrick 2015 //
 // return(x[1]+x[7]+3/4*(x[3]+x[5])+x[8]*0.5)                       //
 // r_xy                                                             //
 //////////////////////////////////////////////////////////////////////
 double rab = (td_out->pars[8] + td_out->pars[2]) +
   (0.75 * (td_out->pars[6] + td_out->pars[4])) +
   td_out->pars[1] * 0.5;
 snprintf(retbuf+strlen(retbuf),4096, "\t%f", rab);

 /////////////////////////////////////////////////
 // Fa - inbreeding coefficient of individual a //
 // sum(x[1]+x[2],x[3],x[4])                    //
 /////////////////////////////////////////////////
 double Fa = td_out->pars[8] + td_out->pars[7] + td_out->pars[6] +
   td_out->pars[5];
 snprintf(retbuf+strlen(retbuf),4096, "\t%f", Fa);
 
 /////////////////////////////////////////////////
 // Fb - inbreeding coefficient of individual b //
 // sum(x[1],x[2],x[5],x[6])                    //
 /////////////////////////////////////////////////
 double Fb = td_out->pars[8] + td_out->pars[7] + td_out->pars[4] +
   td_out->pars[3];
 snprintf(retbuf+strlen(retbuf),4096, "\t%f", Fb);
 
 /////////////////////////////////////////////////////
 // theta / coancestry / kinship coefficient / f_XY //
 /////////////////////////////////////////////////////
 
 double theta =
   td_out->pars[8] +
   0.5 * (td_out->pars[6] + td_out->pars[4] + td_out->pars[2]) +
   0.25 * td_out->pars[1];
 snprintf(retbuf+strlen(retbuf),4096, "\t%f", theta);

 /////////////////////////////////////////////
 // summary statistics from ackerman et al. //
 /////////////////////////////////////////////
 double inbred_relatedness_1_2 = td_out->pars[8] + (td_out->pars[6] / 2.0);
 double inbred_relatedness_2_1 = td_out->pars[8] + (td_out->pars[4] / 2.0);
 double fraternity = td_out->pars[7] + td_out->pars[2];
 double identity = td_out->pars[8];
 double zygosity = td_out->pars[8] + td_out->pars[7] + td_out->pars[2];
 snprintf(retbuf+strlen(retbuf),4096, "\t%f\t%f\t%f\t%f\t%f", inbred_relatedness_1_2,
	 inbred_relatedness_2_1, fraternity, identity, zygosity);
 
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // summary statistics from Non-identifiability of identity coefficients at biallelic loci by miklos csuros //
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 // at least one pair of IBD alleles among three randomly selected ones:
 double eq_11e = td_out->pars[8] + td_out->pars[7] +
   td_out->pars[6] + td_out->pars[4] + td_out->pars[2] + 0.5 * (td_out->pars[5] + td_out->pars[3] + td_out->pars[1]);
 // 
 double eq_11f = 0.5 * (td_out->pars[5] - td_out->pars[3]);
 snprintf(retbuf+strlen(retbuf),4096, "\t%f\t%f", eq_11e, eq_11f);

 ///////////////////////////////
 // optimization of EM output //
 ///////////////////////////////
 
 if (td_out->best == 9) {
   snprintf(retbuf+strlen(retbuf),4096, "\t%f\t%d\t%d\t%f", td_out->ll, td_out->niter, -1,
	   (1.0 * td_out->nkeep) / total_sites);
 } else {
   snprintf(retbuf+strlen(retbuf),4096, "\t%f\t%d\t%f\t%f", td_out->bestll, td_out->niter, td_out->ll,
            (1.0 * td_out->nkeep) / total_sites);
 }

 
 ////////////////////
 // printing 2dsfs //
 ////////////////////
 snprintf(retbuf+strlen(retbuf),4096, "\t%e", td_out->pars_2dsfs[0]);
 for (int j = 1; j < 9; j++) {
   snprintf(retbuf+strlen(retbuf),4096, ",%e", td_out->pars_2dsfs[j]);
 }

 //////////////////////////////////////
 // computing 2dsfs based statistics //
 //////////////////////////////////////
 double r0, r1, king;
 r0 = (td_out->pars_2dsfs[2] + td_out->pars_2dsfs[6]) /
   td_out->pars_2dsfs[4];
 r1 = td_out->pars_2dsfs[4] /
   (td_out->pars_2dsfs[1] + td_out->pars_2dsfs[2] +
    td_out->pars_2dsfs[3] + td_out->pars_2dsfs[5] +
    td_out->pars_2dsfs[6] + td_out->pars_2dsfs[7]);
 king = (td_out->pars_2dsfs[4] -
	 (2 * (td_out->pars_2dsfs[2] + td_out->pars_2dsfs[6]))) /
   (2 * td_out->pars_2dsfs[4] + td_out->pars_2dsfs[1] +
    td_out->pars_2dsfs[3] + td_out->pars_2dsfs[5] +
    td_out->pars_2dsfs[7]);
 
 /////////////////////////////////////
 // printing 2dsfs based statistics //
 /////////////////////////////////////
 
 snprintf(retbuf+strlen(retbuf),4096, "\t%f\t%f\t%f", r0, r1, king);

 snprintf(retbuf+strlen(retbuf),4096, "\t%f\t%d\n", td_out->ll_2dsfs, td_out->niter_2dsfs);
 // fprintf(stderr,"retbuf:%s\n",retbuf);
 return strdup(retbuf);
}


/// go here
void * turbothread(void *threadarg){
  worker_args * td;
  td = ( worker_args * ) threadarg;

  for(int i=td->a;i<td->b;i++){
    bool nosites=false;
    anal1(mp[i].a,mp[i].b, td, minMaf, nosites);
    //collate results
    char buf[4096];
    if (do_inbred){
      if(td->best==0)
	snprintf(buf,4096,"%d\t%f\t%f\t%f\t%d\t%f\t%d\n",mp[i].a, p10[0],p10[1],td->bestll,-1,((double)td->nkeep)/((double)overall_number_of_sites), td->nkeep);
      if(td->best==1)
	snprintf(buf,4096,"%d\t%f\t%f\t%f\t%d\t%f\t%d\n",mp[i].a,p01[0],p01[1],td->bestll,-1,((double)td->nkeep)/((double)overall_number_of_sites), td->nkeep);
      if(td->best==2)
	snprintf(buf,4096,"%d\t%f\t%f\t%f\t%d\t%f\t%d\n",mp[i].a,td->pars[0],td->pars[1],td->bestll,td->niter,((double)td->nkeep)/((double)overall_number_of_sites), td->nkeep);
      mp[i].res=strdup(buf);
    } else if (nosites){
      mp[i].res = formatoutputnosites(mp[i].a, mp[i].b);
    } else{
      mp[i].res =formatoutput(mp[i].a,mp[i].b,td,total_sites);
    }
    fwrite(mp[i].res,sizeof(char),strlen(mp[i].res),spillfiles[td->thread_id]);
    fflush(spillfiles[td->thread_id]);
  }
  pthread_exit(NULL);
}

int main_analysis2(std::vector<double> &freq,double **gls,int num_threads,FILE *output,int total_sites){

  //initalize jobids
  if(do_inbred==0){
    for(int i=0;i<nind;i++){
      if(pair1!=-1)
	i=pair1;
      for(int j=(i+1);j<nind;j++){
	if(pair2!=-1)
	  j=pair2;
	mypair tmp;
	tmp.a=i;tmp.b=j;
	mp.push_back(tmp);
	if(pair2!=-1)
	  break;
      }
      if(pair1!=-1)
	break;
    }
  }else{
    for(int i=0;i<nind;i++){
      mypair tmp;
      tmp.a=tmp.b=i;
      mp.push_back(tmp);
    }
  }
  if(nBootstrap>0 &&mp.size()>1){
    fprintf(stderr,"\t-> You need to specify -a and -b for doing bootstrap\n");
    return 0;
  }
  if(nBootstrap>0){
    for(int i=0;i<nBootstrap;i++)
      mp.push_back(mp[0]);

  }

#ifdef DB_MP
  fprintf(stderr, "combinations: ");
  for(auto &val: mp){
    fprintf(stderr, "%d-%d;", val.a, val.b);
  }
  fprintf(stderr, "\n");
#endif
  std::random_shuffle(mp.begin(),mp.end());
  fprintf(stderr,"\t-> length of joblist:%lu\n",mp.size());

  //initialize threads ids
  pthread_t threads[num_threads];
  worker_args **all_args = new worker_args*[num_threads];
  int block=mp.size()/num_threads;

  for(int i=0;i<num_threads;i++){
    int first = i==0?0:all_args[i-1]->b;
    int second = first+block;
    all_args[i] = new worker_args(first, second, &freq, gls, overall_number_of_sites);
      all_args[i]->thread_id=i;
  }
  all_args[num_threads-1]->b = mp.size();

   for(int i=0;i<num_threads;i++)
     assert(pthread_create(&threads[i],NULL,turbothread,all_args[i])==0);

   for(int i=0;i<num_threads;i++)
     assert(pthread_join(threads[i], NULL)==0);
   
   if(do_inbred)
     fprintf(output,"Ind\tZ=0\tZ=1\tloglh\tnIter\tcoverage\tsites\n");
   else{
     if (ids.size()) {
       fprintf(output,
	       "a\tb\tida\tidb\tnSites\tJ9\tJ8\tJ7\tJ6\tJ5\tJ4\tJ3\tJ2\tJ1\trab\tFa\tFb\ttheta\tinbred_relatedness_1_2\tinbred_relatedness_2_1\tfraternity\tidentity\tzygosity\t2of3_IDB\tF_diff_a_b\tloglh\tnIter\tbestoptimll\tcoverage\t2dsfs\tR0\tR1\tKING\t2dsfs_loglike\t2dsfsf_niter\n");
     } else {
       fprintf(output, "a\tb\tnSites\tJ9\tJ8\tJ7\tJ6\tJ5\tJ4\tJ3\tJ2\tJ1\trab\tFa\tFb\ttheta\tinbred_relatedness_1_2\tinbred_relatedness_2_1\tfraternity\tidentity\tzygosity\t2of3_IDB\tFDiff\tloglh\tnIter\tbestoptimll\tcoverage\t2dsfs\tR0\tR1\tKING\t2dsfs_loglike\t2dsfsf_niter\n");
     }
   }
   
   std::sort(mp.begin(),mp.end(),myfunction);


   for(int i=0;i<mp.size();i++){
     fwrite(mp[i].res,sizeof(char),strlen(mp[i].res),output);
     free(mp[i].res);
   }
   FINISHED=1;
   for(int i=0;i<num_threads;i++)
     delete all_args[i];
   delete [] all_args;
  return 0;


  
    if(do_simple){
      fprintf(stderr, "\t-> setting Jacquard coefficient related to inbreeding (1-6) to zero\n");
    }
     
}

int my_atoi(char *s){
  // https://stackoverflow.com/a/3850567
  char *endptr = s;
  int value = (int)strtol(s, &endptr, 10);
  if(*endptr == *s ){
    fprintf(stderr, "pair '%s' is not a 0-based index\n", s);
    exit(0);
  }
  return value;
}

int main(int argc, char **argv){
   clock_t t=clock();
   time_t t2=time(NULL);


  if(argc==1)
    print_info(stderr);

  if(strcasecmp(argv[1],"extract_freq_bim")==0)
    return extract_freq_bim(--argc,++argv);
  if(strcasecmp(argv[1],"extract_freq")==0)
    return extract_freq(--argc,++argv);
#ifdef MP_DB
  fprintf(stdout,"#");
  for(int i=0;i<argc;i++)
    fprintf(stdout," %s",argv[i]);
  fprintf(stdout,"\n");
#endif
  
  char *htsfile=NULL;
  char *plinkfile=NULL;
  char *outname=NULL;
  char *region=NULL;
  
  while ((n = getopt(argc, argv, "f:i:t:r:g:G:m:s:F:o:c:e:a:b:n:l:z:p:h:L:T:A:P:O:X:R:B:N:")) >= 0) {
    switch (n) {
    case 'f': freqname = strdup(optarg); break;
    case 'P': plinkfile = strdup(optarg); break;
    case 'O': outname = strdup(optarg); break;
    case 'R': region = strdup(optarg); break;
    case 'i': maxIter = atoi(optarg); break;
    case 'N': ntimes = atoi(optarg); break;
    case 't': tole = atof(optarg); break;
    case 'r': seed = atoi(optarg); break;
    case 'g': gname = strdup(optarg); break;
    case 'G': beaglefile = strdup(optarg); break;      
    case 'm': model = atoi(optarg); break;
    case 'B': nBootstrap = atoi(optarg); break;
    case 's': switchMaf = atoi(optarg); break;
    case 'F': do_inbred = atoi(optarg); break;
    case 'o': do_simple = atoi(optarg); break;
    case 'c': gc = atoi(optarg); break;
    case 'a': pair1 = my_atoi(optarg); break;
    case 'b': pair2 = my_atoi(optarg); break;
    // case 'a': pair1 = atoi(optarg); break;
    // case 'b': pair2 = atoi(optarg); break;
    case 'n': {nind = atoi(optarg); hasDef=1; break;}
    case 'p': num_threads = atoi(optarg);break;
    case 'e': errate = atof(optarg); break;
    case 'l': minMaf = atof(optarg); break;
    case 'h': htsfile = strdup(optarg); break;
    case 'T': free(vcf_format_field);vcf_format_field = strdup(optarg); break;
    case 'A': free(vcf_allele_field);vcf_allele_field = strdup(optarg); break;            
    case 'z': readids(ids,optarg); break;
    case 'L': nsites_nofreqfile = atoi(optarg); break;
    default: {fprintf(stderr,"unknown arg:\n");return 0;}
      print_info(stderr);
    }
  }
  if(outname==NULL){
    fprintf(stderr,"\t-> No output filename declared please supply filename with -O output.res\n");
    return 0;
  }
  std::string plink_fam,plink_bim,plink_bed;
  if(plinkfile){
    fprintf(stderr,"\t-P %s\n",plinkfile);
    std::string p_str =std::string(plinkfile);
    if(p_str.length()>4){
      std::string ext = p_str.substr(p_str.length()-4,p_str.length());
      if (!ext.compare(".bed")||!ext.compare(".bim")||!ext.compare(".fam")){
	std::string front = p_str.substr(0,p_str.length()-4);
	plink_bim = (front+".bim");
	plink_fam = (front+".fam");
	plink_bed = (front+".bed");
      }else{
	plink_bim = (p_str+".bim");
	plink_fam = (p_str+".fam");
	plink_bed = (p_str+".bed");	
      }}else{
      plink_bim = (p_str+".bim");
      plink_fam = (p_str+".fam");
      plink_bed = (p_str+".bed");
    }
    fprintf(stderr,"\t NB make sure plink file only contains autosomes");
    fprintf(stderr,"\t-P %s -> bed:\'%s\' fam:\'%s\' bim:\'%s\'\n",plinkfile,plink_bed.c_str(),plink_fam.c_str(),plink_bim.c_str());

  }
#ifndef __WITH_BCF__
  if(htsfile){
    fprintf(stderr,"\t-> Must compile with -D__WITH_BCF__ for using htsfiles\n");
    return 0;
  }
#endif


  std::vector<double> freq;
  double **gls=NULL;


  
  if (hasDef == 0&&htsfile==NULL &&plinkfile==NULL) {
  fprintf(stderr, "\t-> -n parameter has not been supplied. Will assume that "
    "file contains 2 samples...\n");
  }
  if (ids.size() != 0 && ids.size() != nind) {
    fprintf(stderr, "\t-> Number of names doesnt match the -n parameter\n");
  }

  srand(time(NULL));
  if (seed ==  std::numeric_limits<int>::max()){
    seed=rand();
  }
  fprintf(stderr,"\t-> Seed is: %d\n",seed);
  srand48(seed);
  srand(seed);
  if (htsfile!=NULL){
    fprintf(stderr, "\t-> Will use TAG: '%s' from the VCF file\n", vcf_format_field);
    fprintf(stderr, "\t-> Will use TAG: '%s' in the VCF file as allele frequency if present. Otherwise allele frequencies are estimated from the data\n", vcf_allele_field);
  }

  if (gname != NULL && beaglefile != NULL){
    fprintf(stderr, "cannot provide both a beagle ( -G %s) and glf ( -g %s) file\n", beaglefile, gname);
    return 0;
  }
  

  if ((nind == -1 || (gname == NULL && beaglefile == NULL) )&&htsfile==NULL&&plinkfile==NULL) {
    fprintf(stderr, "\t-> Must supply: \n\t-> -n -g parameters (%d,%s) \n\t-> -n -G parameters (%d,%s) \n\t-> -h file.[vb]cf\n", nind,gname, nind, beaglefile);
    return 0;
  }

  if( nsites_nofreqfile==0 && freqname == NULL && (gname != NULL || beaglefile != NULL)){
    fprintf(stderr, "\t-> Number of genomic sites must be provided (-L <INT>)\n");
    return 0;
  }

  if((freqname!=NULL) && (nColInFile(freqname)!=1)){
    fprintf(stderr, "ERROR: more than one column (%d) in frequency file: %s\n", nColInFile(freqname), freqname);
    return 0;
  }


  
  // if ( freqname == NULL && ( do_simple || do_inbred ) && htsfile==NULL ) {
  //   fprintf(stderr, "\t-> Must supply -f (allele frequency file) if '-o 1' or '-F 1' are enabled\n");
  //   return 0;
  //   }
  
  // if (freqname == NULL && htsfile==NULL &&plinkfile==NULL){
  //   fprintf(stderr, "\t-> Allele frequencies file (-f) is not provided. Only summary statistitics based on 2dsfs will be reported\n");
  //   do_2dsfs_only = 1;
  //   if(!nsites_nofreqfile){
  //     fprintf(stderr, "\t-> Number of genomic sites must be provided (-L <INT>)\n");
  //     return 0;
  //   }
  // }

  if(nBootstrap>0&&(pair1==-1||pair2==-1)){
    fprintf(stderr,"\t-> Must specifiy pair of samples when performing bootstrap replicates\n");
    return 0;
  }
  
  if(plinkfile){
    fprintf(stderr,"\t-> Starting to read plinkfiles\n");
    int famnlines=nlines(plink_fam.c_str());
    int bimnlines=nlines(plink_bim.c_str());
    nind = famnlines;
    overall_number_of_sites = bimnlines;
    fprintf(stderr,"\t-> nlines in .fam:%d\t nlines in .bim:%d\n",famnlines,bimnlines);
    int **imat = bed_to_intMatrix(plink_bed.c_str(),famnlines,bimnlines);//leak
    //imat[ind][site]
    for(int i=0;i<bimnlines;i++){
      double hit =0;
      double asum =0;
      for(int j=0;j<famnlines;j++)
	if(imat[j][i]!=3){
	  hit+=1.0;
	  asum+=imat[j][i];
	}
      freq.push_back(1-asum/hit/2.0);//flip
    }
#if 0 //validated with plink --freq
    for(int i=0;i<bimnlines;i++)
      fprintf(stderr,"freq[%d]:\t%f\n",i,freq[i]);
#endif
    gls = new double*[bimnlines];
    for(int i=0;i<bimnlines;i++){
      gls[i] = new double[3*famnlines];
      for(int j=0;j<famnlines;j++){
	double *pi = gls[i]+j*3;
	if(imat[j][i]==3)
	  pi[0]=pi[1]=pi[2] = 1;
	else if(imat[j][i]==2){
	  pi[0] =1;
	  pi[1]=pi[2] =0;
	}
	else if(imat[j][i]==1){
	  pi[1] =1;
	  pi[0]=pi[2] =0;
	}
	else if(imat[j][i]==0){
	  pi[2] =1;
	  pi[0]=pi[1] =0;
	}
      }
    }
  }else{  
    if(htsfile==NULL){ // && !do_2dsfs_only){

      size_t tmp_freq_size;
      
      if(freqname!=NULL){
        getDouble(freqname,freq);
        tmp_freq_size = freq.size();
      } else {
        tmp_freq_size = nsites_nofreqfile;
      }
        
      if(beaglefile != NULL)
        gls = readBeagle(beaglefile, tmp_freq_size, nind);
      else if(gname != NULL)        
        gls = getGL(gname, tmp_freq_size, nind);

      if(freqname==NULL){
        // estimate allele frequencies
        fprintf(stderr,"\t-> Allele frequencies will be estimated from the data (N=%d)\n",nind);
        for(size_t i=0; i<tmp_freq_size; i++){
          double val = emFrequency_from_data(gls[i], nind, 50, 0.05);
          freq.push_back(val);
        }
      }      
      overall_number_of_sites = freq.size();
    }


#if 0
    for(size_t i=0; i<freq.size(); i++){
      fprintf(stdout, "%f\n", freq[i]);
    }
    
#endif


    
    // if(htsfile==NULL && do_2dsfs_only){
    //   if(beaglefile != NULL)
    //     gls = readBeagle(beaglefile, nsites_nofreqfile, nind);
    //   else if(gname != NULL)        
    //     gls = getGL(gname, nsites_nofreqfile, nind);

    //   overall_number_of_sites = nsites_nofreqfile;
      
    // }
  }

#ifdef __WITH_BCF__

  if(htsfile){
    if(freqname!=NULL){
      fprintf(stderr,"\t-> Allele frequencies provided to -f will be OVERWRITTEN with estimated allele frequencies from the data in the VCF file OR the provided TAG (-A) if present.\n");
      fprintf(stderr,"\t-> Annotate the VCF with the allele frequencies if external frequencies should be used.\n");
      fprintf(stderr, "\t-> Remove -f argument. EXITING\n");
      return 0;      
    }
    gls=readbcfvcf(htsfile,nind,freq,2,minMaf, vcf_format_field, vcf_allele_field,region);
    overall_number_of_sites = freq.size();
  }
  
#endif

  total_sites = overall_number_of_sites * 1.0;
  
  //all data read from either 1) glf/freq 2) hts/vcf/bcf 3)plink
  //now call genotypes if needed
  fprintf(stderr,"\t-> nind:%d overall_number_of_sites:%lu\n",nind,overall_number_of_sites);
  fflush(stderr);
  if(overall_number_of_sites==0)
    return 0;
  if (switchMaf) {
    fprintf(stderr, "\t-> switching frequencies\n");
    for (size_t i = 0; i < overall_number_of_sites; i++)
      freq[i] = 1 - freq[i];
  }
  
  if (gc) {
    if (gc> 1) {
      fprintf(stderr,"\t-> Calling genotypes assuming hwe\n");
      callgenotypesHwe(gls, overall_number_of_sites, nind, freq);
    }
    
    if (gc > 0){
      fprintf(stderr,"\t-> Modelling errors for genotypes (should only be used for called genotypes)\n");    
      callgenotypesEps(gls, overall_number_of_sites, nind, errate);
    }
    fprintf(stderr,"\t-> Done calling genotypes\n");
  }

#ifdef DB_AF //for printout everything
  FILE * af_fp = fopen("db_af.txt", "w");
  for(int i=0;i<freq.size();i++){
    fprintf(af_fp,"%f\n",freq[i]);
  }
  fclose(af_fp);
#endif

  
#ifdef DB_GLS //for printout everything
  FILE * gls_fp = fopen("db_gls.txt", "w");
  for(int i=0;i<freq.size();i++){
    fprintf(gls_fp,"%f",gls[i][0]);    
    for(int ii=1;ii<3*nind;ii++)
      fprintf(gls_fp,"\t%f",gls[i][ii]);
    fprintf(gls_fp, "\n");
  }
  fclose(gls_fp);
#endif
 

  float splittimes[2] = {(float)(clock() - t) / CLOCKS_PER_SEC,(float)(time(NULL) - t2)};
  fprintf(stderr,"\t-> Done reading data from file: %.2f %.2f\n",splittimes[0],splittimes[1]);
  fflush(stderr);
  FILE *output = NULL;
  output = fopen(outname,"wb");
  assert(output);
  for(int i=0;i<num_threads;i++){
    char buf[strlen(outname)+20];
    snprintf(buf,strlen(outname)+20,"%s.spill%d.res",outname,i);
    FILE *fp=NULL;
    fp=fopen(buf,"wb");
    assert(fp!=NULL);
    spillfiles.push_back(fp);
    spillfilesnames.push_back(strdup(buf));
  }
  

  fprintf(stderr,"\t-> Starting analysis now\n");
  main_analysis2(freq,gls,num_threads,output,total_sites);
  
  for (size_t i = 0; i < overall_number_of_sites; i++) {
    delete[] gls[i];
  }
  delete[] gls;
  free(freqname);
  free(gname);
  free(htsfile);
  fclose(output);
  for(int i=0;i<spillfiles.size();i++){
    fclose(spillfiles[i]);
    unlink(spillfilesnames[i]);
    free(spillfilesnames[i]);
  }

  fprintf(stderr, "\t[ALL done] cpu-time used =  %.2f sec (filereading took: %.2f sec)\n", (float)(clock() - t) / CLOCKS_PER_SEC,splittimes[0]);
  fprintf(stderr, "\t[ALL done] walltime used =  %.2f sec (filereading took: %.2f sec)\n", (float)(time(NULL) - t2),splittimes[1]);  
  
  if(outname) free(outname);
  if(region) free(region);
  if(vcf_allele_field) free(vcf_allele_field);
  if(vcf_format_field) free(vcf_format_field);
  return 0;
}
