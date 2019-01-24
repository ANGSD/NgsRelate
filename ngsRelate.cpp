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
#include "filereaders.h"

#ifdef __WITH_BCF__
#include "vcf.h"
#endif

std::vector<char *> posinfo;//<- debug


//this is as close to the bound we will allow
double TINY=1e-8;
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
double tole =1e-8;
int n=-1;

int seed=std::numeric_limits<int>::max();

int model =1;
int gc =0;
double errate = 0.005;
int pair1 =-1;
int pair2 =-1;
int nind =2;
int nsites_2dsfs = 0;
size_t overall_number_of_sites = 0;
int do_2dsfs_only = 0;
int do_inbred=0;
int do_simple=0;
int switchMaf = 0;

double minMaf =0.05;
int hasDef = 0;
double ttol=1e-6;
std::string vcf_format_field = "PL"; // can take PL or GT
std::string vcf_allele_field = "AFngsrelate"; // can take any tag value e.g. AF AF1 etc

std::vector<char *> ids;

float emTole=1e-12;

// https://en.cppreference.com/w/c/numeric/math/isnan
bool is_nan(double x) { return x != x; }

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

double loglike(double *p,double **emis,int len,int dim){
  double ret =0;

  for(int i=0;i<len;i++){
    double tmp = 0;
    for(int j=0;j<dim;j++)
      tmp += p[j]*emis[i][j];
    ret +=log(tmp);
  }
  return ret;
}

void emStep(double *pre,double **emis,double *post,int len,int len2){
  for(int i=0;i<len2;i++){
    if(pre[i]<0||pre[i]>1){
      fprintf(stderr,"Problem with gues in emStep: ");
      for(int j=0;j<len2;j++)
	fprintf(stderr," %f ",pre[0]);
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
    for(int x=0;x<len2;x++){
      post[x] += inner[x];
    }

  }
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


  double F_em1[dim];
  double F_diff1[dim];
  double F_em2[dim];
  double F_diff2[dim];
  double F_diff3[dim];
  double F_tmp[dim];
  niter++;
  emStep(F, emis, F_em1, len,dim);
  // stayin(F_em1);
  
  minus(F_em1, F, F_diff1,dim);

  double sr2 = sumSquare(F_diff1,dim);
  
  if(sqrt(sr2)<ttol){
    //     fprintf(stderr,"sr2 break: %e\n", sqrt(sr2));
    for(int i=0;0&&i<dim;i++)
      F_new[i]  = F_em1[i];
    return 0;
  }
  niter++;
  emStep(F_em1, emis, F_em2, len,dim);
  minus(F_em2, F_em1, F_diff2,dim);

  double sq2 = sumSquare(F_diff2,dim);

  if(sqrt(sq2)<ttol){
    //    fprintf(stderr,"sq2 break: %e\n", sqrt(sq2));
    for(int i=0;0&&i<dim;i++)
      F_new[i]  = F_em2[i];
    return 0;
  }

  minus(F_diff2,F_diff1, F_diff3,dim);

  double sv2 = sumSquare(F_diff3,dim);

  double alpha = sqrt(sr2/sv2);
  alpha = std::max(stepMin,std::min(stepMax,alpha));
  for(size_t i=0;i<dim;i++){
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

int em1(double *sfs,double  **emis, int len, int dim){
  int niter = 0;
  double oldLik,lik;
  oldLik = loglike(sfs,emis,len,dim);

  double tmp[dim];
  int it;
  for(it=0;niter<maxIter;it++) {
    niter++;
    emStep(sfs,emis,tmp,len,dim);
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


int em2(double *sfs,double  **emis, int len,int dim){
  int niter=0;
  double oldLik,lik;
  oldLik = loglike(sfs,emis,len,dim);

  double tmp[dim];
  int it;
  for(it=0;niter<maxIter;it++) {
  emAccel(sfs,emis,tmp,len, niter,dim);
    for(int i=0;i<dim;i++)
      sfs[i]= tmp[i];

  lik = loglike(sfs,emis,len,dim);
    if(isnan(lik)){
      fprintf(stderr,"em2 evaluates like to NaN\n");
      exit(0);
      break;
      //  exit(0);
    }
    if(fabs(lik-oldLik)<tole){
      // fprintf(stderr,"breaking\n");
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
#if 0
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
  fprintf(fp, "   -L <INT>            Number of genomic sites. Must be provided if -f (allele frequency file) is NOT provided \n");
  fprintf(fp, "   -m <INTEGER>        model 0=normalEM 1=acceleratedEM\n");
  fprintf(fp, "   -i <UINTEGER>       Maximum number of EM iterations\n");
  fprintf(fp, "   -t <FLOAT>          Tolerance for breaking EM\n");
  fprintf(fp, "   -r <FLOAT>          Seed for rand\n");
  fprintf(fp, "   -g gfile            Name of genotypellh file\n");
  fprintf(fp, "   -p <INT>            threads (default 4)\n");
  fprintf(fp, "   -c <INT>            Should call genotypes instead?\n");
  fprintf(fp, "   -s <INT>            Should you swich the freq with 1-freq?\n");
  fprintf(fp, "   -F <INT>            Estimate inbreeding instead of estimating the nine jacquard coefficients\n");
  fprintf(fp, "   -o <INT>            estimating the 3 jacquard coefficient, assumming no inbreeding\n");
  fprintf(fp, "   -v <INT>            Verbose. print like per iteration\n");
  fprintf(fp, "   -e <INT>            Errorrates when calling genotypes?\n");
  fprintf(fp, "   -a <INT>            First individual used for analysis? (zero offset)\n");
  fprintf(fp, "   -b <INT>            Second individual used for analysis? (zero offset)\n");
  fprintf(fp, "   -n <INT>            Number of samples in glf.gz\n");
  fprintf(fp, "   -l <INT>            minMaf or 1-Maf filter\n");
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



struct worker_args {
  int thread_id, nkeep, a,b, niter, best, niter_2dsfs;
  double **gls;
  std::vector<double> * freq;
  size_t nsites;
  double ll, bestll, ll_2dsfs;
  double pars[9], pars_2dsfs[9];
  int *keeplist;
  double **emis;
  worker_args(int & id_a, int & id_b, std::vector<double> * f, double ** gls_arg, size_t & s ){
    a=id_a;
    b=id_b;
    nkeep=0;
    best=0;
    bestll=0.0;
    bestll=0.0;    
    freq = f;
    gls=gls_arg;
    nsites = s;
    keeplist = new int[nsites];
    emis = new double*[s];
    for(int i=0;i<s;i++)
      new double[9];//<- will hold, 2dsfs,F and 9jacq emissions depending on context
  }
};


void * do_work(void *threadarg){

  // https://www.tutorialspoint.com/cplusplus/cpp_multithreading.htm
  struct worker_args * td;
  td = ( worker_args * ) threadarg;

  assert(td->nsites>0);
  // init all in each thread
  int *keeplist = td->keeplist;

  for (size_t i = 0; i < td->nsites; i++) {

    if(is_missing(&td->gls[i][3*td->a]))
      continue;
    if(is_missing(&td->gls[i][3*td->b]))
      continue;

    // removing minor allele frequencies
    if ( !do_2dsfs_only && (td->freq->at(i) < minMaf || (1 - td->freq->at(i)) < minMaf))
      continue;

    keeplist[td->nkeep] = i;//dont forget
    td->nkeep++;
  }
  //  fprintf(stderr,"td->nkeep:%d\n",td->nkeep);exit(0);
  if (td->nkeep==0){
    fprintf(stderr, "sites with both %d and %d having data: %d\n", td->a, td->b, td->nkeep);
  }
 

  // fprintf(stderr, "\t-> keeping %d sites for downstream analyses", td->nkeep++);
  double **emis=td->emis;
  if(!do_2dsfs_only){
    emis = new double *[td->nkeep];
    for (int i = 0; i < td->nkeep; i++) {
      emis[i] = new double[9];
    }
  }

  double **emislike_2dsfs = new double *[td->nkeep];
  for (int i = 0; i < td->nkeep; i++) {
    emislike_2dsfs[i] = new double[9];
  }

  if(!do_2dsfs_only){
    emission_ngsrelate9(td->freq, td->gls, emis, keeplist, td->nkeep, td->a, td->b);

    if (model == 0){
      td->niter = em1(td->pars, emis, td->nkeep,9);
    }else if (model == 1){
      td->niter = em2(td->pars, emis, td->nkeep,9);
    }

    double l100000000 = loglike(p100000000, emis, td->nkeep,9);
    double l010000000 = loglike(p010000000, emis, td->nkeep,9);
    double l001000000 = loglike(p001000000, emis, td->nkeep,9);
    double l000100000 = loglike(p000100000, emis, td->nkeep,9);
    double l000010000 = loglike(p000010000, emis, td->nkeep,9);
    double l000001000 = loglike(p000001000, emis, td->nkeep,9);
    double l000000100 = loglike(p000000100, emis, td->nkeep,9);
    double l000000010 = loglike(p000000010, emis, td->nkeep,9);
    double l000000001 = loglike(p000000001, emis, td->nkeep,9);
    double lopt = loglike(td->pars, emis, td->nkeep,9);
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
  }

  emislike_2dsfs_gen(td->gls, emislike_2dsfs, keeplist, td->nkeep, td->a, td->b);

  if (model == 0){
    td->niter_2dsfs = em1(td->pars_2dsfs, emislike_2dsfs, td->nkeep,9);
  }else if (model == 1){
    td->niter_2dsfs = em2(td->pars_2dsfs, emislike_2dsfs, td->nkeep,9);  
  }
  td->ll_2dsfs = loglike(td->pars_2dsfs, emislike_2dsfs, td->nkeep,9);

  // end of work. Teardown

  if(do_2dsfs_only){
    for (int i = 0; i < td->nkeep; i++) {
      delete[] emislike_2dsfs[i];
    }
    delete[] emislike_2dsfs;
  } else {
    for (int i = 0; i < td->nkeep; i++) {
      delete[] emis[i];
      delete[] emislike_2dsfs[i];
    }
    delete[] emis;
    delete[] emislike_2dsfs;
  }
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
  int *keeplist = td->keeplist;
  for (size_t i = 0; i < td->freq->size(); i++) {

    if(is_missing(&td->gls[i][3*td->a]))
      continue;

    // removing minor allele frequencies
    if (td->freq->at(i) < minMaf || (1-td->freq->at(i)) < minMaf)
      continue;

    keeplist[td->nkeep] = i;
    td->nkeep++;
  }
  // fprintf(stderr, "\t-> keeping %d sites for downstream analyses", td->nkeep++);
  double **emis = new double *[td->nkeep];
  for (int i = 0; i < td->nkeep; i++) {
    emis[i] = new double[2];
  }

  emission_ngs_inbred(td->freq,td->gls, emis, keeplist, td->nkeep, td->a);
  if (model == 0)
    td->niter=em1(td->pars,emis,td->nkeep,2);
  else
    td->niter=em2(td->pars,emis,td->nkeep,2);
  td->ll = loglike(td->pars,emis,td->nkeep,2);
  double l01= loglike(p01,emis,td->nkeep,2);
  double l10= loglike(p10,emis,td->nkeep,2);
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

int main_analysis1(std::vector<double> &freq,double **gls,int num_threads,FILE *output,int total_sites){
  pthread_t threads[num_threads];
  if(do_inbred){
    fprintf(output,"Ind\tZ=0\tZ=1\tloglh\tnIter\tcoverage\tsites\n");
    int comparison_ids_inbred = 0;
    std::vector<worker_args> all_args_inbred;
    int fake_person = -1;
    for(int a=0;a<nind;a++){
      if (pair1 != -1)
        a = pair1;
      worker_args td_args_inbred(a, fake_person, &freq, gls, overall_number_of_sites);
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
      if(comparison_ids_inbred-cnt_inbred-num_threads>=0){
        nTimes_inbred = num_threads;
      } else{
        nTimes_inbred = comparison_ids_inbred-cnt_inbred;
      }
#if 0
      fprintf(stderr, "\ngot this far. while %d.\n", cnt_inbred);
      fprintf(stderr, "comparisons: %d\n", comparison_ids_inbred);
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
        worker_args * td_out_inbred = &all_args_inbred[cnt_inbred + i];
        if(td_out_inbred->best==0)
          fprintf(output,"%d\t%f\t%f\t%f\t%d\t%f\t%d\n",td_out_inbred->a, p10[0],p10[1],td_out_inbred->bestll,-1,((double)td_out_inbred->nkeep)/((double)overall_number_of_sites), td_out_inbred->nkeep);
        if(td_out_inbred->best==1)
          fprintf(output,"%d\t%f\t%f\t%f\t%d\t%f\t%d\n",td_out_inbred->a,p01[0],p01[1],td_out_inbred->bestll,-1,((double)td_out_inbred->nkeep)/((double)overall_number_of_sites), td_out_inbred->nkeep);
        if(td_out_inbred->best==2)
          fprintf(output,"%d\t%f\t%f\t%f\t%d\t%f\t%d\n",td_out_inbred->a,td_out_inbred->pars[0],td_out_inbred->pars[1],td_out_inbred->bestll,td_out_inbred->niter,((double)td_out_inbred->nkeep)/((double)overall_number_of_sites), td_out_inbred->nkeep);
      fflush(output);
      }
      cnt_inbred += nTimes_inbred;
      fprintf(stderr, "\t-> Processed %d out of %d\r", cnt_inbred, comparison_ids_inbred);
    }
  } else {
    if(do_simple){
      fprintf(stderr, "\t-> setting Jacquard coefficient related to inbreeding (1-6) to zero\n");
    }
    
    if (ids.size()){
      // fprintf(output,"a\tb\tida\tidb\tnSites\tJ9\tJ8\tJ7\tJ6\tJ5\tJ4\tJ3\tJ2\tJ1\trab\tFa\tFb\ttheta\tloglh\tnIter\tcoverage\n");
      fprintf(output,
              "a\tb\tida\tidb\tnSites\tJ9\tJ8\tJ7\tJ6\tJ5\tJ4\tJ3\tJ2\tJ1\trab\tFa\tFb\ttheta\tinbred_relatedness_1_2\tinbred_relatedness_2_1\tfraternity\tidentity\tzygosity\t2of3_IDB\tF_diff_a_b\tloglh\tnIter\tcoverage\t2dsfs\tR0\tR1\tKING\t2dsfs_loglike\t2dsfsf_niter\n");
    } else {
      fprintf(output, "a\tb\tnSites\tJ9\tJ8\tJ7\tJ6\tJ5\tJ4\tJ3\tJ2\tJ1\trab\tFa\tFb\ttheta\tinbred_relatedness_1_2\tinbred_relatedness_2_1\tfraternity\tidentity\tzygosity\t2of3IDB\tFDiff\tloglh\tnIter\tcoverage\t2dsfs\tR0\tR1\tKING\t2dsfs_loglike\t2dsfsf_niter\n");
    }
    int comparison_ids = 0;
    std::vector<worker_args> all_args;
    for (int a = 0; a < nind; a++) {
      for (int b = a + 1; b < nind; b++) {
        if (pair1 != -1)
          a = pair1;
        if (pair2 != -1)
          b = pair2;

        worker_args td_args(a, b, &freq, gls, overall_number_of_sites);
        
        double parval = 0.0, parsum = 0.0;
        for (int i = 0; i < 9; i++) {
          parval = drand48();
          if(do_2dsfs_only){
            td_args.pars[i] = 0;
          } else if(do_simple && i>=3){
            td_args.pars[i] = 0;
          }else {
            td_args.pars[i] = parval;
            parsum += parval;
          }
        }
        if(parsum>0){
          for (int i = 0; i < 9; i++) {
            td_args.pars[i] /= parsum;
            // fprintf(stderr, "par: %d, %f\n", i, td_args.pars[i]);
          }
        }

        // for 2d sfs parameters
        parval = 0.0, parsum = 0.0;
        for (int i = 0; i < 9; i++) {
          parval = drand48();
          td_args.pars_2dsfs[i] = parval;
          parsum += parval;
        }

        for (int i = 0; i < 9; i++) {
          td_args.pars_2dsfs[i] /= parsum;
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
          fprintf(output, "%d\t%d\t%s\t%s\t%d", td_out->a,
                  td_out->b, ids[td_out->a],
                  ids[td_out->b], td_out->nkeep);
        } else {
          fprintf(output, "%d\t%d\t%d", td_out->a, td_out->b,
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

        for (int j = 0; j < 9; j++) {
          fprintf(output, "\t%f", td_out->pars[j]);
        }

        //////////////////////////////////////////////////////////////////////
        // Measuring Relatedness between Inbred Individuals by Hedrick 2015 //
        // return(x[1]+x[7]+3/4*(x[3]+x[5])+x[8]*0.5)                       //
        // r_xy                                                             //
        //////////////////////////////////////////////////////////////////////
        double rab = (td_out->pars[8] + td_out->pars[2]) +
                     (0.75 * (td_out->pars[6] + td_out->pars[4])) +
                     td_out->pars[1] * 0.5;
        fprintf(output, "\t%f", rab);

        /////////////////////////////////////////////////
        // Fa - inbreeding coefficient of individual a //
        // sum(x[1]+x[2],x[3],x[4])                    //
        /////////////////////////////////////////////////
        double Fa = td_out->pars[8] + td_out->pars[7] + td_out->pars[6] +
                    td_out->pars[5];
        fprintf(output, "\t%f", Fa);
        
        /////////////////////////////////////////////////
        // Fb - inbreeding coefficient of individual b //
        // sum(x[1],x[2],x[5],x[6])                    //
        /////////////////////////////////////////////////
        double Fb = td_out->pars[8] + td_out->pars[7] + td_out->pars[4] +
                    td_out->pars[3];
        fprintf(output, "\t%f", Fb);

        /////////////////////////////////////////////////////
        // theta / coancestry / kinship coefficient / f_XY //
        /////////////////////////////////////////////////////

        double theta =
            td_out->pars[8] +
            0.5 * (td_out->pars[6] + td_out->pars[4] + td_out->pars[2]) +
            0.25 * td_out->pars[1];
        fprintf(output, "\t%f", theta);

        /////////////////////////////////////////////
        // summary statistics from ackerman et al. //
        /////////////////////////////////////////////
        double inbred_relatedness_1_2 = td_out->pars[8] + (td_out->pars[6] / 2.0);
        double inbred_relatedness_2_1 = td_out->pars[8] + (td_out->pars[4] / 2.0);
        double fraternity = td_out->pars[7] + td_out->pars[2];
        double identity = td_out->pars[8];
        double zygosity = td_out->pars[8] + td_out->pars[7] + td_out->pars[2];
        fprintf(output, "\t%f\t%f\t%f\t%f\t%f", inbred_relatedness_1_2,
                inbred_relatedness_2_1, fraternity, identity, zygosity);

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // summary statistics from Non-identifiability of identity coefficients at biallelic loci by miklos csuros //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        // at least one pair of IBD alleles among three randomly selected ones:
        double eq_11e = td_out->pars[8] + td_out->pars[7] +
          td_out->pars[6] + td_out->pars[4] + td_out->pars[2] + 0.5 * (td_out->pars[5] + td_out->pars[3] + td_out->pars[1]);
        // 
        double eq_11f = 0.5 * (td_out->pars[5] - td_out->pars[3]);
        fprintf(output, "\t%f\t%f", eq_11e, eq_11f);
        
        ///////////////////////////////
        // optimization of EM output //
        ///////////////////////////////
        if (td_out->best == 9) {
          fprintf(output, "\t%f\t%d\t%f", td_out->ll, td_out->niter,
                  (1.0 * td_out->nkeep) / total_sites);
        } else {
          fprintf(output, "\t%f;s%d_%f\t%d\t%f", td_out->ll, 9 - td_out->best,
                  td_out->bestll, -1, (1.0 * td_out->nkeep) / total_sites);
        }

        ////////////////////
        // printing 2dsfs //
        ////////////////////
        fprintf(output, "\t%e", td_out->pars_2dsfs[0]);
        for (int j = 1; j < 9; j++) {
          fprintf(output, ",%e", td_out->pars_2dsfs[j]);
        }
        // fprintf(output, "%f\t%f\t%f\n",
        // td_out->pars_2dsfs[0],td_out->pars_2dsfs[1],td_out->pars_2dsfs[2]);
        // fprintf(output, "%f\t%f\t%f\n",
        // td_out->pars_2dsfs[3],td_out->pars_2dsfs[4],td_out->pars_2dsfs[5]);
        // fprintf(output, "%f\t%f\t%f\n",
        // td_out->pars_2dsfs[6],td_out->pars_2dsfs[7],td_out->pars_2dsfs[8]);


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
       
        fprintf(output, "\t%f\t%f\t%f", r0, r1, king);

        fprintf(output, "\t%f\t%d\n", td_out->ll_2dsfs, td_out->niter_2dsfs);
      }
      cnt += nTimes;
      fprintf(stderr, "\t-> Processed     %d out of       %d\r", cnt, comparison_ids);
    }
  }
  fprintf(stderr,"\n");
}

typedef struct{
  int a;
  int b;
}mypair;


int main_analysis2(std::vector<double> &freq,double **gls,int num_threads,FILE *output,int total_sites){
  std::vector<mypair> mp;
  for(int i=0;i<nind;i++)
    for(int j=(i+1);j<nind;j++){
      mypair tmp;tmp.a=i;tmp.b=j;
      mp.push_back(tmp);
    }
  std::random_shuffle(mp.begin(),mp.end());
  fprintf(stderr,"\t-> length of joblist:%lu\n",mp.size());
  return 0;

  pthread_t threads[num_threads];
  if(do_inbred){
    fprintf(output,"Ind\tZ=0\tZ=1\tloglh\tnIter\tcoverage\tsites\n");
    int comparison_ids_inbred = 0;
    std::vector<worker_args> all_args_inbred;
    int fake_person = -1;
    for(int a=0;a<nind;a++){
      if (pair1 != -1)
        a = pair1;
      worker_args td_args_inbred(a, fake_person, &freq, gls, overall_number_of_sites);
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
      if(comparison_ids_inbred-cnt_inbred-num_threads>=0){
        nTimes_inbred = num_threads;
      } else{
        nTimes_inbred = comparison_ids_inbred-cnt_inbred;
      }
#if 0
      fprintf(stderr, "\ngot this far. while %d.\n", cnt_inbred);
      fprintf(stderr, "comparisons: %d\n", comparison_ids_inbred);
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
        worker_args * td_out_inbred = &all_args_inbred[cnt_inbred + i];
        if(td_out_inbred->best==0)
          fprintf(output,"%d\t%f\t%f\t%f\t%d\t%f\t%d\n",td_out_inbred->a, p10[0],p10[1],td_out_inbred->bestll,-1,((double)td_out_inbred->nkeep)/((double)overall_number_of_sites), td_out_inbred->nkeep);
        if(td_out_inbred->best==1)
          fprintf(output,"%d\t%f\t%f\t%f\t%d\t%f\t%d\n",td_out_inbred->a,p01[0],p01[1],td_out_inbred->bestll,-1,((double)td_out_inbred->nkeep)/((double)overall_number_of_sites), td_out_inbred->nkeep);
        if(td_out_inbred->best==2)
          fprintf(output,"%d\t%f\t%f\t%f\t%d\t%f\t%d\n",td_out_inbred->a,td_out_inbred->pars[0],td_out_inbred->pars[1],td_out_inbred->bestll,td_out_inbred->niter,((double)td_out_inbred->nkeep)/((double)overall_number_of_sites), td_out_inbred->nkeep);
      fflush(output);
      }
      cnt_inbred += nTimes_inbred;
      fprintf(stderr, "\t-> Processed %d out of %d\r", cnt_inbred, comparison_ids_inbred);
    }
  } else {
    if(do_simple){
      fprintf(stderr, "\t-> setting Jacquard coefficient related to inbreeding (1-6) to zero\n");
    }
    
    if (ids.size()){
      // fprintf(output,"a\tb\tida\tidb\tnSites\tJ9\tJ8\tJ7\tJ6\tJ5\tJ4\tJ3\tJ2\tJ1\trab\tFa\tFb\ttheta\tloglh\tnIter\tcoverage\n");
      fprintf(output,
              "a\tb\tida\tidb\tnSites\tJ9\tJ8\tJ7\tJ6\tJ5\tJ4\tJ3\tJ2\tJ1\trab\tFa\tFb\ttheta\tinbred_relatedness_1_2\tinbred_relatedness_2_1\tfraternity\tidentity\tzygosity\t2of3_IDB\tF_diff_a_b\tloglh\tnIter\tcoverage\t2dsfs\tR0\tR1\tKING\t2dsfs_loglike\t2dsfsf_niter\n");
    } else {
      fprintf(output, "a\tb\tnSites\tJ9\tJ8\tJ7\tJ6\tJ5\tJ4\tJ3\tJ2\tJ1\trab\tFa\tFb\ttheta\tinbred_relatedness_1_2\tinbred_relatedness_2_1\tfraternity\tidentity\tzygosity\t2of3IDB\tFDiff\tloglh\tnIter\tcoverage\t2dsfs\tR0\tR1\tKING\t2dsfs_loglike\t2dsfsf_niter\n");
    }
    int comparison_ids = 0;
    std::vector<worker_args> all_args;
    for (int a = 0; a < nind; a++) {
      for (int b = a + 1; b < nind; b++) {
        if (pair1 != -1)
          a = pair1;
        if (pair2 != -1)
          b = pair2;

        worker_args td_args(a, b, &freq, gls, overall_number_of_sites);
        
        double parval = 0.0, parsum = 0.0;
        for (int i = 0; i < 9; i++) {
          parval = drand48();
          if(do_2dsfs_only){
            td_args.pars[i] = 0;
          } else if(do_simple && i>=3){
            td_args.pars[i] = 0;
          }else {
            td_args.pars[i] = parval;
            parsum += parval;
          }
        }
        if(parsum>0){
          for (int i = 0; i < 9; i++) {
            td_args.pars[i] /= parsum;
            // fprintf(stderr, "par: %d, %f\n", i, td_args.pars[i]);
          }
        }

        // for 2d sfs parameters
        parval = 0.0, parsum = 0.0;
        for (int i = 0; i < 9; i++) {
          parval = drand48();
          td_args.pars_2dsfs[i] = parval;
          parsum += parval;
        }

        for (int i = 0; i < 9; i++) {
          td_args.pars_2dsfs[i] /= parsum;
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
          fprintf(output, "%d\t%d\t%s\t%s\t%d", td_out->a,
                  td_out->b, ids[td_out->a],
                  ids[td_out->b], td_out->nkeep);
        } else {
          fprintf(output, "%d\t%d\t%d", td_out->a, td_out->b,
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

        for (int j = 0; j < 9; j++) {
          fprintf(output, "\t%f", td_out->pars[j]);
        }

        //////////////////////////////////////////////////////////////////////
        // Measuring Relatedness between Inbred Individuals by Hedrick 2015 //
        // return(x[1]+x[7]+3/4*(x[3]+x[5])+x[8]*0.5)                       //
        // r_xy                                                             //
        //////////////////////////////////////////////////////////////////////
        double rab = (td_out->pars[8] + td_out->pars[2]) +
                     (0.75 * (td_out->pars[6] + td_out->pars[4])) +
                     td_out->pars[1] * 0.5;
        fprintf(output, "\t%f", rab);

        /////////////////////////////////////////////////
        // Fa - inbreeding coefficient of individual a //
        // sum(x[1]+x[2],x[3],x[4])                    //
        /////////////////////////////////////////////////
        double Fa = td_out->pars[8] + td_out->pars[7] + td_out->pars[6] +
                    td_out->pars[5];
        fprintf(output, "\t%f", Fa);
        
        /////////////////////////////////////////////////
        // Fb - inbreeding coefficient of individual b //
        // sum(x[1],x[2],x[5],x[6])                    //
        /////////////////////////////////////////////////
        double Fb = td_out->pars[8] + td_out->pars[7] + td_out->pars[4] +
                    td_out->pars[3];
        fprintf(output, "\t%f", Fb);

        /////////////////////////////////////////////////////
        // theta / coancestry / kinship coefficient / f_XY //
        /////////////////////////////////////////////////////

        double theta =
            td_out->pars[8] +
            0.5 * (td_out->pars[6] + td_out->pars[4] + td_out->pars[2]) +
            0.25 * td_out->pars[1];
        fprintf(output, "\t%f", theta);

        /////////////////////////////////////////////
        // summary statistics from ackerman et al. //
        /////////////////////////////////////////////
        double inbred_relatedness_1_2 = td_out->pars[8] + (td_out->pars[6] / 2.0);
        double inbred_relatedness_2_1 = td_out->pars[8] + (td_out->pars[4] / 2.0);
        double fraternity = td_out->pars[7] + td_out->pars[2];
        double identity = td_out->pars[8];
        double zygosity = td_out->pars[8] + td_out->pars[7] + td_out->pars[2];
        fprintf(output, "\t%f\t%f\t%f\t%f\t%f", inbred_relatedness_1_2,
                inbred_relatedness_2_1, fraternity, identity, zygosity);

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // summary statistics from Non-identifiability of identity coefficients at biallelic loci by miklos csuros //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        // at least one pair of IBD alleles among three randomly selected ones:
        double eq_11e = td_out->pars[8] + td_out->pars[7] +
          td_out->pars[6] + td_out->pars[4] + td_out->pars[2] + 0.5 * (td_out->pars[5] + td_out->pars[3] + td_out->pars[1]);
        // 
        double eq_11f = 0.5 * (td_out->pars[5] - td_out->pars[3]);
        fprintf(output, "\t%f\t%f", eq_11e, eq_11f);
        
        ///////////////////////////////
        // optimization of EM output //
        ///////////////////////////////
        if (td_out->best == 9) {
          fprintf(output, "\t%f\t%d\t%f", td_out->ll, td_out->niter,
                  (1.0 * td_out->nkeep) / total_sites);
        } else {
          fprintf(output, "\t%f;s%d_%f\t%d\t%f", td_out->ll, 9 - td_out->best,
                  td_out->bestll, -1, (1.0 * td_out->nkeep) / total_sites);
        }

        ////////////////////
        // printing 2dsfs //
        ////////////////////
        fprintf(output, "\t%e", td_out->pars_2dsfs[0]);
        for (int j = 1; j < 9; j++) {
          fprintf(output, ",%e", td_out->pars_2dsfs[j]);
        }
        // fprintf(output, "%f\t%f\t%f\n",
        // td_out->pars_2dsfs[0],td_out->pars_2dsfs[1],td_out->pars_2dsfs[2]);
        // fprintf(output, "%f\t%f\t%f\n",
        // td_out->pars_2dsfs[3],td_out->pars_2dsfs[4],td_out->pars_2dsfs[5]);
        // fprintf(output, "%f\t%f\t%f\n",
        // td_out->pars_2dsfs[6],td_out->pars_2dsfs[7],td_out->pars_2dsfs[8]);


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
       
        fprintf(output, "\t%f\t%f\t%f", r0, r1, king);

        fprintf(output, "\t%f\t%d\n", td_out->ll_2dsfs, td_out->niter_2dsfs);
      }
      cnt += nTimes;
      fprintf(stderr, "\t-> Processed     %d out of       %d\r", cnt, comparison_ids);
    }
  }
  fprintf(stderr,"\n");
}

int main(int argc, char **argv){
  if(argc==1)
    print_info(stderr);

  if(strcasecmp(argv[1],"extract_freq_bim")==0)
    return extract_freq_bim(--argc,++argv);
  if(strcasecmp(argv[1],"extract_freq")==0)
    return extract_freq(--argc,++argv);
#if 0
  fprintf(stdout,"#");
  for(int i=0;i<argc;i++)
    fprintf(stdout," %s",argv[i]);
  fprintf(stdout,"\n");
#endif

  char *htsfile=NULL;
  char *plinkfile=NULL;
  const char *outname=NULL;
  int faster =0;
  while ((n = getopt(argc, argv, "f:i:t:r:g:m:s:F:o:c:e:a:b:n:l:z:p:h:L:T:A:P:O:X:")) >= 0) {
    switch (n) {
    case 'f': freqname = strdup(optarg); break;
    case 'P': plinkfile = strdup(optarg); break;
    case 'O': outname = strdup(optarg); break;
    case 'i': maxIter = atoi(optarg); break;
    case 'X': faster = atoi(optarg); break;
    case 't': tole = atof(optarg); break;
    case 'r': seed = atoi(optarg); break;
    case 'g': gname = strdup(optarg); break;
    case 'm': model = atoi(optarg); break;
    case 's': switchMaf = atoi(optarg); break;
    case 'F': do_inbred = atoi(optarg); break;
    case 'o': do_simple = atoi(optarg); break;
    case 'c': gc = atoi(optarg); break;
    case 'a': pair1 = atoi(optarg); break;
    case 'b': pair2 = atoi(optarg); break;
    case 'n': {nind = atoi(optarg); hasDef=1; break;}
    case 'p': num_threads = atoi(optarg);break;
    case 'e': errate = atof(optarg); break;
    case 'l': minMaf = atof(optarg); break;
    case 'h': htsfile = strdup(optarg); break;
    case 'T': vcf_format_field = strdup(optarg); break;
    case 'A': vcf_allele_field = strdup(optarg); break;            
    case 'z': readids(ids,optarg); break;
    case 'L': nsites_2dsfs = atoi(optarg); break;
    default: {fprintf(stderr,"unknown arg:\n");return 0;}
      print_info(stderr);
    }
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

  if (htsfile!=NULL){
    fprintf(stderr, "\t-> Will use TAG: '%s' from the VCF file\n", vcf_format_field.c_str());
    fprintf(stderr, "\t-> Will use TAG: '%s' in the VCF file as allele frequency if present. Otherwise allele frequencies are estimated from the data\n", vcf_allele_field.c_str());
  }
  srand48(seed);

  if ((nind == -1 || gname == NULL)&&htsfile==NULL&&plinkfile==NULL) {
    fprintf(stderr, "\t-> Must supply -n -g parameters (%d,%s) OR -h file.[vb]cf\n", nind,gname);
    return 0;
  }

  if ( freqname == NULL && ( do_simple || do_inbred ) && htsfile==NULL ) {
    fprintf(stderr, "\t-> Must supply -f (allele frequency file) if '-o 1' or '-F 1' are enabled\n");
    return 0;
    }
  
  if (freqname == NULL && htsfile==NULL &&plinkfile==NULL){
    fprintf(stderr, "\t-> Allele frequencies file (-f) is not provided. Only summary statistitics based on 2dsfs will be reported\n");
    do_2dsfs_only = 1;
    if(!nsites_2dsfs){
      fprintf(stderr, "\t-> Number of genomic sites must be provided (-L <INT>)\n");
      return 0;
    }
  }

  
  

  std::vector<double> freq;
  double **gls=NULL;

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
    if(htsfile==NULL && !do_2dsfs_only){
      getDouble(freqname,freq);
      gls = getGL(gname, freq.size(), nind);
      overall_number_of_sites = freq.size();
    }
    if(htsfile==NULL && do_2dsfs_only){
      gls = getGL(gname, nsites_2dsfs, nind);
      overall_number_of_sites = nsites_2dsfs;
    }
  }

#ifdef __WITH_BCF__
  if(htsfile){
    std::vector<double *> tmpgl;
    nind=getgls(htsfile,tmpgl,freq,2,minMaf, vcf_format_field, vcf_allele_field, errate,posinfo);
    gls=new double *[tmpgl.size()];
    for(int i=0;i<tmpgl.size();i++){
      gls[i] = tmpgl[i];
      for(int ii=0;ii<3*nind;ii++)
	gls[i][ii]=exp(gls[i][ii]);
    }
    overall_number_of_sites = freq.size();
  }
#endif

  double total_sites = overall_number_of_sites * 1.0;

  //all data read from either 1) glf/freq 2) hts/vcf/bcf 3)plink
  //now call genotypes if needed

  fprintf(stderr,"\t-> nind:%d overall_number_of_sites:%d\n",nind,overall_number_of_sites);
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
  }

#if 0 //for printout everything
  for(int i=0;i<freq.size();i++){
    fprintf(stdout,"%f",freq[i]);
    for(int ii=0;ii<3*nind;ii++)
      fprintf(stdout,"\t%f",gls[i][ii]);
    fprintf(stdout, "\n");
  }
  return 0;
#endif
  if(outname==NULL){
    fprintf(stderr,"\t-> No output filename declared please supply filename with -O output.res\n");
    return 0;
  }
  FILE *output = NULL;
  output = fopen(outname,"wb");
  assert(output);

  if(faster==0)
    main_analysis1(freq,gls,num_threads,output,total_sites);
  else
    main_analysis2(freq,gls,num_threads,output,total_sites);
  fflush(output);
  for (size_t i = 0; i < overall_number_of_sites; i++) {
    delete[] gls[i];
  }
  delete[] gls;
  free(freqname);
  free(gname);
  free(htsfile);
  fclose(output);
  return 0;
}
