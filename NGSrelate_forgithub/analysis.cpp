#include <cstdio>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include "analysis.h"
#include <cassert>
#include <cstring>
#include <cfloat>
#include "bfgs.h"


double alim[2]={0.001,0.15};
double klim[2]={DBL_EPSILON,1-DBL_EPSILON};
//double klim[2]={0.000001,0.999999};

double myRand(double low,double up){
  assert(up>low);
  double r = (drand48());//r=drand48();
  double range = up-low;
  return r*range+low;
}

void printPars(FILE *fp,para p){
  fprintf(fp,"\tpair=(%d,%d) alpha=%f k0=%f k1=%f k2=%f\n",p.pair[0],p.pair[1],p.a,p.k0,p.k1,p.k2);
}



double calculateA(double k0,double k1, double k2,double phi){
  double ma,mb,xa,xb,m,a,sq,pw;
  if(k2==0){
    ma = 1-log(k1)/log(2);
    mb = 0;
  }
  else{
    pw = pow(k1+2*k2,2);
    sq = sqrt(pw-4*k2);
    xa = (k1+2*k2+sq)/2.0;
    xb = k2/xa;
    ma = 1-log(xa)/log(2);
    mb = 1-log(xb)/log(2);
  }
  m = ma+mb;
  a = -m*log(1-phi);
  
  if(std::isnan(a)){//there is a bug in the bfgs algortihm, that makes the optim algo try outside of parameter space
    fprintf(stderr,"m=%f , ma=%f ,xa=%f, sq=%f  mb=%f  , a=%f  k0=%f ,  k1=%f   , k2=%f \n",m,ma,xa,sq,mb,a,k0,k1,k2);
    fprintf(stderr,"calc.a(k0=%f,k1=%f,k2=%f)\n",k0,k1,k2);
    fprintf(stderr,"std::isnan in calculate.a\n");
    fprintf(stderr,"Will return a very large value\n");
    return DBL_MAX;
  }

  else if(std::isinf(a)){
    fprintf(stderr,"std::isinf in calculate.a -> k0=%f,k1=%f,k2=%f,phi=%f\n",k0,k1,k2,phi);
    if(k2==0)
      fprintf(stderr,"k2=0,\t");
    fprintf(stderr,"ma=%f,mb=%f,m=%f\n",ma,mb,m);
	

  }
  
  return a;
}



void dump(char *fname,double **m,int x,int y){
  FILE *fp = fopen(fname,"w");
  for(int i=0;i<x;i++){
    for(int j=0;j<y;j++)
      fprintf(fp,"%f ",m[i][j]);
    fprintf(fp,"\n");
  }
  fclose(fp);
}
void dump(char *fname,double *m,int x){
  FILE *fp = fopen(fname,"w");
  for(int i=0;i<x;i++)
    fprintf(fp,"%f ",m[i]);
  fprintf(fp,"\n");
  fclose(fp);
}

template <typename T>
void fdump(char *fname,T *m,size_t x){
  FILE *fp = fopen(fname,"w");
  assert(x==fwrite(m,sizeof(T),x,fp));
  fclose(fp);
}
template <typename T>
void fdump(char *fname,T **m,int x,int y){
  FILE *fp = fopen(fname,"w");
  for(int i=0;i<x;i++)
    assert(y==fwrite(m[i],sizeof(T),y,fp));
  
  fclose(fp);
}


double max(double a[3]){
  double m=a[0];
  for(int i=1;i<2;i++)
    if(a[i]>m)
      m=a[i];
  return m;
}

/*
  emission ar put in the third argument.
  This function will allocate this matrix;

 */
void emis(int p1,int p2,const perChr &pc,double **ret){

  //arrays of 3*nSites of GLs
  double *l1 = pc.gl[p1];
  double *l2 = pc.gl[p2];

  //arrays of the frequencies
  double *f1,*f2;
  f1 = pc.qerf;//freqA
  f2 = pc.freq;//freqa
  //  fprintf(stderr,"f1=%f f2=%f\n",exp(f1[0]),exp(f2[0]));
  
  //(AA,AA)
  for(int s=0;s<pc.nSites;s++){
    double c1 = l1[3*s]+l2[3*s];
    //    fprintf(stderr,"%f %f c1=%f\n",l1[3*s],l2[3*s],c1);
    ret[0][s] = 4*f1[s]+c1;
    ret[1][s] = 3*f1[s]+c1;
    ret[2][s] = 2*f1[s]+c1;

  }

  //AA,Aa
  for(int s=0;s<pc.nSites;s++){

    double c1 = l1[3*s]+l2[3*s+1];
    ret[0][s] = addProtect2(ret[0][s],M_LN2+3*f1[s]+f2[s]+c1);
    ret[1][s] = addProtect2(ret[1][s],2*f1[s]+f2[s]+c1);
    //no third
  }
 
  //AA,aa
  for(int s=0;s<pc.nSites;s++) {
    double c1 = l1[3*s]+l2[3*s+2];
    ret[0][s] = addProtect2(ret[0][s],2*(f1[s]+f2[s])+c1);
    //no sec,third
  }
 
  //Aa,AA
  for(int s=0;s<pc.nSites;s++) {
    double c1 = l1[3*s+1]+l2[3*s];
    ret[0][s] = addProtect2(ret[0][s],M_LN2+3*f1[s]+f2[s]+c1);
    ret[1][s] = addProtect2(ret[1][s],2*f1[s]+f2[s]+c1);

   
  }
    //Aa,Aa
  for(int s=0;s<pc.nSites;s++) {
    double c1 = l1[3*s+1]+l2[3*s+1];
    ret[0][s] = addProtect2(ret[0][s],2*M_LN2+2*(f1[s]+f2[s])+c1);
    ret[1][s] = addProtect2(ret[1][s],f1[s]+f2[s]+c1);
    ret[2][s] = addProtect2(ret[2][s],M_LN2+f1[s]+f2[s]+c1);
  }
  //Aa,aa
  for(int s=0;s<pc.nSites;s++) {
     double c1 = l1[3*s+1]+l2[3*s+2];
    ret[0][s] = addProtect2(ret[0][s],M_LN2+f1[s]+3*f2[s]+c1);
    ret[1][s] = addProtect2(ret[1][s],f1[s]+2*f2[s]+c1);
  }
  //aa,AA
  for(int s=0;s<pc.nSites;s++) {
    double c1 = l1[3*s+2]+l2[3*s];
    ret[0][s] = addProtect2(ret[0][s],2*(f1[s]+f2[s])+c1);
    //no sec,third
  }
  //aa,Aa
  for(int s=0;s<pc.nSites;s++) {
    double c1 = l1[3*s+2]+l2[3*s+1];
    ret[0][s] = addProtect2(ret[0][s],M_LN2+f1[s]+3*f2[s]+c1);
    ret[1][s] = addProtect2(ret[1][s],f1[s]+2*f2[s]+c1);
    //no sec,third
  }
  for(int s=0;s<pc.nSites;s++) {
    double c1 = l1[3*s+2]+l2[3*s+2];
    ret[0][s] = addProtect2(ret[0][s],4*f2[s]+c1);
    ret[1][s] = addProtect2(ret[1][s],3*f2[s]+c1);
    ret[2][s] = addProtect2(ret[2][s],2*f2[s]+c1);
    //no sec,third

  } 
  
}


double *diffPos(const int *pos,int len){
  double *ret = new double [len+1];
  ret[0] = .0;
  for(int i=0;i<len-1;i++)
    ret[i+1] = (pos[i+1]-pos[i])/1e6;
  ret[len] = .0;
  return ret;

}

int T00 =0;
int T01 =1;
int T02 =2;
int T10 =3;
int T11 =4;
int T12 =5;
int T20 =6;
int T21 =7;
int T22 =8;

//pars is a,k0,k1,k2
//return matrix is
/*0  , 1 , 2 , 3 , 4 , 5 ,6  , 7 , 8
  T00,T01,T02,T10,T11,T12,T20,T21,T22
*/

void fixunderflow(double *a,int l, double tole){
  for(int i=0;i<l;i++)
    if(a[i]<tole)
      a[i]=0.0;

}
/*
  results will be put in last argument

 */
void trans(const double *pars,double *pos,int l,double **ret) {

  double a = pars[0];
  double k0= pars[1];
  double k1= pars[2];
  double k2= pars[3];
  //  fprintf(stderr,"[%s] a=%f k0=%f k1=%f k=%f\n",__FUNCTION__,a,k0,k1,k2);  
  //fprintf(stderr,"a=%f k0=%f k1=%f k2=%f\n",a,k0,k1,k2);
  //allocate return matrix
  //make a*pos :  t*a
  double *at = new double[l];
  for(int i=0;i<l;i++)
    at[i] = a*pos[i];

  //make T01
  for(int i=0;i<l;i++)
    ret[T01][i] = (1.0-exp(-at[i]))*k1;
  //make T02
  for(int i=0;i<l;i++)
    ret[T02][i] = exp(-at[i]*k1)*k2/(k1-1)+exp(-at[i])*k1+exp(-at[i])*k0*k1/(k1-1)+k2;
  if(k1==1)
    for(int i=0;i<l;i++)
      ret[T02][i] = 0;
  fixunderflow(ret[T02],l,1e-15);
  

  //make T00
  for(int i=0;i<l;i++)
    ret[T00][i] = 1-ret[T01][i]-ret[T02][i];

  //make T10,T12
  if(k1==0){
    memcpy(ret[T10],ret[T01],sizeof(double)*l);
    memcpy(ret[T12],ret[T01],sizeof(double)*l);
  }else{
    for(int i=0;i<l;i++){
      ret[T10][i] = ret[T01][i]*k0/k1;
      ret[T12][i] = ret[T01][i]*k2/k1;
      
    }
  }
  
  //make T00
  for(int i=0;i<l;i++)
    ret[T11][i] = 1-ret[T10][i]-ret[T12][i];
  
  //make T21
  memcpy(ret[T21],ret[T01],sizeof(double)*l);


  //make T20;
  for(int i=0;i<l;i++)
    ret[T20][i] = exp(-at[i]*k1)*k0/(k1-1.0)+exp(-at[i])*k1+exp(-at[i])*k2*k1/(k1-1)-(1-k1)*k0/(k1-1);// # why this trick?

  if(k1==1)
    for(int i=0;i<l;i++)
      ret[T20][i] = 0;
  fixunderflow(ret[T02],l,1e-15);
 //make T22
  for(int i=0;i<l;i++)
    ret[T22][i] = 1-ret[T21][i]-ret[T20][i];
  
  //fix underflow
  for(int i=0;i<9;i++)
    for(int s=0;s<l;s++)
      if(ret[i][s]<1e-15)
	ret[i][s]=0;
  
}

//alpha not used in this function

double fastlike(double *pars,double **eprob,double **tprob,int nSites) {
  //  fprintf(stderr,"[%s] %f %f %f = ",__FUNCTION__,pars[0],pars[1],pars[2]);
  double logres[3]={log(pars[0]),log(pars[1]),log(pars[2])};
  //  fprintf(stderr,"init logres: %f %f %f\n",logres[0],logres[1],logres[2]);
  double lp[3];

  double m=0;
  for(int i=1;i<nSites;i++){
    lp[0] = exp(logres[0]+eprob[0][i-1]-m);
    lp[1] = exp(logres[1]+eprob[1][i-1]-m);
    lp[2] = exp(logres[2]+eprob[2][i-1]-m);
    //    fprintf(stderr,"lps[%d]: %f %f %f\n",i,lp[0],lp[1],lp[2]);
    logres[0] = m+log(tprob[T00][i]*lp[0]+tprob[T10][i]*lp[1]+tprob[T20][i]*lp[2]);
    logres[1] = m+log(tprob[T01][i]*lp[0]+tprob[T11][i]*lp[1]+tprob[T21][i]*lp[2]);
    logres[2] = m+log(tprob[T02][i]*lp[0]+tprob[T12][i]*lp[1]+tprob[T22][i]*lp[2]);
    //    fprintf(stderr,"logres[%d]: %f %f %f\n",i,logres[0],logres[1],logres[2]);
    for(int ii=0;ii<3;ii++)
      lp[ii] = log(lp[ii])+m;
    m = max(lp);
    //    fprintf(stderr,"m[%d]=%f\n",i,m);
    if(std::isnan(m)||std::isnan(logres[0])){
      fprintf(stderr,"Probs in likelihood at site:%d\n",i);
      exit(0);
    }
  }
  
  lp[0] = eprob[0][nSites-1]+logres[0];
  lp[1] = eprob[1][nSites-1]+logres[1];
  lp[2] = eprob[2][nSites-1]+logres[2];

  m=max(lp);
  double res = m+log(exp(lp[0]-m) + exp(lp[1]-m)+ exp(lp[2]-m));
  //  fprintf(stderr,"[%s]\tloglike=%f\n",__FUNCTION__,res);
  return res;
}



double **forward(double *pars,double **eprob,double **tprob,int nSites,double &loglike){
  //fprintf(stderr,"[%s] %f %f %f\n",__FUNCTION__,pars[0],pars[1],pars[2]);
  double logres[3]={log(pars[0]),log(pars[1]),log(pars[2])};
  //  fprintf(stderr,"init logres: %f %f %f\n",logres[0],logres[1],logres[2]);

  double **res = new double*[3];
  for(int i=0;i<3;i++){
    res[i] = new double[nSites];
    res[i][0] = logres[i]+eprob[i][0]; 
    //    fprintf(stderr,"init=%f\n",res[i][0]);
  }
  

  double lp[3];

  double m=0;
  for(int i=1;i<nSites;i++){
    lp[0] = exp(res[0][i-1]-m);
    lp[1] = exp(res[1][i-1]-m);
    lp[2] = exp(res[2][i-1]-m);
    //    fprintf(stderr,"lps[%d]: %f %f %f\n",i,lp[0],lp[1],lp[2]);
    res[0][i] = m+log(tprob[T00][i]*lp[0]+tprob[T10][i]*lp[1]+tprob[T20][i]*lp[2])+eprob[0][i];
    res[1][i] = m+log(tprob[T01][i]*lp[0]+tprob[T11][i]*lp[1]+tprob[T21][i]*lp[2])+eprob[1][i];
    res[2][i] = m+log(tprob[T02][i]*lp[0]+tprob[T12][i]*lp[1]+tprob[T22][i]*lp[2])+eprob[2][i];
    //    fprintf(stderr,"logres[%d]: %f %f %f\n",i,logres[0],logres[1],logres[2]);
    m = res[0][i];
    for(int ii=1;ii<3;ii++)
      if(res[ii][i]>m)
	m=res[ii][i];
    
    // fprintf(stderr,"m[%d]=%f\n",i,m);
    // exit(0);
    if(std::isnan(m)||std::isnan(logres[0])){
      fprintf(stderr,"Probs in forward at site:%d\n",i);
      exit(0);
    }
  }
  loglike = m+log(exp(res[0][nSites-1]-m) + exp(res[1][nSites-1]-m)+ exp(res[2][nSites-1]-m));
  //  fprintf(stderr,"[%s]\tloglike=%f\n",__FUNCTION__,loglike);

  //  double res = m+log(exp(lp[0]-m) + exp(lp[1]-m)+ exp(lp[2]-m));
  return res;
}


double **backward(double *pars,double **eprob,double **tprob,int nSites,double &loglike){
  //fprintf(stderr,"[%s] %f %f %f\n",__FUNCTION__,pars[0],pars[1],pars[2]);
  double logres[3]={log(pars[0]),log(pars[1]),log(pars[2])};
  //  fprintf(stderr,"init logres: %f %f %f\n",logres[0],logres[1],logres[2]);

  double **res = new double*[3];
  for(int i=0;i<3;i++){
    res[i] = new double[nSites];
    res[i][nSites-1] = 0.0 ;
  }
  

  double lp[3];

  double m=0;
  for(int i=nSites-2;i>=0;i--){
    lp[0] = exp(eprob[0][i+1]+res[0][i+1]-m);
    lp[1] = exp(eprob[1][i+1]+res[1][i+1]-m);
    lp[2] = exp(eprob[2][i+1]+res[2][i+1]-m);
    //    fprintf(stderr,"lps[%d]: %f %f %f\n",i,lp[0],lp[1],lp[2]);

    res[0][i] = m+log(tprob[T00][i+1]*lp[0]+tprob[T01][i+1]*lp[1]+tprob[T02][i+1]*lp[2]);
    res[1][i] = m+log(tprob[T10][i+1]*lp[0]+tprob[T11][i+1]*lp[1]+tprob[T12][i+1]*lp[2]);
    res[2][i] = m+log(tprob[T20][i+1]*lp[0]+tprob[T21][i+1]*lp[1]+tprob[T22][i+1]*lp[2]);
    //    fprintf(stderr,"logres[%d]: %f %f %f\n",i,logres[0],logres[1],logres[2]);
    //  exit(0);
    m = res[0][i];
    for(int ii=1;ii<3;ii++)
      if(res[ii][i]>m)
	m=res[ii][i];
    
    //    fprintf(stderr,"m[%d]=%f\n",i,m);
    // exit(0);
    if(std::isnan(m)||std::isnan(logres[0])){
      fprintf(stderr,"Probs in backward at site:%d\n",i);
      exit(0);
    }
  }
  for(int i=0;i<3;i++)
    logres[i] += res[i][0]+eprob[i][0];

  m=max(logres);


  loglike = m+log(exp(logres[0]-m) + exp(logres[1]-m)+ exp(logres[2]-m));
  //fprintf(stderr,"[%s]\tloglike=%f\n",__FUNCTION__,loglike);
  //  double res = m+log(exp(lp[0]-m) + exp(lp[1]-m)+ exp(lp[2]-m));
  return res;
}

double max(double a,double b,double c){
  return std::max(a,std::max(b,c));
}
int whichmax(double a,double b,double c){
  double aa[] = {a,b,c};
  int r =0;
  for(int i=1;i<3;i++)
    if(aa[i]>aa[r])
      r=i;
  return r;
}


char *viterbi(double *pars,double **eprob,double **tprob,int nSites){
  //  fprintf(stderr,"[%s] %f %f %f\n",__FUNCTION__,pars[0],pars[1],pars[2]);
  double logres[3]={log(pars[0]),log(pars[1]),log(pars[2])};
  //  fprintf(stderr,"init logres: %f %f %f\n",logres[0],logres[1],logres[2]);

  double **res = new double*[3];
  char **ptr = new char*[3];
  for(int i=0;i<3;i++){
    ptr[i] = new char[nSites];
    res[i] = new double[nSites];
    res[i][0] = logres[i]+eprob[i][0]; 
  }
  
  for(int i=1;i<nSites;i++){
    res[0][i] = eprob[0][i]+max(log(tprob[T00][i])+res[0][i-1],log(tprob[T10][i])+res[1][i-1],log(tprob[T20][i])+res[2][i-1]);
    res[1][i] = eprob[1][i]+max(log(tprob[T01][i])+res[0][i-1],log(tprob[T11][i])+res[1][i-1],log(tprob[T21][i])+res[2][i-1]);
    res[2][i] = eprob[2][i]+max(log(tprob[T02][i])+res[0][i-1],log(tprob[T12][i])+res[1][i-1],log(tprob[T22][i])+res[2][i-1]);
    ptr[0][i] =        whichmax(log(tprob[T00][i])+res[0][i-1],log(tprob[T10][i])+res[1][i-1],log(tprob[T20][i])+res[2][i-1]);
    ptr[1][i] =        whichmax(log(tprob[T01][i])+res[0][i-1],log(tprob[T11][i])+res[1][i-1],log(tprob[T21][i])+res[2][i-1]);
    ptr[2][i] =        whichmax(log(tprob[T02][i])+res[0][i-1],log(tprob[T12][i])+res[1][i-1],log(tprob[T22][i])+res[2][i-1]);
  }
  
  double loglike = max(res[0][nSites-1],res[1][nSites-1],res[2][nSites-1]);
  //  fprintf(stderr,"[%s] \tloglike=%f\n",__FUNCTION__,loglike);

  
  char *vit = new char[nSites];
  //  fprintf(stderr,"[%s] %f %f %f nsites:%d\n",__FUNCTION__,pars[0],pars[1],pars[2],nSites);
  vit[nSites-1] = whichmax(ptr[0][nSites-1],ptr[1][nSites-1],ptr[2][nSites-1]);  
  for(int i=(nSites-1);i>0;i--)
    if(vit[i]==2)
      vit[i-1] = ptr[2][i];
    else if(vit[i]==1)
      vit[i-1] = ptr[1][i];
    else 
      vit[i-1] = ptr[0][i];
  
  return vit;
}




double **post(double **fw,double **bw,int nSites,double lik){
  double **res = new double*[3];
  for(int i=0;i<3;i++){
    res[i] = new double[nSites];
  }
  for(int i=0;i<3;i++)
    for(int s=0;s<nSites;s++){
      res[i][s] = exp(fw[i][s]+bw[i][s]-lik);
      //      fprintf(stderr,"fw=%f bw=%f lik=%f\n",fw[i][s],bw[i][s],lik);
      //      exit(0);
    }
  return res;
}


void forward_backward_decode_viterbi(double *pars,genome &g){
  assert(pars[0]>0);
  if(fabs(pars[1]+pars[2]+pars[3]-1)>1e-6)
    fprintf(stderr,"alhpa=%f k0=%f k1=%f k2=%f sum:%f\n",pars[0],pars[1],pars[2],pars[3],pars[1]+pars[2]+pars[3]);
  double loglikef,loglikeb;
  loglikef=loglikeb =0;
  for(size_t i=0;i<g.results.size();i++){
    hmmRes &res = g.results[i];
    double tmp =0;
    res.forward = forward(pars+1,res.emis,res.trans,res.nSites,tmp);
    loglikef += tmp;
    tmp=0;
    res.backward = backward(pars+1,res.emis,res.trans,res.nSites,tmp);
    loglikeb += tmp;
    res.post = post(res.forward,res.backward,res.nSites,tmp);
    res.viterbi = viterbi(pars+1,res.emis,res.trans,res.nSites);
    //    fprintf(stderr,"in main:%p\n",res.viterbi);
  }
  fprintf(stderr,"loglikeforward:%f loglikebackward:%f \n",loglikef,loglikeb);
}






typedef struct{
  double **eprob;
  int nSites;
  double *pos;
  double **tprob;
  double *tsk;
}toOptim;


void norm(double *a){
  double s = a[0]+a[1]+a[2];
  for(int i=0;i<3;i++)
    a[i] /=s;

}



typedef struct toOptim2_t{
  const genome &g;
  const std::vector<perChr> &pc;
  const para &p;
  double *tsk;
  toOptim2_t(const genome &gg,const std::vector<perChr> &pcc,const para &pp) : g(gg),pc(pcc), p(pp) {};
}toOptim2;

double calcLike(double *pars,const genome &g){
  //  fprintf(stderr,"alhpa=%f k0=%f k1=%f k2=%f\n",pars[0],pars[1],pars[2],pars[3]);
  assert(pars[0]>0);
  if(fabs(pars[1]+pars[2]+pars[3]-1)>1e-6)
    fprintf(stderr,"alhpa=%f k0=%f k1=%f k2=%f sum:%e\n",pars[0],pars[1],pars[2],pars[3],pars[1]+pars[2]+pars[3]-1);

  //  assert(pars[1]+pars[2]+pars[3]==1);

  double lik[g.results.size()];
  double totLik =0;
  for(size_t i=0;i<g.results.size();i++){
    hmmRes results = g.results[i];
    trans(pars,results.dpos,results.nSites+1,results.trans);
    lik[i] = fastlike(pars+1,results.emis,results.trans,results.nSites);
    //fprintf(stderr,"i=%lu=%f\n",i,lik[i]);
    totLik += lik[i];
  }
  return totLik;
}

//case of CalcA==TRUE and k2==0. So only one free parameter
double bfgs_call_k2zero_calcA2(const double* pars,const void *dats){
  fprintf(stderr,"[%s] %f %f %f\n",__FUNCTION__,pars[0],pars[1],pars[2]);
  toOptim2 *to = (toOptim2*)dats;

  double *inV = to->tsk;
  inV[1]=pars[0];
  inV[2]=1-inV[1];
  inV[3]=0;
  inV[0]= calculateA(inV[1],inV[2],inV[3],PHI);
  if(0&&(inV[0]<alim[0]||inV[0]>alim[1])){
    fprintf(stderr,"CalcA is outside bounds\n");
    return DBL_MAX;
  }
  double lik=-calcLike(inV,to->g);
  //  fprintf(stderr,"lik=%f\n",lik);
  return lik;
}
double run_optim_k2zero_calcA2(const genome &g,const std::vector<perChr>&pc,para&p){
  toOptim2 *to = new toOptim2(g,pc,p);
  double tsk[4];//tmparray used of saving stackpointers
  to->tsk = tsk;
 
  double pars[1]={myRand(0,1)};
  //  fprintf(stderr,"pars=%f\n",pars[0]);
  double lbd[1]={klim[0]};
  double ubd[1]={klim[1]};
  int nbd[1]={2};

  double opt= -findmax_bfgs(1,pars,(void *)to, bfgs_call_k2zero_calcA2,NULL,lbd, ubd,nbd, -1);
  p.k0=pars[0];
  p.k1=1-p.k0;
  p.k2=0;
  p.a=calculateA(p.k0,p.k1,p.k2,PHI);
  return opt;
}


//pars is k0,k1
double bfgs_call_full_calcA2(const double* pars,const void *dats){
  //fprintf(stderr,"[%s] %f %f %f\n",__FUNCTION__,pars[0],pars[1],pars[2]);
  toOptim2 *to = (toOptim2*)dats;
  double *inV = to->tsk;
  if(pars[0]+pars[1]>1)
    return DBL_MAX;

  inV[1]=pars[0];
  inV[2]=pars[1];
  inV[3]=1-inV[1]-inV[2];
  inV[0]= calculateA(inV[1],inV[2],inV[3],PHI);
  if(std::isnan(inV[0]))
    return DBL_MAX;
  double lik=-calcLike(inV,to->g);
  //  fprintf(stderr,"lik=%f\n",lik);
  return lik;
}
double run_optim_full_calcA2(const genome &g,const std::vector<perChr>&pc,para&p){
  toOptim2 *to = new toOptim2(g,pc,p);
  double tsk[4];//tmparray used of saving stackpointers
  to->tsk = tsk;
  double pars[2]; pars[0]=myRand(0,1);pars[1]=myRand(0,1-pars[0]);
  //fix last pars
  double lbd[2]={klim[0],klim[0]};
  double ubd[2]={klim[1],klim[1]};
  int nbd[2]={2,2};
  double opt= -findmax_bfgs(2,pars,(void *)to, bfgs_call_full_calcA2,NULL,lbd, ubd,nbd, -1);
  p.k0=pars[0];
  p.k1=pars[1];
  p.k2=1-pars[0]-pars[1];
  p.a=calculateA(p.k0,p.k1,p.k2,PHI);
  return opt;
}

//pars is a,k0,k1
double bfgs_call_full2(const double* pars,const void *dats){
  fprintf(stderr,"[%s] %f %f %f\n",__FUNCTION__,pars[0],pars[1],pars[2]);
  toOptim2 *to = (toOptim2*)dats;
  double *inV = to->tsk;
  if(pars[1]+pars[2]>1)
    return DBL_MAX;

  inV[0]= pars[0];
  inV[1]=pars[1];
  inV[2]=pars[2];
  inV[3]=1-inV[1]-inV[2];

  if(std::isnan(inV[0]))
    return DBL_MAX;
  double lik=-calcLike(inV,to->g);
  //  fprintf(stderr,"lik=%f\n",lik);
  return lik;
}
double run_optim_full2(const genome &g,const std::vector<perChr>&pc,para&p){
  toOptim2 *to = new toOptim2(g,pc,p);
  double tsk[4];//tmparray used of saving stackpointers
  to->tsk = tsk;
  double pars[3];pars[0]=myRand(alim[0],alim[1]); pars[1]=myRand(0,1);pars[2]=myRand(0,1-pars[1]);
  //double pars[3];pars[0]=myRand(alim[0],alim[1]); pars[1]=myRand(0,1);pars[2]=1-pars[1];pars[2]=0;
  //  double pars[3]={0.045854, 0.823746,0.176254};
  //fix last pars
  double lbd[3]={alim[0],klim[0],klim[0]};
  double ubd[3]={alim[1],klim[1],klim[1]};
  int nbd[3]={2,2,2};
  double opt= -findmax_bfgs(3,pars,(void *)to, bfgs_call_full2,NULL,lbd, ubd,nbd, -1);
  p.a=pars[0];
  p.k0=pars[1];
  p.k1=pars[2];
  p.k2=1-pars[1]-pars[2];

  return opt;
}

//pars is a,k0,k1
double bfgs_call_full3(const double* pars,const void *dats){
  //fprintf(stderr,"[%s] %f %f %f\n",__FUNCTION__,pars[0],pars[1],pars[2],pars[3]);
  toOptim2 *to = (toOptim2*)dats;
  double *inV = to->tsk;
  if(pars[1]+pars[2]>1)
    return DBL_MAX;

  inV[0]= pars[0];
  inV[1]=pars[1];
  inV[2]=pars[2];
  inV[3]=pars[3];
  
  double tsum=inV[1]+inV[2]+inV[3];
  
  inV[1] /= tsum;
  inV[2] /= tsum;
  inV[3] /= tsum;


  if(std::isnan(inV[0]))
    return DBL_MAX;
  double lik=-calcLike(inV,to->g);
  //  fprintf(stderr,"lik=%f\n",lik);
  return lik;
}
double run_optim_full3(const genome &g,const std::vector<perChr>&pc,para&p){
  toOptim2 *to = new toOptim2(g,pc,p);
  double tsk[4];//tmparray used of saving stackpointers
  to->tsk = tsk;
  //  double pars[3];pars[0]=myRand(alim[0],alim[1]); pars[1]=myRand(0,1);pars[2]=myRand(0,1-pars[1]);
  double pars[4];
  pars[0]=myRand(alim[0],alim[1]);
  pars[1]=myRand(0,1);
  pars[2]=myRand(0,pars[1]);
  pars[3]=1-pars[1]-pars[2];
  //  double pars[3]={0.045854, 0.823746,0.176254};
  //fix last pars
  double lbd[4]={alim[0],klim[0],klim[0],klim[0]};
  double ubd[4]={alim[1],klim[1],klim[1],klim[1]};
  int nbd[4]={2,2,2,2};
  double opt= -findmax_bfgs(4,pars,(void *)to, bfgs_call_full3,NULL,lbd, ubd,nbd, -1);
  p.a=pars[0];
  double tsum = pars[1]+pars[2]+pars[3];
  p.k0=pars[1]/tsum;
  p.k1=pars[2]/tsum;
  p.k2=pars[3]/tsum;

  return opt;
}


//this function allocates and precomputes the data structures that doesn't depend on the a,k0,k1,k2
genome mkGenome(const std::vector<perChr> &pc,const para &p){
  genome g;
  //prep datastructures
  for(size_t i=0;i<pc.size();i++){
    hmmRes res;
    res.pos = pc[i].pos;
    res.nSites = pc[i].nSites;
    res.trans = new double*[9];
    for(int j=0;j<9;j++)
      res.trans[j] = new double[pc[i].nSites+1];
    res.emis = new double*[3];
    for(int j=0;j<3;j++)
      res.emis[j] = new double[pc[i].nSites];
    emis(p.pair[0],p.pair[1],pc[i],res.emis);

    res.dpos = diffPos(pc[i].pos,pc[i].nSites);

    g.results.push_back(res);

  }
  return g;
}

double doOptim(para &p,const genome &g,const std::vector<perChr>&pc,int calcA){
  double lik;
  srand48(0);
  
  if(calcA==1&&p.k0==-1&&p.k1==-1&&p.k2==0&&p.a==-1){
    fprintf(stderr,"1) [%s] Optimizing k0, with (k2=zero,k1=1-k0, calcA=TRUE)\n",__FUNCTION__);
    lik = run_optim_k2zero_calcA2(g,pc,p);
  }else if(calcA==1&&p.k0==-1&&p.k1==-1&&p.k2==-1&&p.a==-1){
    fprintf(stderr,"2) [%s] Optimizing k0,k1 (k2=1-k0-k1, calcA=TRUE)\n",__FUNCTION__);
    lik = run_optim_full_calcA2(g,pc,p);
  }else if(calcA==0&&p.k0==-1&&p.k1==-1&&p.k2==-1&&p.a==-1){
      fprintf(stderr,"3) [%s] Optimizing k0,k1 (k2=1-k0-k1, calcA=FALSE)\n",__FUNCTION__);
      lik = run_optim_full3(g,pc,p);
  }else{
    fprintf(stderr,"[%s] Optimization not implemented for this combination\n",__FUNCTION__);
    exit(0);
  }
  return lik;
} 





hmmRes analysis(const perChr &pc,double *freq,para p,int calcA){
  srand48(0);
  fprintf(stderr,"[%s](p0=%d,p1=%d) a=%f k0=%f k1=%f k2=%f\n",__FUNCTION__,p.pair[0],p.pair[1],p.a,p.k0,p.k1,p.k2);
  hmmRes res;
  res.pos = pc.pos;
  res.nSites = pc.nSites;
  res.trans = new double*[9];
  for(int i=0;i<9;i++)
    res.trans[i] = new double[pc.nSites+1];
  /*
  emis(p.pair[0],p.pair[1],pc,freq,&res.emis);
  double *dpos = diffPos(pc.pos,pc.nSites);
  //calc unrelated likelihood
  double pars[] = {0.08,1,0,0};
  trans(pars,dpos,pc.nSites+1,res.trans);
  res.ulike = fastlike(pars+1,res.emis,res.trans,pc.nSites);

  //calc parent offsprint
  pars[1] = 0;pars[2]=1;

  trans(pars,dpos,pc.nSites+1,res.trans);
  res.polike = fastlike(pars+1,res.emis,res.trans,pc.nSites);
  

  double lik;
  
  
  //  Check if we should perform optimizaiton

  if(p.k0==-1||p.k1==-1||p.k2==-1||p.a==-1){
    fprintf(stderr,"Performing optimation\n");
    if(calcA==0&&p.k0!=-1&&p.k1!=-1&&p.k2!=-1){
      fprintf(stderr,"optimizing alpha only\n");
      lik = run_optim_afree(dpos,res.emis,res.trans,pc.nSites,pars);
    }else if(calcA==0&&p.k0==-1&&p.k1==-1&&p.k2==-1&&p.a==-1){
      fprintf(stderr,"optimizing all parameters \n");
      lik = run_optim_full(dpos,res.emis,res.trans,pc.nSites,pars);
    }else if(calcA==0&&p.k0==-1&&p.k1==-1&&p.k2==0&&p.a==-1){
      fprintf(stderr,"optimizing a,k0,k1 with k2=zero\n");
      lik = run_optim_k2zero(dpos,res.emis,res.trans,pc.nSites,pars);
    }else if(calcA==1&&p.k0==-1&&p.k1==-1&&p.k2==0&&p.a==-1){
      fprintf(stderr,"optimizing k0, with (k2=zero,k1=1-k0, calcA=TRUE\n");
      lik = run_optim_k2zero_calcA(dpos,res.emis,res.trans,pc.nSites,pars);
    }else
      fprintf(stderr,"Optimizatino not implemented for this combination");
  }else{
    pars[0] = p.a;
    pars[1] = p.k0;
    pars[2] = p.k1;
    pars[3] = p.k2;
    fprintf(stderr,"Will perform point llh estimation\n");
  }
  fprintf(stderr,"Will do posterior decode of pars: a=%f k0=%f k1=%f k2=%f\n",pars[0],pars[1],pars[2],pars[3]);

  trans(pars,dpos,pc.nSites+1,res.trans);
  memcpy(res.pars,pars,4*sizeof(double));
  res.like = fastlike(pars+1,res.emis,res.trans,pc.nSites);
  res.forward = forward(pars+1,res.emis,res.trans,pc.nSites);
  res.backward = backward(pars+1,res.emis,res.trans,pc.nSites);
  res.post = post(res.forward,res.backward,pc.nSites,res.like);
  res.viterbi = viterbi(pars+1,res.emis,res.trans,pc.nSites);
  


  
  
  fprintf(stderr,"lik=%f a=%f k0=%f k1=%f k2=%f\n",lik,pars[0],pars[1],pars[2],pars[3]);
  if(0){
    fdump("pp.bin",res.post,3,pc.nSites);
    fdump("tprob.bin",res.trans,9,pc.nSites+1);   
    fdump("eprob.bin",res.emis,3,pc.nSites);
    fdump("forward.bin",res.forward,3,pc.nSites);
  
    fprintf(stderr,"basis lik=%f a=%f k0=%f k1=%f k2=%f\n",res.like,pars[0],pars[1],pars[2],pars[3]);

  }
  */
  return res;
}

