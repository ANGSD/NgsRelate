#include <vector>
#include <cstring>
#include <zlib.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <signal.h>
int SIG_COND=1;
int VERBOSE=1;
int really_kill=3;

//this is as close to the bound we will allow
float emTole=1e-12;


void handler(int s) {

  if(VERBOSE)
    fprintf(stderr,"\n\t-> Caught SIGNAL: Will try to exit nicely (no more threads are created.\n\t\t\t  We will wait for the current threads to finish)\n");
  
  if(--really_kill!=3)
  fprintf(stderr,"\n\t-> If you really want angsd to exit uncleanly ctrl+c: %d more times\n",really_kill+1);
  fflush(stderr);
  if(!really_kill)
    exit(0);
  VERBOSE=0;
  SIG_COND=0;
}

 //we are threading so we want make a nice signal handler for ctrl+c
void catchkill(){
  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;
  sa.sa_handler = handler;
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);  

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
  //  fprintf(stderr,"%f %f %f\n",p[0],p[1],p[2]);
  double ret =0;
  
  for(int i=0;i<len;i++){
    double tmp = 0;
    for(int j=0;j<3;j++)
      tmp += p[j]*emis[i][j];
    //    fprintf(stderr,"tmp:%f\n",log(tmp));
    //exit(0);
    ret +=log(tmp);

  }
  return ret;
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


double sumSquare(double mat[3]){
  double tmp=0;
  for(size_t i=0;i<3;i++){
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

int em1(double *sfs,double  **emis,double tole,int maxIter,int len){
  double oldLik,lik;
  oldLik = loglike(sfs,emis,len);
  fprintf(stderr,"startlik=%f %f %f %f\n",oldLik,sfs[0],sfs[1],sfs[2]);
  fflush(stderr);

  double tmp[3];
  int it;
  for(it=0;SIG_COND&&it<maxIter;it++) {
    emStep1(sfs,emis,tmp,len);
    for(int i=0;i<3;i++)
      sfs[i]= tmp[i];
    lik = loglike(sfs,emis,len);

    fprintf(stderr,"[%d] lik=%f diff=%g\n",it,lik,fabs(lik-oldLik));

    if(fabs(lik-oldLik)<tole){
     
      oldLik=lik;
      break;
    }
    oldLik=lik;
  }
  return it;
}


int em2(double *sfs,double  **emis,double tole,int maxIter,int len){
  double oldLik,lik;
  oldLik = loglike(sfs,emis,len);
  fprintf(stderr,"startlik=%f\n",oldLik);
  fflush(stderr);

  double tmp[3];
  int it;
  for(it=0;SIG_COND&&it<maxIter;it++) {
    emAccel(sfs,emis,tmp,len);

    for(int i=0;i<3;i++)
      sfs[i]= tmp[i];
    lik = loglike(sfs,emis,len);

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

int em3(double *sfs,double  **emis,double tole,int maxIter,int len){
  //  exit(0);
  double oldLik,lik;
  oldLik = loglike(sfs,emis,len);
  fprintf(stderr,"em3startlik=%f\n",oldLik);
  fflush(stderr);

  double tmp[3];
  int it;
  int speedy=1;
  for(it=0;SIG_COND&&it<maxIter;it++) {
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





void emission_ngsrelate(std::vector<double> &freq,double **l1,double **l2,double **emis){

  for(int i=0;i<freq.size();i++){
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

      

      //      fprintf(stderr,"%f %f %f\n",emis[0][0],emis[0][1],emis[0][2]);
      //      exit(0);
    }

}


 double **getGL(const char *fname,int sites, int ngl){
   gzFile gz=Z_NULL;
   if(((gz=gzopen(fname,"rb")))==Z_NULL){
     fprintf(stderr,"Problem opening file:%s\n",fname);
     exit(0);
   }
   

   double **ret=new double*[sites];
   
   for(int i=0;i<sites;i++){
     ret[i] = new double[6];
     gzread(gz,ret[i],sizeof(double)*6);
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
       ret.push_back(atof(tok));

     }

   }
   fprintf(stderr,"frequency file:%s contain %lu number of sites\n",fname,ret.size());
   gzclose(gz);
   delete [] buf;
   return ret;
 }

void print_info(FILE *fp){
  fprintf(fp, "\n");
  fprintf(fp, "Usage: da  [options] \n");
  fprintf(fp, "Options:\n");
  fprintf(fp, "   -p <filename>       plink prefix filename\n");
  fprintf(fp, "   -o <filename>       outputfilename\n");
  fprintf(fp, "   -f <filename>       freqs\n");
  fprintf(fp, "   -m <INTEGER>        model 0=EMnormal 1=accelerated em\n");
  fprintf(fp, "   -b <filename>       file containing the start NI\n");
  fprintf(fp, "   -i <UINTEGER>       maxIter\n");
  fprintf(fp, "   -t <FLOAT>          tolerance for breaking EM\n");
  fprintf(fp, "   -r <FLOAT>          seed for rand\n");
  fprintf(fp, "   -g gfile            genotypellh file\n");
  fprintf(fp, "   -P <INT>            nThreads\n");
  fprintf(fp, "   -c <INT>            should call genotypes instead?");
  fprintf(fp, "   -e <INT>            errorrates when calling genotypes?");
  fprintf(fp, "\n");
  exit(0);
}

void callgenotypes(double **gls,int len,double eps){
  fprintf(stderr,"callgentypes\n");
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


int main(int argc, char *argv[]){
  if(argc==1)
    print_info(stderr);
  char *pname = NULL;
  char *outname = NULL;
  char *freqname=NULL;
  char *gname=NULL;
  char *startk = NULL;
  int maxIter =10;
  double tole =1e-3;
  int n;
  int seed=100;
  int nThreads = 1;
  int model =0;
  int gc =0;
  double errate = 0.005;

  while ((n = getopt(argc, argv, "p:o:f:i:t:r:P:g:m:c:e:")) >= 0) {
    switch (n) {
    case 'p': pname = strdup(optarg); break; 
    case 'o': outname = strdup(optarg); break;
    case 'f': freqname = strdup(optarg); break;
    case 'i': maxIter = atoi(optarg); break;
    case 't': tole = atof(optarg); break;
    case 'r': seed = atoi(optarg); break;
    case 'P': nThreads = atoi(optarg); break;
    case 'g': gname = strdup(optarg); break;
    case 'b': startk = strdup(optarg); break;      
    case 'm': model = atoi(optarg); break;      
    case 'c': gc = atoi(optarg); break;      
    case 'e': errate = atof(optarg); break;      
    default: {fprintf(stderr,"unknown arg:\n");return 0;}
      print_info(stderr);
    }
  }
  srand48(seed);
  catchkill();

  std::vector<double> freq = getDouble(freqname);

  double **gls = getGL(gname,freq.size(),6);
  // print(stdout,freq.size(),6,gls);
  double **l1,**l2;l1=l2=NULL;
  l1=new double*[freq.size()];
  l2=new double*[freq.size()];
  char *keep=new char[freq.size()];
  memset(keep,1,freq.size());
  double **emis=new double*[freq.size()];
  int nkeep=0;
  for(int i=0;i<freq.size();i++){
    l1[i] = new double[3];
    l2[i] = new double[3];
    emis[i] = new double[3];
    for(int j=0;j<3;j++){
      l1[i][j] = gls[i][j];
      l2[i][j] = gls[i][j+3];
    }
    
    //skip of both samples are missing
    if(l1[i][0]==l1[i][1]&&l1[i][0]!=l1[i][2])
      keep[i]=0;
    if(l2[i][0]==l2[i][1]&&l2[i][0]!=l2[i][2])
      keep[i]=0;
    nkeep += keep[i];
  }
  //print(stdout,freq.size(),6,gls);
  //return 0;
  fprintf(stderr,"fraction of missing:%f\n",nkeep/(1.0*freq.size()));
  if(gc){
    callgenotypes(l1,freq.size(),errate);
    callgenotypes(l2,freq.size(),errate);
  }
  emission_ngsrelate(freq,l1,l2,emis);
  //  print(stdout,freq.size(),3,emis);return 0;
  double pars[3];
  pars[0] = drand48();
  pars[1] = drand48()*(1-pars[0]);
  pars[2] = 1-pars[0]-pars[1];
  fprintf(stderr,"loglike:%f\n",loglike(pars,emis,freq.size()));
  int niter;
  if(model==0)
    niter=em1(pars,emis,tole,maxIter,freq.size());
  else if(model==1)
    niter=em2(pars,emis,tole,maxIter,freq.size());
  else//below might not work
    niter=em3(pars,emis,tole,maxIter,freq.size());

  double p100[3]={1,0,0};
  double p010[3]={0,1,0};
  double p001[3]={0,0,1};
  double l100=loglike(p100,emis,freq.size());
  double l010=loglike(p010,emis,freq.size());
  double l001=loglike(p001,emis,freq.size());
  double lopt= loglike(pars,emis,freq.size());
  fprintf(stderr,"%f %f %f\n",l100,l010,l001);

  double likes[4] = {l100,l010,l001,lopt};
  int best = 0;
  for(int i=0;i<4;i++)
    if(likes[i]>likes[best])
      best=i;
  
  if(best==3)
    fprintf(stdout,"%f\t%f\t%f\t%f\t%d\n",pars[0],pars[1],pars[2],lopt,niter);
  if(best==0)
    fprintf(stdout,"%f\t%f\t%f\t%f\t%d\n",p100[0],p100[1],p100[2],l100,-1);
  if(best==1)
    fprintf(stdout,"%f\t%f\t%f\t%f\t%d\n",p010[0],p010[1],p010[2],l010,-1);
  if(best==2)
    fprintf(stdout,"%f\t%f\t%f\t%f\t%d\n",p001[0],p001[1],p001[2],l001,-1);
  
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
