/*

  write stuff

*/


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>
#include <limits>
#include <zlib.h>
#include <vector>
#include <pthread.h>
#include <signal.h>
#include <vector>
#include <sys/stat.h>
#include <cassert>

#include "analysis.h"

#define MAF_START 0.3
#define MAF_ITER 20

double misTol=0.05;



#define LENS 100000 //this is the max number of bytes perline, should make bigger
int SIG_COND =1;//if we catch signal then quit program nicely
 
int seed =0;//<- do more clever


size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}


int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}

template <typename T>
T max(T* a,size_t l){
  if(l==0)
    return 0;
  size_t m=0;
  for(size_t i=1;i<l;i++)
    if(a[i]>a[m])
      m=i;
  return a[m];
}
template <typename T>
T min(T* a,size_t l){
  if(l==0)
    return 0;
  size_t m=0;
  for(size_t i=1;i<l;i++)
    if(a[i]<a[m])
      m=i;
  return a[m];
}




std::vector<char *> dumpedFiles;
FILE *openFile(const char* a,const char* b){
  if(0)
    fprintf(stderr,"[%s] %s %s",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  //  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  if(0&&fexists(c)){//ANDERS DAEMON DRAGON HATES THIS
    fprintf(stderr,"File: %s exists will exist\n",c);
    fflush(stderr);
    exit(0);
  }
  dumpedFiles.push_back(strdup(c));
  FILE *fp = fopen(c,"w");
  delete [] c;
  return fp;
}

gzFile getGz(const char*fname,const char* mode){

  //  fprintf(stderr,"\t-> opening: %s\n",fname);
  gzFile fp=Z_NULL;
  if(NULL==(fp=gzopen(fname,mode))){
    fprintf(stderr,"\t-> Error opening gzFile handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}

gzFile openFileGz(const char* a,const char* b,const char *mode){
  if(0)
    fprintf(stderr,"[%s] %s%s\n",__FUNCTION__,a,b);
  char *c = new char[strlen(a)+strlen(b)+1];
  strcpy(c,a);
  strncat(c,b,strlen(b));
  //  fprintf(stderr,"\t-> Dumping file: %s\n",c);
  dumpedFiles.push_back(strdup(c));
  gzFile fp = getGz(c,mode);
  delete [] c;
  return fp;
}


//some struct will all the data from the beagle file
typedef struct{
  double **genos;
  char **chr;
  int *pos;
  int nSites;
  int nInd;
}bgl;

//utility function for cleaning up out datastruct
void dalloc(bgl &b){
  for(int i=0;i<b.nSites;i++){
    delete [] b.genos[i];
    free(b.chr[i]);
  }
  delete [] b.genos;
  delete [] b.chr;
}


double addProtect2(double a,double b){
  //function does: log(exp(a)+exp(b)) while protecting for underflow
  if(std::isinf(a)&&std::isinf(b))
    return log(0);
  if(std::isinf(a))
    return b;
  if(std::isinf(b))
    return a;
  double maxVal;// = std::max(a,b));
  if(a>b)
    maxVal=a;

  else
    maxVal=b;
  double sumVal = exp(a-maxVal)+exp(b-maxVal);
  return log(sumVal) + maxVal;
}


double addProtect3(double a,double b, double c){
  //function does: log(exp(a)+exp(b)+exp(c)) while protecting for underflow
  double maxVal;// = std::max(a,std::max(b,c));
  if(a>b&&a>c)
    maxVal=a;
  else if(b>c)
    maxVal=b;
  else
    maxVal=c;
  double sumVal = exp(a-maxVal)+exp(b-maxVal)+exp(c-maxVal);
  return log(sumVal) + maxVal;
}




double emFrequency(double *loglike,int numInds, int iter,double start,char *keep,int keepInd){

  if(keepInd == 0)
    return 0.0;
  
  double W0;
  double W1;
  double W2;
  // fprintf(stderr,"start=%f\n",start);
  double p=(double)start;
  double temp_p=(double)start;
  double accu=0.00001;
  double accu2=0;
  double sum;


  int it=0;
  
  for(it=0;it<iter;it++){
    sum= 0;
    double pl = log(p);
    double mpl = log(1-p);
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=loglike[i*3+0]+2*mpl;
      W1=loglike[i*3+1]+M_LN2+pl+mpl;
      W2=loglike[i*3+2]+2*pl;
      sum+= exp(addProtect2(W1,M_LN2+W2)-M_LN2-addProtect3(W0,W1,W2));
      //  fprintf(stderr,"%f %f %f\n",W0,W1,W2);
      if(0&&isnan(sum)){
	//fprintf(stderr,"PRE[%d]: W %f\t%f\t%f sum=%f\n",i,W0,W1,W2,sum);
	exit(0);
      }
    }

    p=sum/keepInd;
    // fprintf(stderr,"it=%d\tp=%f\tsum=%f\tkeepInd=%d\n",it,p,log(sum),keepInd);
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
    temp_p=p;
  }



  if(isnan(p)){
    fprintf(stderr,"[emFrequency] caught nan will not exit\n");
    fprintf(stderr,"logLike (3*nInd). nInd=%d\n",numInds);
    //print_array(stderr,loglike,3*numInds);
    fprintf(stderr,"keepList (nInd)\n");
    //print_array(stderr,keep,numInds);
    fprintf(stderr,"used logLike (3*length(keep))=%d\n",keepInd);

    for(int ii=0;0&&ii<numInds;ii++){
      if(keep!=NULL && keep[ii]==1)
	    fprintf(stderr,"1\t");
	for(int gg=0;gg<3;gg++)
	  fprintf(stderr,"%f\t",loglike[ii*3+gg]);
      fprintf(stderr,"\n");
    }
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=(loglike[i*3+0])*pow(1-p,2);
      W1=(loglike[i*3+1])*2*p*(1-p);
      W2=(loglike[i*3+2])*(pow(p,2));
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      //fprintf(stderr,"p=%f W %f\t%f\t%f sum=%f loglike: %f\n",p,W0,W1,W2,sum,exp(loglike[i*3+2])*pow(1-p,2));
    }
    p=-999;
    // exit(0);
  }
  
  return(p);
}



/*
  Returns the bgl struct containing all data from a beagle file.

  It find the nsamples from counting the header
  It finds the number of sites by queing every line in a std::vector
  After the file has been read intotal it reloops over the lines in the vector and parses data
 */

bgl readBeagle(const char* fname) {
  const char *delims = "\t \n";
  gzFile fp = NULL;
  if(Z_NULL==(fp=gzopen(fname,"r"))){
    fprintf(stderr,"Error opening file: %s\n",fname);
    exit(0);
  }
  
  bgl ret;
  char buf[LENS];



  //find number of columns
  gzgets(fp,buf,LENS);
  strtok(buf,delims);
  int ncols=1;
  while(strtok(NULL,delims))
    ncols++;
  if(0!=(ncols-3) %3 ){
    fprintf(stderr,"Problem parsing beagle file: ncols=%d will exit\n",ncols);
    exit(0);
  }
  ret.nInd = (ncols-3)/3;//this is the number of samples
  
  //read every line into a vector
  std::vector<char*> tmp;
  while(gzgets(fp,buf,LENS))
    tmp.push_back(strdup(buf));
  
  fprintf(stderr,"[%s] NumSites from beagle file=%lu will roll back and parse data\n",__FUNCTION__,tmp.size());
  //now we now the number of sites
  ret.nSites=tmp.size();
  ret.chr = new char*[ret.nSites];
  ret.pos = new int[ret.nSites];
  ret.genos= new double*[ret.nSites];

  //then loop over the vector and parsing every line
  for(int s=0;SIG_COND&& (s<ret.nSites);s++){
    ret.chr[s] = strdup(strtok(tmp[s],"_"));
    ret.pos[s] = atoi(strtok(NULL,delims));
    strtok(NULL,delims);//major
    strtok(NULL,delims);//minor
    ret.genos[s]= new double[3*ret.nInd];
    for(int i=0;i<ret.nInd*3;i++){
      ret.genos[s][i] = atof(strtok(NULL,delims));
      if(ret.genos[s][i]<0){
	fprintf(stderr,"Likelihoods must be positive\n");
	fprintf(stderr,"site %d ind %d geno %d has value %f\n",s,int(i*1.0/3),i%3,ret.genos[s][i]);
	exit(0);
      }
    }
    for(int i=0;i<ret.nInd;i++){
      double tmpS = 0.0;
      for(int g=0;g<3;g++)
	tmpS += ret.genos[s][i*3+g];
      if(!(tmpS>0)){
	fprintf(stderr,"The sum of likelihoods for a genotypes must be positive\n");
	fprintf(stderr,"individual %d site %d has sum %f\n",i,s,tmpS);
	exit(0);
      } 
      for(int g=0;1&&g<3;g++)
	ret.genos[s][i*3+g] = log(ret.genos[s][i*3+g]);
    }
    free(tmp[s]);

  }


  //  keeps=ret.keeps;
  gzclose(fp); //clean up filepointer
  return ret;
}

FILE *getFILE(const char*fname,const char* mode){

  FILE *fp;
  if(NULL==(fp=fopen(fname,mode))){
    fprintf(stderr,"\t->Error opening FILE handle for file:%s exiting\n",fname);
    exit(0);
  }
  return fp;
}

double **allocDouble(size_t x,size_t y){
  double **ret= new double*[x];
  for(size_t i=0;i<x;i++)
    ret[i] = new double[y];
  return ret;
}


void dalloc(size_t x,double **ret){
  for(size_t i=0;i<x;i++)
    delete [] ret[i] ;
  delete [] ret;
}
//same as above but swapped.... to lazy to change code
void dalloc(double **ret,size_t x){
  for(size_t i=0;i<x;i++)
    delete [] ret[i] ;
  delete [] ret;
}

void printDouble(double **ret,size_t x,size_t y,FILE *fp){
  for(size_t i=0;i<x;i++){
    for(size_t j=0;j<y;j++)
      fprintf(fp,"%.20f ",ret[i][j]);
    fprintf(fp,"\n");
  }
}

int printer =0;


void info(){
  fprintf(stderr,"-beagle/-bin -freqfile -outnames -a -k0 -k1 -k2 -calcA -outnames -pair1 pair2 -selectChr\n");
  exit(0);
}

typedef struct{
  const char *beagle;
  const char *bin;
  const char *freqfile;
  const char *outnames;
  para p;
  int calcA;
  int selectChr;
}cArg;

void printArg(FILE *fp,const cArg &ca){
  fprintf(fp,"-----------------\n%s\n\t-beagle\t\t%s\n",__FILE__,ca.beagle);
  fprintf(fp,"\t-bin\t\t%s\n",ca.bin);
  fprintf(fp,"\t-freqfile\t%s\n",ca.freqfile);
  fprintf(fp,"\t-outnames\t%s\n",ca.outnames);
  fprintf(fp,"\t-calcA\t\t%d\n",ca.calcA);
  printPars(fp,ca.p);
}



int VERBOSE =1;
void handler(int s) {
  if(VERBOSE)
    fprintf(stderr,"Caught SIGNAL: %d will try to exit nicely (no more threads are created, we will wait for the current threads to finish)\n",s);
  VERBOSE=0;
  SIG_COND=0;
}

void marshall_dump(gzFile fp,const perChr& d){
  int strL = strlen(d.name);
  assert(sizeof(int)==gzwrite(fp,&strL,sizeof(int)));
  assert(strL==gzwrite(fp,d.name,strL));
  assert(sizeof(int)==gzwrite(fp,&d.nInd,sizeof(int)));
  assert(sizeof(int)==gzwrite(fp,&d.nSites,sizeof(int)));
  assert((int)(d.nSites*sizeof(int))==gzwrite(fp,d.pos,sizeof(int)*d.nSites));
  for(int i=0;i<d.nInd;i++)
    assert((int)(3*d.nSites*sizeof(double))==gzwrite(fp,d.gl[i],3*sizeof(double)*d.nSites));
}
/*
void fdump(gzFile fp, const hmm&res,const char *chr){
  int strL = strlen(chr);
  assert(sizeof(int)==gzwrite(fp,&strL,sizeof(int)));
  assert(strL==gzwrite(fp,chr,strL));
  gzputc(fp,'\0');
  assert(sizeof(int)==gzwrite(fp,&res.nSites,sizeof(int)));
  assert(4*sizeof(double)==gzwrite(fp,res.pars,4*sizeof(double)));
  assert(sizeof(double)==gzwrite(fp,&res.like,sizeof(double)));
  assert(sizeof(double)==gzwrite(fp,&res.ulike,sizeof(double)));
  assert(sizeof(double)==gzwrite(fp,&res.polike,sizeof(double)));
  fprintf(stderr,"res.nSites=%d\n",res.nSites);
  assert(res.nSites*sizeof(int)==gzwrite(fp,res.pos,sizeof(int)*res.nSites));
  for(int i=0;i<3;i++)
    assert(res.nSites*sizeof(double)==gzwrite(fp,res.post[i],sizeof(double)*res.nSites));
  for(int i=0;i<3;i++)
    assert(res.nSites*sizeof(double)==gzwrite(fp,res.forward[i],sizeof(double)*res.nSites));
  for(int i=0;i<3;i++)
    assert(res.nSites*sizeof(double)==gzwrite(fp,res.backward[i],sizeof(double)*res.nSites));
  for(int i=0;i<3;i++)
    assert(res.nSites*sizeof(double)==gzwrite(fp,res.emis[i],sizeof(double)*res.nSites));
  for(int i=0;i<9;i++)
    assert((1+res.nSites)*sizeof(double)==gzwrite(fp,res.trans[i],sizeof(double)*(1+res.nSites)));
  assert(res.nSites==gzwrite(fp,res.viterbi,res.nSites));
}

*/
std::vector<perChr> marshall_read(gzFile fp){
  fprintf(stderr,"Reading binary:\n");
  std::vector<perChr> ret;
 int tmp;
  while(gzread(fp,&tmp,sizeof(int))){
    perChr pc;
    pc.name =(char*) malloc(tmp+1);
    assert(tmp==gzread(fp,pc.name,tmp));
    pc.name[tmp]='\0';
    assert(sizeof(int)==gzread(fp,&pc.nInd,sizeof(int)));
    assert(sizeof(int)==gzread(fp,&pc.nSites,sizeof(int)));
    //    fprintf(stderr,"nind=%d nsites=%d\n",pc.nInd,pc.nSites);
    pc.pos = new int[pc.nSites];
    pc.gl = new double*[pc.nSites];
    
    assert(pc.nSites*sizeof(int)==gzread(fp,pc.pos,sizeof(int)*pc.nSites));
    for(int i=0;i<pc.nInd;i++){
      pc.gl[i] = new double [3*pc.nSites];
      assert(3*pc.nSites*sizeof(double)==gzread(fp,pc.gl[i],3*sizeof(double)*pc.nSites));
    }
    ret.push_back(pc);
  }
  
  return ret;
}

/*
 char *keep = new char[in.nSites];
  memset(keep,0,in.nSites);
  for(int i=0;i<in.nSites;i++){
    if(in.mafs[i]>minMaf)
      keep[i] =1;
  }
*/

void printStuff(const std::vector<perChr>& pc ){
  for(uint i=0;i<pc.size();i++)
    fprintf(stderr,"i=%d chr=%s nSites=%d (%d,%d) (%f,%f)\n",i,pc[i].name,pc[i].nSites,pc[i].pos[0],pc[i].pos[pc[i].nSites-1],min<double>(pc[i].freq,pc[i].nSites),max<double>(pc[i].freq,pc[i].nSites));


}

template <typename T>
void transpose(T***in,size_t in_x,size_t in_y){
  fprintf(stderr,"inx=%lu iny=%lu\n",in_x,in_y);
  T **tmp = new T*[in_y];
  for(int i=0;i<in_y;i++)
    tmp[i] = new T [in_x];
  
  for(int i=0;i<in_y;i++)
    for(int j=0;j<in_x;j++)
      tmp[i][j] = (*in)[j][i];
  
  (*in) = tmp;
}


template <typename T>
void transpose2(T***in,size_t in_x,size_t in_y){
  int nInd = in_y;
  int nSites = in_x;

  T **tmp = new T*[nInd];
  for(size_t i=0;i<in_y;i++)
    tmp[i] = new T [3*nSites];
  
  for(int i=0;i<nInd;i++)
    for(int j=0;j<nSites;j++)
      for(int o=0;o<3;o++)
	tmp[i][j*3+o] = (*in)[j][i*3+o];
  
  (*in) = tmp;
}

 
void setFreq(perChr &pc){

  pc.keeps=new char*[pc.nSites]; // array nSites x nInd 0 if missing info
  pc.keepInd = new int[pc.nSites];
  
  for(int s=0;s<pc.nSites;s++) {
    pc.keeps[s] = new char[pc.nInd];
    int nKeep =0;
    for(int i=0;i<pc.nInd;i++){
      //check if llh has enough information
      double mmin=std::min(pc.gl[i][s*3],std::min(pc.gl[i][s*3+1],pc.gl[i][s*3+2]));
      double mmax=std::max(pc.gl[i][s*3],std::max(pc.gl[i][s*3+1],pc.gl[i][s*3+2]));
      if(fabs(mmax-mmin)<misTol){
	//	fprintf(stderr,"site skipped from analysis\n");
	pc.keeps[s][i] =0;
      }else{
	pc.keeps[s][i] =1;
	nKeep++;
      }
    }
    pc.keepInd[s] = nKeep;
  }
 
  double tmp[3*pc.nInd];
  pc.freq = new double[ pc.nSites];
  pc.qerf = new double [pc.nSites];
  for(int s=0;s<pc.nSites;s++){
    for(int i=0;i<pc.nInd;i++)
      for(int o=0;o<3;o++)
	tmp[i*3+o] = pc.gl[i][s*3+o];
    double af = emFrequency(tmp,pc.nInd,MAF_ITER,MAF_START,pc.keeps[s],pc.keepInd[s]);
    assert(af<=1 && af>=0);
    pc.freq[s] = log(af);
    pc.qerf[s] = log(1-af);
  }
  
}



std::vector<perChr> makeDat(const bgl& in){
  std::vector<int> splits;
  int start=0;
  int stop=1;
  for(int i=1;i<in.nSites;i++){
    if(strcmp(in.chr[i],in.chr[i-1])!=0)
      splits.push_back(i);
  }
  splits.push_back(in.nSites);
  fprintf(stderr,"[%s] We observe: %lu number of chromosomes\n",__FUNCTION__,splits.size());
  std::vector<perChr> ret;
  
  for(size_t i=0;i<splits.size();i++){
    //    fprintf(stderr,"i=%d val=%d\n",i,splits[i]);
    perChr pc;
    pc.name = strdup(in.chr[splits[i]-1]);
    pc.nInd=in.nInd;
    pc.nSites = i==0?splits[i]:splits[i]-splits[i-1];
    //    fprintf(stderr,"[%s] nsites=%d\n",__FUNCTION__,pc.nSites);
    pc.pos = new int[pc.nSites];
    pc.gl = new double*[pc.nSites];
    pc.keeps = new char*[pc.nSites];
    pc.keepInd = new int [pc.nSites];
    start = i!=0?splits[i-1]:0;
    stop = splits[i];
    int cnt =0;
    for(int j=start;j<stop;j++){
      pc.pos[cnt] = in.pos[j];
      pc.gl[cnt] = in.genos[j];
      cnt++;
    }

    transpose2<double>(&pc.gl,pc.nSites,pc.nInd);
    //    fprintf(stderr,"LLLH: %f %f %f\n",pc.gl[0][0],pc.gl[1][0],pc.gl[2][0]);
    //    exit(0);
    //    transpose<char>(&pc.keeps,pc.nSites,pc.nInd);
    ret.push_back(pc);
  }
  
  return ret;
}


double *readDouble(const char*fname,size_t hint){
  FILE *fp = NULL;
    fp = getFILE(fname,"r");
  assert(fp!=NULL);
  char buf[fsize(fname)+1];
  if(fsize(fname)!=fread(buf,sizeof(char),fsize(fname),fp)){
    fprintf(stderr,"Problems reading file: %s\n will exit\n",fname);
    exit(0);
  }
  buf[fsize(fname)]='\0';
  std::vector<double> res;
  res.push_back(atof(strtok(buf,"\t\n ")));
  char *tok=NULL;
  while((tok=strtok(NULL,"\t\n "))) {  
    //fprintf(stderr,"%s\n",tok);
    res.push_back(atof(tok));

  }
  //  fprintf(stderr,"size of prior=%lu\n",res.size());
  if(hint!=res.size()){
    fprintf(stderr,"problem with size of dimension of prior %lu vs %lu\n",hint,res.size());
    for(uint i=0;0&&i<res.size();i++)
      fprintf(stderr,"%d=%f\n",i,res[i]);
    exit(0);
  }
  double *ret = new double[res.size()];
  for(uint i=0;i<res.size();i++)
    ret[i] = log(res[i]);
  fclose(fp);
  return ret;
}


int main(int argc, char **argv){

  fprintf(stderr,"NGSrelate buildtime: (%s:%s)\n",__DATE__,__TIME__);
  if(argc==1){// if no arguments, print info on program
    info();
    return 0;
  }
  //below for catching ctrl+c, and dumping files
  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = 0;
  sa.sa_handler = handler;
  sigaction(SIGPIPE, &sa, 0);
  sigaction(SIGINT, &sa, 0);  

  //set default parameters
  cArg ca;
  ca.beagle=NULL;
  ca.bin =NULL;
  ca.freqfile=NULL;
  ca.outnames=NULL;
  ca.p.a=ca.p.k0=ca.p.k1=ca.p.k2=-1;
  ca.calcA=1;
  ca.p.pair[0]=0;
  ca.p.pair[1]=1;
  ca.selectChr =-1;


  // reading arguments
  argv++;
  while(*argv){
    if(strcmp(*argv,"-beagle")==0 )  ca.beagle=*++argv; 
    else if(strcmp(*argv,"-bin")==0 )  ca.bin=*++argv; 
    else if(strcmp(*argv,"-freqfile")==0 )  ca.freqfile=*++argv; 
    else if(strcmp(*argv,"-outnames")==0 )  ca.outnames=*++argv; 
    //else if(strcmp(*argv,"-minMaf")==0 )  minMaf = atoi(*++argv); 
    else if(strcmp(*argv,"-a")==0 )  ca.p.a = atof(*++argv); 
    else if(strcmp(*argv,"-k0")==0 )  ca.p.k0 = atof(*++argv); 
    else if(strcmp(*argv,"-k1")==0 )  ca.p.k1 = atof(*++argv); 
    else if(strcmp(*argv,"-k2")==0 )  ca.p.k2 = atof(*++argv); 
    else if(strcmp(*argv,"-calcA")==0 )  ca.calcA = atoi(*++argv);
    else if(strcmp(*argv,"-pair1")==0 )  ca.p.pair[0] = atoi(*++argv); 
    else if(strcmp(*argv,"-pair2")==0 )  ca.p.pair[1] = atoi(*++argv); 
    else if(strcmp(*argv,"-selectChr")==0 )  ca.selectChr = atoi(*++argv); 
    else{
      fprintf(stderr,"Unknown arg:%s\n",*argv);
      info();
      return 0;
    }
    ++argv;
  }
  if(ca.beagle==NULL&&ca.bin==NULL){
    fprintf(stderr,"Please supply input data file: -beagle OR -bin");
    info();
  }else if(ca.beagle!=NULL&&ca.bin!=NULL){
    fprintf(stderr,"Please supply input data file: -beagle OR -bin");
    info();
  }
  if(ca.outnames==NULL){
    if(ca.beagle!=NULL)
      ca.outnames=ca.beagle;
    else
      ca.outnames = ca.bin;
    fprintf(stderr,"Will use: %s as prefix for output\n",ca.outnames);
   
  }
  FILE *flog=openFile(ca.outnames,".log");
  for(int i=0;i<argc;i++){
    //    fprintf(flog,"%s ",argv[i]);
    //fprintf(stderr,"%s ",argv[i]);
  }
  fprintf(stderr,"\n");fprintf(flog,"\n");
  
  clock_t t=clock();//how long time does the run take
  time_t t2=time(NULL);
  
  std::vector<perChr> pd;
  if(ca.beagle!=NULL){
    bgl d=readBeagle(ca.beagle);
    fprintf(stderr,"Input beaglefile has dim: nsites=%d nind=%d\n",d.nSites,d.nInd);
    pd = makeDat(d);
    gzFile dfile = openFileGz(ca.outnames,".bin.gz","wb");
    for(uint i=0;i<pd.size();i++)
      marshall_dump(dfile,pd[i]);
    gzclose(dfile);
  }else{
    gzFile dfile = getGz(ca.bin,"rb");
    pd = marshall_read(dfile);

    gzclose(dfile);
  }
  fprintf(stderr,"number of chrs:%lu\n",pd.size());
  //calculate GL frequencies
  int totSites = 0;
  for(uint i=0;i<pd.size();i++){
    totSites += pd[i].nSites;
    setFreq(pd[i]);
  }

  fprintf(stderr,"totSites=%d\n",totSites);
  //  printStuff(pd);
  //set seed
  srand(seed);
  
  //lets print the freqs;
  FILE *fp =openFile(ca.outnames,".GL.freq");
  for(size_t j=0;j<pd.size();j++)
    for(int i=0;i<pd[j].nSites;i++)
      fprintf(fp,"%f ",exp(pd[j].freq[i]));
  fclose(fp);
  
  
  double *freq = NULL;
  if(ca.freqfile!=NULL){
    freq = readDouble(ca.freqfile,totSites);
    //    fprintf(stderr,"freq=%f\n",freq[0]);
    int cnt =0;
    for(size_t j=0;j<pd.size();j++)
      for(int i=0;i<pd[j].nSites;i++){
	pd[j].freq[i] = freq[cnt++];
	pd[j].qerf[i] =log(1-exp(pd[j].freq[i]));
      }
    delete [] freq ;
    freq=NULL;
  }
  //now over vector containing the input data is done,
  //lets now allocate the datastructures we need.


  genome mkGenome(const std::vector<perChr> &pd,const para &p);
  double calcLike(double *pars,const genome &g);
  double doOptim(para &p,const genome &g,const std::vector<perChr>&pc,int calcA);
  void forward_backward_decode_viterbi(double *pars,genome &g);

  //if we want to run analysis on a single chr.
  if(ca.selectChr!=-1){
    assert(ca.selectChr>=-1 && ca.selectChr<(int)pd.size());
    perChr tmp = pd[ca.selectChr];
    pd.clear();
    pd.push_back(tmp);
  }
  printArg(stderr,ca);
  printArg(flog,ca);
  
  genome g = mkGenome(pd,ca.p);
  double lik;
  if(ca.p.k0!=-1&&ca.p.k1!=-1&&ca.p.k2!=-1){
    if(ca.calcA==1 && ca.p.a==-1){
      fprintf(stderr,"Calculating point estimate of llh\n");
      fprintf(stderr,"p.a=%f\n",ca.p.a);
      ca.p.a=calculateA(ca.p.k0,ca.p.k1,ca.p.k2,PHI);
      fprintf(stderr,"p.a=%f\n",ca.p.a);
      double pars[] = {ca.p.a,ca.p.k0,ca.p.k1,ca.p.k2};
      lik = calcLike(pars,g);
    }
  }else{
    fprintf(stderr,"Will now run optimizaton\n");
    lik = doOptim(ca.p,g,pd,ca.calcA);
  }
  fprintf(stderr,"%f %f %f %f %f\n",ca.p.a,ca.p.k0,ca.p.k1,ca.p.k2,lik);
  FILE *fres=openFile(ca.outnames,".res");
  fprintf(fres,"%f %f %f %f %f\n",ca.p.a,ca.p.k0,ca.p.k1,ca.p.k2,lik);
  fclose(fres);

  //now do the forward,backward, decode and viterbi
  double pars[] = {ca.p.a,ca.p.k0,ca.p.k1,ca.p.k2};
  forward_backward_decode_viterbi(pars,g);
  pars[1] = 1; pars[2] = pars[3] = 0;
  g.ulike = calcLike(pars,g);
  pars[1] = 0; pars[2] = 1;
  g.polike = calcLike(pars,g);
  fprintf(stderr,"ulike:%f polike:%f\n",g.ulike,g.polike);
  
  gzFile add = openFileGz(ca.outnames,".all.gz","wb");
  for(int i=0;i<pd.size();i++){
    perChr en = pd[i];
    hmmRes to = g.results[i];
    assert(en.nSites==to.nSites);
    for(int j=0;j<en.nSites;j++){
      gzprintf(add,"%s\t%d\t%f\t",en.name,to.pos[j],en.freq[j]);
      gzprintf(add,"%d\t",to.viterbi[j]);
      gzprintf(add,"%f\t%f\t%f\n",to.post[0][j],to.post[1][j],to.post[2][j]);
    }
  }
  gzclose(add);
  /*

  std::vector<hmm> all_res;
  for(size_t i =0;i<pd.size();i++){
    hmm res = analysis(pd[i],freq,p,calcA);
    all_res.push_back(res);
    char cname[1024];
    sprintf(cname,"%s.%lu.bres.gz",outfiles,i);
    gzFile bo = openFileGz(cname,".bres.gz","wb");
    //    fdump(bo,res,pd[i].name);
    gzclose(bo);
  }
  
  FILE *fres=openFile(outfiles,".res");
  for(size_t i =0;i<all_res.size();i++){
    double *p = all_res[i].pars;
    fprintf(fres,"%f %f %f %f %f\n",p[0],p[1],p[2],p[3],all_res[i].like);
  }
  fclose(fres);

  */
  // print to log file
  
  
  for(int i=0;1&&i<dumpedFiles.size();i++){
    fprintf(stderr,"dumpedfiles are: %s\n",dumpedFiles[i]);
    free(dumpedFiles[i]);
  }



  fprintf(flog,"\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(flog,"\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fprintf(stderr,"\t[ALL done] cpu-time used =  %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
  fprintf(stderr,"\t[ALL done] walltime used =  %.2f sec\n", (float)(time(NULL) - t2));  
  fclose(flog); 
  
  return 0;

}
