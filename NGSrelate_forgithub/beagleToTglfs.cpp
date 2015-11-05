/*
  program to convert beagle style files to binary 10llh with a position.
  for -tglf input in angsd. This program could be improved quite alot....

 */

#include <sys/stat.h>
#include <zlib.h>
#include <cassert>
#include <vector>
#include "kstring.h"
#include <cstdio>
#include <cmath>

static int majorminor[4][4] = 
  {
    {0,1,2,3},
    {1,4,5,6},
    {2,5,7,8},
    {3,6,8,9}
  };




char *inname=NULL;
char *outname=NULL;
int intName = 1;
int fexists(const char* str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str, &buffer )==0 ); /// @return Function returns 1 if file exists.
}


gzFile setArgs(int argc,char **argv,int &nInd){
  if(argc<3){
    fprintf(stderr,"\t-> Supply in.beagle prefix.out [intName]");
    exit(0);
  }
  inname = strdup(argv[1]);
  outname = argv[2];
  if(argc==4)
    intName = atoi(argv[3]);
  
  gzFile gz = Z_NULL;
  gz=gzopen(inname,"rb");
 
  size_t lens = 1000000;
  char *buffer = new char[lens];
  const char *delims = "\t \n";
  int nCol=1;
  nInd =0;
  gzgets(gz,buffer,lens); 
  strtok(buffer,delims);
  while(strtok(NULL,delims))
    nCol++;
  if(nCol % 3 ){
    fprintf(stderr,"\t-> Number of columns should be a multiple of 3, nCol=%d\n",nCol);
    exit(0);
  } 
  nInd=nCol/3-1;
  delete[] buffer;
  fprintf(stderr,"\t-> inname: %s \t outname: %s\t intName: %d\t nInd=%d\n",inname,outname,intName,nInd);
  return gz;
}


int main(int argc,char **argv){

  char refToChar[256] = {
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
  
  int nInd;
  gzFile gz = setArgs(argc,argv,nInd);
  
  std::vector<char *> chromo;
  std::vector<int> pos;
  kstring_t major,minor;
  major.s=minor.s=NULL;
  major.l=major.m=minor.l=minor.m =0;
  
  
  double *perSiteLlh = new double[10*nInd];


  size_t lens = 1000000;
  char *buffer = new char[lens];
  const char *delims = "\t \n";
  const char *delims2 = "_\t \n";

  
  int nSites=0;
  int positions =0;//every site is a new position across different chunks
  
  char oNam[1024];
  snprintf(oNam,1024,"%s.glf.gz",outname);
  fprintf(stderr,"\t-> Ofile:%s\n",oNam);
  
  gzFile gz_o = Z_NULL;
  gz_o=gzopen(oNam,"w6h");
  if(gz_o==Z_NULL){
    fprintf(stderr,"problem opening file: %s\n",oNam);
  }
  
  while(gzgets(gz,buffer,lens)){
    for(int i=0;i<10*nInd;i++)
      perSiteLlh[i] =log(0);
    assert(buffer!=NULL);
    
    if(intName){
      chromo.push_back(strdup(strtok(buffer,delims2)));
      pos.push_back(atoi(strtok(NULL,delims2))-1);
    }else{
      chromo.push_back(strdup(strtok(buffer,delims)));
      pos.push_back(positions++);
    }
    char IMaj = refToChar[ strtok(NULL,delims)[1]];
    char IMin = refToChar[strtok(NULL,delims)[0]];
    ksprintf(&major,"%c",IMaj);
    ksprintf(&minor,"%c",IMin);

    //    fprintf(stderr,"majmaj:%d\tmajmin:%d\tminmin:%d\n",majorminor[IMaj][IMaj],majorminor[IMaj][IMin],majorminor[IMin][IMin]);
    
    for(int i=0;i<nInd;i++){
      char *tsk = strtok(NULL,delims);assert(tsk!=NULL);
      perSiteLlh[i*10+majorminor[IMaj][IMaj]] = log(atof(tsk));
      tsk = strtok(NULL,delims);assert(tsk!=NULL);
      perSiteLlh[i*10+majorminor[IMaj][IMin]] = log(atof(tsk));
      tsk = strtok(NULL,delims);assert(tsk!=NULL);
      perSiteLlh[i*10+majorminor[IMin][IMin]] = log(atof(tsk));
    }
    
    nSites++;
    gzwrite(gz_o,perSiteLlh,sizeof(double)*10*nInd);
      //    fwrite(perSiteLlh,sizeof(double),10*nInd,gz_o); 
  }
  gzclose(gz_o);
  snprintf(oNam,1024,"%s.pos.gz",outname);
  fprintf(stderr,"\t-> Ofile:%s\n",oNam);
  gz_o=gzopen(oNam,"w6h");
  if(gz_o==Z_NULL){
    fprintf(stderr,"problem opening file: %s\n",oNam);
    oNam;
  }
  assert(chromo.size()==pos.size());
  assert(pos.size()==major.l);
  assert(major.l==minor.l);

  for(uint s=0;s<chromo.size();s++)
    gzprintf(gz_o,"%s\t%d\t%d\t%d\n",chromo[s],pos[s],major.s[s],minor.s[s]);
  
  gzclose(gz_o);
  return 0;
}
