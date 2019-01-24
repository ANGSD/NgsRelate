#include <cstring>
#include <map>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <zlib.h>
#include <vector>
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

size_t nlines(const char *fname){
  FILE *fp = NULL;
  fp=fopen(fname,"rb");
  assert(fp);
  size_t nlines = 0;
  while(!feof(fp)){
    char ch = fgetc(fp);
    if(ch == '\n')
	nlines++;
  }
  return nlines;
}



int **bed_to_intMatrix(const char* file, int nrow,int ncol) {

  //  const unsigned char recode[4] = {'\x01', '\x00', '\x02', '\x03'};
  //0,1,2: 3 is missing 
  const unsigned char recode[4] = { '\x02','\x01', '\x03', '\x00'};
  const unsigned char mask = '\x03';


  FILE *in = fopen(file, "r");
  if (!in){
    printf("Couln't open input file: %s", file);
    exit(0);
  }
  unsigned char start[3];
  if (fread(start, 1, 3, in)!=3){
    printf("Failed to read first 3 bytes");
    exit(0);
  }
  if (start[0]!='\x6C' || start[1]!='\x1B'){
    printf("Input file does not appear to be a .bed file (%X, %X)", 
	   start[0], start[1]);
    exit(0);
  }
  /* Create output object */
  
  
  int **returnMat =new int*[nrow];
  for(int i=0;i<nrow;i++){
    returnMat[i] = new int[ncol];
    for(int j=0;j<ncol;j++)
      returnMat[i][j] = -1;
  }
  int ncell = nrow*ncol;
  unsigned char *result = new unsigned char[nrow*ncol]; 
  memset(result, 0x00, ncell);

  /* Read in data */

  int snp_major = start[2];
  int part=0, ij=0, i=0, j=0;
  while (1) {
    unsigned char byte;
    if (!part) {
      if (feof(in) || !fread(&byte, 1, 1, in)) {
	printf("Unexpected end of file reached");
	exit(0);
      }
      part = 4;
    }
    unsigned char code = byte & mask;
    byte = byte >> 2;
    part--;
    result[ij] = recode[code];
    returnMat[i][j] = result[ij];
    if(returnMat[i][j]==3)
      returnMat[i][j]=1;
    else if(returnMat[i][j]==1)
      returnMat[i][j]=3;
    else if(returnMat[i][j]<0 || returnMat[i][j]>3){
      printf("Problem in bed file at position=(%d,%d)=%d\n",i,j,returnMat[i][j]);
      exit(0);
    }
    // printf("(%d,%d)=%d ",i,j,result[ij]);
    if (snp_major) {
      ij++;
      i++;
      if (i==nrow) {
	i = part = 0;
	j++;
	if (j==ncol)
	  break;
      }
    }	
    else {
      ij += nrow;
      j++;
      if (j==ncol){
	j = part = 0;
	i++;
	if (i==nrow)
	  break;
	ij = i;
      }
    }
  }
  fclose(in);
  delete [] result;
  return returnMat;
}


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
