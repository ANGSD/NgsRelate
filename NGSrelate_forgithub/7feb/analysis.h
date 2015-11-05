
typedef struct{
  char *name;
  int nInd;
  int nSites;
  int *pos;
  double *freq;//in log
  double *qerf;//in log 1-freq
  double **gl;//in log
  char **keeps;
  int *keepInd;
}perChr;

typedef struct{
  double a;
  double k0;
  double k1;
  double k2;
  int pair[2];
}para;


typedef struct{
  double pars[4];
  double like;
  int nSites;
  int *pos;
  double **post;
  double **forward;
  double **backward;
  double **emis;
  double **trans;
  char *viterbi;
  double polike;
  double ulike;
}hmm;


hmm analysis(const perChr &pc,double *freq,para p,int calcA);
double addProtect2(double a,double b);
double addProtect3(double a,double b, double c);
