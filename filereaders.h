int extract_freq(int argc,char **argv);
int extract_freq_bim(int argc,char **argv);
void readids(std::vector<char *> &ids,char *fname);
size_t nlines(const char *fname);
int **bed_to_intMatrix(const char* file, int nrow,int ncol);
size_t getDouble(const char *fname,std::vector<double> &ret);
double **getGL(const char *fname, int sites, int nInd);
