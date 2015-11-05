//alloc by allocIntMatrix
//deall by killMatrix

#include <string>


typedef struct  {
  int x;
  int y;
  int** matrix;
}iMatrix ;

//alloc by allocDoubleMatrix
//deall by killMatrix
typedef struct{
  int x;
  int y;
  double** matrix;
} dMatrix ;

//alloc by allocDoubleMatrix
//deall by killMatrix
typedef struct{
  int x;
  int * array;
}iArray;

typedef struct{
  int x; //length
  int numTrue; //number of items;
  bool * array;
}bArray;

//alloc by allocDoubleArray
//deall by killArray
typedef struct {
  int x;
  double * array;
}dArray;


typedef struct {
  /*
    matrix[d][s]
    d := depth
    s := snp
  */
  dMatrix *dprime;
  dMatrix *pba;
  dMatrix *pbA;
  dMatrix *pBa;
  dMatrix *pBA;
  dMatrix *rmisc;
  dMatrix *D;
  dMatrix *lod;
}snpMatrix;


typedef struct  {
  iMatrix *data;
  int doPrune;
  int LD;
  int back;
  double prune_val;
  std::string data_filename;
  std::string usedSnps_filename;
}toCargs;


typedef struct {
  iMatrix *data;
  iArray *usedSnps;
}fromCres;
