#include <iostream>
#ifndef _types_h
#define _types_h
#include "types.h"
#endif

#include "alloc.h"


iMatrix *allocIntMatrix(int x, int y){
  iMatrix *tmp = new iMatrix();
  int **ppi = new int*[x];
  int *curPtr = new int [x * y];
  
  for( int i = 0; i < x; ++i) {
      *(ppi + i) = curPtr;
      curPtr += y;
    }
  for (int i=0;i<x;i++)
    for(int n=0;n<y;n++)
      ppi[i][n]=0;
  tmp->matrix=ppi;
  tmp->x=x;
  tmp->y=y;
  return tmp;
}


void killArray(dArray *var){
  delete [] var->array;
  delete var;
}

void killArray(iArray *var){
  delete [] var->array;
  delete var;
}


void killDoubleMatrix(double **var){
  delete [] *var;
  delete [] var;
}


void killIntMatrix(int **var){
  delete [] *var;
  delete [] var;
}

void killMatrix(dMatrix *var){
  killDoubleMatrix(var->matrix);
  delete var;
}


void killMatrix(iMatrix *var){
  killIntMatrix(var->matrix);
  delete var;
}




void killSnpMatrix(snpMatrix *mat){
  killMatrix(mat->pba);
  killMatrix(mat->dprime);
  killMatrix(mat->pbA);
  killMatrix(mat->pBa);
  killMatrix(mat->pBA);
  killMatrix(mat->rmisc);
  killMatrix(mat->D);
  killMatrix(mat->lod);
  delete mat;
}







dMatrix *allocDoubleMatrix(int x, int y){
  dMatrix *tmp = new dMatrix();
  double **ppi = new double*[x];
  double *curPtr = new double [x * y];
  
  for( int i = 0; i < x; ++i) {
      *(ppi + i) = curPtr;
      curPtr += y;
    }
  for (int i=0;i<x;i++)
    for(int n=0;n<y;n++)
      ppi[i][n]=0;
  tmp->matrix=ppi;
  tmp->x=x;
  tmp->y=y;
  return tmp;
}



iArray *allocIntArray(int i){
  int *tmp= new int[i];
  iArray *intArray = new iArray();
  for (int j=0;j<i;j++)
    tmp[j]=0;
  intArray->x=i;
  intArray->array=tmp;
  return intArray;
}

dArray *allocDoubleArray(int i){
  double *tmp= new double[i];
  dArray *darray = new dArray();
  for (int j=0;j<i;j++)
    tmp[j]=0.0;
  darray->x=i;
  darray->array=tmp;
  return darray;
}
//copy constructor
dArray *allocDoubleArray(dArray *other){
  dArray *darray = new dArray();
  double *tmp= new double[other->x];
  
  for (int j=0;j<other->x;j++)
    tmp[j]=other->array[j];
  darray->x=other->x;
  darray->array=tmp;
  return darray;
}

