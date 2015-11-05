
#ifndef _types_h
#define _types_h
#include "types.h"
#endif




void killArray(dArray *var);
void killArray(iArray *var);
void killDoubleMatrix(double **var);
void killIntMatrix(int **var);
void killMatrix(dMatrix *var);
void killMatrix(iMatrix *var);

dMatrix *allocDoubleMatrix(int x, int y);
iArray *allocIntArray(int i);
dArray *allocDoubleArray(int i);
dArray *allocDoubleArray(dArray *other);
iMatrix *allocIntMatrix(int x, int y);
void killSnpMatrix(snpMatrix *mat);
