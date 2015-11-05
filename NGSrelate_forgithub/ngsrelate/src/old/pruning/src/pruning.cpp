
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <queue>
#include <deque>
#include <map>         
#include <string>
#include "prune.h"
#ifndef _types_h
#define _types_h
#include "types.h"
#endif

#ifndef _alloc_h
#define _alloc_h
#include "alloc.h"
#endif

using namespace std;

int r_int_to_c_int(SEXP i){
  int retVal;
  PROTECT(i = coerceVector(i, INTSXP));
  if(TYPEOF(i) != INTSXP)
      Rprintf(" expected input is int, you gave wrong type\n");
  retVal = INTEGER(i)[0];
  UNPROTECT(1);
  return retVal;
}

double r_real_to_c_double(SEXP i){
  double retVal;
  PROTECT(i = coerceVector(i, REALSXP));
  if(TYPEOF(i) != REALSXP)
      Rprintf(" expected input is int, you gave wrong type\n");
  retVal = REAL(i)[0];
  UNPROTECT(1);
  return retVal;
}


int r_bool_to_c_int(SEXP i){
  int retVal =0;
  PROTECT(i = coerceVector(i, LGLSXP));
  if(TYPEOF(i) != LGLSXP)
    Rprintf(" expected input is bool, you gave wrong type\n");
  if(LOGICAL(i)[0])
    retVal = 1;
  UNPROTECT(1);
  return retVal;
}


SEXP c_int_to_r_bool(int i){
  SEXP retVal = R_NilValue;
  PROTECT(retVal =allocVector(LGLSXP,i));
  if (i!=0)
    LOGICAL(retVal)[0]=TRUE;
  else
    LOGICAL(retVal)[0]=FALSE;
  UNPROTECT(1);
  return retVal;
}

SEXP c_int_to_r_int(int i){
  SEXP retVal = R_NilValue;
  PROTECT(retVal =allocVector(INTSXP,1));
  INTEGER(retVal)[0]=i;

  UNPROTECT(1);
  return retVal;
}





SEXP c_double_to_r_real(double i){
  SEXP retVal = R_NilValue;
  PROTECT(retVal =allocVector(REALSXP,1));
  REAL(retVal)[0] =i;
  UNPROTECT(1);
  return retVal;
}



dArray *r_array_to_c_dArray(SEXP vec){
  PROTECT(vec);
  int len = length(vec);
  dArray *retVal = allocDoubleArray(len);

  for(int i=0;i<len;i++)
    retVal->array[i] = REAL(vec)[i];
  UNPROTECT(1);
  return retVal;
}

iArray *r_array_to_c_iArray(SEXP vec){
  PROTECT(vec);
  int len = length(vec);
  //  printf("len of vec is:%d\n",len);
  iArray *retVal = allocIntArray(len);
  
  for(int i=0;i<len;i++)
    retVal->array[i] = INTEGER(vec)[i];
  UNPROTECT(1);
  return retVal;
}

SEXP c_iArray_to_r_array(iArray *vec){
  SEXP retVal= R_NilValue;
  if(vec==NULL){
    return retVal;
  }
  int len = vec->x;

  PROTECT(retVal=allocVector(INTSXP,len));
  for(int i=0;i<len;i++)
    INTEGER(retVal)[i] = vec->array[i];
  UNPROTECT(1);
  return retVal;
}


SEXP c_dArray_to_r_array(dArray *vec){
  SEXP retVal= R_NilValue;
  if(vec==NULL){
    return retVal;
  }


  int len = vec->x;

  PROTECT(retVal=allocVector(REALSXP,len));
  for(int i=0;i<len;i++){
    REAL(retVal)[i] = vec->array[i];
  }
  UNPROTECT(1);
  return retVal;
}


iMatrix *r_matrix_to_c_iMatrix(SEXP v){
  if(TYPEOF(v) != INTSXP)
    Rprintf(" input v wrong type\n");
  
  int rows,cols;
  SEXP dim = R_NilValue;
  PROTECT(v);
  PROTECT(dim = getAttrib(v, R_DimSymbol));
  if (length(dim) == 2){
    rows = INTEGER(dim)[0];
    cols = INTEGER(dim)[1];
    //Rprintf("Information: The input contains %i samples with %i snps\n", rows, cols);
  } else
    Rprintf("wrong size of dimension of matrix \n");
  iMatrix *retVal = allocIntMatrix(rows,cols);
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++)
      retVal->matrix[i][j] =INTEGER(v)[j*rows+i];
  
  UNPROTECT(2);
  return retVal;
}


SEXP c_iMatrix_to_r_matrix(iMatrix *v){
  SEXP retVal = R_NilValue;
  if(v==NULL){
    return retVal;
  }

  int rows = v->x;
  int cols = v->y;

  PROTECT(retVal=allocMatrix(INTSXP,rows,cols));
  
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++)
      INTEGER(retVal)[j*rows+i] = v->matrix[i][j];
  
  UNPROTECT(1);
  return retVal;
}


SEXP c_dMatrix_to_r_matrix(dMatrix *v){
  SEXP retVal = R_NilValue;
  if(v==NULL){
    return retVal;
  }

  int rows = v->x;
  int cols = v->y;


  PROTECT(retVal=allocMatrix(REALSXP,rows,cols));
  
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++){
      REAL(retVal)[j*rows+i] = v->matrix[i][j];
      //      printf("conving %f\t",v->matrix[i][j]);
    }
  UNPROTECT(1);
  return retVal;
}


SEXP getListElement(SEXP list,const char *str)     {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  
  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}




extern "C" {
  //this is for getting the data through command args
  SEXP interface(SEXP Rdata, SEXP Rprune, SEXP RLD, SEXP Rback) { 
    
    //init returnValues
    //    printf("start of interface\n");
    SEXP ans = R_NilValue;
    SEXP ans_name = R_NilValue;
    SEXP class_name = R_NilValue;

    //init structure cointain all values needed to run program
    toCargs *pars = new toCargs();
    
    //now input values
    //these are required
    if(Rdata==R_NilValue){
      Rprintf("Must supply arguments\n");
      return ans;
    }
    if(Rprune!=R_NilValue){
      pars->prune_val = r_real_to_c_double(Rprune);
      pars->doPrune =1;
    }
    pars->data = r_matrix_to_c_iMatrix(Rdata);
    //these requires only one change to struct
    pars->back = r_int_to_c_int(Rback);
    pars->LD = r_bool_to_c_int(RLD);
    

  //now do calculation and optimization
  prune object;
  //    print_functionPars(pars);
  fromCres *fromObject = object.main_run(pars);
    
    //now prepare results for R

    //now the data

  int itemNumber  = 0 ;
  PROTECT(ans    = allocVector(VECSXP, 2)); 
  
  SET_VECTOR_ELT(ans, itemNumber++, c_iMatrix_to_r_matrix(fromObject->data));
  SET_VECTOR_ELT(ans, itemNumber++, c_iArray_to_r_array(fromObject->usedSnps));
    
    
    //now input the main names 

  itemNumber  = 0 ;
  PROTECT(ans_name = allocVector(STRSXP, 2));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("data"));
  SET_STRING_ELT(ans_name,itemNumber++ , mkChar("usedSnps"));
  
    //now install the names
  setAttrib(ans, R_NamesSymbol, ans_name);

    //now set up class name
  PROTECT(class_name = allocVector(STRSXP, 1));
  SET_STRING_ELT(class_name, 0, mkChar("pruning.class"));
  classgets(ans, class_name);  
  
  setAttrib(ans, R_NamesSymbol, ans_name);
  UNPROTECT(3);
  
  //clean up
  
  killArray(fromObject->usedSnps);
  killMatrix(fromObject->data);
  killMatrix(pars->data);
  return(ans) ;
  }
}



