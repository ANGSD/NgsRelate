/*
 This is pruning

*/

#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>



using namespace std;


#ifndef _types_h
#define _types_h
#include "types.h"
#endif

#ifndef _alloc_h
#define _alloc_h
#include "alloc.h"
#endif


#include "ld.h"
#include "prune.h"

dMatrix *myAbs(dMatrix *matr){
  dMatrix *retVal = allocDoubleMatrix(matr->x,matr->y);
  for (int i=0;i<retVal->x;i++)
    for (int j=0;j<retVal->y;j++)
      retVal->matrix[i][j]=fabs(matr->matrix[i][j]);
  
  killMatrix(matr);
  return retVal;
    
}

iMatrix *transpose(iMatrix *mat){
  iMatrix *retMat = allocIntMatrix(mat->y,mat->x);
  for (int i=0;i<retMat->x;i++)
    for (int j=0;j<retMat->y;j++)
      retMat->matrix[i][j]=mat->matrix[j][i];
  return retMat;
}

iArray *generateIndices(iArray *keepList){
  int numTrue =0;
  for(int i=0;i<keepList->x;i++)
    if(keepList->array[i]==1)
      numTrue++;

  iArray *returnArray = allocIntArray(numTrue);
  int atPos = 0;
  for(int i=0;i<keepList->x;i++)
    if(keepList->array[i]==1){
      returnArray->array[atPos] = i;
      atPos++;
    }
  return returnArray;
	
}

dArray *extract(dMatrix *var,iArray *choose){
  dArray *retVal = allocDoubleArray(choose->x-1);
   for (int i=1;i<var->y;i++){
    retVal->array[i-1]=var->matrix[choose->array[i]][i];
  }
  return retVal;
}


iArray *extract(iArray *var,iArray *choose){
  iArray *retVal = allocIntArray(var->x-1);
  for (int i=1;i<var->x;i++){
    retVal->array[i-1]=var->array[i-1-choose->array[i]];
  }
  return retVal;
}



iMatrix *revCols(iMatrix *matr,int offset){
  iMatrix *retVal = allocIntMatrix(matr->x,matr->y);
  
  for (int x=0;x<matr->x;x++){
    for (int i=0;i<matr->y;i++){
      retVal->matrix[x][i]=matr->matrix[x][matr->y-1-offset-i];
    }
  }
  return retVal;
}


dMatrix *revCols(dMatrix *matr,int offset){
  dMatrix *retVal = allocDoubleMatrix(matr->x,matr->y);
  
  for (int x=0;x<matr->x;x++){
    for (int i=0;i<matr->y;i++){
      retVal->matrix[x][i]=matr->matrix[x][matr->y-1-offset-i];
    }
  }
  return retVal;
}


snpMatrix *snp_pair_range(iMatrix *v, int  depth){
  snpMatrix *retVal = new snpMatrix();
  int need_signed_r = 0;
  
  int rows = v->x;
  //  int cols = v->y;

  int width =v->y-1; //number of snps to process

  dMatrix *dprime = allocDoubleMatrix( depth,width);
  dMatrix *rmisc  = allocDoubleMatrix( depth,width);
  dMatrix *lod    = allocDoubleMatrix( depth,width);
  dMatrix *D      = allocDoubleMatrix( depth,width);
  dMatrix *pBA      = allocDoubleMatrix(  depth,width);
  dMatrix *pBa      = allocDoubleMatrix( depth,width);
  dMatrix *pbA      = allocDoubleMatrix(  depth,width);
  dMatrix *pba      = allocDoubleMatrix(  depth,width);
  iMatrix *transposed=transpose(v);
  for(int idx_j = 0; idx_j < depth; idx_j++){
    for (int idx_i = 0  ; idx_i < v->y - 1 - idx_j; idx_i++) {
      geno_cptr res = get_geno_count(transposed->matrix[idx_i],transposed->matrix[idx_i+1+idx_j],rows);
      dprime->matrix[idx_j][idx_i] = res->dprime;
      D->matrix[idx_j][idx_i] = res->bigD/((res->total * 2.0)*(res->total * 2.0));
      pBA->matrix[idx_j][idx_i]  = res->u/(res->total * 2.0);
      pBa->matrix[idx_j][idx_i] = res->v/(res->total * 2.0);
      pbA->matrix[idx_j][idx_i] = res->w/(res->total * 2.0);
      pba->matrix[idx_j][idx_i] = res->x/(res->total * 2.0);
      
      if (need_signed_r) {
	if (res->rsq2 > 0) {
	  rmisc->matrix[idx_j][idx_i] =  res->sign_of_r * sqrt(res->rsq2);
	  //	  REAL(rmisc)[offset + idx_j * width] = res->sign_of_r * sqrt(res->rsq2);
	}else {
	  rmisc->matrix[idx_j][idx_i] = -2;
	  //REAL(rmisc)[offset + idx_j * width] = -2;
	}
      } else {
	//	cout << res->rsq2 <<endl;
	rmisc->matrix[idx_j][idx_i] =  res->rsq2;
	//REAL(rmisc)[offset + idx_j * width] = res->rsq2;
      }
      //  REAL(lod)[offset + idx_j * width] = res->lod;
      lod->matrix[idx_j][idx_i] = res->lod;
      free(res->expt);
      free(res);
  
    }
  }
  
  killMatrix(transposed);
  //collect result and return in a struct
  retVal->dprime=dprime;
  retVal->D=D;
  retVal->pBA=pBA;
  retVal->pBa=pBa;
  retVal->pbA=pbA;
  retVal->pba=pba;
  retVal->rmisc = rmisc;
  retVal->lod = lod;
  return retVal;
}



iArray *pruning(snpMatrix *ld,int ld_choose,int back,double prune_val){
  dMatrix *mea;
  if(ld_choose)
    mea = revCols(ld->rmisc,0);
  else
    mea = revCols(ld->D,0);
  if(ld_choose==0)
    mea = myAbs(mea); 
  
  /*
    the procedure is strange,
    if column sum is greater the pruneValue,
    then remove the entire SNP locus.
    then do a SHIFTUP operation on the next BACK-1 loci.
    SHIFT defined as removing the diagonal, and moving the rest up and appending zero
    EXAMPLE:
     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
[1,]    1    6   11   16   21   26   31   36   41    46
[2,]    2    7   12   17   22   27   32   37   42    47
[3,]    3    8   13   18   23   28   33   38   43    48
[4,]    4    9   14   19   24   29   34   39   44    49
[5,]    5   10   15   20   25   30   35   40   45    50
   removing SNP 3 will turn the above to

     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
[1,]    1    6   16   21   26   31   36   41   46
[2,]    2    7   18   22   27   32   37   42   47
[3,]    3    8   19   24   28   33   38   43   48
[4,]    4    9   20   25   30   34   39   44   49
[5,]    5   10    0    0    0    0   40   45   50


  */
	      
  
  int snp = mea->y+1; //number of snps before process;
  iArray *retVal = allocIntArray(ld->rmisc->y+1); //snp length plus one
  retVal->array[0] = 1; // always keep first
  if(back>1) {
    for(int nsnp=1;nsnp<snp;nsnp++) {//iterate thourgh SNP's
      int doStripping=0;
      for(int j=0;j<back;j++)
	if( mea->matrix[j][nsnp-1] >prune_val){
	
	  doStripping=1;
	  break;//jump out of inner loop
	}
      if(doStripping) {
	//should update datastructure.
	retVal->array[nsnp] = 0; //col will be excluded from the keeplist
	if(back>2){
	  for (int p=0 ; p < back-1 ; p++){//number of colums to update
	    if(p+nsnp>=mea->y) //if no more colums exists, exit
	      break;
			    
	    for(int q=p+1;q<back-1;q++){
	      
	      mea->matrix[q][nsnp+p] = mea->matrix[q+1][p+nsnp];
	    }
	    
	    mea->matrix[back-1][nsnp+p] = 0;//input zero at bottom    
	  }
	}else //back must =2 so just remove 
	  
	  mea->matrix[1][nsnp] = 0;
      }
      else {
	
	//we shouldn't remove this snp so set 1 in keeplist
	retVal->array[nsnp] = 1;
      }
    }
  }
	

  else{
    //back is 1 so just check if less than prunevalue
    for (int i=1;i<retVal->x;i++)
      if(mea->matrix[0][i-1] <= prune_val)
	retVal->array[i] = 1;
  }
  killMatrix(mea);
  return retVal;
}

// keep<-! apply(is.na(rbind(data[pair[1],],data[pair[2],])),2,any) & maf> min&maf<1-min
iMatrix *extractOK(iArray *okList,iMatrix *matr){
  
  iMatrix *returnMatrix = allocIntMatrix(matr->x,okList->x);
  for(int j=0;j<matr->x;j++)
    for(int i=0;i<okList->x;i++){
      returnMatrix->matrix[j][i] = matr->matrix [j] [okList->array[i]];
    }
  
  return returnMatrix;

}

iArray *extractOK(iArray *okList,iArray *array){
  
  iArray *returnArray = allocIntArray(okList->x);
  for(int j=0;j<okList->x;j++)
    returnArray->array[j] = array->array[okList->array[j]];
  return returnArray;
}


fromCres *getPars(iMatrix *data ,int back,int ld_choose,int doPrune, double prune_val){
  iMatrix *tmp = revCols(data,0);
  snpMatrix *ld = snp_pair_range(tmp,back);
  iArray *iArrayTmp = pruning(ld,ld_choose,back,prune_val);
  iArray *okList = generateIndices(iArrayTmp);
      
  iMatrix *iMatrixTmp = extractOK(okList,data);

  fromCres *returnVal = new fromCres();
  returnVal->data = iMatrixTmp;
  returnVal->usedSnps =  iArrayTmp;
  //now clean up.
  killSnpMatrix(ld);
  killArray(okList);
  killMatrix(tmp);
  return returnVal;
}


fromCres *prune::main_run(toCargs *pars) {
 
  int back=pars->back;
  int ld_choose =pars->LD; //this is rsq2 if ld_choose=0 then snpMatrix->D
  int doPrune = pars->doPrune;
  double prune_val = pars->prune_val;
 
  fromCres *res =  getPars(pars->data,back,ld_choose,doPrune,prune_val);
 
  return res;
}

