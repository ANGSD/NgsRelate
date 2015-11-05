#include <iostream>
#include <fstream>

#include <sstream>
#include <vector>
#include <sys/stat.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#ifndef _types_h
#define _types_h
#include "types.h"
#endif


#include "alloc.h"
#include "prune.h"



/// Checking for file existence, using stat.h.
int fexists(std::string str){///@param str Filename given as a string.
  struct stat buffer ;
  return (stat(str.c_str(), &buffer )==0 ); /// @return Function returns 1 if file exists.
}




void write_iMatrix_to_file(std::string str, iMatrix *mat){
  printf("\t-> Will write stripped genotype file to file: %s\n",str.c_str());
  std::ofstream myfile;
  myfile.open (str.c_str());

  for (int i=0;i < mat->x ;i++){//loop through the keeplist
    for( int j=0; j <mat->y ; j++)
      myfile << mat->matrix[i][j]<< " ";
    myfile << std::endl;
  }
  myfile.close();
}

void write_iArray_to_file(std::string str, iArray *mat){
  printf("\t-> Wil write usedSnps list to file:%s ",str.c_str());
  
  std::ofstream myfile;
  myfile.open (str.c_str());
  
  for(int i=0 ; i<mat->x ; i++)
    myfile << mat->array[i] << " ";
  myfile << std::endl;
  
  myfile.close();
  
  int remains=0;
  for(int i=0 ; i<mat->x ; i++)
    if(mat->array[i])
      remains++;
  printf(" %d of %d SNP's remain\n",remains,mat->x);
  
}



///convert string to double
double to_double(std::string a){///@param a is the string to convert.
  double str;
  std::stringstream ss;
  ss<<a;
  if(!(ss>>str)){
    printf("Error converting :\"%s\" will exit\n",a.c_str());
    exit(0);
  }
  ss>>str;
  return str;///@return The double value.
}


///Converts a string to integer type.
int to_int(std::string a,int lowerlimit, int upperlimit){///@param a a String to convert.
  int str;
  float check;
  std::stringstream ss;
  std::stringstream ss2;
  if(!a.compare("NA")||!a.compare("Na")||!a.compare("na")){
    return 0;
  }
  ss<<a;
  if(!(ss>>str)){
    printf("Error converting :\"%s\" will exit\n",a.c_str());
    exit(0);
  }

  //now check if numberis really an int
  ss2<<a;
  if(!(ss2>>check)){
    printf("Error converting :\"%s\" will exit\n",a.c_str());
    exit(0);
  }

  if(check!=str){
    printf("String you want as in integer is not really an integer: \"%s\"\n",a.c_str());
    puts("Just a warning");
  }
  //check values
  if(lowerlimit !=upperlimit)
    if(str<lowerlimit||str>upperlimit){
      printf("Value in genotype is out of bounds, should be {0...3} value is:%d\n",str);
      exit(0);
    }
  return str;///@return Function returns the integer value.
}


void typecast_stringarray_to_int_matrix(const std::vector<std::string> &tokens,iMatrix *mat){
  int start=0;
  for (int i=0;i<mat->x;i++)
    for(int j=0;j<mat->y;j++)
      mat->matrix[i][j] = to_int( tokens.at(start++),0,3);
}


/*
  Mother-of-all string tokenizers
  1. args = string to be splittet
  2. args = a vector<string> reference that will contain the splittet string
  3. args = delimiter list in a string like 
  
  get_lexemes(STRING_TO_SPLIT,RESULTS_STRING_VECTOR,DELIM=" \t:,;")
  this will split on whitespace,tab,colons,comma and semicolons
  return value will be number of elems tokenized
*/
///Mother-of-all string tokenizers! Function will split a string from a given set of delimiters.
int get_lexemes(std::string str,std::vector<std::string> &lexems,const std::string delims){
  ///@param str The string to be splittet.
  ///@param lexems A reference to a list of strings in which the splittet string will be put as seperate elements. A queue.
  ///@param delims A string containing all the chars that the function will split by.
  char *token = strtok( const_cast<char*>( str.c_str() ),delims.c_str());
  int numToks = 0;
  while ( token != NULL ){
    numToks++;
    lexems.push_back( token);
    token = strtok( NULL,delims.c_str() );
  }
  return numToks;///@return The number of strings that has been splittet.
}



iMatrix *readmatrix(std::string filename,const std::string delim=",;: \t"){
  ///@param filename A filename to read.@param delim A string of delimiters.
  std::vector<std::string> tokens;

  const int SIZE=500000;
  char buffer[SIZE];
 
  std::ifstream pFile (filename.c_str(),std::ios::in);
 
  
  if(!pFile){
    std::cout <<"Problems opening file" <<filename<<std::endl;
    exit(0);
  }

  
  std::string tmp_string;
  int doFirstRow =1;
  int itemsInFirstRow=0;
  int numRows =0;
  while(!pFile.eof()){
    pFile.getline(buffer,SIZE);
    tmp_string = std::string(buffer);
    
    if(doFirstRow){
      //if file has a emptystart line
      itemsInFirstRow = get_lexemes(tmp_string,tokens,delim);
      
      if (itemsInFirstRow==0)
	continue;
      // printf("items in first rwo:%d\n",itemsInFirstRow);
      doFirstRow=0;
      numRows++;
    }
    else{
      int nItems = get_lexemes(tmp_string,tokens,delim);
	//if line is empty
	if(nItems==0)
	  continue;
	numRows++;
	if(nItems!=itemsInFirstRow){
	  printf("row length mismatch at line:%d numitems is:%d shouldn't be:%d\t will exit\n",numRows,itemsInFirstRow,nItems);
	  exit(0);
	}
      }
  }

  iMatrix *data_ = allocIntMatrix(numRows,itemsInFirstRow);
  //now we have a token array of string coerce the types now
  typecast_stringarray_to_int_matrix(tokens,data_);
  //copy(tokens.begin(), tokens.end(), ostream_iterator<string>(cout, ", "));
  
  printf("\t-> Dimension of genotype datafile is (%d,%d)\n",data_->x,data_->y);
  return data_; 
}



toCargs *readArgs(int argc, char *argv[]){
  toCargs *pars = new toCargs();
  
  int start = 1;
  while (start<argc){
    if(!strcmp(argv[start],"-g")){
      pars->data = readmatrix(argv[start+1]);
    }
    else if(!strcmp(argv[start],"-LD")){
      pars->LD = atoi( argv[start+1]);
    }
    else if(!strcmp(argv[start],"-back")){
      pars->back =to_int ( argv[start+1],0,0);
    }
    else if(!strcmp(argv[start],"-threshold")){
      pars->prune_val =to_double(std::string(argv[start+1]));
    }
    else if(!strcmp(argv[start],"-filename")){
      pars->data_filename = std::string(argv[start+1]);
    }
    else if(!strcmp(argv[start],"-snps")){
      pars->usedSnps_filename = std::string(argv[start+1]);
    }
    else{
      printf("\t->Command line option not recognized: \t%s\n",argv[start]);
      exit(0);
    }
    start+=2;
  }
  return pars;
}

void print_functionPars(toCargs *pars){
  printf("\t-> function args are\n");
  printf("\t-> back=%d\t prune_val=%f\t LD=%d\n",pars->back,pars->prune_val,pars->LD);
}



int main(int argc, char *argv[]){
  toCargs *pars;
  if(argc<12){
    std::cout << "\t-> Please input -g GENOTYPEFILE -LD {0,1} -back BACK -threshold prune_val -filename OUTPUTDATA -snps KEEPLIST\n";
    return 0;
  }
  //get commandline args
  pars = readArgs(argc,argv);
  pars->doPrune = 1;
  
  //check if fileexists
  if(fexists(pars->data_filename)){
    printf("\t-> outputfile:%s exists, will exit\n",pars->data_filename.c_str());
    return 0;
  }
  if(fexists(pars->usedSnps_filename)){
    printf("\t-> outputfile:%s exists will exit\n",pars->usedSnps_filename.c_str());
    return 0;
  }

  //run prune
  prune object;
  print_functionPars(pars);
  fromCres *fromObject = object.main_run(pars);


  //write results
  write_iMatrix_to_file(pars->data_filename,fromObject->data);
  write_iArray_to_file(pars->usedSnps_filename,fromObject->usedSnps);

  killArray(fromObject->usedSnps);
  killMatrix(fromObject->data);
  killMatrix(pars->data);

  return 0;
}
