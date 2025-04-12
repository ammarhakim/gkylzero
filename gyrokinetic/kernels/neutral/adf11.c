#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adf11.h"

char* adas_data_dir = ""; // /Users/bernard/gkylzero/kernels/neutral/";

/* Read adas adf11 data file.
 * Input - filename, either the full path to an adf11 file or
 *      the filename within the adas_data_dir
 * Output - atomic_data structure containing log10 of Te, ne, and
 *      the rate coefficients.
 */
atomic_data loadadf11(char* filename){
  char* filepath, *lineptr, *endptr;
  FILE *ptr;
  ssize_t status;
  size_t n=0;
  atomic_data data;
  int j,ind;

  /* if(filename[5]=='r'){ */
  /*   fprintf(stderr,"Error. Only metastable unresolved files supported"); */
  /*   exit(EXIT_FAILURE); */
  /* } */
  
  filepath = get_adas_file_loc(filename);
  data.reaction_type = (char*) malloc(sizeof(char)*4);
  sprintf(data.reaction_type, "%.3s", filename);
  
  ptr = fopen(filepath,"r");
  if(!ptr){
    fprintf(stderr,"Error opening file '%s'\n",filepath);
    exit(EXIT_FAILURE);
  }

  lineptr = NULL;
  status = getline(&lineptr,&n,ptr);
  data.name = (char*)malloc(sizeof(char)*30);
  data.symbol = (char*)malloc(sizeof(char)*3);
  //read header, name is shifted as it is /ELEMENT NAME
  status = sscanf(lineptr,"%5d%5d%5d%*d%*d%s",&data.Z,&data.ne_intervals,&data.te_intervals,data.name);
  data.name = data.name+1;
  if(status == -1){
    fprintf(stderr,"Failed to read header of ADAS file");
  }
  data.logNe = (double*)malloc(sizeof(double)*data.ne_intervals);
  data.logT = (double*)malloc(sizeof(double)*data.te_intervals);
  data.logData = (double*)malloc(sizeof(double)*data.ne_intervals*data.te_intervals*data.Z);
  data.charge_states = (int*) malloc(sizeof(int)*data.Z);

  status = getline(&lineptr,&n,ptr);//skip line
  //Read density abscissa
  //Data is in 8 columns. Partial last row if not perfectly divisible 
  for(int i=0;i<ceil((data.ne_intervals-.01)/8.0);i++){
    status = getline(&lineptr,&n,ptr);
    if(i<data.ne_intervals/8){
      sscanf(lineptr,"%10lf%10lf%10lf%10lf%10lf%10lf%10lf%10lf",&data.logNe[i*8],
	     &data.logNe[i*8+1],&data.logNe[i*8+2],&data.logNe[i*8+3],&data.logNe[i*8+4],
	     &data.logNe[i*8+5],&data.logNe[i*8+6],&data.logNe[i*8+7]);
    }else{
      j = i*8;
      data.logNe[j] = strtod(lineptr,&endptr);
      while(j<data.ne_intervals-1){
	j++;
	data.logNe[j] = strtod(endptr,&endptr);
	printf("j=%d",j);
      }
    }
  }

  //Read temperature abscissa
  for(int i=0;i<ceil((data.te_intervals-.01)/8.0);i++){
    status = getline(&lineptr,&n,ptr);
    if(i<data.te_intervals/8){
      sscanf(lineptr,"%10lf%10lf%10lf%10lf%10lf%10lf%10lf%10lf",&data.logT[i*8],
	     &data.logT[i*8+1],&data.logT[i*8+2],&data.logT[i*8+3],&data.logT[i*8+4],
	     &data.logT[i*8+5],&data.logT[i*8+6],&data.logT[i*8+7]);
    }else{
      j = i*8;
      data.logT[j] = strtod(lineptr,&endptr);
      while(j<data.te_intervals-1){
	j++;
	data.logT[j] = strtod(endptr,&endptr);
      }
    }
  }

  //Read rate coefficient - loop through densities for each Te
  for(int k=0;k<data.Z;k++){
    status = getline(&lineptr,&n,ptr);
    sscanf(lineptr,"%*s%*s%*d%*s%*s%*d%*s%*s%d",&data.charge_states[k]);
    for(j=0;j<data.te_intervals;j++){
      for(int i=0;i<ceil((data.ne_intervals-.01)/8);i++){
	ind = k*data.ne_intervals*data.te_intervals + j*data.ne_intervals+i*8;
	status = getline(&lineptr,&n,ptr);
	if(i<data.ne_intervals/8){
	  sscanf(lineptr,"%10lf%10lf%10lf%10lf%10lf%10lf%10lf%10lf",&data.logData[ind],
		 &data.logData[ind+1],&data.logData[ind+2],&data.logData[ind+3],&data.logData[ind+4],
		 &data.logData[ind+5],&data.logData[ind+6],&data.logData[ind+7]);
	}else{
	  status = getline(&lineptr,&n,ptr);
	  data.logData[ind] = strtod(lineptr,&endptr);
	  while(ind<data.ne_intervals*data.te_intervals*data.Z-1){
	    ind++;
	    data.logData[ind] = strtod(endptr,&endptr);
	  }
	}
      }
    }
  }
  free(lineptr);
  
  fclose(ptr);
  
  //convert units to SI
  for(int i=0;i<data.ne_intervals;i++){
    data.logNe[i] = data.logNe[i]+6;
  }
  
  for(int k=0;k<data.Z;k++){
    for(j=0;j<data.te_intervals;j++){
      for(int i=0;i<data.ne_intervals;i++){
	ind = k*data.ne_intervals*data.te_intervals + j*data.ne_intervals+i;
	data.logData[ind]=data.logData[ind]-6;
      }
    }
  }
  printf("\n");
  return data;
}

/* Get location of adas file
 *
 */
char* get_adas_file_loc(char* filename){
  if(strcasecmp(filename,"none")==0){
    return NULL;
  }else if(file_isreg(concat(adas_data_dir,filename))==1){
    return concat(adas_data_dir,filename);
  }else if(file_isreg(filename)==1){
    return filename;
  }else{
    printf("\nFile is not on system. Not able to directly download");
    return "Error";
  }
}

/*
 * Concatenate two strings 
 */
char* concat(const char *s1, const char *s2)
{
    char *result = malloc(strlen(s1) + strlen(s2) + 1); // +1 for the null-terminator
    if(result == NULL){
      fprintf(stderr,"Error in allocating memory for string in adf11.c");
      exit(EXIT_FAILURE);
    }
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}
  
/* Check if path to file is a regular file
 */  
int file_isreg(const char *path) {
  struct stat st;
  if (stat(path, &st) < 0){
    return -1;
  }
  return S_ISREG(st.st_mode);
}
