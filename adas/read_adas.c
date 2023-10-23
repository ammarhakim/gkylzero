#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gkyl_array.h>
#include <gkyl_range.h>
#include <gkyl_dg_iz.h>
#include <gkyl_dg_recomb.h>
#include <read_adas.h>

void
read_adas_field_iz(enum gkyl_dg_iz_type type_ion, struct adas_field *data, const char *base) {
  char *fname[4000];
  if (type_ion == GKYL_IZ_H) {
    strcpy(fname, base);
    data->NT = 29, data->NN = 24;
    strcpy(fname, base);  
    strcat(fname, "/adas/ioniz_h.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);
    strcat(fname, "/adas/logT_h.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_h.npy");
    data->logN = fopen(fname, "rb");    
    data->Zmax = 1;
    data->Eiz = malloc(sizeof(double)*data->Zmax);
    static const double Eiz_loc[] = {13.6};
    memcpy(data->Eiz, Eiz_loc, sizeof(Eiz_loc));
  }
  else if (type_ion == GKYL_IZ_HE) {
    data->NT = 30;
    data->NN = 24;
    strcpy(fname, base);  
    strcat(fname, "/adas/ioniz_he.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_he.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_he.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 2;
    data->Eiz = malloc(sizeof(double)*data->Zmax);
    static const double Eiz_loc[] = {24.6, 54.4};
    memcpy(data->Eiz, Eiz_loc, sizeof(Eiz_loc));
  }
  else if (type_ion == GKYL_IZ_LI) {
    data->NT = 25;
    data->NN = 16;
    strcpy(fname, base);  
    strcat(fname, "/adas/ioniz_li.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_li.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_li.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 3;
    data->Eiz = malloc(sizeof(double)*data->Zmax);
    static const double Eiz_loc[] = {5.4, 75.6, 122.4};
    memcpy(data->Eiz, Eiz_loc, sizeof(Eiz_loc));
  }
  else if (type_ion == GKYL_IZ_BE) {
    data->NT = 25;
    data->NN = 24;
    strcpy(fname, base);  
    strcat(fname, "/adas/ioniz_be.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_be.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_be.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 4;
    data->Eiz = malloc(sizeof(double)*data->Zmax);
    static const double Eiz_loc[] = {9.3, 18.2, 153.9, 217.7};
    memcpy(data->Eiz, Eiz_loc, sizeof(Eiz_loc));
  }
  else if (type_ion == GKYL_IZ_B) {
    data->NT = 48;
    data->NN = 26;
    strcpy(fname, base);  
    strcat(fname, "/adas/ioniz_b.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_b.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_b.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 5;
    data->Eiz = malloc(sizeof(double)*data->Zmax);
    static const double Eiz_loc[] = {8.3, 25.2, 37.9, 259.4, 340.2};
    memcpy(data->Eiz, Eiz_loc, sizeof(Eiz_loc));
  }
  else if (type_ion == GKYL_IZ_C) {
    data->NT = 30;
    data->NN = 24;
    strcpy(fname, base);  
    strcat(fname, "/adas/ioniz_c.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_c.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_c.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 6;
    data->Eiz = malloc(sizeof(double)*data->Zmax);
    static const double Eiz_loc[] = {11.3, 24.4, 47.9, 64.5, 392.1, 490.0};
    memcpy(data->Eiz, Eiz_loc, sizeof(Eiz_loc));
  }
  else if (type_ion == GKYL_IZ_N) {
    data->NT = 30;
    data->NN = 24;
    strcpy(fname, base);  
    strcat(fname, "/adas/ioniz_n.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_n.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_n.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 7;
    data->Eiz = malloc(sizeof(double)*data->Zmax);
    static const double Eiz_loc[] = {14.5, 29.6, 47.5, 77.5, 97.9, 552.1, 667.0};
    memcpy(data->Eiz, Eiz_loc, sizeof(Eiz_loc));
  }
  else if (type_ion == GKYL_IZ_O) {
    data->NT = 30;
    data->NN = 24;
    strcpy(fname, base);  
    strcat(fname, "/adas/ioniz_o.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_o.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_o.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 8;
    data->Eiz = malloc(sizeof(double)*data->Zmax);
    static const double Eiz_loc[] = {13.6, 35.1, 54.9, 77.4, 113.9, 138.1, 739.3, 871.4};
    memcpy(data->Eiz, Eiz_loc, sizeof(Eiz_loc));
  }
  else fprintf(stderr, "Incorrect ion type for ionization.");
}
 
void
read_adas_field_recomb(enum gkyl_dg_recomb_type type_ion, struct adas_field *data, const char *base) {
  char fname[4000];
  if (type_ion == GKYL_RECOMB_H) {
    strcpy(fname, base);
    data->NT = 29, data->NN = 24;
    strcpy(fname, base);  
    strcat(fname, "/adas/recomb_h.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_h.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_h.npy");
    data->logN = fopen(fname, "rb");    
    data->Zmax = 1;
   }
  else if (type_ion == GKYL_RECOMB_HE) {
    data->NT = 30;
    data->NN = 24;
    strcpy(fname, base);  
    strcat(fname, "/adas/recomb_he.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_he.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_he.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 2;
   }
  else if (type_ion == GKYL_RECOMB_LI) {
    data->NT = 25;
    data->NN = 16;
    strcpy(fname, base);  
    strcat(fname, "/adas/recomb_li.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_li.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_li.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 3;
  }
  else if (type_ion == GKYL_RECOMB_BE) {
    data->NT = 25;
    data->NN = 24;
    strcpy(fname, base);  
    strcat(fname, "/adas/recomb_be.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_be.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_be.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 4;
  }
  else if (type_ion == GKYL_RECOMB_B) {
    data->NT = 48;
    data->NN = 26;
    strcpy(fname, base);  
    strcat(fname, "/adas/recomb_b.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_b.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_b.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 5;
  }
  else if (type_ion == GKYL_RECOMB_C) {
    data->NT = 30;
    data->NN = 24;
    strcpy(fname, base);  
    strcat(fname, "/adas/recomb_c.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_c.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_c.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 6;
  }
  else if (type_ion == GKYL_RECOMB_N) {
    data->NT = 30;
    data->NN = 24;
    strcpy(fname, base);  
    strcat(fname, "/adas/recomb_n.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_n.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_n.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 7;
  }
  else if (type_ion == GKYL_RECOMB_O) {
    data->NT = 30;
    data->NN = 24;
    strcpy(fname, base);  
    strcat(fname, "/adas/recomb_o.npy");
    data->logData = fopen(fname,"rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logT_o.npy");
    data->logT = fopen(fname, "rb");
    strcpy(fname, base);  
    strcat(fname, "/adas/logN_o.npy");
    data->logN = fopen(fname, "rb"); 
    data->Zmax = 8;
  }
  else fprintf(stderr, "Incorrect ion type for recombination.");
}
