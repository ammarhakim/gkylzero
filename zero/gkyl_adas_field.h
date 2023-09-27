#pragma once

typedef struct adas_field {
  FILE *logData;
  FILE *logT;
  FILE *logN;
  long NT;
  long NN; 
  int Zmax;
  struct gkyl_array fld;
  double *Eiz;
} adas_field;
