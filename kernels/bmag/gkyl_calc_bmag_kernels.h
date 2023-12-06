#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void bmag_2x_Ser_p1( const double **psibyr, const double *psibyr2, const double *bphi, double *bmagout, double scale_factorR, double scale_factorZ); 
GKYL_CU_DH void bmag_2x_Ser_p2( const double **psibyr, const double *psibyr2, const double *bphi, double *bmagout, double scale_factorR, double scale_factorZ); 
EXTERN_C_END 
