#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void cart_bmag_2x_Ser_p1( const double **psibyr, const double *psibyr2, double *bmag_out, double *br_out, double *bz_out, double scale_factorR, double scale_factorZ); 
GKYL_CU_DH void cart_bmag_2x_Ser_p2( const double **psibyr, const double *psibyr2, double *bmag_out, double *br_out, double *bz_out, double scale_factorR, double scale_factorZ); 
GKYL_CU_DH void cart_bmag_2x_Tensor_p1( const double **psibyr, const double *psibyr2, double *bmag_out, double *br_out, double *bz_out, double scale_factorR, double scale_factorZ); 
GKYL_CU_DH void cart_bmag_2x_Tensor_p2( const double **psibyr, const double *psibyr2, double *bmag_out, double *br_out, double *bz_out, double scale_factorR, double scale_factorZ); 
EXTERN_C_END 
