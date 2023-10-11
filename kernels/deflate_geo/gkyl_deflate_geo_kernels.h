#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void deflate_geo_1x_Ser_p1(const double *fld, double *fld_deflated); 
GKYL_CU_DH void deflate_geo_1x_Ser_p2(const double *fld, double *fld_deflated); 
GKYL_CU_DH void deflate_geo_2x_Ser_p1(const double *fld, double *fld_deflated); 
GKYL_CU_DH void deflate_geo_2x_Ser_p2(const double *fld, double *fld_deflated); 
EXTERN_C_END 
