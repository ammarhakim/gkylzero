#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void deflate_zsurf_lo_2x_ser_p1_remy(const double *fld, double *fld_deflated); 
GKYL_CU_DH void deflate_zsurf_up_2x_ser_p1_remy(const double *fld, double *fld_deflated); 
GKYL_CU_DH void deflate_zsurf_lo_2x_ser_p2_remy(const double *fld, double *fld_deflated); 
GKYL_CU_DH void deflate_zsurf_up_2x_ser_p2_remy(const double *fld, double *fld_deflated); 
EXTERN_C_END 
