#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
 
EXTERN_C_BEG 
 
GKYL_CU_DH void prim_lbo_copy_sol(const struct gkyl_mat *rhs, const int nc, const int vdim,
  double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq); 

EXTERN_C_END 
