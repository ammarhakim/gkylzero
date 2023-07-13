#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
 
EXTERN_C_BEG 

GKYL_CU_DH void vlasov_prim_vars_1x1v_ser_p1(const double *moms, double* prim_vars); 

GKYL_CU_DH void vlasov_prim_vars_1x2v_ser_p1(const double *moms, double* prim_vars); 

GKYL_CU_DH void vlasov_prim_vars_1x3v_ser_p1(const double *moms, double* prim_vars); 

GKYL_CU_DH void vlasov_prim_vars_2x2v_ser_p1(const double *moms, double* prim_vars); 

GKYL_CU_DH void vlasov_prim_vars_2x3v_ser_p1(const double *moms, double* prim_vars); 

GKYL_CU_DH void vlasov_prim_vars_3x3v_ser_p1(const double *moms, double* prim_vars); 


GKYL_CU_DH void gyrokinetic_prim_vars_1x1v_ser_p1(const double *moms, double* prim_vars); 

GKYL_CU_DH void gyrokinetic_prim_vars_1x2v_ser_p1(const double *moms, double* prim_vars); 

GKYL_CU_DH void gyrokinetic_prim_vars_1x3v_ser_p1(const double *moms, double* prim_vars); 

GKYL_CU_DH void gyrokinetic_prim_vars_1x2v_ser_p1(const double *moms, double* prim_vars); 

GKYL_CU_DH void gyrokinetic_prim_vars_2x2v_ser_p1(const double *moms, double* prim_vars); 

GKYL_CU_DH void gyrokinetic_prim_vars_3x2v_ser_p1(const double *moms, double* prim_vars); 

EXTERN_C_END 
