#pragma once 
#include <math.h> 
#include <gkyl_util.h>
#include <gkyl_mat.h>

EXTERN_C_BEG

GKYL_CU_DH void gyrokinetic_self_prim_moments_1x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *moms, const double *boundary_corrections);
GKYL_CU_DH void gyrokinetic_self_prim_moments_1x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *moms, const double *boundary_corrections);
GKYL_CU_DH void gyrokinetic_self_prim_moments_2x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *moms, const double *boundary_corrections);
GKYL_CU_DH void gyrokinetic_self_prim_moments_3x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *moms, const double *boundary_corrections);

EXTERN_C_END 
