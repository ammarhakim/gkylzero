#pragma once 
#include <math.h> 
#include <gkyl_util.h>
#include <gkyl_mat.h>

EXTERN_C_BEG

//GKYL_CU_DH void vlasov_self_prim_moments_1x1v_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *m0, const double *m1,
//const double *m2, const double *cM, const double *cE, double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq);
GKYL_CU_DH void vlasov_self_prim_moments_1x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *m0, const double *m1,
  const double *m2, const double *cM, const double *cE, double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq);

//GKYL_CU_DH void vlasov_self_prim_moments_1x2v_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *m0, const double *m1,
//const double *m2, const double *cM, const double *cE, double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq);
GKYL_CU_DH void vlasov_self_prim_moments_1x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *m0, const double *m1,
  const double *m2, const double *cM, const double *cE, double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq);

//GKYL_CU_DH void vlasov_self_prim_moments_1x3v_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *m0, const double *m1,
//const double *m2, const double *cM, const double *cE, double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq);
GKYL_CU_DH void vlasov_self_prim_moments_1x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *m0, const double *m1,
  const double *m2, const double *cM, const double *cE, double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq);

//GKYL_CU_DH void vlasov_self_prim_moments_2x2v_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *m0, const double *m1,
//const double *m2, const double *cM, const double *cE, double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq);
GKYL_CU_DH void vlasov_self_prim_moments_2x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *m0, const double *m1,
  const double *m2, const double *cM, const double *cE, double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq);

//GKYL_CU_DH void vlasov_self_prim_moments_2x3v_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *m0, const double *m1,
//const double *m2, const double *cM, const double *cE, double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq);

//GKYL_CU_DH void vlasov_self_prim_moments_3x3v_ser_p1(struct gkyl_mat *A, struct gkyl_mat *rhs, const double *m0, const double *m1,
//const double *m2, const double *cM, const double *cE, double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq);


EXTERN_C_END 
