#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
 
EXTERN_C_BEG 
 
GKYL_CU_DH void prim_lbo_copy_sol(const struct gkyl_mat *rhs, const int nc, const int vdim, double* GKYL_RESTRICT u, double* GKYL_RESTRICT vtSq); 
 
GKYL_CU_DH void vlasov_self_prim_moments_1x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x1v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_1x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_1x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_2x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x2v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_2x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_3x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_3x3v_ser_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 


GKYL_CU_DH void vlasov_self_prim_moments_1x1v_max_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x1v_max_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_1x2v_max_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x2v_max_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_1x3v_max_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x3v_max_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_2x2v_max_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x2v_max_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_2x3v_max_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x3v_max_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_3x3v_max_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_3x3v_max_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 


GKYL_CU_DH void vlasov_self_prim_moments_1x1v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x1v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_1x2v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x2v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_1x3v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_1x3v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_2x2v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x2v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_2x3v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_2x3v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 

GKYL_CU_DH void vlasov_self_prim_moments_3x3v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs,
  const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 
GKYL_CU_DH void vlasov_cross_prim_moments_3x3v_tensor_p2(struct gkyl_mat *A, struct gkyl_mat *rhs, const double betaGreenep1, const double mSelf,
  const double nuSelf, const double *uSelf, const double *vtSqSelf, const double mOther, const double nuOther, const double *uOther,
  const double *vtSqOther, const double *m0, const double *m1, const double *m2, const double *cM, const double *cE); 


EXTERN_C_END 
