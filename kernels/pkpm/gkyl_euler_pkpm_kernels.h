#pragma once 
#include <math.h> 
#include <gkyl_mat.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void pkpm_vars_pressure_1x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH int pkpm_vars_set_1x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_1x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim); 
GKYL_CU_DH void pkpm_vars_integrated_1x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* pkpm_int_vars); 
GKYL_CU_DH void euler_pkpm_source_1x_ser_p1(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* out); 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p1(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_1x_ser_p1(const double *dxv, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *priml, const double *primc, const double *primr, 
      const double *nu, double* GKYL_RESTRICT pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfx_1x_ser_p1(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *priml, const double *primc, const double *primr,
      const double *p_ijl, const double *p_ijc, const double *p_ijr,
      const double *euler_pkpml, const double *euler_pkpmc, const double *euler_pkpmr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_1x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH int pkpm_vars_set_1x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_1x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim); 
GKYL_CU_DH void pkpm_vars_integrated_1x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* pkpm_int_vars); 
GKYL_CU_DH void euler_pkpm_source_1x_ser_p2(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* out); 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p2(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_1x_ser_p2(const double *dxv, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *priml, const double *primc, const double *primr, 
      const double *nu, double* GKYL_RESTRICT pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfx_1x_ser_p2(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *priml, const double *primc, const double *primr,
      const double *p_ijl, const double *p_ijc, const double *p_ijr,
      const double *euler_pkpml, const double *euler_pkpmc, const double *euler_pkpmr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_2x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH int pkpm_vars_set_2x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_2x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim); 
GKYL_CU_DH void pkpm_vars_integrated_2x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* pkpm_int_vars); 
GKYL_CU_DH void euler_pkpm_source_2x_ser_p1(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* out); 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p1(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_2x_ser_p1(const double *dxv, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *priml, const double *primc, const double *primr, 
      const double *nu, double* GKYL_RESTRICT pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfx_2x_ser_p1(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *priml, const double *primc, const double *primr,
      const double *p_ijl, const double *p_ijc, const double *p_ijr,
      const double *euler_pkpml, const double *euler_pkpmc, const double *euler_pkpmr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_2x_ser_p1(const double *dxv, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *priml, const double *primc, const double *primr, 
      const double *nu, double* GKYL_RESTRICT pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfy_2x_ser_p1(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *priml, const double *primc, const double *primr,
      const double *p_ijl, const double *p_ijc, const double *p_ijr,
      const double *euler_pkpml, const double *euler_pkpmc, const double *euler_pkpmr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_2x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH int pkpm_vars_set_2x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_2x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim); 
GKYL_CU_DH void pkpm_vars_integrated_2x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* pkpm_int_vars); 
GKYL_CU_DH void euler_pkpm_source_2x_ser_p2(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* out); 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p2(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_2x_ser_p2(const double *dxv, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *priml, const double *primc, const double *primr, 
      const double *nu, double* GKYL_RESTRICT pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfx_2x_ser_p2(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *priml, const double *primc, const double *primr,
      const double *p_ijl, const double *p_ijc, const double *p_ijr,
      const double *euler_pkpml, const double *euler_pkpmc, const double *euler_pkpmr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_2x_ser_p2(const double *dxv, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *priml, const double *primc, const double *primr, 
      const double *nu, double* GKYL_RESTRICT pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfy_2x_ser_p2(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *priml, const double *primc, const double *primr,
      const double *p_ijl, const double *p_ijc, const double *p_ijr,
      const double *euler_pkpml, const double *euler_pkpmc, const double *euler_pkpmr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_3x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH int pkpm_vars_set_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_3x_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim); 
GKYL_CU_DH void pkpm_vars_integrated_3x_ser_p1(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* pkpm_int_vars); 
GKYL_CU_DH void euler_pkpm_source_3x_ser_p1(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* out); 
GKYL_CU_DH double euler_pkpm_vol_3x_ser_p1(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_3x_ser_p1(const double *dxv, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *priml, const double *primc, const double *primr, 
      const double *nu, double* GKYL_RESTRICT pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfx_3x_ser_p1(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *priml, const double *primc, const double *primr,
      const double *p_ijl, const double *p_ijc, const double *p_ijr,
      const double *euler_pkpml, const double *euler_pkpmc, const double *euler_pkpmr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_3x_ser_p1(const double *dxv, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *priml, const double *primc, const double *primr, 
      const double *nu, double* GKYL_RESTRICT pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfy_3x_ser_p1(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *priml, const double *primc, const double *primr,
      const double *p_ijl, const double *p_ijc, const double *p_ijr,
      const double *euler_pkpml, const double *euler_pkpmc, const double *euler_pkpmr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_z_3x_ser_p1(const double *dxv, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *priml, const double *primc, const double *primr, 
      const double *nu, double* GKYL_RESTRICT pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfz_3x_ser_p1(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *priml, const double *primc, const double *primr,
      const double *p_ijl, const double *p_ijc, const double *p_ijr,
      const double *euler_pkpml, const double *euler_pkpmc, const double *euler_pkpmr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_3x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH int pkpm_vars_set_3x_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_3x_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim); 
GKYL_CU_DH void pkpm_vars_integrated_3x_ser_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* pkpm_int_vars); 
GKYL_CU_DH void euler_pkpm_source_3x_ser_p2(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* out); 
GKYL_CU_DH double euler_pkpm_vol_3x_ser_p2(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_3x_ser_p2(const double *dxv, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *priml, const double *primc, const double *primr, 
      const double *nu, double* GKYL_RESTRICT pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfx_3x_ser_p2(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *priml, const double *primc, const double *primr,
      const double *p_ijl, const double *p_ijc, const double *p_ijr,
      const double *euler_pkpml, const double *euler_pkpmc, const double *euler_pkpmr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_3x_ser_p2(const double *dxv, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *priml, const double *primc, const double *primr, 
      const double *nu, double* GKYL_RESTRICT pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfy_3x_ser_p2(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *priml, const double *primc, const double *primr,
      const double *p_ijl, const double *p_ijc, const double *p_ijr,
      const double *euler_pkpml, const double *euler_pkpmc, const double *euler_pkpmr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_z_3x_ser_p2(const double *dxv, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *priml, const double *primc, const double *primr, 
      const double *nu, double* GKYL_RESTRICT pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfz_3x_ser_p2(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *priml, const double *primc, const double *primr,
      const double *p_ijl, const double *p_ijc, const double *p_ijr,
      const double *euler_pkpml, const double *euler_pkpmc, const double *euler_pkpmr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_vars_pressure_2x_tensor_p2(const double *bvar, const double *vlasov_pkpm_moms, double* GKYL_RESTRICT p_ij); 
GKYL_CU_DH int pkpm_vars_set_2x_tensor_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
  const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *p_ij, const double *pkpm_div_ppar); 
GKYL_CU_DH void pkpm_vars_copy_2x_tensor_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT prim); 
GKYL_CU_DH void pkpm_vars_integrated_2x_tensor_p2(const double *vlasov_pkpm_moms, const double *euler_pkpm, 
  const double* prim, double* pkpm_int_vars); 
GKYL_CU_DH void euler_pkpm_source_2x_tensor_p2(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* out); 
GKYL_CU_DH double euler_pkpm_vol_2x_tensor_p2(const double *w, const double *dxv, const double *prim, const double *p_ij, const double *euler_pkpm, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_x_2x_tensor_p2(const double *dxv, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *priml, const double *primc, const double *primr, 
      const double *nu, double* GKYL_RESTRICT pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfx_2x_tensor_p2(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *priml, const double *primc, const double *primr,
      const double *p_ijl, const double *p_ijc, const double *p_ijr,
      const double *euler_pkpml, const double *euler_pkpmc, const double *euler_pkpmr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_vars_accel_y_2x_tensor_p2(const double *dxv, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *priml, const double *primc, const double *primr, 
      const double *nu, double* GKYL_RESTRICT pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfy_2x_tensor_p2(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *priml, const double *primc, const double *primr,
      const double *p_ijl, const double *p_ijc, const double *p_ijr,
      const double *euler_pkpml, const double *euler_pkpmc, const double *euler_pkpmr, double* GKYL_RESTRICT out); 

EXTERN_C_END 
