#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void euler_pkpm_prim_vars_1x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, 
  double* u_i, double* p_ij, double* T_ij, double* rho_inv, double* T_perp_over_m, double* T_perp_over_m_inv); 
GKYL_CU_DH void euler_pkpm_source_1x_ser_p1(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* out); 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p1(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_mirror_force_1x1v_ser_p1(const double *w, const double *dxv, 
  const double* T_perp_over_m, const double* T_perp_over_m_inv, 
  const double *nu_vthsq, const double* pkpm_accel_vars, 
  const double* f, const double* F_k_p_1, 
  double* g_dist_source, double* F_k_m_1); 
GKYL_CU_DH void euler_pkpm_recovery_x_1x_ser_p1(const double *dxv, double nuHyp, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *p_ijl, const double *p_ijc, const double *p_ijr, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr, 
      const double *statevecl, const double *statevecc, const double *statevecr, 
      const double *pkpm_div_ppar, const double *rho_inv, const double *T_perp_over_m, 
      const double *T_perp_over_m_inv, const double *nu, 
      double* div_p, double* pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfx_1x_ser_p1(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *T_ijl, const double *T_ijc, const double *T_ijr, 
      const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pkpm_prim_vars_1x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, 
  double* u_i, double* p_ij, double* T_ij, double* rho_inv, double* T_perp_over_m, double* T_perp_over_m_inv); 
GKYL_CU_DH void euler_pkpm_source_1x_ser_p2(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* out); 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p2(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_mirror_force_1x1v_ser_p2(const double *w, const double *dxv, 
  const double* T_perp_over_m, const double* T_perp_over_m_inv, 
  const double *nu_vthsq, const double* pkpm_accel_vars, 
  const double* f, const double* F_k_p_1, 
  double* g_dist_source, double* F_k_m_1); 
GKYL_CU_DH void euler_pkpm_recovery_x_1x_ser_p2(const double *dxv, double nuHyp, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *p_ijl, const double *p_ijc, const double *p_ijr, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr, 
      const double *statevecl, const double *statevecc, const double *statevecr, 
      const double *pkpm_div_ppar, const double *rho_inv, const double *T_perp_over_m, 
      const double *T_perp_over_m_inv, const double *nu, 
      double* div_p, double* pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfx_1x_ser_p2(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *T_ijl, const double *T_ijc, const double *T_ijr, 
      const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void pkpm_dist_mirror_force_1x1v_tensor_p2(const double *w, const double *dxv, 
  const double* T_perp_over_m, const double* T_perp_over_m_inv, 
  const double *nu_vthsq, const double* pkpm_accel_vars, 
  const double* f, const double* F_k_p_1, 
  double* g_dist_source, double* F_k_m_1); 

GKYL_CU_DH void euler_pkpm_prim_vars_2x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, 
  double* u_i, double* p_ij, double* T_ij, double* rho_inv, double* T_perp_over_m, double* T_perp_over_m_inv); 
GKYL_CU_DH void euler_pkpm_source_2x_ser_p1(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* out); 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p1(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_mirror_force_2x1v_ser_p1(const double *w, const double *dxv, 
  const double* T_perp_over_m, const double* T_perp_over_m_inv, 
  const double *nu_vthsq, const double* pkpm_accel_vars, 
  const double* f, const double* F_k_p_1, 
  double* g_dist_source, double* F_k_m_1); 
GKYL_CU_DH void euler_pkpm_recovery_x_2x_ser_p1(const double *dxv, double nuHyp, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *p_ijl, const double *p_ijc, const double *p_ijr, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr, 
      const double *statevecl, const double *statevecc, const double *statevecr, 
      const double *pkpm_div_ppar, const double *rho_inv, const double *T_perp_over_m, 
      const double *T_perp_over_m_inv, const double *nu, 
      double* div_p, double* pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfx_2x_ser_p1(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *T_ijl, const double *T_ijc, const double *T_ijr, 
      const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_y_2x_ser_p1(const double *dxv, double nuHyp, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *p_ijl, const double *p_ijc, const double *p_ijr, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr, 
      const double *statevecl, const double *statevecc, const double *statevecr, 
      const double *pkpm_div_ppar, const double *rho_inv, const double *T_perp_over_m, 
      const double *T_perp_over_m_inv, const double *nu, 
      double* div_p, double* pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfy_2x_ser_p1(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *T_ijl, const double *T_ijc, const double *T_ijr, 
      const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pkpm_prim_vars_2x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, 
  double* u_i, double* p_ij, double* T_ij, double* rho_inv, double* T_perp_over_m, double* T_perp_over_m_inv); 
GKYL_CU_DH void euler_pkpm_source_2x_ser_p2(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* out); 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p2(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_mirror_force_2x1v_ser_p2(const double *w, const double *dxv, 
  const double* T_perp_over_m, const double* T_perp_over_m_inv, 
  const double *nu_vthsq, const double* pkpm_accel_vars, 
  const double* f, const double* F_k_p_1, 
  double* g_dist_source, double* F_k_m_1); 
GKYL_CU_DH void euler_pkpm_recovery_x_2x_ser_p2(const double *dxv, double nuHyp, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *p_ijl, const double *p_ijc, const double *p_ijr, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr, 
      const double *statevecl, const double *statevecc, const double *statevecr, 
      const double *pkpm_div_ppar, const double *rho_inv, const double *T_perp_over_m, 
      const double *T_perp_over_m_inv, const double *nu, 
      double* div_p, double* pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfx_2x_ser_p2(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *T_ijl, const double *T_ijc, const double *T_ijr, 
      const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_y_2x_ser_p2(const double *dxv, double nuHyp, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *p_ijl, const double *p_ijc, const double *p_ijr, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr, 
      const double *statevecl, const double *statevecc, const double *statevecr, 
      const double *pkpm_div_ppar, const double *rho_inv, const double *T_perp_over_m, 
      const double *T_perp_over_m_inv, const double *nu, 
      double* div_p, double* pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfy_2x_ser_p2(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *T_ijl, const double *T_ijc, const double *T_ijr, 
      const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pkpm_prim_vars_3x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, 
  double* u_i, double* p_ij, double* T_ij, double* rho_inv, double* T_perp_over_m, double* T_perp_over_m_inv); 
GKYL_CU_DH void euler_pkpm_source_3x_ser_p1(const double *qmem, const double *vlasov_pkpm_moms, const double *euler_pkpm, double* out); 
GKYL_CU_DH double euler_pkpm_vol_3x_ser_p1(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void pkpm_dist_mirror_force_3x1v_ser_p1(const double *w, const double *dxv, 
  const double* T_perp_over_m, const double* T_perp_over_m_inv, 
  const double *nu_vthsq, const double* pkpm_accel_vars, 
  const double* f, const double* F_k_p_1, 
  double* g_dist_source, double* F_k_m_1); 
GKYL_CU_DH void euler_pkpm_recovery_x_3x_ser_p1(const double *dxv, double nuHyp, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *p_ijl, const double *p_ijc, const double *p_ijr, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr, 
      const double *statevecl, const double *statevecc, const double *statevecr, 
      const double *pkpm_div_ppar, const double *rho_inv, const double *T_perp_over_m, 
      const double *T_perp_over_m_inv, const double *nu, 
      double* div_p, double* pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfx_3x_ser_p1(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *T_ijl, const double *T_ijc, const double *T_ijr, 
      const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_y_3x_ser_p1(const double *dxv, double nuHyp, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *p_ijl, const double *p_ijc, const double *p_ijr, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr, 
      const double *statevecl, const double *statevecc, const double *statevecr, 
      const double *pkpm_div_ppar, const double *rho_inv, const double *T_perp_over_m, 
      const double *T_perp_over_m_inv, const double *nu, 
      double* div_p, double* pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfy_3x_ser_p1(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *T_ijl, const double *T_ijc, const double *T_ijr, 
      const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_z_3x_ser_p1(const double *dxv, double nuHyp, 
      const double *bvarl, const double *bvarc, const double *bvarr, 
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *p_ijl, const double *p_ijc, const double *p_ijr, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr, 
      const double *statevecl, const double *statevecc, const double *statevecr, 
      const double *pkpm_div_ppar, const double *rho_inv, const double *T_perp_over_m, 
      const double *T_perp_over_m_inv, const double *nu, 
      double* div_p, double* pkpm_accel_vars); 
GKYL_CU_DH void euler_pkpm_surfz_3x_ser_p1(const double *w, const double *dxv, 
      const double *vlasov_pkpm_momsl, const double *vlasov_pkpm_momsc, const double *vlasov_pkpm_momsr,
      const double *u_il, const double *u_ic, const double *u_ir, 
      const double *T_ijl, const double *T_ijc, const double *T_ijr, 
      const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

EXTERN_C_END 
