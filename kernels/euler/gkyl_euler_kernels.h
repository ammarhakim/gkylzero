#pragma once 
#include <math.h> 
#include <gkyl_util.h> 
EXTERN_C_BEG 

GKYL_CU_DH void euler_pressure_1x_ser_p1(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_1x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_1x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_1x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfx_1x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pkpm_prim_vars_1x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, double* u_i, double* u_perp_i, double* rhou_perp_i, double* p_perp, double* p_ij); 
GKYL_CU_DH void euler_pkpm_source_1x_ser_p1(const double *qmem, const double *nu, const double *nu_vth_sq, const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *rhou_perp_i, const double *p_perp, double* out); 
GKYL_CU_DH void euler_pkpm_p_force_1x_ser_p1(const double *bvar, const double *div_p, const double *vlasov_pkpm_moms, 
  const double *p_perp, const double *div_b, double* p_force); 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p1(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *u_perp_i, const double *p_perp, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_x_1x_ser_p1(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_x_1x1v_ser_p1(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfx_1x_ser_p1(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pressure_1x_ser_p2(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_1x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_1x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_1x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfx_1x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pkpm_prim_vars_1x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, double* u_i, double* u_perp_i, double* rhou_perp_i, double* p_perp, double* p_ij); 
GKYL_CU_DH void euler_pkpm_source_1x_ser_p2(const double *qmem, const double *nu, const double *nu_vth_sq, const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *rhou_perp_i, const double *p_perp, double* out); 
GKYL_CU_DH void euler_pkpm_p_force_1x_ser_p2(const double *bvar, const double *div_p, const double *vlasov_pkpm_moms, 
  const double *p_perp, const double *div_b, double* p_force); 
GKYL_CU_DH double euler_pkpm_vol_1x_ser_p2(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *u_perp_i, const double *p_perp, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_x_1x_ser_p2(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_x_1x1v_ser_p2(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfx_1x_ser_p2(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pressure_2x_ser_p1(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_2x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_2x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfy_2x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_2x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfx_2x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfy_2x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pkpm_prim_vars_2x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, double* u_i, double* u_perp_i, double* rhou_perp_i, double* p_perp, double* p_ij); 
GKYL_CU_DH void euler_pkpm_source_2x_ser_p1(const double *qmem, const double *nu, const double *nu_vth_sq, const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *rhou_perp_i, const double *p_perp, double* out); 
GKYL_CU_DH void euler_pkpm_p_force_2x_ser_p1(const double *bvar, const double *div_p, const double *vlasov_pkpm_moms, 
  const double *p_perp, const double *div_b, double* p_force); 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p1(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *u_perp_i, const double *p_perp, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_x_2x_ser_p1(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_x_2x1v_ser_p1(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfx_2x_ser_p1(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_y_2x_ser_p1(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_y_2x1v_ser_p1(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfy_2x_ser_p1(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pressure_2x_ser_p2(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_2x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_2x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfy_2x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_2x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfx_2x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfy_2x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pkpm_prim_vars_2x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, double* u_i, double* u_perp_i, double* rhou_perp_i, double* p_perp, double* p_ij); 
GKYL_CU_DH void euler_pkpm_source_2x_ser_p2(const double *qmem, const double *nu, const double *nu_vth_sq, const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *rhou_perp_i, const double *p_perp, double* out); 
GKYL_CU_DH void euler_pkpm_p_force_2x_ser_p2(const double *bvar, const double *div_p, const double *vlasov_pkpm_moms, 
  const double *p_perp, const double *div_b, double* p_force); 
GKYL_CU_DH double euler_pkpm_vol_2x_ser_p2(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *u_perp_i, const double *p_perp, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_x_2x_ser_p2(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_x_2x1v_ser_p2(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfx_2x_ser_p2(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_y_2x_ser_p2(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_y_2x1v_ser_p2(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfy_2x_ser_p2(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pressure_3x_ser_p1(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_3x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_3x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfy_3x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfz_3x_ser_p1(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_3x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfx_3x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfy_3x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfz_3x_ser_p1(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pkpm_prim_vars_3x_ser_p1(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, double* u_i, double* u_perp_i, double* rhou_perp_i, double* p_perp, double* p_ij); 
GKYL_CU_DH void euler_pkpm_source_3x_ser_p1(const double *qmem, const double *nu, const double *nu_vth_sq, const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *rhou_perp_i, const double *p_perp, double* out); 
GKYL_CU_DH void euler_pkpm_p_force_3x_ser_p1(const double *bvar, const double *div_p, const double *vlasov_pkpm_moms, 
  const double *p_perp, const double *div_b, double* p_force); 
GKYL_CU_DH double euler_pkpm_vol_3x_ser_p1(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *u_perp_i, const double *p_perp, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_x_3x_ser_p1(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_x_3x1v_ser_p1(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfx_3x_ser_p1(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_y_3x_ser_p1(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_y_3x1v_ser_p1(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfy_3x_ser_p1(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_z_3x_ser_p1(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_z_3x1v_ser_p1(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfz_3x_ser_p1(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pressure_3x_ser_p2(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_3x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_3x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfy_3x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfz_3x_ser_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_3x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfx_3x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfy_3x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfz_3x_ser_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pkpm_prim_vars_3x_ser_p2(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, double* u_i, double* u_perp_i, double* rhou_perp_i, double* p_perp, double* p_ij); 
GKYL_CU_DH void euler_pkpm_source_3x_ser_p2(const double *qmem, const double *nu, const double *nu_vth_sq, const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *rhou_perp_i, const double *p_perp, double* out); 
GKYL_CU_DH void euler_pkpm_p_force_3x_ser_p2(const double *bvar, const double *div_p, const double *vlasov_pkpm_moms, 
  const double *p_perp, const double *div_b, double* p_force); 
GKYL_CU_DH double euler_pkpm_vol_3x_ser_p2(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *u_perp_i, const double *p_perp, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_x_3x_ser_p2(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_x_3x1v_ser_p2(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfx_3x_ser_p2(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_y_3x_ser_p2(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_y_3x1v_ser_p2(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfy_3x_ser_p2(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_z_3x_ser_p2(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_z_3x1v_ser_p2(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfz_3x_ser_p2(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pressure_1x_tensor_p2(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_1x_tensor_p2(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_1x_tensor_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_1x_tensor_p2(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfx_1x_tensor_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pkpm_prim_vars_1x_tensor_p2(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, double* u_i, double* u_perp_i, double* rhou_perp_i, double* p_perp, double* p_ij); 
GKYL_CU_DH void euler_pkpm_source_1x_tensor_p2(const double *qmem, const double *nu, const double *nu_vth_sq, const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *rhou_perp_i, const double *p_perp, double* out); 
GKYL_CU_DH void euler_pkpm_p_force_1x_tensor_p2(const double *bvar, const double *div_p, const double *vlasov_pkpm_moms, 
  const double *p_perp, const double *div_b, double* p_force); 
GKYL_CU_DH double euler_pkpm_vol_1x_tensor_p2(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *u_perp_i, const double *p_perp, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_x_1x_tensor_p2(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_x_1x1v_tensor_p2(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfx_1x_tensor_p2(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pressure_2x_tensor_p2(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_2x_tensor_p2(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_2x_tensor_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfy_2x_tensor_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_2x_tensor_p2(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfx_2x_tensor_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfy_2x_tensor_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pkpm_prim_vars_2x_tensor_p2(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, double* u_i, double* u_perp_i, double* rhou_perp_i, double* p_perp, double* p_ij); 
GKYL_CU_DH void euler_pkpm_source_2x_tensor_p2(const double *qmem, const double *nu, const double *nu_vth_sq, const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *rhou_perp_i, const double *p_perp, double* out); 
GKYL_CU_DH void euler_pkpm_p_force_2x_tensor_p2(const double *bvar, const double *div_p, const double *vlasov_pkpm_moms, 
  const double *p_perp, const double *div_b, double* p_force); 
GKYL_CU_DH double euler_pkpm_vol_2x_tensor_p2(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *u_perp_i, const double *p_perp, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_x_2x_tensor_p2(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_x_2x1v_tensor_p2(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfx_2x_tensor_p2(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_y_2x_tensor_p2(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_y_2x1v_tensor_p2(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfy_2x_tensor_p2(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pressure_3x_tensor_p2(const double gas_gamma, const double *uvar, const double *statevec, double* pressure); 
GKYL_CU_DH double euler_vol_3x_tensor_p2(const double *w, const double *dxv, const double gas_gamma, const double *uvar, const double *pvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfx_3x_tensor_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfy_3x_tensor_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_surfz_3x_tensor_p2(const double *w, const double *dxv, const double gas_gamma, const double *ul, const double *uc, const double *ur, const double *pl, const double *pc, const double *pr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH double euler_iso_vol_3x_tensor_p2(const double *w, const double *dxv, const double vth, const double *uvar, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfx_3x_tensor_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfy_3x_tensor_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_iso_surfz_3x_tensor_p2(const double *w, const double *dxv, const double vth, const double *uvarl, const double *uvarc, const double *uvarr, const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

GKYL_CU_DH void euler_pkpm_prim_vars_3x_tensor_p2(const double *bvar, const double *vlasov_pkpm_moms, const double *statevec, double* u_i, double* u_perp_i, double* rhou_perp_i, double* p_perp, double* p_ij); 
GKYL_CU_DH void euler_pkpm_source_3x_tensor_p2(const double *qmem, const double *nu, const double *nu_vth_sq, const double *vlasov_pkpm_moms, const double *euler_pkpm, const double *rhou_perp_i, const double *p_perp, double* out); 
GKYL_CU_DH void euler_pkpm_p_force_3x_tensor_p2(const double *bvar, const double *div_p, const double *vlasov_pkpm_moms, 
  const double *p_perp, const double *div_b, double* p_force); 
GKYL_CU_DH double euler_pkpm_vol_3x_tensor_p2(const double *w, const double *dxv, const double *u_i, const double *div_p, const double *u_perp_i, const double *p_perp, const double *statevec, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_x_3x_tensor_p2(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_x_3x1v_tensor_p2(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfx_3x_tensor_p2(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_y_3x_tensor_p2(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_y_3x1v_tensor_p2(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfy_3x_tensor_p2(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 
GKYL_CU_DH void euler_pkpm_recovery_z_3x_tensor_p2(const double *dxv, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *p_ijl, const double *p_ijc, const double *p_ijr, 
    double* div_b, double* bb_grad_u, double* div_p) ; 
GKYL_CU_DH void vlasov_pkpm_upwind_p_force_z_3x1v_tensor_p2(const double *w, const double *dxv, double mass, 
    const double *bvarl, const double *bvarc, const double *bvarr, 
    const double *fl, const double *fc, const double *fr, double* out) ; 
GKYL_CU_DH void euler_pkpm_surfz_3x_tensor_p2(const double *w, const double *dxv, 
    const double *u_il, const double *u_ic, const double *u_ir, 
    const double *u_perp_il, const double *u_perp_ic, const double *u_perp_ir, 
    const double *p_perpl, const double *p_perpc, const double *p_perpr, 
    const double *statevecl, const double *statevecc, const double *statevecr, double* GKYL_RESTRICT out); 

EXTERN_C_END 
