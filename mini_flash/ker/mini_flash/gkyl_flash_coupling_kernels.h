#pragma once

#include <gkyl_util.h>

void flash_coupling_mom_1x1v_ser_p1(const double *dxv, 
  const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
  const double *pkpm_div_ppar, double *flux_div_p);

void flash_coupling_ener_1x1v_ser_p1(const double *dxv, const double *vlasov_pkpm_moms, const double *pkpm_prim, 
  const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r,
  double *flux_p_gradu);

void flash_coupling_div_q_1x1v_ser_p1(const double *w, const double *dxv, double mass,
  const double *fl, const double *fc, const double *fr,
  const double *bvar_c, const double *pkpm_max_b,
  const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
  double *flux_div_q);
