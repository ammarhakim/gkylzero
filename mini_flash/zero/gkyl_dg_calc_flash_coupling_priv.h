#pragma once

#include <gkyl_array_ops.h>
#include <gkyl_dg_calc_flash_coupling.h>
#include <gkyl_flash_coupling_kernels.h>

// Function type for computing momentum flux term, div(P) - grad(p)
typedef void (*flash_coupling_mom_t)(const double *dxv, 
  const double *vlasov_pkpm_moms_l, const double *vlasov_pkpm_moms_c, const double *vlasov_pkpm_moms_r,
  const double *pkpm_div_ppar, double *flux_div_p);

// Function type for computing energy flux term, P:grad(u)
typedef void (*flash_coupling_ener_t)(const double *dxv, const double *vlasov_pkpm_moms, const double *pkpm_prim, 
  const double *pkpm_prim_surf_l, const double *pkpm_prim_surf_c, const double *pkpm_prim_surf_r,
  double *flux_p_gradu);

// Function type for computing heat flux energy transfer term, div(q)
typedef void (*flash_coupling_div_q_t)(const double *w, const double *dxv, double mass,
  const double *fl, const double *fc, const double *fr,
  const double *bvar_c, const double *pkpm_max_b,
  const double *bvar_surf_l, const double *bvar_surf_c, const double *bvar_surf_r, 
  double *flux_div_q);

struct gkyl_flash_coupling_kernels {
  flash_coupling_mom_t mom_flux_ker[3];
  flash_coupling_ener_t ener_flux_ker[3];
  flash_coupling_div_q_t div_q_flux_ker[3];
};

// Updater structure for computing the contributions
// to the Euler fluxes from the kinetic equation, div(P).
struct gkyl_dg_calc_flash_coupling {
  int cdim, pdim;
  struct gkyl_rect_grid *phase_grid, *conf_grid_kinetic;

  int num_conf_basis;

  struct gkyl_flash_coupling_kernels ker; // Kernels for computing flux contributions
};

// Kernel lists
static const flash_coupling_mom_t flash_coupling_mom_kernels[3][3] = {
  {flash_coupling_mom_1x1v_ser_p1, NULL, NULL},
  {NULL, NULL, NULL},
  {NULL, NULL, NULL}
};

static const flash_coupling_ener_t flash_coupling_ener_kernels[3][3] = {
  {flash_coupling_ener_1x1v_ser_p1, NULL, NULL},
  {NULL, NULL, NULL},
  {NULL, NULL, NULL}
};

static const flash_coupling_div_q_t flash_coupling_div_q_kernels[3][3] = {
  {flash_coupling_div_q_1x1v_ser_p1, NULL, NULL},
  {NULL, NULL, NULL},
  {NULL, NULL, NULL}
};

static const flash_coupling_mom_t
choose_flash_coupling_mom_kernel(int dir, int cdim) {
  return flash_coupling_mom_kernels[cdim-1][dir];
}

static flash_coupling_ener_t
choose_flash_coupling_ener_kernel(int dir, int cdim) {
  return flash_coupling_ener_kernels[cdim-1][dir];
}

static flash_coupling_div_q_t
choose_flash_coupling_div_q_kernel(int dir, int cdim) {
  return flash_coupling_div_q_kernels[cdim-1][dir];
}