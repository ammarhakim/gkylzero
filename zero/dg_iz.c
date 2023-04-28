#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_iz.h>
#include <gkyl_dg_iz_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

gkyl_dg_iz*
gkyl_dg_iz_new(const struct gkyl_basis* cbasis, 
  double elem_charge, double mass_elc, enum gkyl_dg_iz_type type_ion, 
  bool use_gpu)
{
  gkyl_dg_iz *up = gkyl_malloc(sizeof(gkyl_dg_iz));

  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  up->cdim = cdim;
  up->poly_order = poly_order;

  up->elem_charge = elem_charge;
  up->mass_elc = mass_elc;
  if (type_ion == GKYL_H) {
    up->E = 13.6;
    up->P = 0.0;
    up->A = 0.291e-7;
    up->K = 0.39;
    up->X = 0.232;
  }
  else if (type_ion == GKYL_AR) {
    up->E = 15.8;
    up->P = 1.0;
    up->A = 0.599e-7;
    up->K = 0.26;
    up->X = 0.136;
  }

  const gkyl_iz_temp_kern_list *iz_temp_kernels;
  const gkyl_iz_react_rate_kern_list *iz_react_rate_kernels;

  up->iz_temp = CK(iz_temp_kernels, cdim, poly_order);
  up->react_rate = CK(iz_react_rate_kernels, cdim, poly_order);
  
  return up;
}

void gkyl_dg_iz_temp(const struct gkyl_dg_iz *iz,
  const struct gkyl_range *update_rng, const struct gkyl_array *vth_sq_elc,
  struct gkyl_array *vth_sq_iz)
{
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_rng);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(update_rng, iter.idx);
    const double *vth_sq_elc_d = gkyl_array_cfetch(vth_sq_elc, loc);
    double *vth_sq_iz_d = gkyl_array_fetch(vth_sq_iz, loc);

    iz->iz_temp(iz->elem_charge, iz->mass_elc, iz->E, vth_sq_elc_d, vth_sq_iz_d);
  }
}

void gkyl_dg_iz_react_rate(const struct gkyl_dg_iz *iz,
  const struct gkyl_range *update_rng, const struct gkyl_range *phase_rng, 
  const struct gkyl_array *n_neut, const struct gkyl_array *vth_sq_neut, const struct gkyl_array *vth_sq_elc,
  struct gkyl_array *cflrate, struct gkyl_array *coef_iz)
{
  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<update_rng->ndim; ++d) rem_dir[d] = 1;

  gkyl_range_iter_init(&conf_iter, update_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(update_rng, conf_iter.idx);

    // Calculate these fields from M0, M1, M2
    // Reference proj_maxwellian_on_basis
    const double *n_neut_d = gkyl_array_cfetch(n_neut, loc);
    const double *vth_sq_neut_d = gkyl_array_cfetch(vth_sq_neut, loc);
    const double *vth_sq_elc_d = gkyl_array_cfetch(vth_sq_elc, loc);
    double *coef_iz_d = gkyl_array_fetch(coef_iz, loc);

    double cflr = iz->react_rate(iz->elem_charge, iz->mass_elc, 
      iz->E, iz->A, iz->K, iz->P, iz->X, 
      n_neut_d, vth_sq_neut_d, vth_sq_elc_d, 
      coef_iz_d);

    gkyl_range_deflate(&vel_rng, phase_rng, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
    // cfl associated with reaction is a *phase space* cfl
    // Need to loop over velocity space for each configuration space cell
    // to get total cfl rate in each phase space cell
    while (gkyl_range_iter_next(&vel_iter)) {
      long cfl_idx = gkyl_range_idx(&vel_rng, vel_iter.idx);
      double *cflrate_d = gkyl_array_fetch(cflrate, cfl_idx);
      cflrate_d[0] += cflr; // frequencies are additive
    }
  }
}

void
gkyl_dg_iz_release(gkyl_dg_iz* iz)
{
  free(iz);
}
