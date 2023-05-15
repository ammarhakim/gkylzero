#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_iz.h>
#include <gkyl_dg_iz_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

gkyl_dg_iz*
gkyl_dg_iz_new(struct gkyl_basis* cbasis, struct gkyl_basis* pbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng, 
  double elem_charge, double mass_elc, enum gkyl_dg_iz_type type_ion, 
  bool is_gk, bool use_gpu)
{
  gkyl_dg_iz *up = gkyl_malloc(sizeof(gkyl_dg_iz));

  int cdim = cbasis->ndim;
  int pdim = pbasis->ndim;
  int poly_order = cbasis->poly_order;
  up->cdim = cdim;
  up->poly_order = poly_order;
  up->basis = cbasis;
  up->conf_rng = conf_rng;
  up->phase_rng = phase_rng; 

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

  // Establish vdim for gk and vlasov species
  // will be either 1 & 1 or 2 & 3
  // for gk and vlasov, respectively
  int vdim_gk, vdim_vl; 
  int vdim = pdim - up->cdim;
  if (is_gk) {
    if (vdim == 1) {
      vdim_gk = vdim;
      vdim_vl = vdim;
    }
    else {
      vdim_gk = vdim;
      vdim_vl = vdim+1;
    }
  }
  else if (vdim == 1) {
    vdim_gk = vdim;
    vdim_vl = vdim;
  }
  else {
    vdim_gk = vdim-1;
    vdim_vl = vdim;
  }
  up->vdim_gk = vdim_gk;
  up->vdim_vl = vdim_vl; 
    
  // allocate fields for prim mom calculation
  up->u_sq_temp = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->conf_rng->volume);
  up->m2_temp = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->conf_rng->volume);
  up->vth_sq_neut = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->conf_rng->volume);
  up->vth_sq_elc = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->conf_rng->volume);

  gkyl_dg_bin_op_mem *mem;
  up->mem = gkyl_dg_bin_op_mem_new(up->vth_sq_elc->size, up->basis->num_basis);

  up->react_rate = CK(ser_iz_react_rate_kernels, cdim, poly_order);
  assert(up->react_rate);
  
  return up;
}

void gkyl_dg_iz_react_rate(const struct gkyl_dg_iz *iz,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_neut, struct gkyl_array *vth_sq_iz,
  struct gkyl_array *cflrate, struct gkyl_array *coef_iz)
{
  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int pidx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<iz->conf_rng->ndim; ++d) rem_dir[d] = 1;
  
  // Calculate vtsq from moms using bin ops
  // neutral vtsq
  for (int d=1; d<iz->vdim_vl+1;) {
    gkyl_dg_mul_op(*iz->basis, 0, iz->u_sq_temp, d, moms_neut, d, moms_neut);
    gkyl_array_accumulate(iz->m2_temp, -1.0, iz->u_sq_temp); // -u.u
  }
  gkyl_dg_mul_op(*iz->basis, 0, iz->m2_temp, 0, iz->m2_temp, 0, moms_neut); // -u.u*m0
  gkyl_array_accumulate_offset(iz->m2_temp, 1.0, moms_neut, iz->vdim_vl+1); // m2-u.u*m0
  gkyl_dg_div_op(iz->mem, *iz->basis, 0, iz->vth_sq_neut, 0, iz->m2_temp, 0, moms_neut); // (m2-u.u*m0)/m0
  gkyl_array_scale(iz->vth_sq_neut, 1/iz->vdim_vl); // (m2-u.u*m0)/(vdim*m0)

  // elc vtsq
  gkyl_dg_mul_op(*iz->basis, 0, iz->m2_temp, 1, moms_elc, 1, moms_elc); // upar*upar
  gkyl_dg_mul_op(*iz->basis, 0, iz->m2_temp, 0, iz->m2_temp, 0, moms_elc); // uparSq*m0
  gkyl_array_accumulate_offset(iz->m2_temp, 1.0, moms_elc, 2); // uparSq*m0 - m2
  gkyl_dg_div_op(iz->mem, *iz->basis, 0, iz->vth_sq_elc, 0, iz->m2_temp, 0, moms_elc); // (uparSq*m0 - m2)/m0
  gkyl_array_scale(iz->vth_sq_elc, -1/iz->vdim_vl); // (m2 - uparSq*m0) / (vdim*m0)
    
  // Calculate vt_sq_iz
  gkyl_array_copy_range(vth_sq_iz, iz->vth_sq_elc, *iz->conf_rng); 
  gkyl_array_scale(vth_sq_iz, 1/2.0);
  gkyl_array_shiftc0(vth_sq_iz, -iz->E*iz->elem_charge/(3*iz->mass_elc)); 
  
  gkyl_range_iter_init(&conf_iter, iz->conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(iz->conf_rng, conf_iter.idx);
   
    // Reference proj_maxwellian_on_basis
    const double *m0_neut_d = gkyl_array_cfetch(moms_neut, loc);
    const double *vth_sq_neut_d = gkyl_array_cfetch(iz->vth_sq_neut, loc);
    const double *vth_sq_elc_d = gkyl_array_cfetch(iz->vth_sq_elc, loc);
    double *coef_iz_d = gkyl_array_fetch(coef_iz, loc);

    double cflr = iz->react_rate(iz->elem_charge, iz->mass_elc,
      iz->E, iz->A, iz->K, iz->P, iz->X,
      m0_neut_d, vth_sq_neut_d, vth_sq_elc_d,
      coef_iz_d);

    gkyl_range_deflate(&vel_rng, iz->phase_rng, rem_dir, conf_iter.idx);
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
