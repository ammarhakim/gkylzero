#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_cx.h>
#include <gkyl_dg_cx_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

gkyl_dg_cx*
gkyl_dg_cx_new(const struct gkyl_rect_grid *grid,
  struct gkyl_basis *cbasis, struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng,
  double mass_ion, enum gkyl_dg_cx_type type_ion, 
  bool is_gk, bool use_gpu)
{
  gkyl_dg_cx *up = gkyl_malloc(sizeof(gkyl_dg_cx));
  
  int num_basis = pbasis->num_basis;
  int pdim = pbasis->ndim;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  up->basis = cbasis;
  up->grid = *grid; 
  up->cdim = cdim; 
  up->poly_order = poly_order;
  up->conf_rng = conf_rng;
  up->phase_rng = phase_rng;

  // Establish vdim for gk and vlasov species
  // will be either 1 & 1 or 2 & 3
  // for gk and vlasov, respectively
  int vdim_vl; 
  int vdim = pdim - up->cdim;
  if (is_gk) {
    if (vdim == 1) {
      vdim_vl = vdim;
    }
    else {
      vdim_vl = vdim+1;
    }
  }
  else  {
    vdim_vl = vdim;
  }
  up->vdim_vl = vdim_vl;

  /* this assumes ion mass = neut mass */
  up->mass_ion = mass_ion;

  if (type_ion == GKYL_H) {
    up->a = 1.12e-18;
    up->b = 7.15e-20;
  }
  else if (type_ion == GKYL_D) {
    up->a = 1.09e-18;
    up->b = 7.15e-20;
  }
  else if (type_ion == GKYL_NE) {
    up->a = 7.95e-19;
    up->b = 5.65e-20;
  }
  else if (type_ion == GKYL_HE) {
    up->a = 6.484e-19;
    up->b = 4.350e-20;
  }  

  // allocate fields for prim mom calculation
  up->m2_temp = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->conf_rng->volume);
  up->udrift_neut = gkyl_array_new(GKYL_DOUBLE, vdim_vl*up->basis->num_basis, up->conf_rng->volume);
  up->udrift_ion = gkyl_array_new(GKYL_DOUBLE, vdim_vl*up->basis->num_basis, up->conf_rng->volume);
  up->vth_sq_neut = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->conf_rng->volume);
  up->vth_sq_ion = gkyl_array_new(GKYL_DOUBLE, up->basis->num_basis, up->conf_rng->volume);

  gkyl_dg_bin_op_mem *mem;
  up->mem = gkyl_dg_bin_op_mem_new(up->vth_sq_neut->size, up->basis->num_basis);
  
  up->react_rate = CK(ser_cx_react_rate_kernels, cdim, poly_order);
  assert(up->react_rate);
  
  return up;
}

void gkyl_dg_cx_react_rate(const struct gkyl_dg_cx *cx, const struct gkyl_array *moms_ion,
  const struct gkyl_array *moms_neut, const struct gkyl_array *b_hat,
  struct gkyl_array *cflrate, struct gkyl_array *coef_cx)
{
  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int rem_dir[GKYL_MAX_DIM] = { 0 };
  for (int d=0; d<cx->conf_rng->ndim; ++d) rem_dir[d] = 1;

  // Calculate neutral udrift and vt_sq
  for (int d=0; d<cx->vdim_vl; d+=1) {
    gkyl_dg_div_op_range(cx->mem, *cx->basis, d, cx->udrift_neut, d+1, moms_neut, 0, moms_neut, cx->conf_rng);
  }
  
  gkyl_dg_dot_product_op_range(*cx->basis, cx->m2_temp, cx->udrift_neut, cx->udrift_neut, cx->conf_rng); // u.u
  gkyl_dg_mul_op_range(*cx->basis, 0, cx->m2_temp, 0, cx->m2_temp, 0, moms_neut, cx->conf_rng); //u.u*m0
  gkyl_array_accumulate_offset_range(cx->m2_temp, -1.0, moms_neut, (cx->vdim_vl+1)*cx->basis->num_basis, *cx->conf_rng); // u.u*m0 - m2
  gkyl_dg_div_op_range(cx->mem, *cx->basis, 0, cx->vth_sq_neut, 0, cx->m2_temp, 0, moms_neut, cx->conf_rng); // (u.u*m0 - m2)/m0
  gkyl_array_scale_range(cx->vth_sq_neut, -1/cx->vdim_vl, *cx->conf_rng); // (m2-u.u*m0)/(vdim*m0)

  // Calculate ion udrift and vt_sq
  for (int d=0; d<cx->vdim_vl; d+=1) {
    gkyl_dg_div_op_range(cx->mem, *cx->basis, d, cx->udrift_ion, 1, moms_ion, 0, moms_ion, cx->conf_rng);
    gkyl_dg_mul_op_range(*cx->basis, d, cx->udrift_ion, d, cx->udrift_ion, d, b_hat, cx->conf_rng);
  }

  gkyl_dg_dot_product_op_range(*cx->basis, cx->m2_temp, cx->udrift_ion, cx->udrift_ion, cx->conf_rng);
  gkyl_dg_mul_op_range(*cx->basis, 0, cx->m2_temp, 0, cx->m2_temp, 0, moms_ion, cx->conf_rng);
  gkyl_array_accumulate_offset_range(cx->m2_temp, -1.0, moms_ion, (cx->vdim_vl+1)*cx->basis->num_basis, *cx->conf_rng); // u.u*m0 - m2
  gkyl_dg_div_op_range(cx->mem, *cx->basis, 0, cx->vth_sq_ion, 0, cx->m2_temp, 0, moms_ion, cx->conf_rng); // (u.u*m0 - m2)/m0
  gkyl_array_scale_range(cx->vth_sq_ion, -1/cx->vdim_vl, *cx->conf_rng); // (m2-u.u*m0)/(vdim*m0)
  
  gkyl_range_iter_init(&conf_iter, cx->conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(cx->conf_rng, conf_iter.idx);

    const double *m0_neut_d = gkyl_array_cfetch(moms_neut, loc);
    const double *u_neut_d = gkyl_array_cfetch(cx->udrift_neut, loc);
    const double *u_ion_d = gkyl_array_cfetch(cx->udrift_ion, loc);
    const double *vth_sq_neut_d = gkyl_array_cfetch(cx->vth_sq_neut, loc);
    const double *vth_sq_ion_d = gkyl_array_cfetch(cx->vth_sq_ion, loc);

    double *coef_cx_d = gkyl_array_fetch(coef_cx, loc);

    // Calculate vt_sq min for ion, neut (use same for now to test 1x1v)
    double vth_sq_ion_min;
    double vth_sq_neut_min;
    double TempMin = 0.0; 
    for (int d=0; d<cx->vdim_vl; d++) {
      TempMin = TempMin + (1./3.)*(cx->mass_ion/6.)*cx->grid.dx[cx->cdim+d];
    }
    vth_sq_ion_min = TempMin/cx->mass_ion;
    vth_sq_neut_min = TempMin/cx->mass_ion;
    
    double cflr = cx->react_rate(cx->a, cx->b,
      m0_neut_d, u_ion_d, u_neut_d, vth_sq_ion_d,
      vth_sq_ion_min, vth_sq_neut_d, vth_sq_neut_min,
      coef_cx_d);

    gkyl_range_deflate(&vel_rng, cx->phase_rng, rem_dir, conf_iter.idx);
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
gkyl_dg_cx_release(gkyl_dg_cx* cx)
{
  free(cx);
}
