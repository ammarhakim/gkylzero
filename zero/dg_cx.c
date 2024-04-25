#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_prim_vars_vlasov.h>
#include <gkyl_dg_prim_vars_transform.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_cx.h>
#include <gkyl_dg_cx_priv.h>
#include <gkyl_util.h>

gkyl_dg_cx*
gkyl_dg_cx_new(struct gkyl_dg_cx_inp *inp, bool use_gpu)
{
  gkyl_dg_cx *up = gkyl_malloc(sizeof(struct gkyl_dg_cx));

  up->cbasis = inp->cbasis;
  up->pbasis = inp->pbasis;
  up->conf_rng = inp->conf_rng;
  up->conf_rng_ext = inp->conf_rng_ext;
  up->phase_rng = inp->phase_rng;
  up->grid = inp->grid;
  up->mass_ion = inp->mass_ion;
  up->mass_neut = inp->mass_neut;
  up->type_ion = inp->type_ion;
  up->vt_sq_ion_min = inp->vt_sq_ion_min;
  up->vt_sq_neut_min = inp->vt_sq_neut_min;

  int cdim = up->cbasis->ndim;
  int pdim = up->pbasis->ndim;
  int poly_order = up->cbasis->poly_order;
  up->cdim = cdim;
  up->use_gpu = use_gpu;
  up->vdim_gk = 2; 
  up->vdim_vl = 3; // assume 2x3v or 3x3v (will change for true axisym Vlasov)  

  if (up->type_ion == GKYL_ION_H) {
    up->a = 1.12e-18;
    up->b = 7.15e-20;
  }
  else if (up->type_ion == GKYL_ION_D) {
    up->a = 1.09e-18;
    up->b = 7.15e-20;
  }
  else if (up->type_ion == GKYL_ION_HE) {
    up->a = 6.484e-19;
    up->b = 4.350e-20;
  } 
  else if (up->type_ion == GKYL_ION_NE) {
    up->a = 7.95e-19;
    up->b = 5.65e-20;
  }

  up->calc_prim_vars_ion = gkyl_dg_prim_vars_transform_new(up->cbasis, up->pbasis, up->conf_rng, "prim_vlasov", use_gpu);
  up->calc_prim_vars_neut = gkyl_dg_prim_vars_vlasov_new(up->cbasis, up->pbasis, "prim", use_gpu);
  up->calc_prim_vars_neut_gk = gkyl_dg_prim_vars_transform_new(up->cbasis, up->pbasis, up->conf_rng, "prim_gk", use_gpu);
  
  up->on_dev = up; // CPU eqn obj points to itself
  
  up->react_rate = ser_cx_react_rate_kernels[cv_index[cdim].vdim[up->vdim_vl]].kernels[poly_order];
  assert(up->react_rate);
  
  return up;
}

void gkyl_dg_cx_coll(const struct gkyl_dg_cx *up, const double vtsq_min_ion,
  const double vtsq_min_neut, const struct gkyl_array *moms_ion,
  const struct gkyl_array *moms_neut, const struct gkyl_array *b_i,
  struct gkyl_array *prim_vars_ion, struct gkyl_array *prim_vars_neut,
  struct gkyl_array *prim_vars_neut_gk, struct gkyl_array *coef_cx, struct gkyl_array *cflrate)
{
  #ifdef GKYL_HAVE_CUDA
  if(gkyl_array_is_cu_dev(coef_recomb)) {
    return gkyl_dg_cx_coll_cu(up, vtsq_min_ion, vtsq_min_neut, moms_ion, moms_neut, b_i,
      prim_vars_ion, prim_vars_neut, prim_vars_neut_gk, coef_cx, cflrate);
  }
#endif
  gkyl_dg_prim_vars_transform_set_auxfields(up->calc_prim_vars_ion, 
      (struct gkyl_dg_prim_vars_auxfields) {.b_i = b_i});
  gkyl_dg_prim_vars_transform_set_auxfields(up->calc_prim_vars_neut_gk, 
      (struct gkyl_dg_prim_vars_auxfields) {.b_i = b_i});
  
  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  gkyl_range_iter_init(&conf_iter, up->conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(up->conf_rng, conf_iter.idx);

    const double *moms_ion_d = gkyl_array_cfetch(moms_ion, loc);
    const double *moms_neut_d = gkyl_array_cfetch(moms_neut, loc);
    const double *m0_neut_d = &moms_neut_d[0];

    double *prim_vars_ion_d = gkyl_array_fetch(prim_vars_ion, loc);
    double *prim_vars_neut_d = gkyl_array_fetch(prim_vars_neut, loc);
    double *prim_vars_neut_gk_d = gkyl_array_fetch(prim_vars_neut_gk, loc);
    double *coef_cx_d = gkyl_array_fetch(coef_cx, loc);

    up->calc_prim_vars_ion->kernel(up->calc_prim_vars_ion, conf_iter.idx,
				   moms_ion_d, prim_vars_ion_d);
    up->calc_prim_vars_neut->kernel(up->calc_prim_vars_neut, conf_iter.idx,
				    moms_neut_d, prim_vars_neut_d);
    up->calc_prim_vars_neut_gk->kernel(up->calc_prim_vars_neut_gk, conf_iter.idx,
				    moms_neut_d, prim_vars_neut_gk_d);
    
    double cflr = up->react_rate(up->a, up->b, vtsq_min_ion, vtsq_min_neut,
      m0_neut_d, prim_vars_ion_d, prim_vars_neut_d, coef_cx_d);
    
    /* gkyl_range_deflate(&vel_rng, up->phase_rng, rem_dir, conf_iter.idx); */
    /* gkyl_range_iter_no_split_init(&vel_iter, &vel_rng); */
    /* // cfl associated with reaction is a *phase space* cfl */
    /* // Need to loop over velocity space for each configuration space cell */
    /* // to get total cfl rate in each phase space cell */
    /* while (gkyl_range_iter_next(&vel_iter)) { */
    /*   long cfl_idx = gkyl_range_idx(&vel_rng, vel_iter.idx); */
    /*   double *cflrate_d = gkyl_array_fetch(cflrate, cfl_idx); */
    /*   cflrate_d[0] += cflr; // frequencies are additive */
    /* } */
  }
}

void
gkyl_dg_cx_release(gkyl_dg_cx* cx)
{
  free(cx);
}
