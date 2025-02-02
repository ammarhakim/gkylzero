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

gkyl_dg_cx*
gkyl_dg_cx_new(struct gkyl_dg_cx_inp *inp, bool use_gpu)
{
  gkyl_dg_cx *up = gkyl_malloc(sizeof(struct gkyl_dg_cx));

  up->cbasis = inp->cbasis;
  up->pbasis_gk = inp->pbasis_gk;
  up->pbasis_vl = inp->pbasis_vl;
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
  int poly_order = up->cbasis->poly_order;
  up->cdim = cdim;
  up->use_gpu = use_gpu;
  up->vdim_gk = up->pbasis_gk->ndim - cdim; 
  up->vdim_vl = up->pbasis_vl->ndim - cdim; 

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
  
  //up->on_dev = up; // CPU eqn obj points to itself
  
  //up->react_rate = ser_cx_react_rate_kernels[cv_index[cdim].vdim[up->vdim_vl]].kernels[poly_order];

  if (!use_gpu) {
    up->kernels = gkyl_malloc(sizeof(struct gkyl_dg_cx_kernels));
  }
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->kernels = gkyl_cu_malloc(sizeof(struct gkyl_dg_cx_kernels));
  }
#endif

  dg_cx_choose_kernel(up->kernels, *up->pbasis_vl, *up->cbasis, use_gpu);
  //assert(up->react_rate);
  
  return up;
}

void gkyl_dg_cx_coll(const struct gkyl_dg_cx *up, 
  struct gkyl_array *prim_vars_ion, struct gkyl_array *prim_vars_neut,
  struct gkyl_array *upar_b_i, struct gkyl_array *coef_cx, struct gkyl_array *cflrate)
{
#ifdef GKYL_HAVE_CUDA
  if(gkyl_array_is_cu_dev(coef_cx)) {
    return gkyl_dg_cx_coll_cu(up, prim_vars_ion, prim_vars_neut, 
      upar_b_i, coef_cx, cflrate);
  }
#endif
  
  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  gkyl_range_iter_init(&conf_iter, up->conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long loc = gkyl_range_idx(up->conf_rng, conf_iter.idx);

    const double *prim_vars_ion_d = gkyl_array_cfetch(prim_vars_ion, loc);
    const double *prim_vars_neut_d = gkyl_array_cfetch(prim_vars_neut, loc);
    const double *upar_b_i_d = gkyl_array_cfetch(upar_b_i, loc);

    double *coef_cx_d = gkyl_array_fetch(coef_cx, loc);
    
    double cflr = up->kernels->react_rate(up->a, up->b, up->vt_sq_ion_min, up->vt_sq_neut_min,
      prim_vars_ion_d, prim_vars_neut_d, upar_b_i_d, coef_cx_d);
    
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
