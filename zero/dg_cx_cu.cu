/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_prim_vars_vlasov.h>
#include <gkyl_dg_prim_vars_gyrokinetic.h>
#include <gkyl_dg_prim_vars_transform.h>
#include <gkyl_dg_prim_vars_type.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_cx.h>
#include <gkyl_dg_cx_priv.h>
#include <gkyl_util.h>
#include <gkyl_const.h>
}

// CUDA kernel to set device pointers to kernels.
__global__ static void
gkyl_dg_cx_set_cu_ker_ptrs(struct gkyl_dg_cx_kernels *kernels,
  struct gkyl_basis pbasis_vl, int tblidx)
{
  enum gkyl_basis_type b_type = pbasis_vl.b_type;
  int poly_order = pbasis_vl.poly_order;

  switch (b_type) {
    case GKYL_BASIS_MODAL_HYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->react_rate = ser_cx_react_rate_kernels[tblidx].kernels[poly_order];
      break;
    default:
      assert(false);
      break;
  }
};

void
dg_cx_choose_kernel_cu(struct gkyl_dg_cx_kernels *kernels,
  struct gkyl_basis pbasis_vl, struct gkyl_basis cbasis)
{
  int pdim = pbasis_vl.ndim;
  int cdim = cbasis.ndim;
  int vdim = pdim - cdim;

  assert(cv_index[cdim].vdim[vdim] != -1);
  gkyl_dg_cx_set_cu_ker_ptrs<<<1,1>>>(kernels, pbasis_vl, cv_index[cdim].vdim[vdim]);
}

__global__ static void
gkyl_cx_react_rate_cu_ker(struct gkyl_dg_cx_kernels *kernels, const struct gkyl_range conf_rng, 
  const struct gkyl_dg_prim_vars_type *calc_prim_vars_ion, const struct gkyl_dg_prim_vars_type *calc_prim_vars_neut,
  const struct gkyl_dg_prim_vars_type *calc_prim_vars_neut_gk, const struct gkyl_array *moms_ion,
  const struct gkyl_array *moms_neut, struct gkyl_array *prim_vars_ion, struct gkyl_array *prim_vars_neut, 
  struct gkyl_array *prim_vars_neut_gk, double vtsq_min_ion, double vtsq_min_neut, struct gkyl_array *coef_cx,
  double a, double b)
{
  int cidx[GKYL_MAX_CDIM];
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < conf_rng.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&conf_rng, tid, cidx);
    long loc = gkyl_range_idx(&conf_rng, cidx);

    const double *moms_ion_d = (const double*) gkyl_array_cfetch(moms_ion, loc);
    const double *moms_neut_d = (const double*) gkyl_array_cfetch(moms_neut, loc);
    const double *m0_neut_d = &moms_neut_d[0];

    double *prim_vars_ion_d = (double*) gkyl_array_fetch(prim_vars_ion, loc);
    double *prim_vars_neut_d = (double*) gkyl_array_fetch(prim_vars_neut, loc);
    double *prim_vars_neut_gk_d = (double*) gkyl_array_fetch(prim_vars_neut_gk, loc);
    double *coef_cx_d = (double*) gkyl_array_fetch(coef_cx, loc);
    
    calc_prim_vars_ion->kernel(calc_prim_vars_ion, cidx,
				   moms_ion_d, prim_vars_ion_d);
    calc_prim_vars_neut->kernel(calc_prim_vars_neut, cidx,
				    moms_neut_d, prim_vars_neut_d);
    calc_prim_vars_neut_gk->kernel(calc_prim_vars_neut_gk, cidx,
				    moms_neut_d, prim_vars_neut_gk_d);

    // call the cx kernel
    double cflr = kernels->react_rate(a, b, vtsq_min_ion, vtsq_min_neut, m0_neut_d,
      prim_vars_ion_d, prim_vars_neut_d, coef_cx_d);
  }
}

void gkyl_dg_cx_coll_cu(const struct gkyl_dg_cx *up, const double vtsq_min_ion,
  const double vtsq_min_neut, const struct gkyl_array *moms_ion,
  const struct gkyl_array *moms_neut, const struct gkyl_array *b_i,
  struct gkyl_array *prim_vars_ion, struct gkyl_array *prim_vars_neut,
  struct gkyl_array *prim_vars_neut_gk, struct gkyl_array *coef_cx, struct gkyl_array *cflrate)
{

  gkyl_dg_prim_vars_transform_set_auxfields(up->calc_prim_vars_ion, 
    (struct gkyl_dg_prim_vars_auxfields) {.b_i = b_i});
  gkyl_dg_prim_vars_transform_set_auxfields(up->calc_prim_vars_neut_gk, 
    (struct gkyl_dg_prim_vars_auxfields) {.b_i = b_i});
  
  gkyl_cx_react_rate_cu_ker<<<up->conf_rng->nblocks, up->conf_rng->nthreads>>>(up->kernels, *up->conf_rng,
    up->calc_prim_vars_ion->on_dev, up->calc_prim_vars_neut->on_dev, up->calc_prim_vars_neut_gk->on_dev,
    moms_ion->on_dev, moms_neut->on_dev, prim_vars_ion->on_dev, prim_vars_neut->on_dev,
    prim_vars_neut_gk->on_dev,vtsq_min_ion, vtsq_min_neut, coef_cx->on_dev, up->a, up->b);
}
