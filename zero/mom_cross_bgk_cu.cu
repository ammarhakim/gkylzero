/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_cross_bgk.h>
#include <gkyl_mom_cross_bgk_priv.h>
}

__global__ void
gkyl_mom_cross_bgk_gyrokinetic_advance_cu_kernel(gkyl_mom_cross_bgk_gyrokinetic *up,
  struct gkyl_range conf_range, double beta,
  double m_self, const struct gkyl_array *moms_self,
  double m_other, const struct gkyl_array *moms_other,
  const struct gkyl_array *nu_sr, const struct gkyl_array *nu_rs, 
  struct gkyl_array *moms_cross)
{ 
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc_conf = gkyl_range_idx(&conf_range, idx);

    const double *in1_d = (const double*) gkyl_array_cfetch(moms_self, loc_conf);
    const double *in2_d = (const double*) gkyl_array_cfetch(moms_other, loc_conf);
    const double *nu_sr_d = (const double*) gkyl_array_cfetch(nu_sr, loc_conf);
    const double *nu_rs_d = (const double*) gkyl_array_cfetch(nu_rs, loc_conf);
    double *out_d = (double*) gkyl_array_fetch(moms_cross, loc_conf);

    up->mom_cross_calc(beta, m_self, in1_d, m_other, in2_d, nu_sr_d, nu_rs_d, out_d);
  }
}

// Host-side wrapper for cross BGK moments
void
gkyl_mom_cross_bgk_gyrokinetic_advance_cu(gkyl_mom_cross_bgk_gyrokinetic *up,
  const struct gkyl_range *conf_range, double beta,
  double m_self, const struct gkyl_array *moms_self,
  double m_other, const struct gkyl_array *moms_other,
  const struct gkyl_array *nu_sr, const struct gkyl_array *nu_rs, 
  struct gkyl_array *moms_cross)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_mom_cross_bgk_gyrokinetic_advance_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, 
    *conf_range, beta, m_self, moms_self->on_dev, m_other, moms_other->on_dev, 
    nu_sr->on_dev, nu_rs->on_dev, moms_cross->on_dev);
}

__global__
static void
set_mom_cross_bgk_cu_ptrs(struct gkyl_mom_cross_bgk_gyrokinetic *up, int cdim, int vdim, int poly_order)
{
  up->mom_cross_calc = choose_mom_cross_bgk_gyrokinetic_kern(cdim, vdim, poly_order); 
}

gkyl_mom_cross_bgk_gyrokinetic* 
gkyl_mom_cross_bgk_gyrokinetic_cu_dev_new(const struct gkyl_basis *phase_basis, const struct gkyl_basis *conf_basis)
{
  struct gkyl_mom_cross_bgk_gyrokinetic *up = (struct gkyl_mom_cross_bgk_gyrokinetic*) gkyl_malloc(sizeof(*up));
  up->use_gpu = true;
  
  // copy struct to device
  struct gkyl_mom_cross_bgk_gyrokinetic *up_cu = (struct gkyl_mom_cross_bgk_gyrokinetic*) gkyl_cu_malloc(sizeof(*up_cu));
  
  set_mom_cross_bgk_cu_ptrs<<<1,1>>>(up_cu, conf_basis->ndim,
    phase_basis->ndim-conf_basis->ndim, phase_basis->poly_order);

  up->on_dev = up_cu;  
  return up;   
}