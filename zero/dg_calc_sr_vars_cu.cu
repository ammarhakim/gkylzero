/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_calc_sr_vars_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_calc_sr_vars_GammaV2_cu_kernel(struct gkyl_basis cbasis, struct gkyl_basis pbasis, 
  struct gkyl_range range, const struct gkyl_array* V, struct gkyl_array* GammaV2)
{
  int cdim = cbasis.ndim, pdim = pbasis.ndim, vdim = pdim-cdim;
  int poly_order = cbasis.poly_order;
  sr_t sr_GammaV2 = choose_ser_sr_GammaV2_kern(cdim, vdim, poly_order);

  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&range, idx);

    const double *V_d = (const double*) gkyl_array_cfetch(V, start);

    double *GammaV2_d = (double*) gkyl_array_fetch(GammaV2, start);

    sr_GammaV2(V_d, GammaV2_d);
  }
}

// Host-side wrapper for GammaV^2 = 1/(1 - V^2/c^2) operation
void
gkyl_calc_sr_vars_GammaV2_cu(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV2)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_calc_sr_vars_GammaV2_cu_kernel<<<nblocks, nthreads>>>(*cbasis, *pbasis, *range, 
    V->on_dev, GammaV2->on_dev);
}

__global__ void
gkyl_calc_sr_vars_GammaV_inv_cu_kernel(struct gkyl_basis cbasis, struct gkyl_basis pbasis, 
  struct gkyl_range range, const struct gkyl_array* V, struct gkyl_array* GammaV_inv)
{
  int cdim = cbasis.ndim, pdim = pbasis.ndim, vdim = pdim-cdim;
  int poly_order = cbasis.poly_order;
  sr_t sr_GammaV_inv = choose_ser_sr_GammaV_inv_kern(cdim, vdim, poly_order);

  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&range, idx);

    const double *V_d = (const double*) gkyl_array_cfetch(V, start);

    double *GammaV_inv_d = (double*) gkyl_array_fetch(GammaV_inv, start);

    sr_GammaV_inv(V_d, GammaV_inv_d);
  }
}

// Host-side wrapper for GammaV^{-1} = sqrt(1 - V^2/c^2) operation
void
gkyl_calc_sr_vars_GammaV_inv_cu(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV_inv)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_calc_sr_vars_GammaV_inv_cu_kernel<<<nblocks, nthreads>>>(*cbasis, *pbasis, *range, 
    V->on_dev, GammaV_inv->on_dev);
}
