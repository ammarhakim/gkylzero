/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_dg_calc_em_vars_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_calc_em_vars_bvar_cu_kernel(struct gkyl_basis cbasis,
  struct gkyl_range range, 
  const struct gkyl_array* em, struct gkyl_array* bvar)
{
  int cdim = cbasis.ndim;
  int poly_order = cbasis.poly_order;

  em_t em_bvar = choose_ser_em_bvar_kern(cdim, poly_order);

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

    const double *em_d = (const double*) gkyl_array_cfetch(em, start);

    double *bvar_d = (double*) gkyl_array_fetch(bvar, start);

    em_bvar(em_d, bvar_d);
  }
}

// Host-side wrapper for magnetic field unit vector calculation
void
gkyl_calc_em_vars_bvar_cu(const struct gkyl_basis* cbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* em, struct gkyl_array* bvar)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_calc_em_vars_bvar_cu_kernel<<<nblocks, nthreads>>>(*cbasis, *range, em->on_dev, bvar->on_dev);
}

__global__ void
gkyl_calc_em_vars_ExB_cu_kernel(struct gkyl_basis cbasis,
  struct gkyl_range range, 
  const struct gkyl_array* em, struct gkyl_array* ExB)
{
  int cdim = cbasis.ndim;
  int poly_order = cbasis.poly_order;

  em_t em_ExB = choose_ser_em_ExB_kern(cdim, poly_order);

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

    const double *em_d = (const double*) gkyl_array_cfetch(em, start);

    double *ExB_d = (double*) gkyl_array_fetch(ExB, start);

    em_ExB(em_d, ExB_d);
  }
}

// Host-side wrapper for E x B velocity calculation
void
gkyl_calc_em_vars_ExB_cu(const struct gkyl_basis* cbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* em, struct gkyl_array* ExB)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_calc_em_vars_ExB_cu_kernel<<<nblocks, nthreads>>>(*cbasis, *range, em->on_dev, ExB->on_dev);
}

__global__ void
gkyl_calc_em_vars_pkpm_kappa_inv_b_cu_kernel(struct gkyl_basis cbasis,
  struct gkyl_range range, 
  const struct gkyl_array* em, struct gkyl_array* kappa_inv_b)
{
  int cdim = cbasis.ndim;
  int poly_order = cbasis.poly_order;

  em_pkpm_kappa_inv_b_t em_pkpm_kappa_inv_b = choose_ser_em_pkpm_kappa_inv_b_kern(cdim, poly_order);

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

    const double *bvar_d = (const double*) gkyl_array_cfetch(bvar, start);
    const double *ExB_d = (const double*) gkyl_array_cfetch(ExB, start);

    double *kappa_inv_b_d = (double*) gkyl_array_fetch(kappa_inv_b, start);

    em_pkpm_kappa_inv_b(bvar_d, ExB_d, kappa_inv_b_d);
  }
}

// Host-side wrapper for b/kappa (kappa Lorentz boost factor for E x B velocity) calculation
void
gkyl_calc_em_vars_kappa_inv_b_cu(const struct gkyl_basis* cbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* bvar, const struct gkyl_array* ExB, 
  struct gkyl_array* kappa_inv_b)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_calc_em_vars_pkpm_kappa_inv_b_cu_kernel<<<nblocks, nthreads>>>(*cbasis, *range, 
    bvar->on_dev, ExB->on_dev, kappa_inv_b->on_dev);
}
