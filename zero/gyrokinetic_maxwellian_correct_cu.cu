/* -*- c++ -*- */

extern "C" {
#include <gkyl_gyrokinetic_maxwellian_correct.h>
#include <gkyl_gyrokinetic_maxwellian_correct_priv.h>
#include <gkyl_range.h>
}

__global__ static void
gkyl_gyrokinetic_maxwellian_correct_all_moments_abs_diff_cu_ker(struct gkyl_range conf_range, 
  int nc, const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter, 
  struct gkyl_array *abs_diff_moms)
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
    long loc = gkyl_range_idx(&conf_range, idx);

    const double *moms_local = (const double*) gkyl_array_cfetch(moms_iter, loc);
    const double *moms_target_local = (const double*) gkyl_array_cfetch(moms_target, loc);
    double *abs_diff_moms_local = (double*) gkyl_array_fetch(abs_diff_moms, loc);
    // Compute the absolute value of the difference of cell averages 
    // Note: max error found by follow-up thread-safe reduction operation
    // Also: for density and temperature, this error is a relative error compared to the target moment value
    // so that we can converge to the correct target moments in SI units and minimize finite precision issues.
    abs_diff_moms_local[0] = fabs(moms_local[0*nc] - moms_target_local[0*nc])/moms_target_local[0*nc];
    abs_diff_moms_local[2] = fabs(moms_local[2*nc] - moms_target_local[2*nc])/moms_target_local[2*nc];
    // However, u_par may be ~ 0 and if it is, we need to use absolute error. We can converge safely using
    // absolute error if u_par ~ O(1). Otherwise, we use relative error for u_par. 
    if (fabs(moms_target_local[1*nc]) < 1.0) {
      abs_diff_moms_local[1] = fabs(moms_local[1*nc] - moms_target_local[1*nc]);
    }
    else {
      abs_diff_moms_local[1] = fabs(moms_local[1*nc] - moms_target_local[1*nc])/moms_target_local[1*nc];
    }
  }
}

void
gkyl_gyrokinetic_maxwellian_correct_all_moments_abs_diff_cu(const struct gkyl_range *conf_range, 
  int nc, const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter, 
  struct gkyl_array *moms_abs_diff)
{
  int nblocks = conf_range->nblocks, nthreads = conf_range->nthreads;
  gkyl_gyrokinetic_maxwellian_correct_all_moments_abs_diff_cu_ker<<<nblocks, nthreads>>>(*conf_range, 
    nc, moms_target->on_dev, moms_iter->on_dev, moms_abs_diff->on_dev);
}

__global__ static void
gkyl_gyrokinetic_bimaxwellian_correct_all_moments_abs_diff_cu_ker(struct gkyl_range conf_range, 
  int nc, const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter, 
  struct gkyl_array *abs_diff_moms)
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
    long loc = gkyl_range_idx(&conf_range, idx);

    const double *moms_local = (const double*) gkyl_array_cfetch(moms_iter, loc);
    const double *moms_target_local = (const double*) gkyl_array_cfetch(moms_target, loc);
    double *abs_diff_moms_local = (double*) gkyl_array_fetch(abs_diff_moms, loc);
    // Compute the absolute value of the difference of cell averages 
    // Note: max error found by follow-up thread-safe reduction operation
    // Also: for density, T_par/m, T_perp/m this error is a relative error compared to the target moment value
    // so that we can converge to the correct target moments in SI units and minimize finite precision issues.
    abs_diff_moms_local[0] = fabs(moms_local[0*nc] - moms_target_local[0*nc])/moms_target_local[0*nc];
    abs_diff_moms_local[2] = fabs(moms_local[2*nc] - moms_target_local[2*nc])/moms_target_local[2*nc];
    abs_diff_moms_local[3] = fabs(moms_local[3*nc] - moms_target_local[3*nc])/moms_target_local[3*nc];
    // However, u_par may be ~ 0 and if it is, we need to use absolute error. We can converge safely using
    // absolute error if u_par ~ O(1). Otherwise, we use relative error for u_par. 
    if (fabs(moms_target_local[1*nc]) < 1.0) {
      abs_diff_moms_local[1] = fabs(moms_local[1*nc] - moms_target_local[1*nc]);
    }
    else {
      abs_diff_moms_local[1] = fabs(moms_local[1*nc] - moms_target_local[1*nc])/moms_target_local[1*nc];
    }
  }
}

void
gkyl_gyrokinetic_bimaxwellian_correct_all_moments_abs_diff_cu(const struct gkyl_range *conf_range, 
  int nc, const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter, 
  struct gkyl_array *moms_abs_diff)
{
  int nblocks = conf_range->nblocks, nthreads = conf_range->nthreads;
  gkyl_gyrokinetic_bimaxwellian_correct_all_moments_abs_diff_cu_ker<<<nblocks, nthreads>>>(*conf_range, 
    nc, moms_target->on_dev, moms_iter->on_dev, moms_abs_diff->on_dev);
}
