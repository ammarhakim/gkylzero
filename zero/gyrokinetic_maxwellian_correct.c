#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_array_reduce.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_gyrokinetic_maxwellian_correct.h>
#include <gkyl_gyrokinetic_maxwellian_correct_priv.h>
#include <gkyl_gyrokinetic_maxwellian_moments.h>
#include <gkyl_proj_maxwellian_on_basis.h>

#include <assert.h>

struct gkyl_gyrokinetic_maxwellian_correct*
gkyl_gyrokinetic_maxwellian_correct_inew(const struct gkyl_gyrokinetic_maxwellian_correct_inp *inp)
{
  gkyl_gyrokinetic_maxwellian_correct *up = gkyl_malloc(sizeof(*up));
  up->eps = inp->eps;
  up->max_iter = inp->max_iter;
  up->use_gpu = inp->use_gpu;

  up->conf_basis = *inp->conf_basis;
  up->phase_basis = *inp->phase_basis;
  up->num_conf_basis = up->conf_basis.num_basis;
  up->gk_geom = gkyl_gk_geometry_acquire(inp->gk_geom);
  up->divide_jacobgeo = inp->divide_jacobgeo;
  up->use_last_converged = inp->use_last_converged;
  up->mass = inp->mass;
  
  long conf_range_ncells = inp->conf_range->volume;
  long conf_range_ext_ncells = inp->conf_range_ext->volume;

  // Individual moment memory: the iteration of the moments, the differences (d) and differences of differences (dd)
  // We are correcting three moments: n, u_par, T/m in the gyrokinetic system
  if (up->use_gpu) {
    up->moms_iter = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*inp->conf_basis->num_basis, conf_range_ext_ncells);
    up->d_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*inp->conf_basis->num_basis, conf_range_ext_ncells);
    up->dd_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3*inp->conf_basis->num_basis, conf_range_ext_ncells);
    // Two additional GPU-specific allocations for iterating over the grid to find the absolute value of 
    // the difference between the target and iterative moments, and the GPU-side array for performing the
    // thread-safe reduction to find the maximum error on the grid.
    up->abs_diff_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, 3, conf_range_ext_ncells);
    up->error_cu = gkyl_cu_malloc(sizeof(double[3]));
    // Arrays for correcting the density 
    up->num_ratio = gkyl_array_cu_dev_new(GKYL_DOUBLE, inp->conf_basis->num_basis, conf_range_ext_ncells);
    up->mem = gkyl_dg_bin_op_mem_cu_dev_new(conf_range_ncells, up->num_conf_basis);
  }
  else {
    up->moms_iter = gkyl_array_new(GKYL_DOUBLE, 3*inp->conf_basis->num_basis, conf_range_ext_ncells);
    up->d_moms = gkyl_array_new(GKYL_DOUBLE, 3*inp->conf_basis->num_basis, conf_range_ext_ncells);
    up->dd_moms = gkyl_array_new(GKYL_DOUBLE, 3*inp->conf_basis->num_basis, conf_range_ext_ncells);
    // Arrays for correcting the density 
    up->num_ratio = gkyl_array_new(GKYL_DOUBLE, inp->conf_basis->num_basis, conf_range_ext_ncells);
    up->mem = gkyl_dg_bin_op_mem_new(conf_range_ncells, up->num_conf_basis);
  }

  // Moments structure 
  struct gkyl_gyrokinetic_maxwellian_moments_inp inp_mom = {
    .phase_grid = inp->phase_grid,
    .conf_basis = inp->conf_basis,
    .phase_basis = inp->phase_basis,
    .conf_range =  inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .vel_range = inp->vel_range,
    .mass = inp->mass, 
    .gk_geom = inp->gk_geom, 
    .divide_jacobgeo = inp->divide_jacobgeo, 
    .use_gpu = inp->use_gpu,
  };
  up->moments_up = gkyl_gyrokinetic_maxwellian_moments_inew( &inp_mom );

  // Create a projection updater for projecting the gyrokinetic Maxwellian
  up->proj_max_prim = gkyl_proj_maxwellian_on_basis_new(inp->phase_grid,
    inp->conf_basis, inp->phase_basis, inp->phase_basis->poly_order+1, inp->use_gpu);

  return up;
}

void
gkyl_gyrokinetic_maxwellian_correct_density_moment(gkyl_gyrokinetic_maxwellian_correct *up,
  struct gkyl_array *f_max, const struct gkyl_array *moms_target, 
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range)
{
  // Correct the density of the projected gyrokinetic Maxwellian distribution function through rescaling. 
  gkyl_gyrokinetic_maxwellian_density_moment_advance(up->moments_up, phase_range, conf_range, f_max, up->num_ratio);

  // compute number density ratio: num_ratio = n/n0
  // 0th component of moms_target is the target density
  gkyl_dg_div_op_range(up->mem, up->conf_basis, 0, up->num_ratio,
    0, moms_target, 0, up->num_ratio, conf_range);

  // rescale distribution function
  gkyl_dg_mul_conf_phase_op_range(&up->conf_basis, &up->phase_basis,
    f_max, up->num_ratio, f_max, conf_range, phase_range);  
}

struct gkyl_gyrokinetic_maxwellian_correct_status
gkyl_gyrokinetic_maxwellian_correct_all_moments(gkyl_gyrokinetic_maxwellian_correct *up,
  struct gkyl_array *f_max, const struct gkyl_array *moms_target, 
  const struct gkyl_range *phase_range, const struct gkyl_range *conf_range)
{
  int nc = up->num_conf_basis;
  double tol = up->eps;  // tolerance of the iterative scheme
  int max_iter = up->max_iter;

  int niter = 0;
  bool corr_status = true;
  int ispositive_f_max = true;

  // Set initial max error to start the iteration.
  up->error[0] = 1.0;
  up->error[1] = 1.0;
  up->error[2] = 1.0;
  // Copy the initial max error to GPU so initial error is set correctly (no uninitialized values).
  if (up->use_gpu) {
    gkyl_cu_memcpy(up->error_cu, up->error, sizeof(double[3]), GKYL_CU_MEMCPY_H2D);
  }

  // Clear the differences prior to iteration
  gkyl_array_clear(up->d_moms, 0.0);
  gkyl_array_clear(up->dd_moms, 0.0);

  // Iteration loop, max_iter iterations is usually sufficient for machine precision moments
  while ((ispositive_f_max) &&  ((niter < max_iter) && ((fabs(up->error[0]) > tol) || (fabs(up->error[1]) > tol) ||
    (fabs(up->error[2]) > tol))))
  {
    // 1. Calculate the Maxwellian moments (n, u_par, T/m) from the projected Maxwellian distribution
    gkyl_gyrokinetic_maxwellian_moments_advance(up->moments_up, phase_range, conf_range, f_max, up->moms_iter);

    // a. Calculate  ddMi^(k+1) =  Mi_corr - Mi_new
    // ddn = n_target - n;
    // Compute out = out + a*inp. Returns out.
    gkyl_array_set(up->dd_moms, -1.0, up->moms_iter);
    gkyl_array_accumulate(up->dd_moms, 1.0, moms_target);

    // b. Calculate  dMi^(k+1) = dn^k + ddMi^(k+1) | where dn^0 = 0
    // dm_new = dm_old + ddn;
    gkyl_array_accumulate(up->d_moms, 1.0, up->dd_moms);

    // End the iteration early if all moments converge
    if ((niter % 1) == 0) {
      if (up->use_gpu) {
        // We insure the reduction to find the maximum error is thread-safe on GPUs
        // by first calling a specialized kernel for computing the absolute value 
        // of the difference of the cell averages, then calling reduce_range.
        gkyl_gyrokinetic_maxwellian_correct_all_moments_abs_diff_cu(conf_range, 
          nc, moms_target, up->moms_iter, up->abs_diff_moms);
        gkyl_array_reduce_range(up->error_cu, up->abs_diff_moms, GKYL_MAX, conf_range);
        gkyl_cu_memcpy(up->error, up->error_cu, sizeof(double[3]), GKYL_CU_MEMCPY_D2H);
      }
      else {
        struct gkyl_range_iter biter;

        // Reset the maximum error
        up->error[0] = 0.0;
        up->error[1] = 0.0;
        up->error[2] = 0.0;
        // Iterate over the input configuration-space range to find the maximum error
        gkyl_range_iter_init(&biter, conf_range);
        while (gkyl_range_iter_next(&biter)){
          long midx = gkyl_range_idx(conf_range, biter.idx);
          const double *moms_local = gkyl_array_cfetch(up->moms_iter, midx);
          const double *moms_target_local = gkyl_array_cfetch(moms_target, midx);
          // Check the error in the absolute value of the cell average
          // Note: for density and temperature, this error is a relative error compared to the target moment value
          // so that we can converge to the correct target moments in SI units and minimize finite precision issues.
          up->error[0] = fmax(fabs(moms_local[0*nc] - moms_target_local[0*nc])/moms_target_local[0*nc],fabs(up->error[0]));
          up->error[2] = fmax(fabs(moms_local[2*nc] - moms_target_local[2*nc])/moms_target_local[2*nc],fabs(up->error[2]));
          // However, u_par may be ~ 0 and if it is, we need to use absolute error. We can converge safely using
          // absolute error if u_par ~ O(1). Otherwise, we use relative error for u_par. 
          if (fabs(moms_target_local[1*nc]) < 1.0) {
            up->error[1] = fmax(fabs(moms_local[1*nc] - moms_target_local[1*nc]),fabs(up->error[1]));
          }
          else {
            up->error[1] = fmax(fabs(moms_local[1*nc] - moms_target_local[1*nc])/moms_target_local[1*nc],fabs(up->error[1]));
          }
          // Check if density and temperature are positive, if they aren't we will break out of the iteration
          ispositive_f_max = (moms_local[0*nc]>0.0) && ispositive_f_max;
          ispositive_f_max = (moms_local[2*nc]>0.0) && ispositive_f_max;
        }
      }
    }

    // c. Calculate  n^(k+1) = M^k + dM^(k+1)
    // n = n_target + dm_new;
    gkyl_array_set(up->moms_iter, 1.0, moms_target);
    gkyl_array_accumulate(up->moms_iter, 1.0, up->d_moms);

    // 2. Update the gyrokinetic Maxwellian distribution function using the corrected moments.
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(up->proj_max_prim,
      phase_range, conf_range, up->moms_iter, 
      up->gk_geom->bmag, up->gk_geom->bmag, up->mass, f_max);
    // Correct the density before the next iteration
    gkyl_gyrokinetic_maxwellian_correct_density_moment(up,
      f_max, up->moms_iter, phase_range, conf_range);

    niter += 1;
  }
  if ((niter < max_iter) && (ispositive_f_max) && ((fabs(up->error[0]) < tol) && (fabs(up->error[1]) < tol) &&
    (fabs(up->error[2]) < tol))) {
    corr_status = 0;
  } 
  else {
    corr_status = 1;
  }

  // If the algorithm fails to converge and we are *not* using the results of the failed convergence,
  // we project the distribution function with the basic target moments.
  // We correct the density and then recompute moments/errors for this new projections
  if (corr_status == 1 && !up->use_last_converged) { 
    gkyl_proj_gkmaxwellian_on_basis_prim_mom(up->proj_max_prim,
      phase_range, conf_range, moms_target, 
      up->gk_geom->bmag, up->gk_geom->bmag, up->mass, f_max);
    gkyl_gyrokinetic_maxwellian_correct_density_moment(up,
      f_max, moms_target, phase_range, conf_range);

    gkyl_gyrokinetic_maxwellian_moments_advance(up->moments_up, phase_range, conf_range, f_max, up->moms_iter);

    if (up->use_gpu) {
      // We insure the reduction to find the maximum error is thread-safe on GPUs
      // by first calling a specialized kernel for computing the absolute value 
      // of the difference of the cell averages, then calling reduce_range.
      gkyl_gyrokinetic_maxwellian_correct_all_moments_abs_diff_cu(conf_range, 
        nc, moms_target, up->moms_iter, up->abs_diff_moms);
      gkyl_array_reduce_range(up->error_cu, up->abs_diff_moms, GKYL_MAX, conf_range);
      gkyl_cu_memcpy(up->error, up->error_cu, sizeof(double[3]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      struct gkyl_range_iter biter;

      // Reset the maximum error
      up->error[0] = 0.0;
      up->error[1] = 0.0;
      up->error[2] = 0.0;
      // Iterate over the input configuration-space range to find the maximum error
      gkyl_range_iter_init(&biter, conf_range);
      while (gkyl_range_iter_next(&biter)){
        long midx = gkyl_range_idx(conf_range, biter.idx);
        const double *moms_local = gkyl_array_cfetch(up->moms_iter, midx);
        const double *moms_target_local = gkyl_array_cfetch(moms_target, midx);
        // Check the error in the absolute value of the cell average
        // Note: for density and temperature, this error is a relative error compared to the target moment value
        // so that we can converge to the correct target moments in SI units and minimize finite precision issues.
        up->error[0] = fmax(fabs(moms_local[0*nc] - moms_target_local[0*nc])/moms_target_local[0*nc],fabs(up->error[0]));
        up->error[2] = fmax(fabs(moms_local[2*nc] - moms_target_local[2*nc])/moms_target_local[2*nc],fabs(up->error[2]));
        // However, u_par may be ~ 0 and if it is, we need to use absolute error. We can converge safely using
        // absolute error if u_par ~ O(1). Otherwise, we use relative error for u_par. 
        if (fabs(moms_target_local[1*nc]) < 1.0) {
          up->error[1] = fmax(fabs(moms_local[1*nc] - moms_target_local[1*nc]),fabs(up->error[1]));
        }
        else {
          up->error[1] = fmax(fabs(moms_local[1*nc] - moms_target_local[1*nc])/moms_target_local[1*nc],fabs(up->error[1]));
        }
      }
    }
  }

  return (struct gkyl_gyrokinetic_maxwellian_correct_status) {
    .iter_converged = corr_status,
    .num_iter = niter,
    .error[0] = up->error[0], 
    .error[1] = up->error[1],
    .error[2] = up->error[2], 
  };  
}

void 
gkyl_gyrokinetic_maxwellian_correct_release(gkyl_gyrokinetic_maxwellian_correct *up)
{
  gkyl_array_release(up->moms_iter);
  gkyl_array_release(up->d_moms);
  gkyl_array_release(up->dd_moms);
  gkyl_array_release(up->num_ratio);
  gkyl_dg_bin_op_mem_release(up->mem);
  if (up->use_gpu) {
    gkyl_array_release(up->abs_diff_moms);
    gkyl_cu_free(up->error_cu);
  }

  gkyl_gyrokinetic_maxwellian_moments_release(up->moments_up);
  gkyl_proj_maxwellian_on_basis_release(up->proj_max_prim);

  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

void 
gkyl_gyrokinetic_maxwellian_correct_all_moments_abs_diff_cu(const struct gkyl_range *conf_range, 
  int nc, const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter, 
  struct gkyl_array *moms_abs_diff)
{
  assert(false);
}

#endif
