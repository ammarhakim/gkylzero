#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_array_reduce.h>
#include <gkyl_vlasov_lte_correct.h>
#include <gkyl_vlasov_lte_correct_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_vlasov_lte_moments.h>
#include <gkyl_vlasov_lte_proj_on_basis.h>

#include <assert.h>

struct gkyl_vlasov_lte_correct*
gkyl_vlasov_lte_correct_inew(const struct gkyl_vlasov_lte_correct_inp *inp)
{
  gkyl_vlasov_lte_correct *up = gkyl_malloc(sizeof(*up));
  up->eps = inp->eps;
  up->max_iter = inp->max_iter;
  up->use_gpu = inp->use_gpu;

  up->model_id = inp->model_id;
  up->conf_basis = *inp->conf_basis;
  up->phase_basis = *inp->phase_basis;
  up->vdim = up->phase_basis.ndim - up->conf_basis.ndim;
  up->num_conf_basis = up->conf_basis.num_basis;

  long conf_local_ncells = inp->conf_range->volume;
  long conf_local_ext_ncells = inp->conf_range_ext->volume;

  // Individual moment memory: the iteration of the moments, the differences (d) and differences of differences (dd)
  if (up->use_gpu) {
    up->moms_iter = gkyl_array_cu_dev_new(GKYL_DOUBLE, (up->vdim+2)*inp->conf_basis->num_basis, conf_local_ext_ncells);
    up->d_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, (up->vdim+2)*inp->conf_basis->num_basis, conf_local_ext_ncells);
    up->dd_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, (up->vdim+2)*inp->conf_basis->num_basis, conf_local_ext_ncells);
    // Two additional GPU-specific allocations for iterating over the grid to find the absolute value of 
    // the difference between the target and iterative moments, and the GPU-side array for performing the
    // thread-safe reduction to find the maximum error on the grid.
    up->abs_diff_moms = gkyl_array_cu_dev_new(GKYL_DOUBLE, (up->vdim+2), conf_local_ext_ncells);
    up->error_cu = gkyl_cu_malloc(sizeof(double[5]));
  }
  else {
    up->moms_iter = gkyl_array_new(GKYL_DOUBLE, (up->vdim+2)*inp->conf_basis->num_basis, conf_local_ext_ncells);
    up->d_moms = gkyl_array_new(GKYL_DOUBLE, (up->vdim+2)*inp->conf_basis->num_basis, conf_local_ext_ncells);
    up->dd_moms = gkyl_array_new(GKYL_DOUBLE, (up->vdim+2)*inp->conf_basis->num_basis, conf_local_ext_ncells);
  }

  // Moments structure 
  struct gkyl_vlasov_lte_moments_inp inp_mom = {
    .phase_grid = inp->phase_grid,
    .conf_basis = inp->conf_basis,
    .phase_basis = inp->phase_basis,
    .conf_range =  inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .vel_range = inp->vel_range,
    .p_over_gamma = inp->p_over_gamma,
    .gamma = inp->gamma,
    .gamma_inv = inp->gamma_inv,
    .h_ij_inv = inp->h_ij_inv,
    .det_h = inp->det_h,
    .model_id = inp->model_id,
    .mass = inp->mass,
    .use_gpu = inp->use_gpu,
  };
  up->moments_up = gkyl_vlasov_lte_moments_inew( &inp_mom );

  // Create a projection updater for projecting the LTE distribution function
  // Projection routine also corrects the density before returning 
  // the LTE distribution function.
  struct gkyl_vlasov_lte_proj_on_basis_inp inp_proj = {
    .phase_grid = inp->phase_grid,
    .conf_basis = inp->conf_basis,
    .phase_basis = inp->phase_basis,
    .conf_range =  inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .vel_range = inp->vel_range,
    .p_over_gamma = inp->p_over_gamma,
    .gamma = inp->gamma,
    .gamma_inv = inp->gamma_inv,
    .h_ij_inv = inp->h_ij_inv,  
    .det_h = inp->det_h,
    .model_id = inp->model_id,
    .mass = inp->mass,
    .use_gpu = inp->use_gpu,
  };
  up->proj_lte = gkyl_vlasov_lte_proj_on_basis_inew( &inp_proj );

  return up;
}

struct gkyl_vlasov_lte_correct_status
gkyl_vlasov_lte_correct_all_moments(gkyl_vlasov_lte_correct *c_corr,
  struct gkyl_array *f_lte, const struct gkyl_array *moms_target, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  int vdim = c_corr->vdim;
  int nc = c_corr->num_conf_basis;
  double tol = c_corr->eps;  // tolerance of the iterative scheme
  int max_iter = c_corr->max_iter;

  int niter = 0;
  bool corr_status = true;
  int ispositive_flte = true;

  // Set initial max error to start the iteration.
  for (int i=0; i<5; ++i) {
    if (i < vdim+2) {
      c_corr->error[i] = 1.0;  
    }
    else {
      c_corr->error[i] = 0.0;
    }
  }
  // Copy the initial max error to GPU so initial error is set correctly (no uninitialized values).
  if (c_corr->use_gpu) {
    gkyl_cu_memcpy(c_corr->error_cu, c_corr->error, sizeof(double[5]), GKYL_CU_MEMCPY_H2D);
  }

  // Clear the differences prior to iteration
  gkyl_array_clear(c_corr->d_moms, 0.0);
  gkyl_array_clear(c_corr->dd_moms, 0.0);

  // Iteration loop, max_iter iterations is usually sufficient (for all vdim) for machine precision moments
  while ((ispositive_flte) &&  ((niter < max_iter) && ((fabs(c_corr->error[0]) > tol) || (fabs(c_corr->error[1]) > tol) ||
    (fabs(c_corr->error[2]) > tol) || (fabs(c_corr->error[3]) > tol) || (fabs(c_corr->error[4]) > tol))))
  {
    // 1. Calculate the LTE moments (n, V_drift, T) from the projected LTE distribution
    gkyl_vlasov_lte_moments_advance(c_corr->moments_up, phase_local, conf_local, f_lte, c_corr->moms_iter);

    // a. Calculate  ddMi^(k+1) =  Mi_corr - Mi_new
    // ddn = n_target - n;
    // Compute out = out + a*inp. Returns out.
    gkyl_array_set(c_corr->dd_moms, -1.0, c_corr->moms_iter);
    gkyl_array_accumulate(c_corr->dd_moms, 1.0, moms_target);

    // b. Calculate  dMi^(k+1) = dn^k + ddMi^(k+1) | where dn^0 = 0
    // dm_new = dm_old + ddn;
    gkyl_array_accumulate(c_corr->d_moms, 1.0, c_corr->dd_moms);

    // End the iteration early if all moments converge
    if ((niter % 1) == 0) {
      if (c_corr->use_gpu) {
        // We insure the reduction to find the maximum error is thread-safe on GPUs
        // by first calling a specialized kernel for computing the absolute value 
        // of the difference of the cell averages, then calling reduce_range.
        gkyl_vlasov_lte_correct_all_moments_abs_diff_cu(conf_local, 
          vdim, nc, moms_target, c_corr->moms_iter, c_corr->abs_diff_moms);
        gkyl_array_reduce_range(c_corr->error_cu, c_corr->abs_diff_moms, GKYL_MAX, conf_local);
        gkyl_cu_memcpy(c_corr->error, c_corr->error_cu, sizeof(double[5]), GKYL_CU_MEMCPY_D2H);
      }
      else {
        struct gkyl_range_iter biter;

        // Reset the maximum error
        for (int i=0; i<vdim+2; ++i) {
          c_corr->error[i] = 0.0;
        }
        // Iterate over the input configuration-space range to find the maximum error
        printf("\n------------------ Per Cell Error (targ tol: %1.3e) -------------------\n",tol);
        gkyl_range_iter_init(&biter, conf_local);
        while (gkyl_range_iter_next(&biter)){
          long midx = gkyl_range_idx(conf_local, biter.idx);
          const double *moms_local = gkyl_array_cfetch(c_corr->moms_iter, midx);
          const double *moms_target_local = gkyl_array_cfetch(moms_target, midx);
          // Check the error in the absolute value of the cell average
          for (int d=0; d<vdim+2; ++d) {
            c_corr->error[d] = fmax(fabs(moms_local[d*nc] - moms_target_local[d*nc]),fabs(c_corr->error[d]));
            if (midx == 35){
              printf("niter[%d] at index %d = %ld, error: %1.16e",niter, d, midx,c_corr->error[d]);  // Printing error for each dimension
              printf(" (M): %1.16e (M_targ.): %1.16e \n",moms_local[d*nc],moms_target_local[d*nc]);
            }
            if (d == 0 || d == vdim+1) // Temp or density
              ispositive_flte = ((moms_local[d*nc]>0)) && ispositive_flte ;
          }
        }
      }
    }

    // c. Calculate  n^(k+1) = M^k + dM^(k+1)
    // n = n_target + dm_new;
    gkyl_array_set(c_corr->moms_iter, 1.0, moms_target);
    gkyl_array_accumulate(c_corr->moms_iter, 1.0, c_corr->d_moms);

    // 2. Update the LTE distribution function using the corrected moments.
    // Projection routine also corrects the density before the next iteration.
    gkyl_vlasov_lte_proj_on_basis_advance(c_corr->proj_lte, 
      phase_local, conf_local, c_corr->moms_iter, f_lte);

    niter += 1;
  }
  if ((niter < max_iter) && (ispositive_flte) && ((fabs(c_corr->error[0]) < tol) && (fabs(c_corr->error[1]) < tol) &&
    (fabs(c_corr->error[2]) < tol) && (fabs(c_corr->error[3]) < tol) && (fabs(c_corr->error[4]) < tol))) {
    corr_status = 0;
  } 
  else {
    corr_status = 1;
  }

  // If the algorithm fails (density fails to converge)!
  // Project the distribution function with the basic moments.
  // Projection routine internally corrects the density.
  // Recomputes moments/errors for this new projections
  if (corr_status == 1) {
    gkyl_vlasov_lte_proj_on_basis_advance(c_corr->proj_lte, 
      phase_local, conf_local, moms_target, f_lte);

    printf("Correction Algorithm failed!\n");

    gkyl_vlasov_lte_moments_advance(c_corr->moments_up, phase_local, conf_local, f_lte, c_corr->moms_iter);

    if (c_corr->use_gpu) {
      // We insure the reduction to find the maximum error is thread-safe on GPUs
      // by first calling a specialized kernel for computing the absolute value 
      // of the difference of the cell averages, then calling reduce_range.
      gkyl_vlasov_lte_correct_all_moments_abs_diff_cu(conf_local, 
        vdim, nc, moms_target, c_corr->moms_iter, c_corr->abs_diff_moms);
      gkyl_array_reduce_range(c_corr->error_cu, c_corr->abs_diff_moms, GKYL_MAX, conf_local);
      gkyl_cu_memcpy(c_corr->error, c_corr->error_cu, sizeof(double[5]), GKYL_CU_MEMCPY_D2H);
    }
    else {
      struct gkyl_range_iter biter;

      // Reset the maximum error
      for (int i=0; i<vdim+2; ++i) {
        c_corr->error[i] = 0.0;
      }
        // Iterate over the input configuration-space range to find the maximum error
        gkyl_range_iter_init(&biter, conf_local);
        while (gkyl_range_iter_next(&biter)){
          long midx = gkyl_range_idx(conf_local, biter.idx);
          const double *moms_local = gkyl_array_cfetch(c_corr->moms_iter, midx);
          const double *moms_target_local = gkyl_array_cfetch(moms_target, midx);
          // Check the error in the absolute value of the cell average
          for (int d=0; d<vdim+2; ++d) {
            c_corr->error[d] = fmax(fabs(moms_local[d*nc] - moms_target_local[d*nc]),fabs(c_corr->error[d]));
            if (d == 0 || d == vdim+1) // Temp or density
              ispositive_flte = ((moms_local[d*nc]>0)) && ispositive_flte ;
          }
        }
    }
  }

  return (struct gkyl_vlasov_lte_correct_status) {
    .iter_converged = corr_status,
    .num_iter = niter,
    .error[0] = c_corr->error[0], 
    .error[1] = c_corr->error[1],
    .error[2] = c_corr->error[2], 
    .error[3] = c_corr->error[3], 
    .error[4] = c_corr->error[4], 
  };  
}

void 
gkyl_vlasov_lte_correct_release(gkyl_vlasov_lte_correct *c_corr)
{
  gkyl_array_release(c_corr->moms_iter);
  gkyl_array_release(c_corr->d_moms);
  gkyl_array_release(c_corr->dd_moms);
  if (c_corr->use_gpu) {
    gkyl_array_release(c_corr->abs_diff_moms);
    gkyl_cu_free(c_corr->error_cu);
  }

  gkyl_vlasov_lte_moments_release(c_corr->moments_up);
  gkyl_vlasov_lte_proj_on_basis_release(c_corr->proj_lte);

  gkyl_free(c_corr);
}

#ifndef GKYL_HAVE_CUDA

void 
gkyl_vlasov_lte_correct_all_moments_abs_diff_cu(const struct gkyl_range *conf_range, 
  int vdim, int nc, 
  const struct gkyl_array *moms_target, const struct gkyl_array *moms_iter, 
  struct gkyl_array *moms_abs_diff)
{
  assert(false);
}

#endif
