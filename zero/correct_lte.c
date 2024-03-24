#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_correct_lte.h>
#include <gkyl_correct_lte_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_lte_moments.h>
#include <gkyl_proj_vlasov_lte_on_basis.h>

gkyl_correct_vlasov_lte *
gkyl_correct_vlasov_lte_inew(const struct gkyl_correct_vlasov_lte_inp *inp)
{
  gkyl_correct_vlasov_lte *up = gkyl_malloc(sizeof(*up));
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
    up->error_cu = gkyl_cu_malloc(sizeof(double[5]));
  }
  else {
    up->moms_iter = gkyl_array_new(GKYL_DOUBLE, (up->vdim+2)*inp->conf_basis->num_basis, conf_local_ext_ncells);
    up->d_moms = gkyl_array_new(GKYL_DOUBLE, (up->vdim+2)*inp->conf_basis->num_basis, conf_local_ext_ncells);
    up->dd_moms = gkyl_array_new(GKYL_DOUBLE, (up->vdim+2)*inp->conf_basis->num_basis, conf_local_ext_ncells);
  }

  // Moments structure 
  struct gkyl_lte_moments_vlasov_inp inp_mom = {
    .phase_grid = inp->phase_grid,
    .conf_basis = inp->conf_basis,
    .phase_basis = inp->phase_basis,
    .conf_range =  inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .vel_range = inp->vel_range,
    .p_over_gamma = inp->p_over_gamma,
    .gamma = inp->gamma,
    .gamma_inv = inp->gamma_inv,
    .model_id = inp->model_id,
    .mass = inp->mass,
    .use_gpu = inp->use_gpu,
  };
  up->moments_up = gkyl_lte_moments_inew( &inp_mom );

  // Create a projection updater for projecting the LTE distribution function
  // Projection routine also corrects the density before returning 
  // the LTE distribution function.
  struct gkyl_proj_vlasov_lte_inp inp_proj = {
    .phase_grid = inp->phase_grid,
    .conf_basis = inp->conf_basis,
    .phase_basis = inp->phase_basis,
    .conf_range =  inp->conf_range,
    .conf_range_ext = inp->conf_range_ext,
    .vel_range = inp->vel_range,
    .p_over_gamma = inp->p_over_gamma,
    .gamma = inp->gamma,
    .gamma_inv = inp->gamma_inv,
    .model_id = inp->model_id,
    .mass = inp->mass,
    .use_gpu = inp->use_gpu,
  };
  up->proj_lte = gkyl_proj_vlasov_lte_on_basis_inew( &inp_proj );

  return up;
}

struct gkyl_correct_vlasov_lte_status
gkyl_correct_all_moments_vlasov_lte(gkyl_correct_vlasov_lte *c_corr,
  struct gkyl_array *f_lte, const struct gkyl_array *moms_target, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  int vdim = c_corr->vdim;
  int nc = c_corr->num_conf_basis;
  double tol = c_corr->eps;  // tolerance of the iterative scheme
  int max_iter = c_corr->max_iter;

  int niter = 0;
  bool corr_status = true;

  // Set initial max error 
  for (int i=0; i<5; ++i) {
    if (i < vdim+2) {
      c_corr->error[i] = 1.0;  
    }
    else {
      c_corr->error[i] = 0.0;
    }
  }

  // Clear the differences prior to iteration
  gkyl_array_clear(c_corr->d_moms, 0.0);
  gkyl_array_clear(c_corr->dd_moms, 0.0);

  // Iteration loop, max_iter iterations is usually sufficient (for all vdim) for machine precision moments
  while ((niter < max_iter) && ((fabs(c_corr->error[0]) > tol) || (fabs(c_corr->error[1]) > tol) ||
    (fabs(c_corr->error[2]) > tol) || (fabs(c_corr->error[3]) > tol) || (fabs(c_corr->error[4]) > tol)))
  {
    // 1. Calculate the LTE moments (n, V_drift, T) from the projected LTE distribution
    gkyl_lte_moments_advance(c_corr->moments_up, phase_local, conf_local, f_lte, c_corr->moms_iter);

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
        // Call specialized kernel for finding the error and copying back to CPU

      }
      else {
        struct gkyl_range_iter biter;

        // Reset the maximum error
        for (int i=0; i<vdim+2; ++i) {
          c_corr->error[i] = 0.0;
        }
        // Iterate over the input configuration space range to find the maximum error
        gkyl_range_iter_init(&biter, conf_local);
        while (gkyl_range_iter_next(&biter)){
          long midx = gkyl_range_idx(conf_local, biter.idx);
          const double *moms_local = gkyl_array_cfetch(c_corr->moms_iter, midx);
          const double *moms_target_local = gkyl_array_cfetch(moms_target, midx);
          // Check the error in the absolute value of the cell average
          for (int d=0; d<vdim+2; ++d) {
            c_corr->error[d] = fmax(fabs(moms_local[d*nc] - moms_target_local[d*nc]),fabs(c_corr->error[d]));
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
    gkyl_proj_vlasov_lte_on_basis_advance(c_corr->proj_lte, 
      phase_local, conf_local, c_corr->moms_iter, f_lte);

    niter += 1;
  }
  if ((niter < max_iter) && ((fabs(c_corr->error[0]) < tol) && (fabs(c_corr->error[1]) < tol) &&
    (fabs(c_corr->error[2]) < tol) && (fabs(c_corr->error[3]) < tol) && (fabs(c_corr->error[4]) < tol))) {
    corr_status = 0;
  } 
  else {
    corr_status = 1;
  }

  // If the algorithm fails (density fails to converge)!
  // Project the distribution function with the basic moments.
  // Projection routine internally corrects the density.
  if (corr_status == 1) {
    gkyl_proj_vlasov_lte_on_basis_advance(c_corr->proj_lte, 
      phase_local, conf_local, moms_target, f_lte);
  }

  return (struct gkyl_correct_vlasov_lte_status) {
    .iter_converged = corr_status,
    .num_iter = niter
  };  
}

void 
gkyl_correct_vlasov_lte_release(gkyl_correct_vlasov_lte *c_corr)
{
  gkyl_array_release(c_corr->moms_iter);
  gkyl_array_release(c_corr->d_moms);
  gkyl_array_release(c_corr->dd_moms);

  gkyl_lte_moments_release(c_corr->moments_up);
  gkyl_proj_vlasov_lte_on_basis_release(c_corr->proj_lte);

  gkyl_free(c_corr);
}
