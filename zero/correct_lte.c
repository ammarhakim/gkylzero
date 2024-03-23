#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_correct_lte.h>
#include <gkyl_correct_lte_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_lte_moments.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov_sr.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_mj_on_basis.h>
#include <gkyl_proj_on_basis.h>

gkyl_correct_vlasov_lte *
gkyl_correct_vlasov_lte_inew(const struct gkyl_correct_vlasov_lte_inp *inp)
{
  gkyl_correct_vlasov_lte *up = gkyl_malloc(sizeof(*up));
  up->eps = inp->eps;
  up->max_iter = inp->max_iter;

  up->model_id = inp->model_id;
  up->conf_basis = *inp->conf_basis;
  up->phase_basis = *inp->phase_basis;
  up->vdim = up->phase_basis.ndim - up->conf_basis.ndim;
  up->num_conf_basis = up->conf_basis.num_basis;

  long conf_local_ncells = inp->conf_range->volume;
  long conf_local_ext_ncells = inp->conf_range_ext->volume;

  // Number density ratio: num_ratio = n_target/n0 and bin_op memory to compute ratio
  // Used for fixing the density with simple rescaling
  up->moms_dens_corr = gkyl_array_new(GKYL_DOUBLE, (up->vdim+2)*inp->conf_basis->num_basis, conf_local_ext_ncells);
  up->num_ratio = gkyl_array_new(GKYL_DOUBLE, inp->conf_basis->num_basis, conf_local_ext_ncells);
  up->mem = gkyl_dg_bin_op_mem_new(conf_local_ncells, inp->conf_basis->num_basis);

  // Individual moment memory: the iteration of the moments, the differences (d) and differences of differences (dd)
  up->moms_iter = gkyl_array_new(GKYL_DOUBLE, (up->vdim+2)*inp->conf_basis->num_basis, conf_local_ext_ncells);
  up->d_moms = gkyl_array_new(GKYL_DOUBLE, (up->vdim+2)*inp->conf_basis->num_basis, conf_local_ext_ncells);
  up->dd_moms = gkyl_array_new(GKYL_DOUBLE, (up->vdim+2)*inp->conf_basis->num_basis, conf_local_ext_ncells);

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

  // For relativistic/nonrelativistic models:
  if (up->model_id == GKYL_MODEL_SR) {
    up->proj_mj = gkyl_proj_mj_on_basis_new(inp->phase_grid, inp->conf_basis, inp->phase_basis, 
      inp->conf_basis->poly_order+1, inp->use_gpu);
  } 
  else {
    up->proj_max = gkyl_proj_maxwellian_on_basis_new(inp->phase_grid, inp->conf_basis, inp->phase_basis, 
      inp->conf_basis->poly_order+1, inp->use_gpu);
  }

  return up;
}

struct gkyl_correct_vlasov_lte_status
gkyl_correct_density_moment_vlasov_lte(gkyl_correct_vlasov_lte *c_corr, 
  struct gkyl_array *f_lte, const struct gkyl_array *moms_target, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  int nc = c_corr->num_conf_basis;
  
  // compute the moments
  gkyl_lte_moments_advance(c_corr->moments_up, phase_local, conf_local, f_lte, c_corr->moms_dens_corr);

  // Fetch the density moment from the input LTE moments (0th component)
  gkyl_array_set_offset_range(c_corr->num_ratio, 1.0, c_corr->moms_dens_corr, 0*nc, conf_local);

  // compute number density ratio: num_ratio = n/n0
  // 0th component of moms_target is the target density
  gkyl_dg_div_op_range(c_corr->mem, c_corr->conf_basis, 0, c_corr->num_ratio,
    0, moms_target, 0, c_corr->num_ratio, conf_local);

  // rescale distribution function
  gkyl_dg_mul_conf_phase_op_range(&c_corr->conf_basis, &c_corr->phase_basis,
    f_lte, c_corr->num_ratio, f_lte, conf_local, phase_local);

  return (struct gkyl_correct_vlasov_lte_status) {
    .iter_converged = true,
    .num_iter = 1
  };
}

struct gkyl_correct_vlasov_lte_status
gkyl_correct_all_moments_vlasov_lte(gkyl_correct_vlasov_lte *c_corr,
  struct gkyl_array *f_lte, const struct gkyl_array *moms_target, 
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  int vdim = c_corr->vdim;
  int nc = c_corr->num_conf_basis;
  double eps = c_corr->eps;
  int max_iter = c_corr->max_iter;

  // tolerance of the iterative scheme
  double tol = eps;
  int niter = 0;
  bool corr_status = true;
  
  c_corr->error_n = 1.0;
  c_corr->error_vb[0] = 1.0;
  c_corr->error_vb[1] = 0.0;
  if (c_corr->vdim > 1)
    c_corr->error_vb[1] = 1.0;
  c_corr->error_vb[2] = 0.0;
  if (c_corr->vdim > 2)
    c_corr->error_vb[2] = 1.0;
  c_corr->error_T = 1.0;

  // Clear the differences prior to iteration
  gkyl_array_clear(c_corr->d_moms, 0.0);
  gkyl_array_clear(c_corr->dd_moms, 0.0);

  // Iteration loop, max_iter iterations is usually sufficient (for all vdim) for machine precision moments
  while ((niter < max_iter) && ((fabs(c_corr->error_n) > tol) || (fabs(c_corr->error_vb[0]) > tol) ||
    (fabs(c_corr->error_vb[1]) > tol) || (fabs(c_corr->error_vb[2]) > tol) || (fabs(c_corr->error_T) > tol)))
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
    if ((niter % 1) == 0){
      struct gkyl_range_iter biter;

      // Reset the maximum error
      c_corr->error_n = 0; c_corr->error_T = 0;
      c_corr->error_vb[0] = 0; c_corr->error_vb[1] = 0; c_corr->error_vb[2] = 0;

      // Iterate over the input configuration space range to find the maximum error
      gkyl_range_iter_init(&biter, conf_local);
      while (gkyl_range_iter_next(&biter)){
        long midx = gkyl_range_idx(conf_local, biter.idx);
        const double *moms_local = gkyl_array_cfetch(c_corr->moms_iter, midx);
        const double *moms_target_local = gkyl_array_cfetch(moms_target, midx);

        c_corr->error_n = fmax(fabs(moms_local[0] - moms_target_local[0]),fabs(c_corr->error_n));
        for (int d=0; d<vdim; ++d) {
          c_corr->error_vb[d] = fmax(fabs(moms_local[(d+1)*nc] - moms_target_local[(d+1)*nc]), fabs(c_corr->error_vb[d]));
        }
        c_corr->error_T = fmax(fabs(moms_local[(vdim+1)*nc] - moms_target_local[(vdim+1)*nc]), fabs(c_corr->error_T));
      }
    }

    // c. Calculate  n^(k+1) = M^k + dM^(k+1)
    // n = n_target + dm_new;
    gkyl_array_set(c_corr->moms_iter, 1.0, moms_target);
    gkyl_array_accumulate(c_corr->moms_iter, 1.0, c_corr->d_moms);

    // 2. Update the LTE distribution function using the corrected moments
    if (c_corr->model_id == GKYL_MODEL_SR) {
      gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(c_corr->proj_mj, 
        phase_local, conf_local, c_corr->moms_iter, f_lte);
    }
    else {
      gkyl_proj_maxwellian_on_basis_prim_mom(c_corr->proj_max, 
        phase_local, conf_local, c_corr->moms_iter, f_lte);
    }
    // 3. Correct the n moment to fix the asymptotically approximated MJ function
    gkyl_correct_density_moment_vlasov_lte(c_corr, f_lte, c_corr->moms_iter, phase_local, conf_local);

    niter += 1;
  }
  if ((niter < max_iter) && ((fabs(c_corr->error_n) < tol) && (fabs(c_corr->error_vb[0]) < tol) &&
    (fabs(c_corr->error_vb[1]) < tol) && (fabs(c_corr->error_vb[2]) < tol) && (fabs(c_corr->error_T) < tol))) {
    corr_status = 0;
  } 
  else {
    corr_status = 1;
  }

  // If the algorithm fails (density fails to converge)!
  // Project the distribution function with the basic moments and correct n
  if (corr_status == 1) {
    if (c_corr->model_id == GKYL_MODEL_SR) {
      gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(c_corr->proj_mj, 
        phase_local, conf_local, moms_target, f_lte);
    }
    else {
      gkyl_proj_maxwellian_on_basis_prim_mom(c_corr->proj_max, 
        phase_local, conf_local, moms_target, f_lte);
    }
    gkyl_correct_density_moment_vlasov_lte(c_corr, f_lte, moms_target, phase_local, conf_local);
  }

  return (struct gkyl_correct_vlasov_lte_status) {
    .iter_converged = corr_status,
    .num_iter = niter
  };  
}

void 
gkyl_correct_vlasov_lte_release(gkyl_correct_vlasov_lte *c_corr)
{
  gkyl_array_release(c_corr->num_ratio);
  gkyl_dg_bin_op_mem_release(c_corr->mem);

  gkyl_array_release(c_corr->moms_iter);
  gkyl_array_release(c_corr->d_moms);
  gkyl_array_release(c_corr->dd_moms);

  gkyl_lte_moments_release(c_corr->moments_up);
  if (c_corr->model_id == GKYL_MODEL_SR) {
    gkyl_proj_mj_on_basis_release(c_corr->proj_mj);
  }
  else {
    gkyl_proj_maxwellian_on_basis_release(c_corr->proj_max);
  }

  gkyl_free(c_corr);
}
