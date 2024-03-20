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
#include <gkyl_maxwellian_moments.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov_sr.h>
#include <gkyl_proj_maxwellian_on_basis.h>
#include <gkyl_proj_mj_on_basis.h>
#include <gkyl_proj_on_basis.h>

gkyl_correct_vlasov_lte *
gkyl_correct_vlasov_lte_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_range_ext, const struct gkyl_range *vel_range, 
  const struct gkyl_array *p_over_gamma, const struct gkyl_array *gamma, const struct gkyl_array *gamma_inv,
  enum gkyl_model_id model_id, double mass, bool use_gpu)
{
  gkyl_correct_vlasov_lte *up = gkyl_malloc(sizeof(*up));

  up->model_id = model_id;
  up->grid = *grid;
  up->conf_basis = *conf_basis;
  up->phase_basis = *phase_basis;
  int vdim = up->phase_basis.ndim - up->conf_basis.ndim;
  int poly_order = conf_basis->poly_order;

  long conf_local_ncells = conf_range->volume;
  long conf_local_ext_ncells = conf_range_ext->volume;

  up->num_ratio = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->num_vb = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->V_drift = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->mem = gkyl_dg_bin_op_mem_new(conf_local_ncells, conf_basis->num_basis);

  // Moment memory
  up->n_stationary = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->vbi = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->T_stationary = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);

  // Create a copy for differences (d) and differences of differences (dd)
  up->dn = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->dvbi = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->dT = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);

  up->ddn = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->ddvbi = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->ddT = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);

  // Make the general moments array
  up->moms = gkyl_array_new(GKYL_DOUBLE, (vdim + 2)*conf_basis->num_basis, conf_local_ext_ncells);

  // Moments structure 
  up->moments_up = gkyl_maxwellian_moments_new(grid, conf_basis, phase_basis, conf_range, 
    conf_range_ext, vel_range, p_over_gamma, gamma, gamma_inv,
    up->model_id, mass, use_gpu);

  // For relativistic/nonrelativistic models:
  if (up->model_id == GKYL_MODEL_SR) {

    // Set auxiliary fields for moment updates. 
    struct gkyl_mom_vlasov_sr_auxfields sr_inp = {.p_over_gamma = p_over_gamma, 
      .gamma = 0, .gamma_inv = 0, .V_drift = 0, 
      .GammaV2 = 0, .GammaV_inv = 0};  

    // Make the projection routine for the new MJ distributions
    up->proj_mj = gkyl_proj_mj_on_basis_new(grid, conf_basis, phase_basis, poly_order + 1);

  } else {

    // Maxwellian Projection routine
    up->proj_max = gkyl_proj_maxwellian_on_basis_new(grid, conf_basis, phase_basis, poly_order + 1, use_gpu);
  }

  return up;
}

void 
gkyl_correct_density_moment_vlasov_lte(gkyl_correct_vlasov_lte *c_corr, 
  struct gkyl_array *fin,
  const struct gkyl_array *n_target,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  // vdim
  int vdim = c_corr->phase_basis.ndim - c_corr->conf_basis.ndim;

  // compute the moments
  gkyl_maxwellian_moments_advance(c_corr->moments_up, phase_local, conf_local, fin, c_corr->moms);

  // Pull apart density moment:
  gkyl_array_set_offset_range(c_corr->num_ratio, 1.0, c_corr->moms, 0*c_corr->conf_basis.num_basis, conf_local);

  // compute number density ratio: num_ratio = n/n0
  gkyl_dg_div_op_range(c_corr->mem, c_corr->conf_basis, 0, c_corr->num_ratio,
    0, n_target, 0, c_corr->num_ratio, conf_local);

  // rescale distribution function
  gkyl_dg_mul_conf_phase_op_range(&c_corr->conf_basis, &c_corr->phase_basis,
    fin, c_corr->num_ratio, fin, conf_local, phase_local);

  // Hand back rescaled n:
  gkyl_array_set(c_corr->n_stationary, 1.0, c_corr->num_ratio);
}

void 
gkyl_correct_all_moments_vlasov_lte(gkyl_correct_vlasov_lte *c_corr,
  struct gkyl_array *distf,
  const struct gkyl_array *n_target, const struct gkyl_array *vbi_target, const struct gkyl_array *T_target,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local,
  int poly_order)
{
  int vdim = c_corr->phase_basis.ndim - c_corr->conf_basis.ndim;

  // Copy the intial moments for m*_corr -> m*
  gkyl_array_set(c_corr->n_stationary, 1.0, n_target);
  gkyl_array_set(c_corr->vbi, 1.0, vbi_target);
  gkyl_array_set(c_corr->T_stationary, 1.0, T_target);

  // 0. Project the MJ with the intially correct moments
  if (c_corr->model_id == GKYL_MODEL_SR) {
    gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(c_corr->proj_mj, phase_local, 
      conf_local, n_target, vbi_target, T_target, distf);
  }
  else {
    // Set the moms_target
    gkyl_array_set_offset_range( c_corr->moms, 1.0, n_target, 0*c_corr->conf_basis.num_basis, conf_local);
    gkyl_array_set_offset_range( c_corr->moms, 1.0, vbi_target, 1*c_corr->conf_basis.num_basis, conf_local);
    gkyl_array_set_offset_range( c_corr->moms, 1.0, T_target, (vdim+1)*c_corr->conf_basis.num_basis, conf_local);
    gkyl_proj_maxwellian_on_basis_lab_mom(c_corr->proj_max, phase_local, 
      conf_local, c_corr->moms_target, distf);
  }

  // tolerance of the iterative scheme
  double tol = 1e-12;
  c_corr->niter = 0;
  c_corr->error_n = 1.0;
  c_corr->error_vb[0] = 1.0;
  c_corr->error_vb[1] = 0.0;
  if (vdim > 1)
    c_corr->error_vb[1] = 1.0;
  c_corr->error_vb[2] = 0.0;
  if (vdim > 2)
    c_corr->error_vb[2] = 1.0;
  c_corr->error_T = 1.0;

  // Iteration loop, 100 iterations is usually sufficient (for all vdim) for machine precision moments
  while ((c_corr->niter < 100) && ((fabs(c_corr->error_n) > tol) || (fabs(c_corr->error_vb[0]) > tol) ||
    (fabs(c_corr->error_vb[1]) > tol) || (fabs(c_corr->error_vb[2]) > tol) || (fabs(c_corr->error_T) > tol)))
  {

    // 1. Calculate the new moments
    // calculate the moments of the dist (n, vb, T -> n, vbi, T)
    gkyl_maxwellian_moments_advance(c_corr->moments_up, phase_local, conf_local, distf, c_corr->moms);
    gkyl_array_set_offset_range(c_corr->n_stationary, 1.0, c_corr->moms, 0*c_corr->conf_basis.num_basis, conf_local);
    gkyl_array_set_offset_range(c_corr->vbi, 1.0, c_corr->moms, 1*c_corr->conf_basis.num_basis, conf_local);
    gkyl_array_set_offset_range(c_corr->T_stationary, 1.0, c_corr->moms, (vdim+1)*c_corr->conf_basis.num_basis, conf_local);
    //gkyl_mj_moments_advance(c_corr->mj_moms, distf_mj, c_corr->n_stationary, c_corr->vbi, c_corr->T_stationary, phase_local, conf_local);

    // a. Calculate  ddMi^(k+1) =  Mi_corr - Mi_new
    // ddn = n_target - n;
    //  Compute out = out + a*inp. Returns out.
    gkyl_array_set(c_corr->ddn, -1.0, c_corr->n_stationary);
    gkyl_array_accumulate(c_corr->ddn, 1.0, n_target);
    gkyl_array_set(c_corr->ddvbi, -1.0, c_corr->vbi);
    gkyl_array_accumulate(c_corr->ddvbi, 1.0, vbi_target);
    gkyl_array_set(c_corr->ddT, -1.0, c_corr->T_stationary);
    gkyl_array_accumulate(c_corr->ddT, 1.0, T_target);

    // b. Calculate  dMi^(k+1) = dn^k + ddMi^(k+1) | where dn^0 = 0
    // dm_new = dm_old + ddn;
    gkyl_array_accumulate(c_corr->dn, 1.0, c_corr->ddn);
    gkyl_array_accumulate(c_corr->dvbi, 1.0, c_corr->ddvbi);
    gkyl_array_accumulate(c_corr->dT, 1.0, c_corr->ddT);

    // End the iteration early if all moments converge
    if ((c_corr->niter % 1) == 0){
      struct gkyl_range_iter biter;

      // Reset the maximum error
      c_corr->error_n = 0; c_corr->error_T = 0;
      c_corr->error_vb[0] = 0; c_corr->error_vb[1] = 0; c_corr->error_vb[2] = 0;

      // Iterate over the grid to find the maximum error
      gkyl_range_iter_init(&biter, conf_local);
      while (gkyl_range_iter_next(&biter)){
        long midx = gkyl_range_idx(conf_local, biter.idx);
        const double *n_local = gkyl_array_cfetch(c_corr->n_stationary, midx);
        const double *vbi_local = gkyl_array_cfetch(c_corr->vbi, midx);
        const double *T_local = gkyl_array_cfetch(c_corr->T_stationary, midx);
        const double *n_original_local = gkyl_array_cfetch(n_target, midx);
        const double *vbi_original_local = gkyl_array_cfetch(vbi_target, midx);
        const double *T_original_local = gkyl_array_cfetch(T_target, midx);
        c_corr->error_n = fmax(fabs(n_local[0] - n_original_local[0]),fabs(c_corr->error_n));
        c_corr->error_vb[0] = fmax(fabs(vbi_local[0] - vbi_original_local[0]),fabs(c_corr->error_vb[0]));
        c_corr->error_T = fmax(fabs(T_local[0] - T_original_local[0]),fabs(c_corr->error_T));
        if (vdim > 1)
          c_corr->error_vb[1] = fmax(fabs(vbi_local[poly_order + 1] - vbi_original_local[poly_order + 1]),fabs(c_corr->error_vb[1]));
        if (vdim > 2)
          c_corr->error_vb[2] = fmax(fabs(vbi_local[2 * (poly_order + 1)] - vbi_original_local[2 * (poly_order + 1)]),fabs(c_corr->error_vb[2]));
      }
    }

    // c. Calculate  n^(k+1) = M^k + dM^(k+1)
    // n = n_target + dm_new;
    gkyl_array_set(c_corr->n_stationary, 1.0, n_target);
    gkyl_array_accumulate(c_corr->n_stationary, 1.0, c_corr->dn);
    gkyl_array_set(c_corr->vbi, 1.0, vbi_target);
    gkyl_array_accumulate(c_corr->vbi, 1.0, c_corr->dvbi);
    gkyl_array_set(c_corr->T_stationary, 1.0, T_target);
    gkyl_array_accumulate(c_corr->T_stationary, 1.0, c_corr->dT);

    // 2. Update the dist_mj using the corrected moments
    if (c_corr->model_id == GKYL_MODEL_SR) {
      gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(c_corr->proj_mj, phase_local, 
        conf_local, c_corr->n_stationary, c_corr->vbi, c_corr->T_stationary, distf);
    }
    else {
      // Set the moms_target
      gkyl_array_set_offset_range( c_corr->moms, 1.0, c_corr->n_stationary, 0*c_corr->conf_basis.num_basis, conf_local);
      gkyl_array_set_offset_range( c_corr->moms, 1.0, c_corr->vbi, 1*c_corr->conf_basis.num_basis, conf_local);
      gkyl_array_set_offset_range( c_corr->moms, 1.0, c_corr->T_stationary, (vdim+1)*c_corr->conf_basis.num_basis, conf_local);
      gkyl_proj_maxwellian_on_basis_lab_mom(c_corr->proj_max, phase_local, 
        conf_local, c_corr->moms_target, distf);
    }

    // 3. Correct the n moment to fix the asymptotically approximated MJ function
    gkyl_correct_density_moment_vlasov_lte(c_corr, distf, c_corr->n_stationary, phase_local, conf_local);

    c_corr->niter += 1;
  }
  if ((c_corr->niter < 100) && ((fabs(c_corr->error_n) < tol) && (fabs(c_corr->error_vb[0]) < tol) &&
    (fabs(c_corr->error_vb[1]) < tol) && (fabs(c_corr->error_vb[2]) < tol) && (fabs(c_corr->error_T) < tol)))
    c_corr->status = 0;
  else
    c_corr->status = 1;

  // If the algorithm fails (density fails to converge)!
  // Project the distribution function with the basic moments and correct n
  if (c_corr->status == 1){
    if (c_corr->model_id == GKYL_MODEL_SR) {
      gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(c_corr->proj_mj, phase_local, 
        conf_local, n_target, vbi_target, T_target, distf);
    }
    else {
      gkyl_proj_maxwellian_on_basis_lab_mom(c_corr->proj_max, phase_local, 
        conf_local, c_corr->moms_target, distf);
    }
    gkyl_correct_density_moment_vlasov_lte(c_corr, distf, n_target, phase_local, conf_local);
  }
}

void 
gkyl_correct_vlasov_lte_release(gkyl_correct_vlasov_lte *c_corr)
{

  // gkyl_mom_calc_release(c_corr->ncalc);
  gkyl_array_release(c_corr->num_ratio);
  gkyl_array_release(c_corr->num_vb);
  gkyl_array_release(c_corr->V_drift);
  gkyl_dg_bin_op_mem_release(c_corr->mem);

  gkyl_array_release(c_corr->n_stationary);
  gkyl_array_release(c_corr->vbi);
  gkyl_array_release(c_corr->T_stationary);
  gkyl_array_release(c_corr->dn);
  gkyl_array_release(c_corr->dvbi);
  gkyl_array_release(c_corr->dT);
  gkyl_array_release(c_corr->ddn);
  gkyl_array_release(c_corr->ddvbi);
  gkyl_array_release(c_corr->ddT);

  gkyl_maxwellian_moments_release(c_corr->moments_up);
  gkyl_array_release(c_corr->moms);
  gkyl_array_release(c_corr->moms_target);

  if (c_corr->model_id == GKYL_MODEL_SR) {
    gkyl_proj_mj_on_basis_release(c_corr->proj_mj);
  }
  else {
    gkyl_proj_maxwellian_on_basis_release(c_corr->proj_max);
  }

  gkyl_free(c_corr);
}