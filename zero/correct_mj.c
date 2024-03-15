#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_correct_mj.h>
#include <gkyl_correct_mj_priv.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_mj_moments.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov_sr.h>
#include <gkyl_proj_mj_on_basis.h>
#include <gkyl_proj_on_basis.h>

gkyl_correct_mj *
gkyl_correct_mj_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_range_ext, const struct gkyl_range *vel_range, 
  const struct gkyl_array *p_over_gamma, const struct gkyl_array *gamma, const struct gkyl_array *gamma_inv, 
  bool use_gpu)
{
  gkyl_correct_mj *up = gkyl_malloc(sizeof(*up));

  up->grid = *grid;
  up->conf_basis = *conf_basis;
  up->phase_basis = *phase_basis;
  int vdim = up->phase_basis.ndim - up->conf_basis.ndim;
  int poly_order = conf_basis->poly_order;

  long conf_local_ncells = conf_range->volume;
  long conf_local_ext_ncells = conf_range_ext->volume;

  // Set auxiliary fields for moment updates. 
  struct gkyl_mom_vlasov_sr_auxfields sr_inp = {.p_over_gamma = p_over_gamma, 
    .gamma = 0, .gamma_inv = 0, .V_drift = 0, 
    .GammaV2 = 0, .GammaV_inv = 0};  

  // updated moment calculator for sr N and N*vb moments
  up->ncalc = gkyl_dg_updater_moment_new(grid, conf_basis,
    phase_basis, conf_range, vel_range, GKYL_MODEL_SR, &sr_inp, "M0", 0, 1, use_gpu);
  up->vbicalc = gkyl_dg_updater_moment_new(grid, conf_basis,
    phase_basis, conf_range, vel_range, GKYL_MODEL_SR, &sr_inp, "M1i", 0, 1, use_gpu);

  up->num_ratio = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->num_vb = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->V_drift = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->vb_dot_nvb = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->n_minus_vb_dot_nvb = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->gamma = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->mem = gkyl_dg_bin_op_mem_new(conf_local_ncells, conf_basis->num_basis);

  // Moment memory
  up->n = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->vbi = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->T = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);

  // Create a copy for differences (d) and differences of differences (dd)
  up->dn = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->dvbi = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->dT = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);

  up->ddn = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->ddvbi = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->ddT = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);

  // Make mj moms stucture
  up->mj_moms = gkyl_mj_moments_new(grid, conf_basis, phase_basis, 
    conf_range, vel_range, conf_local_ncells, conf_local_ext_ncells, 
    p_over_gamma, gamma, gamma_inv, false);

  // Make the projection routine for the new MJ distributions
  up->proj_mj = gkyl_proj_mj_on_basis_new(grid,
     conf_basis, phase_basis, poly_order + 1);

  return up;
}

void 
gkyl_correct_mj_fix_n_stationary(gkyl_correct_mj *cmj, 
  struct gkyl_array *fout,
  const struct gkyl_array *n, const struct gkyl_array *vbi,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  // vdim
  int vdim = cmj->phase_basis.ndim - cmj->conf_basis.ndim;

  // compute the sr moments
  gkyl_dg_updater_moment_advance(cmj->ncalc, phase_local, conf_local, fout, cmj->num_ratio);
  gkyl_dg_updater_moment_advance(cmj->vbicalc, phase_local, conf_local, fout, cmj->num_vb);

  // isolate vb by dividing N*vb by N
  for (int d = 0; d < vdim; ++d)
    gkyl_dg_div_op_range(cmj->mem, cmj->conf_basis, d, cmj->V_drift,
    d, cmj->num_vb, 0, cmj->num_ratio, conf_local);

  // calculate gamma from cmj->num_vb
  gkyl_calc_sr_vars_Gamma(&cmj->conf_basis, &cmj->phase_basis,
    conf_local, cmj->V_drift, cmj->gamma);

  // calculate the stationary density, n0 = gamma*(N - vb dot Nvb)
  gkyl_array_clear_range(cmj->vb_dot_nvb, 0.0, conf_local);
  gkyl_array_set(cmj->n_minus_vb_dot_nvb, 1.0, cmj->num_ratio);
  gkyl_dg_dot_product_op_range(cmj->conf_basis,cmj->vb_dot_nvb,cmj->V_drift,cmj->num_vb, conf_local);
  gkyl_array_accumulate_range(cmj->n_minus_vb_dot_nvb, -1.0, cmj->vb_dot_nvb, conf_local);
  gkyl_dg_mul_op_range(cmj->conf_basis,0,cmj->num_ratio,0,cmj->gamma,0,cmj->n_minus_vb_dot_nvb, conf_local);

  // compute number density ratio: num_ratio = n/n0
  gkyl_dg_div_op_range(cmj->mem, cmj->conf_basis, 0, cmj->num_ratio,
    0, n, 0, cmj->num_ratio, conf_local);

  // rescale distribution function
  gkyl_dg_mul_conf_phase_op_range(&cmj->conf_basis, &cmj->phase_basis,
    fout, cmj->num_ratio, fout, conf_local, phase_local);

  // Hand back rescaled n:
  gkyl_array_set(cmj->n, 1.0, cmj->num_ratio);
}

void 
gkyl_correct_mj_fix(gkyl_correct_mj *cmj,
  struct gkyl_array *distf_mj,
  const struct gkyl_array *n_target, const struct gkyl_array *vbi_target, const struct gkyl_array *T_target,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local,
  int poly_order)
{
  int vdim = cmj->phase_basis.ndim - cmj->conf_basis.ndim;

  // Copy the intial moments for m*_corr -> m*
  gkyl_array_set(cmj->n, 1.0, n_target);
  gkyl_array_set(cmj->vbi, 1.0, vbi_target);
  gkyl_array_set(cmj->T, 1.0, T_target);

  // 0. Project the MJ with the intially correct moments
  gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(cmj->proj_mj, phase_local, 
    conf_local, n_target, vbi_target, T_target, distf_mj);

  // tolerance of the iterative scheme
  double tol = 1e-12;
  cmj->niter = 0;
  cmj->error_n = 1.0;
  cmj->error_vb[0] = 1.0;
  cmj->error_vb[1] = 0.0;
  if (vdim > 1)
    cmj->error_vb[1] = 1.0;
  cmj->error_vb[2] = 0.0;
  if (vdim > 2)
    cmj->error_vb[2] = 1.0;
  cmj->error_T = 1.0;

  // Iteration loop, 100 iterations is usually sufficient (for all vdim) for machine precision moments
  while ((cmj->niter < 100) && ((fabs(cmj->error_n) > tol) || (fabs(cmj->error_vb[0]) > tol) ||
    (fabs(cmj->error_vb[1]) > tol) || (fabs(cmj->error_vb[2]) > tol) || (fabs(cmj->error_T) > tol)))
  {

    // 1. Calculate the new moments
    // calculate the moments of the dist (n, vb, T -> n, vbi, T)
    gkyl_mj_moments_advance(cmj->mj_moms, distf_mj, cmj->n, cmj->vbi, cmj->T, phase_local, conf_local);

    // a. Calculate  ddMi^(k+1) =  Mi_corr - Mi_new
    // ddn = n_target - n;
    //  Compute out = out + a*inp. Returns out.
    gkyl_array_set(cmj->ddn, -1.0, cmj->n);
    gkyl_array_accumulate(cmj->ddn, 1.0, n_target);
    gkyl_array_set(cmj->ddvbi, -1.0, cmj->vbi);
    gkyl_array_accumulate(cmj->ddvbi, 1.0, vbi_target);
    gkyl_array_set(cmj->ddT, -1.0, cmj->T);
    gkyl_array_accumulate(cmj->ddT, 1.0, T_target);

    // b. Calculate  dMi^(k+1) = dn^k + ddMi^(k+1) | where dn^0 = 0
    // dm_new = dm_old + ddn;
    gkyl_array_accumulate(cmj->dn, 1.0, cmj->ddn);
    gkyl_array_accumulate(cmj->dvbi, 1.0, cmj->ddvbi);
    gkyl_array_accumulate(cmj->dT, 1.0, cmj->ddT);

    // End the iteration early if all moments converge
    if ((cmj->niter % 1) == 0){
      struct gkyl_range_iter biter;

      // Reset the maximum error
      cmj->error_n = 0; cmj->error_T = 0;
      cmj->error_vb[0] = 0; cmj->error_vb[1] = 0; cmj->error_vb[2] = 0;

      // Iterate over the grid to find the maximum error
      gkyl_range_iter_init(&biter, conf_local);
      while (gkyl_range_iter_next(&biter)){
        long midx = gkyl_range_idx(conf_local, biter.idx);
        const double *n_local = gkyl_array_cfetch(cmj->n, midx);
        const double *vbi_local = gkyl_array_cfetch(cmj->vbi, midx);
        const double *T_local = gkyl_array_cfetch(cmj->T, midx);
        const double *n_original_local = gkyl_array_cfetch(n_target, midx);
        const double *vbi_original_local = gkyl_array_cfetch(vbi_target, midx);
        const double *T_original_local = gkyl_array_cfetch(T_target, midx);
        cmj->error_n = fmax(fabs(n_local[0] - n_original_local[0]),fabs(cmj->error_n));
        cmj->error_vb[0] = fmax(fabs(vbi_local[0] - vbi_original_local[0]),fabs(cmj->error_vb[0]));
        cmj->error_T = fmax(fabs(T_local[0] - T_original_local[0]),fabs(cmj->error_T));
        if (vdim > 1)
          cmj->error_vb[1] = fmax(fabs(vbi_local[poly_order + 1] - vbi_original_local[poly_order + 1]),fabs(cmj->error_vb[1]));
        if (vdim > 2)
          cmj->error_vb[2] = fmax(fabs(vbi_local[2 * (poly_order + 1)] - vbi_original_local[2 * (poly_order + 1)]),fabs(cmj->error_vb[2]));
      }
    }

    // c. Calculate  n^(k+1) = M^k + dM^(k+1)
    // n = n_target + dm_new;
    gkyl_array_set(cmj->n, 1.0, n_target);
    gkyl_array_accumulate(cmj->n, 1.0, cmj->dn);
    gkyl_array_set(cmj->vbi, 1.0, vbi_target);
    gkyl_array_accumulate(cmj->vbi, 1.0, cmj->dvbi);
    gkyl_array_set(cmj->T, 1.0, T_target);
    gkyl_array_accumulate(cmj->T, 1.0, cmj->dT);

    // 2. Update the dist_mj using the corrected moments
    gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(cmj->proj_mj, phase_local, 
      conf_local, cmj->n, cmj->vbi, cmj->T, distf_mj);

    // 3. Correct the n moment to fix the asymptotically approximated MJ function
    gkyl_correct_mj_fix_n_stationary(cmj, distf_mj, cmj->n, cmj->vbi, phase_local, conf_local);

    cmj->niter += 1;
  }
  if ((cmj->niter < 100) && ((fabs(cmj->error_n) < tol) && (fabs(cmj->error_vb[0]) < tol) &&
    (fabs(cmj->error_vb[1]) < tol) && (fabs(cmj->error_vb[2]) < tol) && (fabs(cmj->error_T) < tol)))
    cmj->status = 0;
  else
    cmj->status = 1;

  // If the algorithm fails (density fails to converge)!
  // Project the distribution function with the basic moments and correct n
  if (cmj->status == 1){
    gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(cmj->proj_mj, phase_local, 
      conf_local, n_target, vbi_target, T_target, distf_mj);
    gkyl_correct_mj_fix_n_stationary(cmj, distf_mj, n_target, vbi_target, phase_local, conf_local);
  }
}

void 
gkyl_correct_mj_release(gkyl_correct_mj *cmj)
{

  // gkyl_mom_calc_release(cmj->ncalc);
  gkyl_dg_updater_moment_release(cmj->ncalc);
  gkyl_dg_updater_moment_release(cmj->vbicalc);
  gkyl_array_release(cmj->num_ratio);
  gkyl_array_release(cmj->vb_dot_nvb);
  gkyl_array_release(cmj->n_minus_vb_dot_nvb);
  gkyl_array_release(cmj->num_vb);
  gkyl_array_release(cmj->V_drift);
  gkyl_array_release(cmj->gamma);
  gkyl_dg_bin_op_mem_release(cmj->mem);

  gkyl_array_release(cmj->n);
  gkyl_array_release(cmj->vbi);
  gkyl_array_release(cmj->T);
  gkyl_array_release(cmj->dn);
  gkyl_array_release(cmj->dvbi);
  gkyl_array_release(cmj->dT);
  gkyl_array_release(cmj->ddn);
  gkyl_array_release(cmj->ddvbi);
  gkyl_array_release(cmj->ddT);

  gkyl_mj_moments_release(cmj->mj_moms);
  gkyl_proj_mj_on_basis_release(cmj->proj_mj);

  gkyl_free(cmj);
}