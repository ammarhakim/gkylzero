#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_correct_mj.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_mj_moments.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov_sr.h>
#include <gkyl_proj_mj_on_basis.h>
#include <gkyl_proj_on_basis.h>

struct gkyl_correct_mj
{
  struct gkyl_rect_grid grid;
  struct gkyl_basis conf_basis, phase_basis;

  struct gkyl_dg_updater_moment *m0calc;  
  struct gkyl_dg_updater_moment *m1icalc; 
  struct gkyl_array *num_ratio;  
  struct gkyl_array *num_vb;   
  struct gkyl_array *V_drift;  
  struct gkyl_array *gamma;   
  struct gkyl_array *vb_dot_nvb; 
  struct gkyl_array *n_minus_vb_dot_nvb; 

  struct gkyl_dg_bin_op_mem *mem;     
  struct gkyl_array *m0, *m1i, *m2;
  struct gkyl_array *dm0, *dm1i, *dm2;
  struct gkyl_array *ddm0, *ddm1i, *ddm2;

  struct gkyl_mj_moments *mj_moms;
  struct gkyl_proj_mj_on_basis *proj_mj;
};

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
  up->m0calc = gkyl_dg_updater_moment_new(grid, conf_basis,
    phase_basis, conf_range, vel_range, GKYL_MODEL_SR, &sr_inp, "M0", 0, 1, use_gpu);
  up->m1icalc = gkyl_dg_updater_moment_new(grid, conf_basis,
    phase_basis, conf_range, vel_range, GKYL_MODEL_SR, &sr_inp, "M1i", 0, 1, use_gpu);

  up->num_ratio = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->num_vb = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->V_drift = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->vb_dot_nvb = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->n_minus_vb_dot_nvb = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->gamma = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->mem = gkyl_dg_bin_op_mem_new(conf_local_ncells, conf_basis->num_basis);

  // Moment memory
  up->m0 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->m1i = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->m2 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);

  // Create a copy for differences (d) and differences of differences (dd)
  up->dm0 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->dm1i = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->dm2 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);

  up->ddm0 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->ddm1i = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
  up->ddm2 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);

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
gkyl_correct_mj_fix_m0(gkyl_correct_mj *cmj, 
  struct gkyl_array *fout,
  const struct gkyl_array *m0, const struct gkyl_array *m1i,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  // vdim
  int vdim = cmj->phase_basis.ndim - cmj->conf_basis.ndim;

  // compute the sr moments
  gkyl_dg_updater_moment_advance(cmj->m0calc, phase_local, conf_local, fout, cmj->num_ratio);
  gkyl_dg_updater_moment_advance(cmj->m1icalc, phase_local, conf_local, fout, cmj->num_vb);

  // isolate vb by dividing N*vb by N
  for (int d = 0; d < vdim; ++d)
    gkyl_dg_div_op_range(cmj->mem, cmj->conf_basis, d, cmj->V_drift,
    d, cmj->num_vb, 0, cmj->num_ratio, conf_local);

  // calculate gamma from cmj->num_vb
  gkyl_calc_sr_vars_Gamma(&cmj->conf_basis, &cmj->phase_basis,
    conf_local, cmj->V_drift, cmj->gamma);

  // calculate the stationary density, n0 = gamma*(N - vb dot Nvb)
  gkyl_array_clear_range(cmj->vb_dot_nvb, 0.0, conf_local);
  gkyl_array_clear_range(cmj->n_minus_vb_dot_nvb, 0.0, conf_local);
  gkyl_array_accumulate_range(cmj->n_minus_vb_dot_nvb, 1.0, cmj->num_ratio, conf_local);
  gkyl_dg_dot_product_op_range(cmj->conf_basis,cmj->vb_dot_nvb,cmj->V_drift,cmj->num_vb, conf_local);
  gkyl_array_accumulate_range(cmj->n_minus_vb_dot_nvb, -1.0, cmj->vb_dot_nvb, conf_local);
  gkyl_dg_mul_op_range(cmj->conf_basis,0,cmj->num_ratio,0,cmj->gamma,0,cmj->n_minus_vb_dot_nvb, conf_local);

  // compute number density ratio: num_ratio = n/n0
  gkyl_dg_div_op_range(cmj->mem, cmj->conf_basis, 0, cmj->num_ratio,
    0, m0, 0, cmj->num_ratio, conf_local);

  // rescale distribution function
  gkyl_dg_mul_conf_phase_op_range(&cmj->conf_basis, &cmj->phase_basis,
    fout, cmj->num_ratio, fout, conf_local, phase_local);

  // Hand back rescaled m0:
  gkyl_array_clear_range(cmj->m0, 0.0, conf_local);
  gkyl_array_accumulate_range(cmj->m0, 1.0, cmj->num_ratio, conf_local);
}

void 
gkyl_correct_mj_fix(gkyl_correct_mj *cmj,
  struct gkyl_array *distf_mj,
  const struct gkyl_array *m0_corr, const struct gkyl_array *m1i_corr, const struct gkyl_array *m2_corr,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local,
  int poly_order)
{
  int vdim = cmj->phase_basis.ndim - cmj->conf_basis.ndim;

  // Copy the intial moments for m*_corr -> m*
  gkyl_array_clear_range(cmj->m0, 0.0, conf_local);
  gkyl_array_accumulate_range(cmj->m0, 1.0, m0_corr, conf_local);
  gkyl_array_clear_range(cmj->m1i, 0.0, conf_local);
  gkyl_array_accumulate_range(cmj->m1i, 1.0, m1i_corr, conf_local);
  gkyl_array_clear_range(cmj->m2, 0.0, conf_local);
  gkyl_array_accumulate_range(cmj->m2, 1.0, m2_corr, conf_local);

  // 0. Project the MJ with the intially correct moments
  gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(cmj->proj_mj, phase_local, 
    conf_local, m0_corr, m1i_corr, m2_corr, distf_mj);

  // tolerance of the iterative scheme
  double tol = 1e-12;
  int i = 0;
  double error_n = 1.0;
  double error_vbx = 1.0;
  double error_vby = 0.0;
  double error_vbz = 0.0;
  double error_T = 1.0;

  // Iteration loop, 100 iterations is usually sufficient (for all vdim) for machine precision moments
  while ((i < 100) && ((fabs(error_n) > tol) || (fabs(error_vbx) > tol) ||
    (fabs(error_vby) > tol) || (fabs(error_vbz) > tol) || (fabs(error_T) > tol)))
  {

    // 1. Calculate the new moments
    // calculate the moments of the dist (n, vb, T -> m0, m1i, m2)
    gkyl_mj_moments_advance(cmj->mj_moms, distf_mj, cmj->m0, cmj->m1i, cmj->m2, phase_local, conf_local);

    // a. Calculate  ddMi^(k+1) =  Mi_corr - Mi_new
    // ddm0 = m0_corr - m0;
    //  Compute out = out + a*inp. Returns out.
    gkyl_array_clear_range(cmj->ddm0, 0.0, conf_local);
    gkyl_array_accumulate_range(cmj->ddm0, -1.0, cmj->m0, conf_local);
    gkyl_array_accumulate_range(cmj->ddm0, 1.0, m0_corr, conf_local);
    gkyl_array_clear_range(cmj->ddm1i, 0.0, conf_local);
    gkyl_array_accumulate_range(cmj->ddm1i, -1.0, cmj->m1i, conf_local);
    gkyl_array_accumulate_range(cmj->ddm1i, 1.0, m1i_corr, conf_local);
    gkyl_array_clear_range(cmj->ddm2, 0.0, conf_local);
    gkyl_array_accumulate_range(cmj->ddm2, -1.0, cmj->m2, conf_local);
    gkyl_array_accumulate_range(cmj->ddm2, 1.0, m2_corr, conf_local);

    // b. Calculate  dMi^(k+1) = dM0^k + ddMi^(k+1) | where dM0^0 = 0
    // dm_new = dm_old + ddm0;
    gkyl_array_accumulate_range(cmj->dm0, 1.0, cmj->ddm0, conf_local);
    gkyl_array_accumulate_range(cmj->dm1i, 1.0, cmj->ddm1i, conf_local);
    gkyl_array_accumulate_range(cmj->dm2, 1.0, cmj->ddm2, conf_local);

    // End the iteration early if all moments converge
    if ((i % 5) == 0){
      struct gkyl_range_iter biter;
      gkyl_range_iter_init(&biter, conf_local);
      while (gkyl_range_iter_next(&biter)){
        long midx = gkyl_range_idx(conf_local, biter.idx);
        const double *m0_local = gkyl_array_cfetch(cmj->m0, midx);
        const double *m1i_local = gkyl_array_cfetch(cmj->m1i, midx);
        const double *m2_local = gkyl_array_cfetch(cmj->m2, midx);
        const double *m0_original_local = gkyl_array_cfetch(m0_corr, midx);
        const double *m1i_original_local = gkyl_array_cfetch(m1i_corr, midx);
        const double *m2_original_local = gkyl_array_cfetch(m2_corr, midx);
        error_n = m0_local[0] - m0_original_local[0];
        error_vbx = m1i_local[0] - m1i_original_local[0];
        error_T = m2_local[0] - m2_original_local[0];
        if (vdim > 1)
          error_vby = m1i_local[poly_order + 1] - m1i_original_local[poly_order + 1];
        if (vdim > 2)
          error_vbz = m1i_local[2 * (poly_order + 1)] - m1i_original_local[2 * (poly_order + 1)];
      }
    }

    // c. Calculate  M0^(k+1) = M^k + dM^(k+1)
    // m0 = m0_corr + dm_new;
    gkyl_array_clear_range(cmj->m0, 0.0, conf_local);
    gkyl_array_accumulate_range(cmj->m0, 1.0, m0_corr, conf_local);
    gkyl_array_accumulate_range(cmj->m0, 1.0, cmj->dm0, conf_local);
    gkyl_array_clear_range(cmj->m1i, 0.0, conf_local);
    gkyl_array_accumulate_range(cmj->m1i, 1.0, m1i_corr, conf_local);
    gkyl_array_accumulate_range(cmj->m1i, 1.0, cmj->dm1i, conf_local);
    gkyl_array_clear_range(cmj->m2, 0.0, conf_local);
    gkyl_array_accumulate_range(cmj->m2, 1.0, m2_corr, conf_local);
    gkyl_array_accumulate_range(cmj->m2, 1.0, cmj->dm2, conf_local);

    // 2. Update the dist_mj using the corrected moments
    gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(cmj->proj_mj, phase_local, 
      conf_local, cmj->m0, cmj->m1i, cmj->m2, distf_mj);

    // 3. Correct the M0 moment to fix the asymptotically approximated MJ function
    gkyl_correct_mj_fix_m0(cmj, distf_mj, cmj->m0, cmj->m1i, phase_local, conf_local);

    i += 1;
  }

  // If the algorithm fails (density fails to converge)!
  // Project the distribution function with the basic moments and correct m0
  gkyl_mj_moments_advance(cmj->mj_moms, distf_mj, cmj->m0, cmj->m1i, cmj->m2, phase_local, conf_local);
  double diff = 0.0;
  struct gkyl_range_iter biter;
  gkyl_range_iter_init(&biter, conf_local);
  while (gkyl_range_iter_next(&biter)){
    long midx = gkyl_range_idx(conf_local, biter.idx);
    const double *m0_corr_local = gkyl_array_cfetch(m0_corr, midx);
    const double *m0_local = gkyl_array_cfetch(cmj->m0, midx);
    if (diff < (m0_local[0] - m0_corr_local[0]))
      diff = fabs(m0_local[0] - m0_corr_local[0]);
  }
  if (fabs(diff) > 1e-10){
    gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(cmj->proj_mj, phase_local, 
      conf_local, m0_corr, m1i_corr, m2_corr, distf_mj);
    gkyl_correct_mj_fix_m0(cmj, distf_mj, cmj->m0, cmj->m1i, phase_local, conf_local);
  }
}

void 
gkyl_correct_mj_release(gkyl_correct_mj *cmj)
{

  // gkyl_mom_calc_release(cmj->m0calc);
  gkyl_dg_updater_moment_release(cmj->m0calc);
  gkyl_dg_updater_moment_release(cmj->m1icalc);
  gkyl_array_release(cmj->num_ratio);
  gkyl_array_release(cmj->vb_dot_nvb);
  gkyl_array_release(cmj->n_minus_vb_dot_nvb);
  gkyl_array_release(cmj->num_vb);
  gkyl_array_release(cmj->V_drift);
  gkyl_array_release(cmj->gamma);
  gkyl_dg_bin_op_mem_release(cmj->mem);

  gkyl_array_release(cmj->m0);
  gkyl_array_release(cmj->m1i);
  gkyl_array_release(cmj->m2);
  gkyl_array_release(cmj->dm0);
  gkyl_array_release(cmj->dm1i);
  gkyl_array_release(cmj->dm2);
  gkyl_array_release(cmj->ddm0);
  gkyl_array_release(cmj->ddm1i);
  gkyl_array_release(cmj->ddm2);

  gkyl_mj_moments_release(cmj->mj_moms);
  gkyl_proj_mj_on_basis_release(cmj->proj_mj);

  gkyl_free(cmj);
}