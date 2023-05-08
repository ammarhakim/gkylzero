#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_correct_mj.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_proj_on_basis.h>

#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov_sr.h>

// Temporary
#include <gkyl_array_ops_priv.h>
#include <gkyl_correct_mj.h>
#include <gkyl_array_ops.h>
#include <gkyl_mj_moments.h>
#include <gkyl_proj_mj_on_basis.h>
#include <gkyl_proj_on_basis.h>

struct gkyl_correct_mj
{
  struct gkyl_rect_grid grid;
  struct gkyl_basis conf_basis, phase_basis;

  // struct gkyl_mom_calc *m0calc;
  gkyl_dg_updater_moment *m0calc;  // moment calculator
  gkyl_dg_updater_moment *m1icalc; // moment calculator
  struct gkyl_array *num_ratio;    // number density ratio
  struct gkyl_array *num_vb;       // number density times vb
  struct gkyl_array *gamma;        // dg represented gamma of the shifted frame

  gkyl_dg_bin_op_mem *mem;          // memory for division operator
  struct gkyl_array *m0, *m1i, *m2; // memory for corrective mj scheme
  struct gkyl_array *dm0, *dm1i, *dm2;
  struct gkyl_array *ddm0, *ddm1i, *ddm2;
};

gkyl_correct_mj *
gkyl_correct_mj_new(const struct gkyl_rect_grid *grid,
                    const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
                    const struct gkyl_range *conf_range, const struct gkyl_range *vel_range,
                    long conf_local_ncells, long conf_local_ext_ncells, bool use_gpu)
{
  gkyl_correct_mj *up = gkyl_malloc(sizeof(*up));

  up->grid = *grid;
  up->conf_basis = *conf_basis;
  up->phase_basis = *phase_basis;
  int vdim = up->phase_basis.ndim - up->conf_basis.ndim;

  // updated moment calculator for sr N and N*vb moments
  up->m0calc = gkyl_dg_updater_moment_new(grid, conf_basis,
                                          phase_basis, conf_range, vel_range, GKYL_MODEL_SR, "M0", 0, 1, use_gpu);
  up->m1icalc = gkyl_dg_updater_moment_new(grid, conf_basis,
                                           phase_basis, conf_range, vel_range, GKYL_MODEL_SR, "M1i", 0, 1, use_gpu);

  up->num_ratio = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->num_vb = gkyl_array_new(GKYL_DOUBLE, vdim * conf_basis->num_basis, conf_local_ext_ncells);
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

  return up;
}

void gkyl_correct_mj_fix_m0(gkyl_correct_mj *cmj, const struct gkyl_array *p_over_gamma,
                            struct gkyl_array *fout,
                            const struct gkyl_array *m0,
                            const struct gkyl_array *m1i,
                            const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  // vdim
  int vdim = cmj->phase_basis.ndim - cmj->conf_basis.ndim;

  // compute the sr moments
  gkyl_dg_updater_moment_advance(cmj->m0calc, phase_local, conf_local,
                                 0, 0, 0,
                                 0, 0, 0,
                                 fout, cmj->num_ratio);
  gkyl_dg_updater_moment_advance(cmj->m1icalc, phase_local, conf_local,
                                 p_over_gamma, 0, 0,
                                 0, 0, 0,
                                 fout, cmj->num_vb);

  // isolate vb by dividing N*vb by N
  for (int d = 0; d < vdim; ++d)
  {
    gkyl_dg_div_op_range(cmj->mem, cmj->conf_basis, d, cmj->num_vb,
                         d, cmj->num_vb, 0, cmj->num_ratio, conf_local);
  }

  // compute number density ratio
  gkyl_dg_div_op_range(cmj->mem, cmj->conf_basis, 0, cmj->num_ratio,
                       0, m0, 0, cmj->num_ratio, conf_local);

  // calculate gamma from cmj->num_vb
  gkyl_calc_sr_vars_Gamma(&cmj->conf_basis, &cmj->phase_basis,
                          conf_local, cmj->num_vb, cmj->gamma);

  // multiply the number density ratio by gamma, to account for the frame trans.
  gkyl_dg_mul_op_range(cmj->conf_basis, 0, cmj->num_ratio,
                       0, cmj->num_ratio, 0, cmj->gamma, conf_local);

  // rescale distribution function
  gkyl_dg_mul_conf_phase_op_range(&cmj->conf_basis, &cmj->phase_basis,
                                  fout, cmj->num_ratio, fout, conf_local, phase_local);
}

void gkyl_correct_mj_fix(gkyl_correct_mj *cmj,
                         struct gkyl_array *distf_mj,
                         const struct gkyl_array *m0_corr,
                         const struct gkyl_array *m1i_corr,
                         const struct gkyl_array *m2_corr,
                         const struct gkyl_range *phase_local, const struct gkyl_range *conf_local,
                         int poly_order, const struct gkyl_range *conf_local_ext, const struct gkyl_range *velLocal,
                         const struct gkyl_basis *vel_basis, const struct gkyl_rect_grid *vel_grid)
// Add:: (1) polyorder (2) conf_local_ext (3) velLocal
// ADD ::  (1) (2) vel_basis, (3) vel_grid
// Will have to update alot of other codes dependent on this routine
{
  // vdim
  int vdim = cmj->phase_basis.ndim - cmj->conf_basis.ndim;

  struct gkyl_array *p_over_gamma;
  struct gkyl_array *gamma;
  struct gkyl_array *gamma_inv;
  p_over_gamma = gkyl_array_new(GKYL_DOUBLE, vdim * vel_basis->num_basis, velLocal->volume);
  gamma = gkyl_array_new(GKYL_DOUBLE, vel_basis->num_basis, velLocal->volume);
  gamma_inv = gkyl_array_new(GKYL_DOUBLE, vel_basis->num_basis, velLocal->volume);

  // Make GammaV2, GammaV, GammaV_inv
  gkyl_calc_sr_vars_init_p_vars(vel_grid, vel_basis, velLocal,
                                p_over_gamma, gamma, gamma_inv);

  // Diagnostic ouputs
  /* LOOKS GOOD
  struct gkyl_range_iter biter;
  gkyl_range_iter_init(&biter, velLocal);
  while (gkyl_range_iter_next(&biter))
  {
    long midx = gkyl_range_idx(velLocal, biter.idx);
    const double *p_over_gamma_local = gkyl_array_cfetch(p_over_gamma, midx);
    const double *gamma_local = gkyl_array_cfetch(gamma, midx);
    const double *gamma_inv_local = gkyl_array_cfetch(gamma_inv, midx);
    printf("\n------- p_over_gamma ------\n");
    printf("p_over_gamma: %g\n", p_over_gamma_local[0]);
    printf("\n------- gamma ------\n");
    printf("gamma: %g\n", gamma_local[0]);
    printf("\n------- gamma_inv ------\n");
    printf("gamma_inv: %g\n", gamma_inv_local[0]);
  }
  */

  // Copy the intial moments for m*_corr -> m*
  gkyl_array_clear(cmj->m0, 0.0);
  gkyl_array_accumulate(cmj->m0, 1.0, m0_corr);
  gkyl_array_clear(cmj->m1i, 0.0);
  gkyl_array_accumulate(cmj->m1i, 1.0, m1i_corr);
  gkyl_array_clear(cmj->m2, 0.0);
  gkyl_array_accumulate(cmj->m2, 1.0, m2_corr);

  // Make mj moms stucture
  gkyl_mj_moments *mj_moms = gkyl_mj_moments_new(&cmj->grid, &cmj->conf_basis, &cmj->phase_basis, conf_local, velLocal, conf_local->volume, conf_local_ext->volume, false);

  // Make the projection routine for the new MJ distributions
  gkyl_proj_mj_on_basis *proj_mj = gkyl_proj_mj_on_basis_new(&cmj->grid,
                                                             &cmj->conf_basis, &cmj->phase_basis, poly_order + 1);

  // 0. Project the MJ with the intially correct moments
  gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(proj_mj, phase_local, conf_local, m0_corr, m1i_corr, m2_corr, distf_mj);

  // Iteration loop, 100 iterations is usually more than sufficient (for all vdim)
  for (int i = 0; i < 300; ++i)
  {

    // 1. Correct the M0 moment to fix the asymptotically approximated MJ function
    // gkyl_correct_mj_fix_m0(cmj, p_over_gamma, distf_mj, cmj->m0, cmj->m1i, phase_local, conf_local);

    // 2. Calculate the new moments
    // calculate the moments of the dist (n, vb, T -> m0, m1i, m2)
    gkyl_mj_moments_advance(mj_moms, p_over_gamma, gamma, gamma_inv, distf_mj, cmj->m0, cmj->m1i, cmj->m2, phase_local, conf_local);

    // a. Calculate  ddMi^(k+1) =  Mi_corr - Mi_new
    // ddm0 = m0_corr - m0;
    //  Compute out = out + a*inp. Returns out.
    gkyl_array_clear(cmj->ddm0, 0.0);
    gkyl_array_accumulate(cmj->ddm0, -1.0, cmj->m0);
    gkyl_array_accumulate(cmj->ddm0, 1.0, m0_corr);
    gkyl_array_clear(cmj->ddm1i, 0.0);
    gkyl_array_accumulate(cmj->ddm1i, -1.0, cmj->m1i);
    gkyl_array_accumulate(cmj->ddm1i, 1.0, m1i_corr);
    gkyl_array_clear(cmj->ddm2, 0.0);
    gkyl_array_accumulate(cmj->ddm2, -1.0, cmj->m2);
    gkyl_array_accumulate(cmj->ddm2, 1.0, m2_corr);

    // b. Calculate  dMi^(k+1) = dM0^k + ddMi^(k+1) | where dM0^0 = 0
    // dm_new = dm_old + ddm0;
    gkyl_array_accumulate(cmj->dm0, 1.0, cmj->ddm0);
    gkyl_array_accumulate(cmj->dm1i, 1.0, cmj->ddm1i);
    gkyl_array_accumulate(cmj->dm2, 1.0, cmj->ddm2);

    // c. Calculate  M0^(k+1) = m0 + dm^(k+1)
    // m0 = m0_corr + dm_new;
    gkyl_array_clear(cmj->m0, 0.0);
    gkyl_array_accumulate(cmj->m0, 1.0, m0_corr);
    gkyl_array_accumulate(cmj->m0, 1.0, cmj->dm0);
    gkyl_array_clear(cmj->m1i, 0.0);
    gkyl_array_accumulate(cmj->m1i, 1.0, m1i_corr);
    gkyl_array_accumulate(cmj->m1i, 1.0, cmj->dm1i);
    gkyl_array_clear(cmj->m2, 0.0);
    gkyl_array_accumulate(cmj->m2, 1.0, m2_corr);
    gkyl_array_accumulate(cmj->m2, 1.0, cmj->dm2);

    // 3. Update the dist_mj using the corrected moments
    gkyl_proj_mj_on_basis_fluid_stationary_frame_mom(proj_mj, phase_local, conf_local, cmj->m0, cmj->m1i, cmj->m2, distf_mj);

  } // Iteration loop

  // Release the memory
  gkyl_mj_moments_release(mj_moms);
  gkyl_array_release(p_over_gamma);
  gkyl_array_release(gamma);
  gkyl_array_release(gamma_inv);
}

void gkyl_correct_mj_release(gkyl_correct_mj *cmj)
{

  // gkyl_mom_calc_release(cmj->m0calc);
  gkyl_dg_updater_moment_release(cmj->m0calc);
  gkyl_dg_updater_moment_release(cmj->m1icalc);
  gkyl_array_release(cmj->num_ratio);
  gkyl_array_release(cmj->num_vb);
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

  gkyl_free(cmj);
}
