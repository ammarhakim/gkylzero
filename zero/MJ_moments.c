#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_mj_moments.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_proj_on_basis.h>

#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov_sr.h>
#include <gkyl_array_ops_priv.h>


struct gkyl_mj_moments {
  struct gkyl_rect_grid grid;
  struct gkyl_basis conf_basis, phase_basis;


  //struct gkyl_mom_calc *m0calc;
  gkyl_dg_updater_moment *m0calc; // moment calculator
  gkyl_dg_updater_moment *m1icalc; // moment calculator
  gkyl_dg_updater_moment *m2calc; // moment calculator
  struct gkyl_array *num_ratio; // number density ratio
  struct gkyl_array *num_vb; // number density times vb
  struct gkyl_array *V_drift; // number density times vb
  struct gkyl_array *gamma; // dg represented gamma of the shifted frame
  struct gkyl_array *GammaV2;
  struct gkyl_array *Gamma_inv;
  struct gkyl_array *pressure;
  struct gkyl_array *temperature;

  gkyl_dg_bin_op_mem *mem; // memory for division operator
};


gkyl_mj_moments *
gkyl_mj_moments_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range,
  long conf_local_ncells, long conf_local_ext_ncells, bool use_gpu)
{
  gkyl_mj_moments *up = gkyl_malloc(sizeof(*up));

  up->grid = *grid;
  up->conf_basis = *conf_basis;
  up->phase_basis = *phase_basis;
  int vdim = up->phase_basis.ndim - up->conf_basis.ndim;

  // updated moment calculator for sr N and N*vb moments
  up->m0calc = gkyl_dg_updater_moment_new(grid, conf_basis,
      phase_basis, conf_range, vel_range, GKYL_MODEL_SR, "M0", 0, 1, use_gpu);
  up->m1icalc = gkyl_dg_updater_moment_new(grid, conf_basis,
      phase_basis, conf_range, vel_range, GKYL_MODEL_SR, "M1i", 0, 1, use_gpu);
  up->m2calc = gkyl_dg_updater_moment_new(grid, conf_basis,
      phase_basis, conf_range, vel_range, GKYL_MODEL_SR, "Pressure", 0, 1, use_gpu);


  up->num_ratio = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->num_vb = gkyl_array_new(GKYL_DOUBLE, vdim*conf_basis->num_basis, conf_local_ext_ncells);
  up->V_drift = gkyl_array_new(GKYL_DOUBLE, vdim*conf_basis->num_basis, conf_local_ext_ncells);
  up->Gamma_inv = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->gamma = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->GammaV2 = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->pressure = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->temperature = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->mem = gkyl_dg_bin_op_mem_new(conf_local_ncells, conf_basis->num_basis);
  return up;
}

void gkyl_mj_moments_advance(gkyl_mj_moments *cmj, const struct gkyl_array *p_over_gamma,
  const struct gkyl_array *gamma, const struct gkyl_array *gamma_inv,
  struct gkyl_array *fout,
  struct gkyl_array *m0,
  struct gkyl_array *m1i,
  struct gkyl_array *m2,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  // vdim
  int vdim = cmj->phase_basis.ndim - cmj->conf_basis.ndim;

  // compute the sr moments, lab frame <N> and <Nvb>
  gkyl_dg_updater_moment_advance(cmj->m0calc, phase_local, conf_local,
    0, 0, 0,
    0, 0, 0,
    fout, cmj->num_ratio);
  gkyl_dg_updater_moment_advance(cmj->m1icalc, phase_local, conf_local,
    p_over_gamma, 0, 0,
    0, 0, 0,
    fout, cmj->num_vb);

  // (vb = <Nvb>/<N>) isolate vb by dividing <N*vb> by <N>
  gkyl_dg_div_op_range(cmj->mem, cmj->conf_basis, 0, cmj->V_drift,
    0, cmj->num_vb, 0, cmj->num_ratio, conf_local);

  // (Gamma = 1/sqrt(1-vb^2)) calculate cmj->gamma from cmj->num_vb
  gkyl_calc_sr_vars_Gamma(&cmj->conf_basis, &cmj->phase_basis,
      conf_local, cmj->V_drift, cmj->gamma);

  // (Gamma^2 = 1/(1-vb^2)) calculate
  gkyl_calc_sr_vars_Gamma2(&cmj->conf_basis, &cmj->phase_basis,
      conf_local, cmj->V_drift, cmj->GammaV2);

  // (Gamma_inv = sqrt(1-vb^2))
  gkyl_calc_sr_vars_Gamma_inv(&cmj->conf_basis, &cmj->phase_basis,
      conf_local, cmj->V_drift, cmj->Gamma_inv);

  // compute the pressure moment (stationary frame)
  gkyl_dg_updater_moment_advance(cmj->m2calc, phase_local, conf_local,
        p_over_gamma, gamma,
        gamma_inv, cmj->V_drift,
        cmj->GammaV2, cmj->Gamma_inv,
        fout, cmj->pressure);

  // (n = N/gamma) divide the number density by gamma, to account for the frame trans.
  gkyl_dg_div_op_range(cmj->mem, cmj->conf_basis, 0, cmj->num_ratio,
      0, cmj->num_ratio, 0, cmj->gamma, conf_local);

  // (T = P/n) Calculate from the restframe Pressure, the rest frame temperature
  gkyl_dg_div_op_range(cmj->mem, cmj->conf_basis, 0, cmj->temperature,
    0, cmj->pressure, 0, cmj->num_ratio, conf_local);

}

void
gkyl_mj_moments_release(gkyl_mj_moments* cmj)
{

  //gkyl_mom_calc_release(cmj->m0calc);
  gkyl_dg_updater_moment_release(cmj->m0calc);
  gkyl_dg_updater_moment_release(cmj->m1icalc);
  gkyl_dg_updater_moment_release(cmj->m2calc);
  gkyl_array_release(cmj->num_ratio);
  gkyl_array_release(cmj->num_vb);
  gkyl_array_release(cmj->V_drift);
  gkyl_array_release(cmj->Gamma_inv);
  gkyl_array_release(cmj->gamma);
  gkyl_array_release(cmj->GammaV2);
  gkyl_array_release(cmj->pressure);
  gkyl_array_release(cmj->temperature);
  gkyl_dg_bin_op_mem_release(cmj->mem);

  gkyl_free(cmj);
}
