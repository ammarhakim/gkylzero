#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_MJ_moments.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_proj_on_basis.h>

#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov_sr.h>


struct gkyl_MJ_moments {
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


gkyl_MJ_moments *
gkyl_MJ_moments_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range,
  long conf_local_ncells, long conf_local_ext_ncells, bool use_gpu)
{
  gkyl_MJ_moments *up = gkyl_malloc(sizeof(*up));

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

void gkyl_MJ_moments_fix(gkyl_MJ_moments *cMJ, const struct gkyl_array *p_over_gamma,
  const struct gkyl_array *gamma, const struct gkyl_array *gamma_inv,
  struct gkyl_array *fout,
  struct gkyl_array *m0,
  struct gkyl_array *m1i,
  struct gkyl_array *m2,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  // vdim
  int vdim = cMJ->phase_basis.ndim - cMJ->conf_basis.ndim;

  // compute the sr moments, lab frame <N> and <Nvb>
  gkyl_dg_updater_moment_advance(cMJ->m0calc, phase_local, conf_local,
    0, 0, 0,
    0, 0, 0,
    fout, cMJ->num_ratio);
  gkyl_dg_updater_moment_advance(cMJ->m1icalc, phase_local, conf_local,
    p_over_gamma, 0, 0,
    0, 0, 0,
    fout, cMJ->num_vb);

  // (vb = <Nvb>/<N>) isolate vb by dividing <N*vb> by <N>
  gkyl_dg_div_op_range(cMJ->mem, cMJ->conf_basis, 0, cMJ->V_drift,
    0, cMJ->num_vb, 0, cMJ->num_ratio, conf_local);

  // (Gamma = 1/sqrt(1-vb^2)) calculate cMJ->gamma from cMJ->num_vb
  gkyl_calc_sr_vars_Gamma(&cMJ->conf_basis, &cMJ->phase_basis,
      conf_local, cMJ->V_drift, cMJ->gamma);

  // (Gamma^2 = 1/(1-vb^2)) calculate
  gkyl_calc_sr_vars_Gamma2(&cMJ->conf_basis, &cMJ->phase_basis,
      conf_local, cMJ->V_drift, cMJ->GammaV2);

  // (Gamma_inv = sqrt(1-vb^2))
  gkyl_calc_sr_vars_Gamma_inv(&cMJ->conf_basis, &cMJ->phase_basis,
      conf_local, cMJ->V_drift, cMJ->Gamma_inv);

  // compute the pressure moment (stationary frame)
  gkyl_dg_updater_moment_advance(cMJ->m2calc, phase_local, conf_local,
        p_over_gamma, gamma,
        gamma_inv, cMJ->V_drift,
        cMJ->GammaV2, cMJ->Gamma_inv,
        fout, cMJ->pressure);

  // (n = N/gamma) divide the number density by gamma, to account for the frame trans.
  gkyl_dg_div_op_range(cMJ->mem, cMJ->conf_basis, 0, cMJ->num_ratio,
      0, cMJ->num_ratio, 0, cMJ->gamma, conf_local);

  // (T = P/n) Calculate from the restframe Pressure, the rest frame temperature
  gkyl_dg_div_op_range(cMJ->mem, cMJ->conf_basis, 0, cMJ->temperature,
    0, cMJ->pressure, 0, cMJ->num_ratio, conf_local);

    // Temporarily print out the n, vb, T: Recall the quantities are *1/sqrt(2)
    // due to them being the coefficients
  struct gkyl_range_iter biter;
  gkyl_range_iter_init(&biter, conf_local);
  while (gkyl_range_iter_next(&biter)) {
      long midx = gkyl_range_idx(conf_local, biter.idx);
      const double *num = gkyl_array_cfetch(cMJ->num_ratio, midx);
      const double *vb = gkyl_array_cfetch(cMJ->V_drift, midx);
      const double *pressure = gkyl_array_cfetch(cMJ->pressure, midx);
      const double *temperature = gkyl_array_cfetch(cMJ->temperature, midx);
      printf("\n----------- Ouptuts Start ---------\n");
      printf("num: %1.16g\n",num[0]);
      printf("vb : %1.16g\n",vb[0]);
      printf("T  : %1.16g\n",temperature[0]);
      printf("P  : %1.16g\n",pressure[0]);
      printf("\n----------- Ouptuts End ---------\n");
  }


}

void
gkyl_MJ_moments_release(gkyl_MJ_moments* cMJ)
{

  //gkyl_mom_calc_release(cMJ->m0calc);
  gkyl_dg_updater_moment_release(cMJ->m0calc);
  gkyl_dg_updater_moment_release(cMJ->m1icalc);
  gkyl_dg_updater_moment_release(cMJ->m2calc);
  gkyl_array_release(cMJ->num_ratio);
  gkyl_array_release(cMJ->num_vb);
  gkyl_array_release(cMJ->V_drift);
  gkyl_array_release(cMJ->Gamma_inv);
  gkyl_array_release(cMJ->gamma);
  gkyl_array_release(cMJ->GammaV2);
  gkyl_array_release(cMJ->pressure);
  gkyl_array_release(cMJ->temperature);
  gkyl_dg_bin_op_mem_release(cMJ->mem);

  gkyl_free(cMJ);
}
