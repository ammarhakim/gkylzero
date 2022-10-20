#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_correct_MJ.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_proj_on_basis.h>


struct gkyl_correct_MJ {
  struct gkyl_rect_grid grid;
  struct gkyl_basis conf_basis, phase_basis;

  gkyl_dg_updater_moment *m0calc; // moment calculator
  gkyl_dg_updater_moment *m1icalc; // moment calculator
  struct gkyl_array *num_ratio; // number density ratio
  struct gkyl_array *num_vb; // number density times vb
  struct gkyl_array *gamma; // dg represented gamma of the shifted frame

  gkyl_dg_bin_op_mem *mem; // memory for division operator
};

// // (Added) projects vb in DG represenation into the function gamma
// // up gives the basis, fun_at_ords in our function (gamma) to be evaluated, f is the array we are updating
// static void
// proj_on_basis(const gkyl_correct_MJ *up, const struct gkyl_array *fun_at_ords, double* f)
// {
//   int num_basis = up->num_phase_basis;
//   int tot_quad = up->tot_quad;
//
//   const double* GKYL_RESTRICT weights = up->weights->data;
//   const double* GKYL_RESTRICT basis_at_ords = up->basis_at_ords->data;
//   const double* GKYL_RESTRICT func_at_ords = fun_at_ords->data;
//
//   for (int k=0; k<num_basis; ++k) f[k] = 0.0;
//
//   for (int imu=0; imu<tot_quad; ++imu) {
//     double tmp = weights[imu]*func_at_ords[imu];
//     for (int k=0; k<num_basis; ++k)
//       f[k] += tmp*basis_at_ords[k+num_basis*imu];
//   }
// }


void gamma_simple_test(double t, const double *xn, double* restrict fout, void *ctx)
{
  double x = xn[0];
  double weight = 0.5;
  fout[0] = weight*sqrt(1.0 - 0.5*0.5);
}



gkyl_correct_MJ *
gkyl_correct_MJ_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range,
  long conf_local_ncells, long conf_local_ext_ncells, bool use_gpu)
{
  gkyl_correct_MJ *up = gkyl_malloc(sizeof(*up));

  up->grid = *grid;
  up->conf_basis = *conf_basis;
  up->phase_basis = *phase_basis;
  int vdim = up->phase_basis.ndim - up->conf_basis.ndim;

  // old moment calculator routine
  //struct gkyl_mom_type *m0 = gkyl_mom_vlasov_sr_new(conf_basis, phase_basis,conf_range, vel_range, "M0",use_gpu);
  //up->m0calc = gkyl_mom_calc_new(grid, m0, use_gpu);
  //gkyl_mom_type_release(m0);

  //struct gkyl_mom_type *m1i = gkyl_mom_vlasov_sr_new(conf_basis, phase_basis,conf_range, vel_range, "M1i",use_gpu);
  //up->m1icalc = gkyl_mom_calc_new(grid, m1i, use_gpu);
  //gkyl_mom_type_release(m1i);

  // updated moment calculator for sr
  up->m0calc = gkyl_dg_updater_moment_new(grid, conf_basis,
    phase_basis, conf_range, vel_range, GKYL_MODEL_SR, "M0", 0, 1, use_gpu);
  up->m1icalc = gkyl_dg_updater_moment_new(grid, conf_basis,
      phase_basis, conf_range, vel_range, GKYL_MODEL_SR, "M1i", 0, 1, use_gpu);


  up->num_ratio = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->num_vb = gkyl_array_new(GKYL_DOUBLE, vdim*conf_basis->num_basis, conf_local_ext_ncells);
  up->gamma = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->mem = gkyl_dg_bin_op_mem_new(conf_local_ncells, conf_basis->num_basis);
  return up;
}

void gkyl_correct_MJ_fix(gkyl_correct_MJ *cMJ, const struct gkyl_array *p_over_gamma,
  struct gkyl_array *fout,
  const struct gkyl_array *m0,
  const struct gkyl_array *m1i,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  // vdim
  int vdim = cMJ->phase_basis.ndim - cMJ->conf_basis.ndim;

  // old moment calculator routine
  // calculate number density and beam velocity
  //gkyl_mom_calc_advance(cMJ->m0calc, phase_local, conf_local, fout, cMJ->num_ratio);
  //gkyl_mom_calc_advance(cMJ->m1icalc, phase_local, conf_local, fout, cMJ->num_vb);
  gkyl_dg_updater_moment_advance(cMJ->m0calc, phase_local, conf_local,
    0, 0, 0,
    0, 0, 0,
    fout, cMJ->num_ratio);
  gkyl_dg_updater_moment_advance(cMJ->m1icalc, phase_local, conf_local,
    p_over_gamma, 0, 0,
    0, 0, 0,
    fout, cMJ->num_vb);

    struct gkyl_range_iter witer;
    gkyl_range_iter_init(&witer, conf_local);
    while (gkyl_range_iter_next(&witer)) {
        long midx = gkyl_range_idx(conf_local, witer.idx);
        const double *gamma_val = gkyl_array_cfetch(cMJ->gamma, midx);
        const double *num = gkyl_array_cfetch(cMJ->num_ratio, midx);
        const double *vb = gkyl_array_cfetch(cMJ->num_vb, midx);
        printf("Gamma: %g\n",gamma_val[0]);
        printf("num: %g\n",num[0]);
        printf("vb: %g\n",vb[0]);
    }

  // (Added) isolate vb by dividing N*vb by N
  gkyl_dg_div_op_range(cMJ->mem, cMJ->conf_basis, 0, cMJ->num_vb,
    0, cMJ->num_vb, 0, cMJ->num_ratio, conf_local);

  // compute number density ratio
  gkyl_dg_div_op_range(cMJ->mem, cMJ->conf_basis, 0, cMJ->num_ratio,
    0, m0, 0, cMJ->num_ratio, conf_local);

  // Test when gamma is a constant value given: vb = 0.5
  //gkyl_proj_on_basis *gamma = gkyl_proj_on_basis_new(&cMJ->grid, &cMJ->conf_basis,
  //  cMJ->conf_basis.poly_order+1, 1, gamma_simple_test, NULL);
  //gkyl_proj_on_basis_advance(gamma, 0.0, conf_local, cMJ->gamma);
  //gkyl_proj_on_basis_release(gamma);

  // calculate cMJ->gamma from cMJ->num_vb
  gkyl_calc_sr_vars_Gamma(&cMJ->conf_basis, &cMJ->phase_basis,
      conf_local, cMJ->num_vb, cMJ->gamma);


  // multiply the number density ratio by gamma, to account for the frame trans.
  gkyl_dg_mul_op_range(cMJ->conf_basis, 0, cMJ->num_ratio,
      0, cMJ->num_ratio, 0, cMJ->gamma, conf_local);

  //Try with div op instead, gamma -> gamma_inv
  //gkyl_calc_sr_vars_Gamma_inv(&cMJ->conf_basis, &cMJ->phase_basis,
  //  conf_local, cMJ->num_vb, cMJ->gamma);
  //gkyl_dg_div_op_range(cMJ->mem, cMJ->conf_basis, 0, cMJ->num_ratio,
  //  0, cMJ->num_ratio, 0, cMJ->gamma, conf_local);


  // rescale distribution function
  gkyl_dg_mul_conf_phase_op_range(&cMJ->conf_basis, &cMJ->phase_basis,
    fout, cMJ->num_ratio, fout, conf_local, phase_local);
}

void
gkyl_correct_MJ_release(gkyl_correct_MJ* cMJ)
{
  gkyl_dg_updater_moment_release(cMJ->m0calc);
  gkyl_dg_updater_moment_release(cMJ->m1icalc);
  gkyl_array_release(cMJ->num_ratio);
  gkyl_array_release(cMJ->num_vb);
  gkyl_array_release(cMJ->gamma);
  gkyl_dg_bin_op_mem_release(cMJ->mem);

  gkyl_free(cMJ);
}
