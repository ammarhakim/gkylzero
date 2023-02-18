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


struct gkyl_correct_mj {
  struct gkyl_rect_grid grid;
  struct gkyl_basis conf_basis, phase_basis;


  //struct gkyl_mom_calc *m0calc;
  gkyl_dg_updater_moment *m0calc; // moment calculator
  gkyl_dg_updater_moment *m1icalc; // moment calculator
  struct gkyl_array *num_ratio; // number density ratio
  struct gkyl_array *num_vb; // number density times vb
  struct gkyl_array *gamma; // dg represented gamma of the shifted frame

  gkyl_dg_bin_op_mem *mem; // memory for division operator
};


gkyl_correct_mj *
gkyl_correct_mj_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range,
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
  up->num_vb = gkyl_array_new(GKYL_DOUBLE, vdim*conf_basis->num_basis, conf_local_ext_ncells);
  up->gamma = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->mem = gkyl_dg_bin_op_mem_new(conf_local_ncells, conf_basis->num_basis);
  return up;
}

void gkyl_correct_mj_fix(gkyl_correct_mj *cmj, const struct gkyl_array *p_over_gamma,
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
  for (int d=0; d<vdim; ++d){
    gkyl_dg_div_op_range(cmj->mem, cmj->conf_basis, d, cmj->num_vb,
      d, cmj->num_vb, 0, cmj->num_ratio, conf_local);
  }

  // compute number density ratio
  gkyl_dg_div_op_range(cmj->mem, cmj->conf_basis, 0, cmj->num_ratio,
    0, m0, 0, cmj->num_ratio, conf_local);

  // calculate cmj->gamma from cmj->num_vb
  gkyl_calc_sr_vars_Gamma(&cmj->conf_basis, &cmj->phase_basis,
      conf_local, cmj->num_vb, cmj->gamma);

  // multiply the number density ratio by gamma, to account for the frame trans.
  gkyl_dg_mul_op_range(cmj->conf_basis, 0, cmj->num_ratio,
      0, cmj->num_ratio, 0, cmj->gamma, conf_local);

  // rescale distribution function
  gkyl_dg_mul_conf_phase_op_range(&cmj->conf_basis, &cmj->phase_basis,
    fout, cmj->num_ratio, fout, conf_local, phase_local);

}

void
gkyl_correct_mj_release(gkyl_correct_mj* cmj)
{

  //gkyl_mom_calc_release(cmj->m0calc);
  gkyl_dg_updater_moment_release(cmj->m0calc);
  gkyl_dg_updater_moment_release(cmj->m1icalc);
  gkyl_array_release(cmj->num_ratio);
  gkyl_array_release(cmj->num_vb);
  gkyl_array_release(cmj->gamma);
  gkyl_dg_bin_op_mem_release(cmj->mem);

  gkyl_free(cmj);
}
