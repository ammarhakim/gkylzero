#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_correct_MJ.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_updater_moment.h>
#include <gkyl_proj_on_basis.h>

#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov_sr.h>

// Temporary
#include <gkyl_array_ops_priv.h>


struct gkyl_correct_MJ {
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

void gkyl_correct_MJ_fix(gkyl_correct_MJ *cMJ, const struct gkyl_array *p_over_gamma,
  struct gkyl_array *fout,
  const struct gkyl_array *m0,
  const struct gkyl_array *m1i,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  // vdim
  int vdim = cMJ->phase_basis.ndim - cMJ->conf_basis.ndim;


  // compute the sr moments
  gkyl_dg_updater_moment_advance(cMJ->m0calc, phase_local, conf_local,
    0, 0, 0,
    0, 0, 0,
    fout, cMJ->num_ratio);
  gkyl_dg_updater_moment_advance(cMJ->m1icalc, phase_local, conf_local,
    p_over_gamma, 0, 0,
    0, 0, 0,
    fout, cMJ->num_vb);

  // isolate vb by dividing N*vb by N
  for (int d=0; d<vdim; ++d){
    gkyl_dg_div_op_range(cMJ->mem, cMJ->conf_basis, d, cMJ->num_vb,
      d, cMJ->num_vb, 0, cMJ->num_ratio, conf_local);
  }

  // compute number density ratio
  gkyl_dg_div_op_range(cMJ->mem, cMJ->conf_basis, 0, cMJ->num_ratio,
    0, m0, 0, cMJ->num_ratio, conf_local);

  // calculate cMJ->gamma from cMJ->num_vb
  gkyl_calc_sr_vars_Gamma(&cMJ->conf_basis, &cMJ->phase_basis,
      conf_local, cMJ->num_vb, cMJ->gamma);

  // multiply the number density ratio by gamma, to account for the frame trans.
  gkyl_dg_mul_op_range(cMJ->conf_basis, 0, cMJ->num_ratio,
      0, cMJ->num_ratio, 0, cMJ->gamma, conf_local);

  // rescale distribution function
  gkyl_dg_mul_conf_phase_op_range(&cMJ->conf_basis, &cMJ->phase_basis,
    fout, cMJ->num_ratio, fout, conf_local, phase_local);

}

void
gkyl_correct_MJ_release(gkyl_correct_MJ* cMJ)
{

  //gkyl_mom_calc_release(cMJ->m0calc);
  gkyl_dg_updater_moment_release(cMJ->m0calc);
  gkyl_dg_updater_moment_release(cMJ->m1icalc);
  gkyl_array_release(cMJ->num_ratio);
  gkyl_array_release(cMJ->num_vb);
  gkyl_array_release(cMJ->gamma);
  gkyl_dg_bin_op_mem_release(cMJ->mem);

  gkyl_free(cMJ);
}
