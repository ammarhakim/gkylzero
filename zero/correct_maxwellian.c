#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_correct_maxwellian.h>
#include <gkyl_dg_bin_ops.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>

struct gkyl_correct_maxwellian {
  struct gkyl_rect_grid grid;
  struct gkyl_basis conf_basis, phase_basis;
  
  gkyl_mom_calc *m0calc; // moment calculator
  struct gkyl_array *num_ratio; // number density ratio

  gkyl_dg_bin_op_mem *mem; // memory for division operator
};

gkyl_correct_maxwellian *
gkyl_correct_maxwellian_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis,
  long conf_local_ncells, long conf_local_ext_ncells)    
{
  gkyl_correct_maxwellian *up = gkyl_malloc(sizeof(*up));

  up->grid = *grid;
  up->conf_basis = *conf_basis;
  up->phase_basis = *phase_basis;
  
  struct gkyl_mom_type *m0 = gkyl_mom_vlasov_new(conf_basis, phase_basis, "M0", false);  
  up->m0calc = gkyl_mom_calc_new(grid, m0, false);
  gkyl_mom_type_release(m0);

  up->num_ratio = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);
  up->mem = gkyl_dg_bin_op_mem_new(conf_local_ncells, conf_basis->num_basis);

  return up;
}

void gkyl_correct_maxwellian_fix(gkyl_correct_maxwellian *cmax,
  struct gkyl_array *fout,
  const struct gkyl_array *m0,
  const struct gkyl_range *phase_local, const struct gkyl_range *conf_local)
{
  // calculate number density
  gkyl_mom_calc_advance(cmax->m0calc, phase_local, conf_local, fout, cmax->num_ratio);
  
  // compute number density ratio
  gkyl_dg_div_op_range(cmax->mem, cmax->conf_basis, 0, cmax->num_ratio,
    0, m0, 0, cmax->num_ratio, conf_local);
  
  // rescale distribution function
  gkyl_dg_mul_conf_phase_op_range(&cmax->conf_basis, &cmax->phase_basis,
    fout, cmax->num_ratio, fout, conf_local, phase_local);
}

void
gkyl_correct_maxwellian_release(gkyl_correct_maxwellian* cmax)
{
  gkyl_mom_calc_release(cmax->m0calc);
  gkyl_array_release(cmax->num_ratio);
  gkyl_dg_bin_op_mem_release(cmax->mem);
  
  gkyl_free(cmax);
}
