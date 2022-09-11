#include "gkyl_elem_type.h"
#include "gkyl_mom_type.h"
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_correct_maxwellian.h>
#include <gkyl_mom_calc.h>
#include <gkyl_mom_vlasov.h>

struct gkyl_correct_maxwellian {
  struct gkyl_rect_grid grid;
  gkyl_mom_calc *m0calc;
  struct gkyl_array *num_ratio; // number density ratio
};

gkyl_correct_maxwellian *
gkyl_correct_maxwellian_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, long conf_local_ext_ncells)
{
  gkyl_correct_maxwellian *up = gkyl_malloc(sizeof(*up));

  up->grid = *grid;
  struct gkyl_mom_type *m0 = gkyl_mom_vlasov_new(conf_basis, phase_basis, "M0");  
  up->m0calc = gkyl_mom_calc_new(grid, m0);
  up->num_ratio = gkyl_array_new(GKYL_DOUBLE, conf_basis->num_basis, conf_local_ext_ncells);

  gkyl_mom_type_release(m0);

  return up;
}

void
gkyl_correct_maxwellian_release(gkyl_correct_maxwellian* cmax)
{
  gkyl_mom_calc_release(cmax->m0calc);
  gkyl_array_release(cmax->num_ratio);
  gkyl_free(cmax);
}
