extern "C" {
#include <gkyl_mom_calc.h>
#include <gkyl_util.h>
}

void
gkyl_mom_calc_advance_cu(const gkyl_mom_calc* mcalc,
  const struct gkyl_range phase_range, const struct gkyl_range conf_range,
  const struct gkyl_array* GKYL_RESTRICT fin, struct gkyl_array* GKYL_RESTRICT mout)
{

//  gkyl_mom_calc_advance_cu_ker<<<nblocks, nthreads>>>(mcalc, phase_range, conf_range, fin, mout)

}
