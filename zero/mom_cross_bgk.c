#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_cross_bgk.h>
#include <gkyl_mom_cross_bgk_priv.h>
#include <math.h>

gkyl_mom_cross_bgk_gyrokinetic * 
gkyl_mom_cross_bgk_gyrokinetic_new(const struct gkyl_basis *phase_basis, const struct gkyl_basis *conf_basis)
{
  gkyl_mom_cross_bgk_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  // select the kernel
  up->mom_cross_calc = choose_mom_cross_bgk_gyrokinetic_kern(conf_basis->ndim, phase_basis->ndim-conf_basis->ndim, phase_basis->poly_order); 
}

void gkyl_mom_cross_bgk_gyrokinetic_advance(gkyl_mom_cross_bgk_gyrokinetic *up,
const struct gkyl_range *conf_rng, const double beta,
const double m_self, const struct gkyl_array *moms_self,
const double m_other, const struct gkyl_array *moms_other,
const double nu_sr, const double nu_rs, struct gkyl_array *moms_cross)
{
  struct gkyl_range_iter conf_iter;  
  
  // loop over configuration space cells
  gkyl_range_iter_init(&conf_iter, conf_rng);
  while (gkyl_range_iter_next(&conf_iter)) {
    long midx = gkyl_range_idx(conf_rng, conf_iter.idx);
    const double *in1_d = gkyl_array_cfetch(moms_self, midx);
    const double *in2_d = gkyl_array_cfetch(moms_other, midx);
    double *out_d = gkyl_array_fetch(moms_cross, midx);

    up->mom_cross_calc(beta, m_self, in1_d, m_other, in2_d, nu_sr, nu_rs, out_d);
  }
}

void gkyl_mom_cross_bgk_gyrokinetic_release(gkyl_mom_cross_bgk_gyrokinetic *up)
{
  gkyl_free(up);
}
