#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_dg_calc_em_vars_priv.h>
#include <gkyl_util.h>

void gkyl_calc_em_vars_bvar(const struct gkyl_basis* cbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* em, struct gkyl_array* bvar)
{
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;

  em_t em_bvar = choose_ser_em_bvar_kern(cdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *em_d = gkyl_array_cfetch(em, loc);

    double *bvar_d = gkyl_array_fetch(bvar, loc);
    em_bvar(em_d, bvar_d);
  }
}

void gkyl_calc_em_vars_ExB(const struct gkyl_basis* cbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* em, struct gkyl_array* ExB)
{
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;

  em_t em_ExB = choose_ser_em_ExB_kern(cdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *em_d = gkyl_array_cfetch(em, loc);

    double *ExB_d = gkyl_array_fetch(ExB, loc);
    em_ExB(em_d, ExB_d);
  }
}

void gkyl_calc_em_vars_pkpm_kappa_inv_b(const struct gkyl_basis* cbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* bvar, const struct gkyl_array* ExB,
  struct gkyl_array* kappa_inv_b)
{
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;

  em_pkpm_kappa_inv_b_t em_pkpm_kappa_inv_b = choose_ser_em_pkpm_kappa_inv_b_kern(cdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *bvar_d = gkyl_array_cfetch(bvar, loc);
    const double *ExB_d = gkyl_array_cfetch(ExB, loc);

    double *kappa_inv_b_d = gkyl_array_fetch(kappa_inv_b, loc);
    em_pkpm_kappa_inv_b(bvar_d, ExB_d, kappa_inv_b_d);
  }
}