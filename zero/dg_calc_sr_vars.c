#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_calc_sr_vars_priv.h>
#include <gkyl_util.h>

void gkyl_calc_sr_vars_init_p_vars(const struct gkyl_rect_grid *vgrid, 
  const struct gkyl_basis *vbasis, const struct gkyl_range *vrange,
  struct gkyl_array* gamma, struct gkyl_array* gamma_inv)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(gamma)) {
    return gkyl_calc_sr_vars_init_p_vars_cu(vgrid, vbasis, vrange, gamma, gamma_inv);
  }
#endif

  int vdim = vbasis->ndim;
  int poly_order = vbasis->poly_order; 

  p_vars_t p_vars = choose_ser_sr_p_vars_kern(vdim, poly_order);

  // Cell center array
  double xc[GKYL_MAX_DIM];  
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, vrange);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_rect_grid_cell_center(vgrid, iter.idx, xc);
    long loc = gkyl_range_idx(vrange, iter.idx);

    double *gamma_d = gkyl_array_fetch(gamma, loc);
    double *gamma_inv_d = gkyl_array_fetch(gamma_inv, loc);
    p_vars(xc, vgrid->dx, gamma_d, gamma_inv_d);
  }
}

void gkyl_calc_sr_vars_GammaV2(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV2)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(GammaV2)) {
    return gkyl_calc_sr_vars_GammaV2_cu(cbasis, pbasis, range, V, GammaV2);
  }
#endif

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  sr_t sr_GammaV2 = choose_ser_sr_GammaV2_kern(cdim, vdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *V_d = gkyl_array_cfetch(V, loc);

    double *GammaV2_d = gkyl_array_fetch(GammaV2, loc);
    sr_GammaV2(V_d, GammaV2_d);
  }
}

void gkyl_calc_sr_vars_GammaV_inv(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV_inv)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(GammaV_inv)) {
    return gkyl_calc_sr_vars_GammaV_inv_cu(cbasis, pbasis, range, V, GammaV_inv);
  }
#endif

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  sr_t sr_GammaV_inv = choose_ser_sr_GammaV_inv_kern(cdim, vdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *V_d = gkyl_array_cfetch(V, loc);

    double *GammaV_inv_d = gkyl_array_fetch(GammaV_inv, loc);
    sr_GammaV_inv(V_d, GammaV_inv_d);
  }
}
