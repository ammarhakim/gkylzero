#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_sr_vars.h>
#include <gkyl_dg_calc_sr_vars_priv.h>
#include <gkyl_util.h>

void gkyl_calc_sr_vars_init_p_vars(const struct gkyl_rect_grid *vgrid, 
  const struct gkyl_basis *vbasis, const struct gkyl_range *vrange,
  struct gkyl_array* p_over_gamma, struct gkyl_array* gamma, struct gkyl_array* gamma_inv)
{
  int vdim = vbasis->ndim;
  int poly_order = vbasis->poly_order;
  // Create projection objects for p/gamma, gamma, and 1/gamma
  gkyl_proj_on_basis *p_over_gamma_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = vgrid,
      .basis = vbasis,
      .qtype = GKYL_GAUSS_LOBATTO_QUAD,
      .num_quad = 8,
      .num_ret_vals = vdim,
      .eval = p_over_gamma_func[vdim-1],
      .ctx = 0
    }
  );  
  gkyl_proj_on_basis *gamma_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = vgrid,
      .basis = vbasis,
      .qtype = GKYL_GAUSS_LOBATTO_QUAD,
      .num_quad = 8,
      .num_ret_vals = 1,
      .eval = gamma_func[vdim-1],
      .ctx = 0
    }
  );  
  gkyl_proj_on_basis *gamma_inv_proj = gkyl_proj_on_basis_inew( &(struct gkyl_proj_on_basis_inp) {
      .grid = vgrid,
      .basis = vbasis,
      .qtype = GKYL_GAUSS_LOBATTO_QUAD,
      .num_quad = 8,
      .num_ret_vals = 1,
      .eval = gamma_inv_func[vdim-1],
      .ctx = 0
    }
  );  
  // Run updater
  gkyl_proj_on_basis_advance(p_over_gamma_proj, 0.0, vrange, p_over_gamma);
  gkyl_proj_on_basis_advance(gamma_proj, 0.0, vrange, gamma);
  gkyl_proj_on_basis_advance(gamma_inv_proj, 0.0, vrange, gamma_inv);
  // Free projection object
  gkyl_proj_on_basis_release(p_over_gamma_proj);
  gkyl_proj_on_basis_release(gamma_proj);
  gkyl_proj_on_basis_release(gamma_inv_proj);
}

void gkyl_calc_sr_vars_Gamma2(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV2)
{
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  sr_t sr_Gamma2 = choose_ser_sr_Gamma2_kern(cdim, vdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *V_d = gkyl_array_cfetch(V, loc);

    double *GammaV2_d = gkyl_array_fetch(GammaV2, loc);
    sr_Gamma2(V_d, GammaV2_d);
  }
}

void gkyl_calc_sr_vars_Gamma(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV)
{
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  sr_t sr_Gamma = choose_ser_sr_Gamma_kern(cdim, vdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *V_d = gkyl_array_cfetch(V, loc);

    double *GammaV_d = gkyl_array_fetch(GammaV, loc);
    sr_Gamma(V_d, GammaV_d);
  }
}

void gkyl_calc_sr_vars_Gamma_inv(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* range, 
  const struct gkyl_array* V, struct gkyl_array* GammaV_inv)
{
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  sr_t sr_Gamma_inv = choose_ser_sr_Gamma_inv_kern(cdim, vdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *V_d = gkyl_array_cfetch(V, loc);

    double *GammaV_inv_d = gkyl_array_fetch(GammaV_inv, loc);
    sr_Gamma_inv(V_d, GammaV_inv_d);
  }
}