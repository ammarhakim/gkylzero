#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_canonical_pb_fluid_vars.h>
#include <gkyl_dg_calc_canonical_pb_fluid_vars_priv.h>
#include <gkyl_wv_canonical_pb_fluid.h>
#include <gkyl_util.h>

gkyl_dg_calc_canonical_pb_fluid_vars*
gkyl_dg_calc_canonical_pb_fluid_vars_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_wv_eqn *wv_eqn, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_calc_canonical_pb_fluid_vars_cu_dev_new(conf_grid, 
      conf_basis, wv_eqn);
  } 
#endif     
  gkyl_dg_calc_canonical_pb_fluid_vars *up = gkyl_malloc(sizeof(gkyl_dg_calc_canonical_pb_fluid_vars));

  up->conf_grid = *conf_grid;
  int cdim = conf_basis->ndim;
  int poly_order = conf_basis->poly_order;
  up->cdim = cdim;
  up->alpha = 0.0;
  up->kappa = 0.0;
  up->is_modified = 0; 

  if (wv_eqn->type == GKYL_EQN_CAN_PB_HASEGAWA_MIMA) {
    up->kappa = gkyl_wv_can_pb_hasegawa_mima_kappa(wv_eqn); 
    up->canonical_pb_fluid_source = choose_canonical_pb_fluid_hasegawa_mima_source_kern(conf_basis->b_type, cdim, poly_order);
  }
  else if (wv_eqn->type == GKYL_EQN_CAN_PB_HASEGAWA_WAKATANI) {
    up->alpha = gkyl_wv_can_pb_hasegawa_wakatani_alpha(wv_eqn); 
    up->kappa = gkyl_wv_can_pb_hasegawa_wakatani_kappa(wv_eqn); 
    up->is_modified = gkyl_wv_can_pb_hasegawa_wakatani_is_modified(wv_eqn); 
    up->canonical_pb_fluid_source = choose_canonical_pb_fluid_hasegawa_wakatani_source_kern(conf_basis->b_type, cdim, poly_order);
  }
  else {
    // Default source kernel; immediately returns and does not do anything. 
    up->canonical_pb_fluid_source = choose_canonical_pb_fluid_default_source_kern(conf_basis->b_type, cdim, poly_order);
  }
  for (int d=0; d<cdim; ++d) {
    up->alpha_surf[d] = choose_canonical_pb_fluid_alpha_surf_kern(conf_basis->b_type, d, cdim, poly_order);
    up->alpha_edge_surf[d] = choose_canonical_pb_fluid_alpha_edge_surf_kern(conf_basis->b_type, d, cdim, poly_order);
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_dg_calc_canonical_pb_fluid_vars_alpha_surf(struct gkyl_dg_calc_canonical_pb_fluid_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *conf_ext_range, 
  const struct gkyl_array* phi,
  struct gkyl_array* alpha_surf, struct gkyl_array* sgn_alpha_surf, struct gkyl_array* const_sgn_alpha)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(alpha_surf)) {
    return gkyl_dg_calc_canonical_pb_fluid_vars_alpha_surf_cu(up, conf_range, conf_ext_range, phi, 
      alpha_surf, sgn_alpha_surf, const_sgn_alpha);
  }
#endif
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_edge[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idx);
    long loc = gkyl_range_idx(conf_range, idx);
    gkyl_rect_grid_cell_center(&up->conf_grid, idx, xc);

    double* alpha_surf_d = gkyl_array_fetch(alpha_surf, loc);
    double* sgn_alpha_surf_d = gkyl_array_fetch(sgn_alpha_surf, loc);
    int* const_sgn_alpha_d = gkyl_array_fetch(const_sgn_alpha, loc);

    // Fill in the configuration space alpha_surf
    for (int dir = 0; dir<cdim; ++dir) {
      const_sgn_alpha_d[dir] = up->alpha_surf[dir](xc, up->conf_grid.dx, 
        (const double*) gkyl_array_cfetch(phi, loc),
        alpha_surf_d, sgn_alpha_surf_d);

      // If the configuration space index is at the local configuration space upper value, we
      // we are at the configuration space upper edge and we also need to evaluate 
      // alpha = +1 to avoid evaluating the potential in the ghost cells 
      // where it is not defined when computing the final surface alpha we need
      // (since the surface alpha array stores only the *lower* surface expansion)
      if (idx[dir] == conf_range->upper[dir]) {
        gkyl_copy_int_arr(cdim, idx, idx_edge);
        idx_edge[dir] = idx_edge[dir]+1;
        long loc_ext = gkyl_range_idx(conf_ext_range, idx_edge);

        double* alpha_surf_ext_d = gkyl_array_fetch(alpha_surf, loc_ext);
        double* sgn_alpha_surf_ext_d = gkyl_array_fetch(sgn_alpha_surf, loc_ext);
        int* const_sgn_alpha_ext_d = gkyl_array_fetch(const_sgn_alpha, loc_ext);
        const_sgn_alpha_ext_d[dir] = up->alpha_edge_surf[dir](xc, up->conf_grid.dx, 
          (const double*) gkyl_array_cfetch(phi, loc),
          alpha_surf_ext_d, sgn_alpha_surf_ext_d);
      }  
    }
  }
}

void gkyl_canonical_pb_fluid_vars_source(struct gkyl_dg_calc_canonical_pb_fluid_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array *background_n_gradient, const struct gkyl_array *phi, 
  const struct gkyl_array *fluid, struct gkyl_array *rhs)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(rhs)) {
    return gkyl_canonical_pb_fluid_vars_source_cu(up, conf_range, phi, fluid, rhs);
  }
#endif
  int cdim = up->cdim; 
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *background_n_gradient_d = gkyl_array_cfetch(background_n_gradient, loc);
    const double *phi_d = gkyl_array_cfetch(phi, loc);
    const double *fluid_d = gkyl_array_cfetch(fluid, loc);

    double* rhs_d = gkyl_array_fetch(rhs, loc);

    up->canonical_pb_fluid_source(up->conf_grid.dx, up->alpha, up->kappa, 
      background_n_gradient_d, phi_d, fluid_d, rhs_d);
  }
}

void gkyl_dg_calc_canonical_pb_fluid_vars_release(gkyl_dg_calc_canonical_pb_fluid_vars *up)
{
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
