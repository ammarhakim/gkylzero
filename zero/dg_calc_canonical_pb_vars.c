#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_canonical_pb_vars.h>
#include <gkyl_dg_calc_canonical_pb_vars_priv.h>
#include <gkyl_util.h>

gkyl_dg_calc_canonical_pb_vars*
gkyl_dg_calc_canonical_pb_vars_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_calc_canonical_pb_vars_cu_dev_new(phase_grid, 
      conf_basis, phase_basis, gk_geom);
  } 
#endif     
  gkyl_dg_calc_canonical_pb_vars *up = gkyl_malloc(sizeof(gkyl_dg_calc_canonical_pb_vars));

  up->phase_grid = *phase_grid;
  int cdim = conf_basis->ndim;
  int pdim = phase_basis->ndim;
  int poly_order = phase_basis->poly_order;
  up->cdim = cdim;
  up->pdim = pdim;

  for (int d=0; d<cdim; ++d) {
    up->alpha_surf[d] = choose_canonical_pb_alpha_surf_kern(d, cdim, poly_order);
    up->alpha_surf[d+cdim] = choose_canonical_pb_alpha_surf_v_kern(d, cdim, poly_order);
    up->alpha_edge_surf[d] = choose_canonical_pb_alpha_edge_surf_kern(d, cdim, poly_order);
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_dg_calc_canonical_pb_vars_alpha_surf(struct gkyl_dg_calc_canonical_pb_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,  const struct gkyl_range *phase_ext_range, 
  struct gkyl_array *hamil,
  struct gkyl_array* alpha_surf, struct gkyl_array* sgn_alpha_surf, struct gkyl_array* const_sgn_alpha)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(alpha_surf)) {
    return gkyl_dg_calc_canonical_pb_vars_alpha_surf_cu(up, conf_range, phase_range, phase_ext_range, 
      alpha_surf, sgn_alpha_surf, const_sgn_alpha);
  }
#endif
  int pdim = up->pdim;
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_edge[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idx);
    long loc_conf = gkyl_range_idx(conf_range, idx);
    long loc_phase = gkyl_range_idx(phase_range, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    double* alpha_surf_d = gkyl_array_fetch(alpha_surf, loc_phase);
    double* sgn_alpha_surf_d = gkyl_array_fetch(sgn_alpha_surf, loc_phase);
    int* const_sgn_alpha_d = gkyl_array_fetch(const_sgn_alpha, loc_phase);
    for (int dir = 0; dir<cdim; ++dir) {
      const double *hamil_local =  gkyl_array_cfetch(hamil, loc_phase);
      
      const_sgn_alpha_d[dir] = up->alpha_surf[dir](xc, up->phase_grid.dx, 
        (const double*) gkyl_array_cfetch(hamil, loc_phase),
        alpha_surf_d, sgn_alpha_surf_d);

      const_sgn_alpha_d[dir+cdim] = up->alpha_surf[dir+cdim](xc, up->phase_grid.dx, 
        (const double*) gkyl_array_cfetch(hamil, loc_phase),
        alpha_surf_d, sgn_alpha_surf_d);

      // If the phase space index is at the local configuration space upper value, we
      // we are at the configuration space upper edge and we also need to evaluate 
      // alpha = +1 to avoid evaluating the geometry information in the ghost cells 
      // where it is not defined when computing the final surface alpha we need
      // (since the surface alpha array stores only the *lower* surface expansion)
      if (idx[dir] == conf_range->upper[dir]) {
        gkyl_copy_int_arr(pdim, idx, idx_edge);
        idx_edge[dir] = idx_edge[dir]+1;
        long loc_phase_ext = gkyl_range_idx(phase_ext_range, idx_edge);

        double* alpha_surf_ext_d = gkyl_array_fetch(alpha_surf, loc_phase_ext);
        double* sgn_alpha_surf_ext_d = gkyl_array_fetch(sgn_alpha_surf, loc_phase_ext);
        int* const_sgn_alpha_ext_d = gkyl_array_fetch(const_sgn_alpha, loc_phase_ext);
        const_sgn_alpha_ext_d[dir] = up->alpha_edge_surf[dir](xc, up->phase_grid.dx, 
          (const double*) gkyl_array_cfetch(hamil, loc_phase),
          alpha_surf_ext_d, sgn_alpha_surf_ext_d);
      }  
    }
  }
}

void gkyl_dg_calc_canonical_pb_vars_release(gkyl_dg_calc_canonical_pb_vars *up)
{
  
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}