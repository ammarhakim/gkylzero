#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_gyrokinetic_vars.h>
#include <gkyl_dg_calc_gyrokinetic_vars_priv.h>
#include <gkyl_util.h>

gkyl_dg_calc_gyrokinetic_vars*
gkyl_dg_calc_gyrokinetic_vars_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const double charge, const double mass, const struct gk_geometry *gk_geom, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_calc_gyrokinetic_vars_cu_dev_new(phase_grid, conf_basis, phase_basis, 
      charge, mass, gk_geom);
  } 
#endif     
  gkyl_dg_calc_gyrokinetic_vars *up = gkyl_malloc(sizeof(gkyl_dg_calc_gyrokinetic_vars));

  up->phase_grid = *phase_grid;
  int cdim = conf_basis->ndim;
  int pdim = phase_basis->ndim;
  int vdim = pdim - cdim;
  int poly_order = phase_basis->poly_order;
  up->cdim = cdim;
  up->pdim = pdim;

  up->charge = charge;
  up->mass = mass;
  up->gk_geom = gkyl_gk_geometry_acquire(gk_geom);

  up->Bstar_Bmag = choose_gyrokinetic_Bstar_Bmag_kern(cdim, poly_order);

  for (int d=0; d<cdim; ++d) {
    up->alpha_surf[d] = choose_gyrokinetic_alpha_surf_conf_kern(d, cdim, vdim, poly_order);
    up->alpha_edge_surf[d] = choose_gyrokinetic_alpha_edge_surf_conf_kern(d, cdim, vdim, poly_order);
  }
  up->alpha_surf[cdim] = choose_gyrokinetic_alpha_surf_vpar_kern(cdim, vdim, poly_order);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_dg_calc_gyrokinetic_vars_Bstar_Bmag(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, struct gkyl_array* Bstar_Bmag)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(Bstar_Bmag)) {
    return gkyl_dg_calc_gyrokinetic_vars_Bstar_Bmag_cu(up, conf_range, phase_range, Bstar_Bmag);
  }
#endif
  int pdim = up->pdim;
  int idx[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idx);
    long loc_conf = gkyl_range_idx(conf_range, idx);
    long loc_phase = gkyl_range_idx(phase_range, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    const double *bmag_d = gkyl_array_cfetch(up->gk_geom->bmag, loc_conf);
    const double *jacobtot_inv_d = gkyl_array_cfetch(up->gk_geom->jacobtot_inv, loc_conf);
    const double *cmag_d = gkyl_array_cfetch(up->gk_geom->cmag, loc_conf);
    const double *b_i_d = gkyl_array_cfetch(up->gk_geom->b_i, loc_conf);

    double* Bstar_Bmag_d = gkyl_array_fetch(Bstar_Bmag, loc_phase);
    up->Bstar_Bmag(xc, up->phase_grid.dx, up->charge, up->mass, 
      bmag_d, jacobtot_inv_d, cmag_d, b_i_d, Bstar_Bmag_d);
  }
}

void gkyl_dg_calc_gyrokinetic_vars_alpha_surf(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,  const struct gkyl_range *phase_ext_range, 
  const struct gkyl_array *phi, const struct gkyl_array *Bstar_Bmag, struct gkyl_array* alpha_surf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(alpha_surf)) {
    return gkyl_dg_calc_gyrokinetic_vars_alpha_surf_cu(up, conf_range, phase_range, phase_ext_range, 
      phi, Bstar_Bmag, alpha_surf);
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

    const double *bmag_d = gkyl_array_cfetch(up->gk_geom->bmag, loc_conf);
    const double *jacobtot_inv_d = gkyl_array_cfetch(up->gk_geom->jacobtot_inv, loc_conf);
    const double *cmag_d = gkyl_array_cfetch(up->gk_geom->cmag, loc_conf);
    const double *b_i_d = gkyl_array_cfetch(up->gk_geom->b_i, loc_conf);
    const double *phi_d = gkyl_array_cfetch(phi, loc_conf);
    const double *Bstar_Bmag_d = gkyl_array_cfetch(Bstar_Bmag, loc_phase);

    double* alpha_surf_d = gkyl_array_fetch(alpha_surf, loc_phase);
    for (int dir = 0; dir<cdim+1; ++dir) {
      up->alpha_surf[dir](xc, up->phase_grid.dx, up->charge, up->mass, 
        bmag_d, jacobtot_inv_d, cmag_d, b_i_d, phi_d,  Bstar_Bmag_d, alpha_surf_d);

      // If the phase space index is at the local configuration space upper value, we
      // we are at the configuration space upper edge and we also need to evaluate 
      // alpha = +1 to avoid evaluating the geometry information in the ghost cells 
      // where it is not defined when computing the final surface alpha we need
      // (since the surface alpha array stores only the *lower* surface expansion)
      if (dir < cdim && idx[dir] == conf_range->upper[dir]) {
        gkyl_copy_int_arr(pdim, idx, idx_edge);
        idx_edge[dir] = idx_edge[dir]+1;
        long loc_phase_ext = gkyl_range_idx(phase_ext_range, idx_edge);

        double* alpha_surf_ext_d = gkyl_array_fetch(alpha_surf, loc_phase_ext);
        up->alpha_edge_surf[dir](xc, up->phase_grid.dx, up->charge, up->mass, 
          bmag_d, jacobtot_inv_d, cmag_d, b_i_d, phi_d, Bstar_Bmag_d, alpha_surf_ext_d); 
      }  
    }
  }
}

void gkyl_dg_calc_gyrokinetic_vars_release(gkyl_dg_calc_gyrokinetic_vars *up)
{
  gkyl_gk_geometry_release(up->gk_geom);
  
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}