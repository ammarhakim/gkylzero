#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_gk_rad_vars.h>
#include <gkyl_dg_calc_gk_rad_vars_priv.h>
#include <gkyl_util.h>

gkyl_dg_calc_gk_rad_vars* 
gkyl_dg_calc_gk_rad_vars_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, double charge,
  double mass, const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map, 
  double a, double alpha, double beta, double gamma, double v0, bool use_gpu) 
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_dg_calc_gk_rad_vars_cu_dev_new(phase_grid, conf_basis, phase_basis, 
      charge, mass, gk_geom, vel_map, a, alpha, beta, gamma, v0);
  } 
#endif 

  gkyl_dg_calc_gk_rad_vars *up = gkyl_malloc(sizeof(*up));

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
  up->vel_map = gkyl_velocity_map_acquire(vel_map);

  // Fitting parameters for a given collision type
  up->a = a;
  up->alpha = alpha;
  up->beta = beta;
  up->gamma = gamma;
  up->v0 = v0;

  up->rad_nu_vpar = choose_rad_gyrokinetic_nu_vpar_kern(cdim, vdim, poly_order);
  up->rad_nu_mu = choose_rad_gyrokinetic_nu_mu_kern(cdim, vdim, poly_order);
  up->rad_nI_nu = choose_rad_gyrokinetic_nI_nu_kern(cdim, vdim, poly_order);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_dg_calc_gk_rad_vars_nu_advance(const struct gkyl_dg_calc_gk_rad_vars *up,
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  struct gkyl_array* vnu_surf, struct gkyl_array* vnu, 
  struct gkyl_array* vsqnu_surf, struct gkyl_array* vsqnu)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(vnu_surf)) {
    return gkyl_dg_calc_gk_rad_vars_nu_advance_cu(up, conf_range, phase_range, 
      vnu_surf, vnu, vsqnu_surf, vsqnu);
  }
#endif
  int pdim = up->pdim;
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_vel[2];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idx);

    for (int d=cdim; d<pdim; d++) idx_vel[d-cdim] = iter.idx[d];

    long loc_conf = gkyl_range_idx(conf_range, idx);
    long loc_vel = gkyl_range_idx(&up->vel_map->local_vel, idx_vel);
    long loc_phase = gkyl_range_idx(phase_range, idx);

    const double *bmag_d = gkyl_array_cfetch(up->gk_geom->bmag, loc_conf);

    double* vnu_surf_d = gkyl_array_fetch(vnu_surf, loc_phase);
    double* vnu_d = gkyl_array_fetch(vnu, loc_phase);
    double* vsqnu_surf_d = gkyl_array_fetch(vsqnu_surf, loc_phase);  
    double* vsqnu_d = gkyl_array_fetch(vsqnu, loc_phase);   
    const double *vmap_d = gkyl_array_cfetch(up->vel_map->vmap, loc_vel);
    const double *vmapSq_d = gkyl_array_cfetch(up->vel_map->vmap_sq, loc_vel);

    up->rad_nu_vpar(vmap_d, vmapSq_d, up->charge, up->mass, 
      up->a, up->alpha, up->beta, up->gamma, up->v0, 
      bmag_d, vnu_surf_d, vnu_d);
    up->rad_nu_mu(vmap_d, vmapSq_d, up->charge, up->mass, 
      up->a, up->alpha, up->beta, up->gamma, up->v0, 
      bmag_d, vsqnu_surf_d, vsqnu_d);
  }
}

void gkyl_dg_calc_gk_rad_vars_nI_nu_advance(const struct gkyl_dg_calc_gk_rad_vars *up[GKYL_MAX_RAD_DENSITIES],
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range, 
  const struct gkyl_array* vnu_surf[GKYL_MAX_RAD_DENSITIES],
  const struct gkyl_array* vnu[GKYL_MAX_RAD_DENSITIES], 
  const struct gkyl_array* vsqnu_surf[GKYL_MAX_RAD_DENSITIES],
  const struct gkyl_array* vsqnu[GKYL_MAX_RAD_DENSITIES],
  const struct gkyl_array* n_elc_rad,
  const struct gkyl_array* n_elc,
  const struct gkyl_array *nI, 
  struct gkyl_array* nvnu_surf, struct gkyl_array* nvnu, 
  struct gkyl_array* nvsqnu_surf, struct gkyl_array* nvsqnu)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(vnu_surf)) {
    return gkyl_dg_calc_gk_rad_vars_nI_nu_advance_cu(up, conf_range, phase_range, 
      vnu_surf, vnu, vsqnu_surf, vsqnu, nI, 
      nvnu_surf, nvnu, nvsqnu_surf, nvsqnu);
  }
#endif
  int pdim = up[0]->pdim; // pdim and cdim are constant across densities
  int cdim = up[0]->cdim;
  int idx[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idx);

    long loc_conf = gkyl_range_idx(conf_range, idx);
    long loc_phase = gkyl_range_idx(phase_range, idx);

    gkyl_rect_grid_cell_center(&up[0]->phase_grid, idx, xc);

    const double* ne = gkyl_array_cfetch(n_elc, loc_conf);
    double ne_cell_avg = ne[0]/pow(2.0, cdim/2.0);
    int ne_idx = gkyl_find_nearest_idx(n_elc_rad, ne_cell_avg);
    
    const double* vnu_surf_d = gkyl_array_cfetch(vnu_surf[ne_idx], loc_phase);
    const double* vnu_d = gkyl_array_cfetch(vnu[ne_idx], loc_phase);
    const double* vsqnu_surf_d = gkyl_array_cfetch(vsqnu_surf[ne_idx], loc_phase);  
    const double* vsqnu_d = gkyl_array_cfetch(vsqnu[ne_idx], loc_phase);   

    const double *nI_d = gkyl_array_cfetch(nI, loc_conf);

    double* nvnu_surf_d = gkyl_array_fetch(nvnu_surf, loc_phase);
    double* nvnu_d = gkyl_array_fetch(nvnu, loc_phase);
    double* nvsqnu_surf_d = gkyl_array_fetch(nvsqnu_surf, loc_phase);  
    double* nvsqnu_d = gkyl_array_fetch(nvsqnu, loc_phase);   

    up[ne_idx]->rad_nI_nu(vnu_surf_d, vnu_d, vsqnu_surf_d, vsqnu_d, nI_d, 
      nvnu_surf_d, nvnu_d, nvsqnu_surf_d, nvsqnu_d);
  }
}

void gkyl_dg_calc_gk_rad_vars_release(gkyl_dg_calc_gk_rad_vars *up)
{
  gkyl_gk_geometry_release(up->gk_geom);
  gkyl_velocity_map_release(up->vel_map);
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
