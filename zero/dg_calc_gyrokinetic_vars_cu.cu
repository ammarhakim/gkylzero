/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_gyrokinetic_vars.h>
#include <gkyl_dg_calc_gyrokinetic_vars_priv.h>
#include <gkyl_util.h>
}

static void
gkyl_parallelize_components_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range range, int ncomp)
{
  // Create a 2D thread grid so we launch ncomp*range.volume number of threads and can parallelize over components too
  dimBlock->y = ncomp;
  dimGrid->y = 1;
  dimBlock->x = gkyl_int_div_up(252, ncomp);
  dimGrid->x = gkyl_int_div_up(range.volume, dimBlock->x);
}

__global__ void
gkyl_dg_calc_gyrokinetic_vars_flux_surf_cu_kernel(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  struct gkyl_range conf_ext_range, struct gkyl_range phase_ext_range, const struct gkyl_array *phi, 
  const struct gkyl_array *fin, struct gkyl_array* flux_surf, struct gkyl_array *cflrate)
{ 
  int pdim = up->pdim;
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_edge[GKYL_MAX_DIM], idx_vel[2];
  int idxL[GKYL_MAX_DIM];
  int idx_velL[2];
  double xc[GKYL_MAX_DIM];

  // 2D thread grid
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    for (int d=cdim; d<pdim; d++) idx_vel[d-cdim] = idx[d];

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc_conf = gkyl_range_idx(&conf_range, idx);
    long loc_vel = gkyl_range_idx(&up->vel_map->local_vel, idx_vel);
    long loc_phase = gkyl_range_idx(&phase_range, idx);

    const double *bmag_d = (const double*) gkyl_array_cfetch(up->gk_geom->geo_int.bmag, loc_conf);
    const double *phi_d = (const double*) gkyl_array_cfetch(phi, loc_conf);
    const double *vmap_d = (const double*) gkyl_array_cfetch(up->vel_map->vmap, loc_vel);
    const double *vmapSq_d = (const double*) gkyl_array_cfetch(up->vel_map->vmap_sq, loc_vel);

    double* flux_surf_d = (double*) gkyl_array_fetch(flux_surf, loc_phase);
    double *cflrate_d = (double*) gkyl_array_fetch(cflrate, loc_phase);

    for (int dir = 0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(pdim, idx, idxL);
      idxL[dir] = idx[dir] - 1;
      long locL = gkyl_range_idx(&phase_range, idxL);
      const double *fL = (const double*) gkyl_array_cfetch(fin, locL);
      const double *fR = (const double*) gkyl_array_cfetch(fin, loc_phase);

      const struct gkyl_dg_surf_geom *dgs = gkyl_dg_geom_get_surf(up->dg_geom, dir, idx);
      const struct gkyl_gk_dg_surf_geom *gkdgs = gkyl_gk_dg_geom_get_surf(up->gk_dg_geom, dir, idx);
      cflrate_d[0] += up->flux_surf[dir](xc, up->phase_grid.dx, 
        vmap_d, vmapSq_d, up->charge, up->mass,
        dgs, gkdgs,
        bmag_d, phi_d,  fL, fR, flux_surf_d);

      // If the phase space index is at the local configuration space upper value, we
      // we are at the configuration space upper edge and we also need to evaluate 
      // alpha = +1 to avoid evaluating the geometry information in the ghost cells 
      // where it is not defined when computing the final surface alpha we need
      // (since the surface alpha array stores only the *lower* surface expansion)
      if (idx[dir] == phase_range.upper[dir]) {
        gkyl_copy_int_arr(pdim, idx, idx_edge);
        idx_edge[dir] = idx_edge[dir]+1;
        long loc_conf_ext = gkyl_range_idx(&conf_ext_range, idx_edge);
        long loc_phase_ext = gkyl_range_idx(&phase_ext_range, idx_edge);

        double *cflrate_ext_d = (double*) gkyl_array_fetch(cflrate, loc_phase_ext);
        const double *fL = (const double*)  gkyl_array_cfetch(fin, loc_phase);
        const double *fR = (const double*)  gkyl_array_cfetch(fin, loc_phase_ext);
        const struct gkyl_dg_surf_geom *dgs = gkyl_dg_geom_get_surf(up->dg_geom, dir, idx_edge);
        const struct gkyl_gk_dg_surf_geom *gkdgs = gkyl_gk_dg_geom_get_surf(up->gk_dg_geom, dir, idx_edge);

        double* flux_surf_ext_d = (double*) gkyl_array_fetch(flux_surf, loc_phase_ext);
        cflrate_ext_d[0] = up->flux_edge_surf[dir](xc, up->phase_grid.dx, 
          vmap_d, vmapSq_d, up->charge, up->mass,
          dgs, gkdgs,
          bmag_d, phi_d, fL, fR, flux_surf_ext_d);
      }  
    }
  }
}

__global__ void
gkyl_dg_calc_gyrokinetic_vars_flux_surfvpar_cu_kernel(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  struct gkyl_range conf_ext_range, struct gkyl_range phase_ext_range, struct gkyl_range vpar_range,
  const struct gkyl_array *phi, const struct gkyl_array *fin, 
  struct gkyl_array* flux_surf, struct gkyl_array *cflrate)
{ 
  int pdim = up->pdim;
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_edge[GKYL_MAX_DIM], idx_vel[2];
  int idxL[GKYL_MAX_DIM];
  int idx_velL[2];
  double xc[GKYL_MAX_DIM];

  // 2D thread grid
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < vpar_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&vpar_range, linc1, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    for (int d=cdim; d<pdim; d++) idx_vel[d-cdim] = idx[d];

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc_conf = gkyl_range_idx(&conf_range, idx);
    long loc_vel = gkyl_range_idx(&up->vel_map->local_vel, idx_vel);
    long loc_phase = gkyl_range_idx(&phase_range, idx);

    const double *bmag_d = (const double*) gkyl_array_cfetch(up->gk_geom->geo_int.bmag, loc_conf);
    const double *phi_d = (const double*) gkyl_array_cfetch(phi, loc_conf);
    const double *vmap_d = (const double*) gkyl_array_cfetch(up->vel_map->vmap, loc_vel);
    const double *vmapSq_d = (const double*) gkyl_array_cfetch(up->vel_map->vmap_sq, loc_vel);

    double* flux_surf_d = (double*) gkyl_array_fetch(flux_surf, loc_phase);
    double *cflrate_d = (double*) gkyl_array_fetch(cflrate, loc_phase);

    int dir = cdim;
    gkyl_copy_int_arr(pdim, idx, idxL);
    idxL[dir] = idx[dir] - 1;
    gkyl_copy_int_arr(pdim-cdim, idx_vel, idx_velL);
    idx_velL[0] = idx_velL[0]-1;
    long locL = gkyl_range_idx(&phase_range, idxL);
    long loc_velL = gkyl_range_idx(&up->vel_map->local_vel, idx_velL);
    const double *fL = (const double*) gkyl_array_cfetch(fin, locL);
    const double *fR = (const double*) gkyl_array_cfetch(fin, loc_phase);

    const double *vpL = (const double*) gkyl_array_cfetch(up->vel_map->vmap_prime, loc_velL);
    const double *vpR = (const double*) gkyl_array_cfetch(up->vel_map->vmap_prime, loc_vel);

    const struct gkyl_dg_vol_geom *dgv = gkyl_dg_geom_get_vol(up->dg_geom, idx);
    const struct gkyl_gk_dg_vol_geom *gkdgv = gkyl_gk_dg_geom_get_vol(up->gk_dg_geom, idx);

    cflrate_d[0] += up->flux_surfvpar[0](xc, up->phase_grid.dx, 
      vpL, vpR,
      vmap_d, vmapSq_d, up->charge, up->mass,
      dgv, gkdgv, bmag_d, phi_d,  fL, fR, flux_surf_d);
  }
}

// Host-side wrapper for gyrokinetic surface alpha calculation
void gkyl_dg_calc_gyrokinetic_vars_flux_surf_cu(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_ext_range, const struct gkyl_range *phase_ext_range, const struct gkyl_array *phi, 
  const struct gkyl_array *fin, struct gkyl_array* flux_surf, struct gkyl_array* cflrate)
{
  gkyl_dg_calc_gyrokinetic_vars_flux_surf_cu_kernel<<<phase_range->volume, GKYL_DEFAULT_NUM_THREADS>>>(up->on_dev, 
    *conf_range, *phase_range, *conf_ext_range, *phase_ext_range, phi->on_dev, fin->on_dev,
    flux_surf->on_dev, cflrate->on_dev);

  struct gkyl_range vpar_range;
  int sublower[GKYL_MAX_DIM];
  int subupper[GKYL_MAX_DIM];
  for(int i = 0; i < up->pdim; i++) {
     sublower[i] = phase_range->lower[i];
     subupper[i] = phase_range->upper[i];
  }
  sublower[up->cdim] += 1;
  gkyl_sub_range_init(&vpar_range, phase_ext_range, sublower, subupper);
  gkyl_dg_calc_gyrokinetic_vars_flux_surfvpar_cu_kernel<<<vpar_range.volume, GKYL_DEFAULT_NUM_THREADS>>>(up->on_dev, 
    *conf_range, *phase_range, *conf_ext_range, *phase_ext_range, vpar_range, phi->on_dev, fin->on_dev,
    flux_surf->on_dev, cflrate->on_dev);

}

// CUDA kernel to set device pointers to gyrokinetic vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_gyrokinetic_vars_set_cu_dev_ptrs(struct gkyl_dg_calc_gyrokinetic_vars *up, 
  int cdim, int vdim, int poly_order, enum gkyl_gkmodel_id gkmodel_id)
{
  if (gkmodel_id == GKYL_GK_MODEL_NO_BY) {
    for (int d=0; d<cdim; ++d) {
      up->flux_surf[d] = choose_gyrokinetic_flux_no_by_surf_conf_kern(d, cdim, vdim, poly_order);
      up->flux_edge_surf[d] = choose_gyrokinetic_flux_no_by_edge_surf_conf_kern(d, cdim, vdim, poly_order);
    }
    up->flux_surfvpar[0] = choose_gyrokinetic_flux_no_by_surf_vpar_kern(cdim, vdim, poly_order);
  }
  else {
    for (int d=0; d<cdim; ++d) {
      up->flux_surf[d] = choose_gyrokinetic_flux_surf_conf_kern(d, cdim, vdim, poly_order);
      up->flux_edge_surf[d] = choose_gyrokinetic_flux_edge_surf_conf_kern(d, cdim, vdim, poly_order);
    }
    up->flux_surfvpar[0] = choose_gyrokinetic_flux_surf_vpar_kern(cdim, vdim, poly_order);
  }
}

gkyl_dg_calc_gyrokinetic_vars*
gkyl_dg_calc_gyrokinetic_vars_cu_dev_new(const struct gkyl_rect_grid *phase_grid, 
  const struct gkyl_basis *conf_basis, const struct gkyl_basis *phase_basis, 
  const double charge, const double mass, enum gkyl_gkmodel_id gkmodel_id, 
  const struct gk_geometry *gk_geom, const struct gkyl_dg_geom *dg_geom, 
  const struct gkyl_gk_dg_geom *gk_dg_geom, const struct gkyl_velocity_map *vel_map)
{
  struct gkyl_dg_calc_gyrokinetic_vars *up = (struct gkyl_dg_calc_gyrokinetic_vars*) gkyl_malloc(sizeof(*up));

  up->phase_grid = *phase_grid;
  int cdim = conf_basis->ndim;
  int pdim = phase_basis->ndim;
  int vdim = pdim - cdim;
  int poly_order = phase_basis->poly_order;
  up->cdim = cdim;
  up->pdim = pdim;

  up->charge = charge;
  up->mass = mass;

  // Acquire pointers to on_dev objects so memcpy below copies those too.
  struct gk_geometry *geom_ho = gkyl_gk_geometry_acquire(gk_geom);
  struct gkyl_dg_geom *dg_geom_ho = gkyl_dg_geom_acquire(dg_geom);
  struct gkyl_gk_dg_geom *gk_dg_geom_ho = gkyl_gk_dg_geom_acquire(gk_dg_geom);
  struct gkyl_velocity_map *vel_map_ho = gkyl_velocity_map_acquire(vel_map);
  up->gk_geom = geom_ho->on_dev;
  up->dg_geom = dg_geom_ho->on_dev;
  up->gk_dg_geom = gk_dg_geom_ho->on_dev;
  up->vel_map = vel_map_ho->on_dev;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_gyrokinetic_vars *up_cu = (struct gkyl_dg_calc_gyrokinetic_vars*) gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_gyrokinetic_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_gyrokinetic_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, cdim, vdim, poly_order, gkmodel_id);

  // set parent on_dev pointer
  up->on_dev = up_cu;

  // Updater should store host pointers.
  up->gk_geom = geom_ho; 
  up->dg_geom = dg_geom_ho; 
  up->gk_dg_geom = gk_dg_geom_ho; 
  up->vel_map = vel_map_ho; 
  
  return up;
}
