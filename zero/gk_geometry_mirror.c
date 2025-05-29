#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mirror.h>
#include <gkyl_position_map.h>
#include <gkyl_mirror_grid_gen.h>
#include <gkyl_mirror_geo_gen.h>
#include <gkyl_mirror_geo_dg.h>

struct gk_geometry*
gk_geometry_mirror_init(struct gkyl_gk_geometry_inp *geometry_inp)
{

  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->basis = geometry_inp->geo_basis;
  up->local = geometry_inp->geo_local;
  up->local_ext = geometry_inp->geo_local_ext;
  up->global = geometry_inp->geo_global;
  up->global_ext = geometry_inp->geo_global_ext;
  up->grid = geometry_inp->geo_grid;

  // Define the nodal range for the geometry
  struct gkyl_range nrange;
  int poly_order = up->basis.poly_order;
  int nodes[GKYL_MAX_DIM];
  if (poly_order == 1) {
    for (int d=0; d<up->grid.ndim; ++d)
      nodes[d] = gkyl_range_shape(&up->local, d) + 1;
  }
  if (poly_order == 2) {
    for (int d=0; d<up->grid.ndim; ++d)
      nodes[d] = 2*gkyl_range_shape(&up->local, d) + 1;
  }
  gkyl_range_init_from_shape(&nrange, up->grid.ndim, nodes);

  up->geqdsk_sign_convention = 0.0; // Hardcoded 0. Means psi increases from axis. Always true for mirror geometry.

  // read psi(R,Z) from file
  struct gkyl_rect_grid psi_grid;
  struct gkyl_array *psi = gkyl_grid_array_new_from_file(&psi_grid, geometry_inp->mirror_grid_info.filename_psi);

  // create mirror geometry
  struct gkyl_mirror_grid_gen *mirror_grid =
    gkyl_mirror_grid_gen_inew(&(struct gkyl_mirror_grid_gen_inp) {
        .comp_grid = &up->grid,
        
        .R = { psi_grid.lower[0], psi_grid.upper[0] },
        .Z = { psi_grid.lower[1], psi_grid.upper[1] },
        
        // psi(R,Z) grid size
        .nrcells = psi_grid.cells[0]-1, // cells and not nodes
        .nzcells = psi_grid.cells[1]-1, // cells and not nodes

        .psiRZ = psi,
        .fl_coord = geometry_inp->mirror_grid_info.fl_coord,
        .include_axis = geometry_inp->mirror_grid_info.include_axis,
        .write_psi_cubic = false,
      }
    );

  // Generate the required geometry objects at each node
  struct gkyl_mirror_geo_gen *mirror_geo = 
    gkyl_mirror_geo_gen_inew(&(struct gkyl_mirror_geo_gen_inp) {
        .comp_grid = &up->grid,
        .mirror_grid = mirror_grid,
        .range = up->global,
        .basis = up->basis,
      }
    );

  // Project the geometry objects from the nodes onto DG arrays
  struct gkyl_mirror_geo_dg *mirror_geo_dg = 
    gkyl_mirror_geo_dg_inew(&(struct gkyl_mirror_geo_dg_inp) {
        .comp_grid = &up->grid,
        .mirror_geo = mirror_geo,
        .range = up->global,
        .range_ext = up->global_ext,
        .basis = up->basis,
      }
    );

  // Global geometry is generated on all cores. Copy the local geometry to each process
  if (up->local.volume == up->global.volume) // If local = global
  {
    up->mc2p = gkyl_array_acquire(mirror_geo_dg->mapc2p);
    up->mc2nu_pos = gkyl_array_acquire(mirror_geo_dg->mc2nu_pos);
    up->dxdz = gkyl_array_acquire(mirror_geo_dg->tang);
    up->dzdx = gkyl_array_acquire(mirror_geo_dg->dual);
    up->dualmag = gkyl_array_acquire(mirror_geo_dg->dualmag);
    up->normals = gkyl_array_acquire(mirror_geo_dg->normals);
    up->g_ij = gkyl_array_acquire(mirror_geo_dg->metric_covar);
    up->g_ij_neut = gkyl_array_acquire(mirror_geo_dg->metric_covar_neut);
    up->gij = gkyl_array_acquire(mirror_geo_dg->metric_contr);
    up->gij_neut = gkyl_array_acquire(mirror_geo_dg->metric_contr_neut);
    up->gxxj = gkyl_array_acquire(mirror_geo_dg->gxxj);
    up->gxyj = gkyl_array_acquire(mirror_geo_dg->gxyj);
    up->gyyj = gkyl_array_acquire(mirror_geo_dg->gyyj);
    up->gxzj = gkyl_array_acquire(mirror_geo_dg->gxzj);
    up->jacobgeo = gkyl_array_acquire(mirror_geo_dg->Jc);
    up->jacobgeo_inv = gkyl_array_acquire(mirror_geo_dg->Jc_inv);
    up->jacobtot = gkyl_array_acquire(mirror_geo_dg->JB);
    up->jacobtot_inv = gkyl_array_acquire(mirror_geo_dg->JB_inv);
    up->b_i = gkyl_array_acquire(mirror_geo_dg->b_covar);
    up->bcart = gkyl_array_acquire(mirror_geo_dg->b_cart);
    up->bmag = gkyl_array_acquire(mirror_geo_dg->Bmag);
    up->bmag_inv = gkyl_array_acquire(mirror_geo_dg->Bmag_inv);
    up->bmag_inv_sq = gkyl_array_acquire(mirror_geo_dg->Bmag_inv_sq);
    up->cmag = gkyl_array_acquire(mirror_geo_dg->C);
    up->eps2 = gkyl_array_acquire(mirror_geo_dg->eps2);
  }
  else // Code is running in parallel and local != global
  {
    // Allocate local arrays for the geometry
    up->mc2p = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
    up->mc2nu_pos = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
    up->bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
    up->g_ij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
    up->dxdz = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
    up->dzdx = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
    up->dualmag = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
    up->normals = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
    up->jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->gij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
    up->gij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
    up->b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
    up->bcart = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
    up->cmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->jacobtot = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->bmag_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->bmag_inv_sq = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->gxxj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->gxyj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->gyyj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->gxzj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
    up->eps2= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);

    struct gkyl_range *global_sub_range;
    int intersect = gkyl_sub_range_intersect(global_sub_range, &up->global, &up->local);
    // Copy the geometry data from the global arrays to the local arrays
    gkyl_array_copy_range_to_range(up->mc2p, mirror_geo_dg->mapc2p, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->mc2nu_pos, mirror_geo_dg->mc2nu_pos, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->bmag, mirror_geo_dg->Bmag, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->g_ij, mirror_geo_dg->metric_covar, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->g_ij_neut, mirror_geo_dg->metric_covar_neut, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->dxdz, mirror_geo_dg->tang, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->dzdx, mirror_geo_dg->dual, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->dualmag, mirror_geo_dg->dualmag, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->normals, mirror_geo_dg->normals, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->jacobgeo, mirror_geo_dg->Jc, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->jacobgeo_inv, mirror_geo_dg->Jc_inv, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->gij, mirror_geo_dg->metric_contr, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->gij_neut, mirror_geo_dg->metric_contr_neut, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->b_i, mirror_geo_dg->b_covar, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->bcart, mirror_geo_dg->b_cart, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->cmag, mirror_geo_dg->C, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->jacobtot, mirror_geo_dg->JB, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->jacobtot_inv, mirror_geo_dg->JB_inv, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->bmag_inv, mirror_geo_dg->Bmag_inv, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->bmag_inv_sq, mirror_geo_dg->Bmag_inv_sq, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->gxxj, mirror_geo_dg->gxxj, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->gxyj, mirror_geo_dg->gxyj, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->gyyj, mirror_geo_dg->gyyj, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->gxzj, mirror_geo_dg->gxzj, &up->local, global_sub_range);
    gkyl_array_copy_range_to_range(up->eps2, mirror_geo_dg->eps2, &up->local, global_sub_range);
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself


  gkyl_mirror_grid_gen_release(mirror_grid);
  gkyl_mirror_geo_gen_release(mirror_geo);
  gkyl_mirror_geo_dg_release(mirror_geo_dg);
  gkyl_array_release(psi);

  return up;
}

struct gk_geometry*
gkyl_gk_geometry_mirror_new(struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry* gk_geom_3d;
  struct gk_geometry* gk_geom;

  if (geometry_inp->position_map == 0){
    geometry_inp->position_map = gkyl_position_map_new((struct gkyl_position_map_inp) {}, \
      geometry_inp->grid, geometry_inp->local, geometry_inp->local_ext, geometry_inp->local, \
      geometry_inp->local_ext, geometry_inp->basis);
    gk_geom_3d = gk_geometry_mirror_init(geometry_inp);
    gkyl_position_map_release(geometry_inp->position_map);
  }
  else {
    // First construct the uniform 3d geometry
    gk_geom_3d = gk_geometry_mirror_init(geometry_inp);
    if (geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_POLYNOMIAL || \
        geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_NUMERIC) {
      // The array mc2nu is computed using the uniform geometry, so we need to deflate it
      // Must deflate the 3D uniform geometry in order for the allgather to work
      if(geometry_inp->grid.ndim < 3)
        gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, geometry_inp);
      else
        gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);

      gkyl_position_map_set_bmag(geometry_inp->position_map, geometry_inp->comm, \
        gk_geom->bmag);

      gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry
      gkyl_gk_geometry_release(gk_geom); // release 3d geometry

      // Construct the non-uniform grid
      gk_geom_3d = gk_geometry_mirror_init(geometry_inp);
    }
  }
  return gk_geom_3d;
}