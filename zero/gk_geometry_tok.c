#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_comm.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_priv.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_math.h>
#include <gkyl_nodal_ops.h>

#include <gkyl_tok_geo.h>
#include <gkyl_tok_calc_derived_geo.h>
#include <gkyl_calc_metric.h>
#include <gkyl_calc_bmag.h>


struct gk_geometry*
gk_geometry_tok_init(struct gkyl_gk_geometry_inp *geometry_inp)
{

  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->basis = geometry_inp->geo_basis;
  up->local = geometry_inp->geo_local;
  up->local_ext = geometry_inp->geo_local_ext;
  up->global = geometry_inp->geo_global;
  up->global_ext = geometry_inp->geo_global_ext;
  up->grid = geometry_inp->geo_grid;

  struct gkyl_range nrange;
  struct gkyl_range nrange_quad;
  struct gkyl_range nrange_quad_surf[3];
  double dzc[3] = {0.0};

  int poly_order = up->basis.poly_order;

  // nodes tensor
  int num_nodes_corners[GKYL_MAX_CDIM];
  if (poly_order == 1) {
    for (int d=0; d<up->grid.ndim; ++d)
      num_nodes_corners[d] = gkyl_range_shape(&up->local, d) + 1;
  }
  if (poly_order == 2) {
    for (int d=0; d<up->grid.ndim; ++d)
      num_nodes_corners[d] = 2*gkyl_range_shape(&up->local, d) + 1;
  }

  int num_quad_points = poly_order+1;

  int num_nodes_quad_interior[GKYL_MAX_CDIM];
  for (int d=0; d<up->grid.ndim; ++d)
    num_nodes_quad_interior[d] = gkyl_range_shape(&up->local, d)*num_quad_points;

  int num_nodes_quad_surf_in_dir[up->grid.ndim][GKYL_MAX_CDIM];
  for (int dir=0; dir<up->grid.ndim; ++dir)
    for (int d=0; d<up->grid.ndim; ++d)
      num_nodes_quad_surf_in_dir[dir][d] = d == dir ? gkyl_range_shape(&up->local, d)+1 : gkyl_range_shape(&up->local, d)*num_quad_points;

  gkyl_range_init_from_shape(&nrange, up->grid.ndim, num_nodes_corners);
  gkyl_range_init_from_shape(&nrange_quad, up->grid.ndim, num_nodes_quad_interior);
  for (int dir=0; dir<up->grid.ndim; ++dir)
    gkyl_range_init_from_shape(&nrange_quad_surf[dir], up->grid.ndim, num_nodes_quad_surf_in_dir[dir]);

  // Initialize surface basis
  gkyl_cart_modal_serendip(&up->surf_basis, up->grid.ndim-1, poly_order);
  up->num_surf_basis = up->surf_basis.num_basis;

  // Initialize tokamak geometry object from EFIT
  const struct gkyl_efit_inp inp = geometry_inp->efit_info;
  struct gkyl_tok_geo_grid_inp ginp = geometry_inp->tok_grid_info;
  ginp.cgrid = up->grid;
  ginp.cbasis = up->basis;
  struct gkyl_tok_geo *geo = gkyl_tok_geo_new(&inp, &ginp);
  up->geqdsk_sign_convention = geo->efit->sibry > geo->efit->simag ? 0 : 1;

  // Allocate nodal and modal arrays for corner, interior, and surface geo
  gk_geometry_corn_alloc_nodal(up, nrange);
  gk_geometry_corn_alloc_expansions(up);
  gk_geometry_int_alloc_nodal(up, nrange_quad);
  gk_geometry_int_alloc_expansions(up);
  for (int dir=0; dir<up->grid.ndim; ++dir) {
    gk_geometry_surf_alloc_nodal(up, dir, nrange_quad_surf[dir]);
    gk_geometry_surf_alloc_expansions(up, dir);
  }


  // calculate mapc2p in cylindrical coords at corner nodes for
  // getting cell coordinates (used only for plotting)
  gkyl_tok_geo_calc(up, &nrange, geo, &ginp, geometry_inp->position_map);
  // calculate mapc2p in cylindrical coords at interior nodes for
  // calculating geo quantity volume expansions 
  gkyl_tok_geo_calc_interior(up, &nrange_quad, dzc, geo, &ginp, geometry_inp->position_map);
  // calculate mapc2p in cylindrical coords at surfaces
  for (int dir = 0; dir <up->grid.ndim; dir++)
    gkyl_tok_geo_calc_surface(up, dir, &nrange_quad_surf[dir], dzc, geo, &ginp, geometry_inp->position_map);

  // calculate bmag at corner nodes 
  gkyl_calc_bmag *bcalculator = gkyl_calc_bmag_new(&up->basis, &geo->rzbasis, &up->grid, &geo->rzgrid, false);
  gkyl_calc_bmag_advance(bcalculator, &up->local, &up->local_ext, &up->global, &geo->rzlocal, &geo->rzlocal_ext, geo->efit->bmagzr, up->geo_corn.bmag, up->geo_corn.mc2p, false);
  // calculate bmag at interior nodes
  gkyl_calc_bmag_advance(bcalculator, &up->local, &up->local_ext, &up->global, &geo->rzlocal, &geo->rzlocal_ext, geo->efit->bmagzr, up->geo_int.bmag, up->geo_int.mc2p, true);
  gkyl_calc_bmag_release(bcalculator);
  // Convert bmag to nodal at interior nodes so we can use it to calculate dphidtheta
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&up->basis, &up->grid, false);
  gkyl_nodal_ops_m2n(n2m, &up->basis, &up->grid, &nrange_quad, &up->local, 1, up->geo_int.bmag_nodal, up->geo_int.bmag, true);
  gkyl_nodal_ops_release(n2m);

  // Now calculate the metrics at interior nodes
  struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(&up->basis, &up->grid, &up->global, &up->global_ext, &up->local, &up->local_ext, false);
  gkyl_calc_metric_advance_rz_interior(mcalc, &nrange_quad, up->geo_int.mc2p_nodal_fd, up->geo_int.ddtheta_nodal, up->geo_int.bmag_nodal, dzc, up->geo_int.g_ij, up->geo_int.dxdz, up->geo_int.dzdx, up->geo_int.dualmag, up->geo_int.normals, up->geo_int.jacobgeo, up->geo_int.bcart, &up->local);
  gkyl_array_copy(up->geo_int.jacobgeo_ghost, up->geo_int.jacobgeo);
  // Calculate neutral metrics at interior nodes
  gkyl_calc_metric_advance_rz_neut_interior(mcalc, &nrange_quad, up->geo_int.mc2p_nodal_fd, up->geo_int.ddtheta_nodal, dzc, up->geo_int.g_ij_neut, up->geo_int.gij_neut, &up->local);
  // calculate the derived geometric quantities at interior nodes
  gkyl_tok_calc_derived_geo *jcalculator = gkyl_tok_calc_derived_geo_new(&up->basis, &up->grid, 1, false);
  gkyl_tok_calc_derived_geo_advance(jcalculator, &up->local, up->geo_int.g_ij, up->geo_int.bmag, 
    up->geo_int.jacobgeo, up->geo_int.jacobgeo_inv, up->geo_int.gij, up->geo_int.b_i, up->geo_int.cmag, up->geo_int.jacobtot, up->geo_int.jacobtot_inv, 
    up->geo_int.bmag_inv, up->geo_int.bmag_inv_sq, up->geo_int.gxxj, up->geo_int.gxyj, up->geo_int.gyyj, up->geo_int.gxzj, up->geo_int.eps2);
  gkyl_tok_calc_derived_geo_release(jcalculator);
  // Calculate metrics/derived geo quantities at surface
  for (int dir = 0; dir <up->grid.ndim; dir++) {
    gkyl_calc_metric_advance_rz_surface(mcalc, dir, &nrange_quad_surf[dir], up->geo_surf[dir].mc2p_nodal_fd, up->geo_surf[dir].ddtheta_nodal, up->geo_surf[dir].bmag_nodal, dzc,
      up->geo_surf[dir].jacobgeo_nodal, up->geo_surf[dir].b_i_nodal, up->geo_surf[dir].cmag_nodal, up->geo_surf[dir].jacobtot_inv_nodal, &up->local);
  }
  gkyl_calc_metric_release(mcalc);
  // Calculate surface expansions
  for (int dir = 0; dir <up->grid.ndim; dir++)
    gk_geometry_surf_calc_expansions(up, dir, nrange_quad_surf[dir]);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself

  gkyl_tok_geo_release(geo);
  // Release nodal data
  gk_geometry_corn_release_nodal(up);
  gk_geometry_int_release_nodal(up);
  for (int dir=0; dir<up->grid.ndim; ++dir)
    gk_geometry_surf_release_nodal(up, dir);

  return up;
}

struct gk_geometry*
gkyl_gk_geometry_tok_new(struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry* gk_geom_3d;
  struct gk_geometry* gk_geom;

  if (geometry_inp->position_map == 0){
    geometry_inp->position_map = gkyl_position_map_new((struct gkyl_position_map_inp) {}, \
      geometry_inp->grid, geometry_inp->local, geometry_inp->local_ext, geometry_inp->local, \
      geometry_inp->local_ext, geometry_inp->basis);
    gk_geom_3d = gk_geometry_tok_init(geometry_inp);
    gkyl_position_map_release(geometry_inp->position_map);
  }
  else {
    // First construct the uniform 3d geometry
    gk_geom_3d = gk_geometry_tok_init(geometry_inp);
    if (geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_POLYNOMIAL || \
        geometry_inp->position_map->id == GKYL_PMAP_CONSTANT_DB_NUMERIC) {
      // The array mc2nu is computed using the uniform geometry, so we need to deflate it
      // Must deflate the 3D uniform geometry in order for the allgather to work
      if(geometry_inp->grid.ndim < 3)
        gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, geometry_inp);
      else
        gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);

      gkyl_position_map_set_bmag(geometry_inp->position_map, geometry_inp->comm, \
        gk_geom->geo_int.bmag);

      gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry
      gkyl_gk_geometry_release(gk_geom); // release 3d geometry

      // Construct the non-uniform grid
      gk_geom_3d = gk_geometry_tok_init(geometry_inp);
    }
  }
  return gk_geom_3d;
}


void
gkyl_gk_geometry_tok_set_grid_extents(struct gkyl_efit_inp efit_info, struct gkyl_tok_geo_grid_inp grid_info, double *theta_lo, double *theta_up) {
  struct gkyl_tok_geo *geo = gkyl_tok_geo_new(&efit_info, &grid_info);
  gkyl_tok_geo_set_extent(&grid_info, geo, theta_lo, theta_up);
  gkyl_tok_geo_release(geo);
}
