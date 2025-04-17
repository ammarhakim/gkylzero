#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_comm.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_tok.h>
#include <gkyl_math.h>
#include <gkyl_nodal_ops.h>

#include <gkyl_tok_geo.h>
#include <gkyl_tok_calc_derived_geo.h>
#include <gkyl_calc_metric.h>
#include <gkyl_calc_bmag.h>

void
gk_geometry_tok_surf_alloc(struct gk_geometry* gk_geom, int dir, struct gkyl_range nrange)
{

  int num_surf_quad = gk_geom->surf_basis.num_basis;
  // mapc2p for calculations of tangents
  int num_fd_nodes = 13;
  gk_geom->geo_surf[dir].mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, gk_geom->grid.ndim*num_fd_nodes, nrange.volume);
  gk_geom->geo_surf[dir].mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, gk_geom->grid.ndim, nrange.volume);
  // bmag.metrics and derived geo quantities
  gk_geom->geo_surf[dir].bmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);
  gk_geom->geo_surf[dir].ddtheta_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
  gk_geom->geo_surf[dir].jacobgeo_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);
  gk_geom->geo_surf[dir].b_i_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
  gk_geom->geo_surf[dir].cmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);
  gk_geom->geo_surf[dir].jacobtot_inv_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);

  gk_geom->geo_surf[dir].bmag = gkyl_array_new(GKYL_DOUBLE, 1*num_surf_quad, gk_geom->local_ext.volume);
  gk_geom->geo_surf[dir].jacobgeo = gkyl_array_new(GKYL_DOUBLE, 1*num_surf_quad, gk_geom->local_ext.volume);
  gk_geom->geo_surf[dir].jacobgeo_sync = gkyl_array_new(GKYL_DOUBLE, 1*num_surf_quad, gk_geom->local_ext.volume);
  gk_geom->geo_surf[dir].b_i = gkyl_array_new(GKYL_DOUBLE, 3*num_surf_quad, gk_geom->local_ext.volume);
  gk_geom->geo_surf[dir].cmag = gkyl_array_new(GKYL_DOUBLE, 1*num_surf_quad, gk_geom->local_ext.volume);
  gk_geom->geo_surf[dir].jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, 1*num_surf_quad, gk_geom->local_ext.volume);
}


void
gk_geometry_tok_surf_release_nodal(struct gk_geometry* gk_geom, int dir)
{
  gkyl_array_release(gk_geom->geo_surf[dir].mc2p_nodal_fd);
  gkyl_array_release(gk_geom->geo_surf[dir].mc2p_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].bmag_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].ddtheta_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].jacobgeo_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].b_i_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].cmag_nodal);
  gkyl_array_release(gk_geom->geo_surf[dir].jacobtot_inv_nodal);
}

void
gk_geometry_surf_calc_expansions(struct gk_geometry* gk_geom, int dir, 
  struct gkyl_range nrange_quad_surf)
{
  struct gk_geom_surf up_surf = gk_geom->geo_surf[dir];
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&gk_geom->basis, &gk_geom->grid, false);

  struct gkyl_range local_ext_in_dir;
  int lower[3] = {gk_geom->local.lower[0], gk_geom->local.lower[1], gk_geom->local.lower[2]};
  int upper[3] = {gk_geom->local.upper[0], gk_geom->local.upper[1], gk_geom->local.upper[2]};
  upper[dir]+=1;
  gkyl_sub_range_init(&local_ext_in_dir, &gk_geom->local_ext, lower, upper);


  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 1, up_surf.bmag_nodal, up_surf.bmag, dir);
  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 1, up_surf.jacobgeo_nodal, up_surf.jacobgeo, dir);
  gkyl_array_copy(up_surf.jacobgeo_sync, up_surf.jacobgeo);
  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 3, up_surf.b_i_nodal, up_surf.b_i, dir);
  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 1, up_surf.cmag_nodal, up_surf.cmag, dir);
  gkyl_nodal_ops_n2m_surface(n2m, &gk_geom->surf_basis, &gk_geom->grid, &nrange_quad_surf, &local_ext_in_dir, 1, up_surf.jacobtot_inv_nodal, up_surf.jacobtot_inv, dir);
  gkyl_nodal_ops_release(n2m);
}

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

  // Initialize surface basis abd allocate surface geo
  gkyl_cart_modal_serendip(&up->surf_basis, up->grid.ndim-1, poly_order);
  for (int dir=0; dir<up->grid.ndim; ++dir) {
    gk_geometry_tok_surf_alloc(up, dir, nrange_quad_surf[dir]);
  }

  int num_fd_nodes = 13;
  struct gkyl_array* mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*num_fd_nodes, nrange_quad.volume);
  struct gkyl_array* mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim, nrange.volume);
  struct gkyl_array* mc2p_nodal_quad = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim, nrange_quad.volume);
  struct gkyl_array *mc2p_quad = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
  up->mc2p = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);

  struct gkyl_array* map_mc2nu_nodal = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim, nrange_quad.volume);
  up->mc2nu_pos = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);

  struct gkyl_array* ddtheta_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange_quad.volume);
  struct gkyl_array* bmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange_quad.volume);

  // bmag, metrics and derived geo quantities
  up->bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->g_ij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->dxdz = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->dzdx = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->dualmag = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->normals = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo_ghost = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
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

  const struct gkyl_efit_inp inp = geometry_inp->efit_info;
  struct gkyl_tok_geo_grid_inp ginp = geometry_inp->tok_grid_info;
  ginp.cgrid = up->grid;
  ginp.cbasis = up->basis;
  struct gkyl_tok_geo *geo = gkyl_tok_geo_new(&inp, &ginp);
  up->geqdsk_sign_convention = geo->efit->sibry > geo->efit->simag ? 0 : 1;


  // calculate mapc2p in cylindrical coords at corner nodes for
  // getting cell coordinates (used only for plotting)
  gkyl_tok_geo_calc(up, &nrange, geo, &ginp, mc2p_nodal, up->mc2p, 
    map_mc2nu_nodal, up->mc2nu_pos, geometry_inp->position_map);
  // calculate mapc2p in cylindrical coords at interior nodes for
  // calculating geo quantity volume expansions 
  gkyl_tok_geo_calc_interior(up, &nrange_quad, dzc, geo, &ginp, mc2p_nodal_quad, mc2p_quad, mc2p_nodal_fd, 
    ddtheta_nodal, geometry_inp->position_map);
  // calculate mapc2p in cylindrical coords at surfaces
  for (int dir = 0; dir <up->grid.ndim; dir++) {
    gkyl_tok_geo_calc_surface(up, dir, &nrange_quad_surf[dir], dzc, geo, &ginp, up->geo_surf[dir].mc2p_nodal,
      up->geo_surf[dir].mc2p_nodal_fd, up->geo_surf[dir].ddtheta_nodal, up->geo_surf[dir].bmag_nodal, geometry_inp->position_map);
  }

  // calculate bmag at interior nodes
  gkyl_calc_bmag *bcalculator = gkyl_calc_bmag_new(&up->basis, &geo->rzbasis, &up->grid, &geo->rzgrid, false);
  gkyl_calc_bmag_advance(bcalculator, &up->local, &up->local_ext, &up->global, &geo->rzlocal, &geo->rzlocal_ext, geo->efit->bmagzr, up->bmag, mc2p_quad, true);
  gkyl_calc_bmag_release(bcalculator);
  // Convert bmag to nodal so we can use it to calculate dphidtheta
  struct gkyl_nodal_ops *n2m = gkyl_nodal_ops_new(&up->basis, &up->grid, false);
  gkyl_nodal_ops_m2n(n2m, &up->basis, &up->grid, &nrange, &up->local, 1, bmag_nodal, up->bmag, true);
  gkyl_nodal_ops_release(n2m);

  // Now calculate the metrics
  struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(&up->basis, &up->grid, &up->global, &up->global_ext, &up->local, &up->local_ext, false);
  gkyl_calc_metric_advance_rz_interior(mcalc, &nrange_quad, mc2p_nodal_fd, ddtheta_nodal, bmag_nodal, dzc, up->g_ij, up->dxdz, up->dzdx, up->dualmag, up->normals, up->jacobgeo, up->bcart, &up->local);
  gkyl_array_copy(up->jacobgeo_ghost, up->jacobgeo);
  // Calculate neutral metrics
  gkyl_calc_metric_advance_rz_neut(mcalc, &nrange, mc2p_nodal_fd, ddtheta_nodal, dzc, up->g_ij_neut, up->gij_neut, &up->local);
  // calculate the derived geometric quantities
  gkyl_tok_calc_derived_geo *jcalculator = gkyl_tok_calc_derived_geo_new(&up->basis, &up->grid, 1, false);
  gkyl_tok_calc_derived_geo_advance(jcalculator, &up->local, up->g_ij, up->bmag, 
    up->jacobgeo, up->jacobgeo_inv, up->gij, up->b_i, up->cmag, up->jacobtot, up->jacobtot_inv, 
    up->bmag_inv, up->bmag_inv_sq, up->gxxj, up->gxyj, up->gyyj, up->gxzj, up->eps2);
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
  // Release interior nodal data
  gkyl_array_release(map_mc2nu_nodal);
  gkyl_array_release(mc2p_nodal_fd);
  gkyl_array_release(mc2p_nodal);
  gkyl_array_release(mc2p_nodal_quad);
  gkyl_array_release(mc2p_quad);
  gkyl_array_release(ddtheta_nodal);
  gkyl_array_release(bmag_nodal);
  // Release surface nodal data
  for (int dir=0; dir<up->grid.ndim; ++dir)
    gk_geometry_tok_surf_release_nodal(up, dir);

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

      geometry_inp->position_map->to_optimize = true;
      gkyl_comm_array_allgather_host(geometry_inp->comm, &geometry_inp->local, \
      &geometry_inp->global, gk_geom->bmag, (struct gkyl_array*) geometry_inp->position_map->bmag_ctx->bmag);

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
