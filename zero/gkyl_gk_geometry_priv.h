#include <gkyl_gk_geometry.h>
#include <gkyl_nodal_ops.h>


static double calc_running_coord(double coord_lo, int i, double dx) {
  double dels[2] = {1.0/sqrt(3), 1.0-1.0/sqrt(3) };
  double coord = coord_lo;
  for(int j = 0; j < i; j++)
    coord+=dels[j%2]*dx;
  return coord;
}


static double calc_running_surf_coord(double coord_lo, int i, double dx) {
  double dels[3] = {(1.0-1.0/sqrt(3))/2.0, 1.0/sqrt(3), (1.0-1.0/sqrt(3))/2.0 };
  double coord = coord_lo;
  for(int j = 0; j < i; j++)
    coord+=dels[j%3]*dx;
  return coord;
}

static void
gk_geometry_surf_alloc_nodal(struct gk_geometry* gk_geom, int dir, struct gkyl_range nrange)
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

}

static void
gk_geometry_surf_alloc_expansions(struct gk_geometry* up, int dir)
{
  up->geo_surf[dir].bmag = gkyl_array_new(GKYL_DOUBLE, 1*up->surf_basis.num_basis, up->local_ext.volume);
  up->geo_surf[dir].jacobgeo = gkyl_array_new(GKYL_DOUBLE, 1*up->surf_basis.num_basis, up->local_ext.volume);
  up->geo_surf[dir].jacobgeo_sync = gkyl_array_new(GKYL_DOUBLE, 1*up->surf_basis.num_basis, up->local_ext.volume);
  up->geo_surf[dir].b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->surf_basis.num_basis, up->local_ext.volume);
  up->geo_surf[dir].cmag = gkyl_array_new(GKYL_DOUBLE, 1*up->surf_basis.num_basis, up->local_ext.volume);
  up->geo_surf[dir].jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, 1*up->surf_basis.num_basis, up->local_ext.volume);
}



static void
gk_geometry_surf_release_nodal(struct gk_geometry* gk_geom, int dir)
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

static void
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


