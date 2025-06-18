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
  // mapc2p for calculations of tangents
  int num_fd_nodes = 13;
  gk_geom->geo_surf[dir].mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, 3*num_fd_nodes, nrange.volume);
  gk_geom->geo_surf[dir].mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
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
  up->geo_surf[dir].bmag = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].jacobgeo = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].jacobgeo_sync = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].cmag = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
  up->geo_surf[dir].jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, 1*up->num_surf_basis, up->local_ext.volume);
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
gk_geometry_int_alloc_nodal(struct gk_geometry* gk_geom, struct gkyl_range nrange)
{
  // mapc2p for calculations of tangents
  int num_fd_nodes = 13;
  gk_geom->geo_int.mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, 3*num_fd_nodes, nrange.volume);
  gk_geom->geo_int.mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
  // bmag.metrics and derived geo quantities
  gk_geom->geo_int.bmag_nodal = gkyl_array_new(GKYL_DOUBLE, 1, nrange.volume);
  gk_geom->geo_int.ddtheta_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
}

static void
gk_geometry_int_alloc_expansions(struct gk_geometry* up)
{
  // mapc2p
  up->geo_int.mc2p = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  // bmag, metrics and derived geo quantities
  up->geo_int.bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.g_ij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.dxdz = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.dzdx = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.dualmag = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.normals = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.jacobgeo_ghost = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.gij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.gij_neut = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.bcart = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->geo_int.cmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.jacobtot = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.jacobtot_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.bmag_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.bmag_inv_sq = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.gxxj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.gxyj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.gyyj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.gxzj= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->geo_int.eps2= gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
}

static void
gk_geometry_int_release_nodal(struct gk_geometry* gk_geom)
{
  gkyl_array_release(gk_geom->geo_int.mc2p_nodal_fd);
  gkyl_array_release(gk_geom->geo_int.mc2p_nodal);
  gkyl_array_release(gk_geom->geo_int.bmag_nodal);
  gkyl_array_release(gk_geom->geo_int.ddtheta_nodal);
}

static void
gk_geometry_corn_alloc_nodal(struct gk_geometry* gk_geom, struct gkyl_range nrange)
{
  // mapc2p
  gk_geom->geo_corn.mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
  gk_geom->geo_corn.mc2nu_pos_nodal = gkyl_array_new(GKYL_DOUBLE, 3, nrange.volume);
}

static void
gk_geometry_corn_alloc_expansions(struct gk_geometry* up)
{
  // mapc2p
  up->geo_corn.mc2p = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  up->geo_corn.mc2nu_pos = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  // bmag
  up->geo_corn.bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  // deflated quantities for plotting
  up->geo_corn.mc2p_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
  up->geo_corn.mc2nu_pos_deflated = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);
}

static void
gk_geometry_corn_release_nodal(struct gk_geometry* gk_geom)
{
  gkyl_array_release(gk_geom->geo_corn.mc2p_nodal);
  gkyl_array_release(gk_geom->geo_corn.mc2nu_pos_nodal);
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


