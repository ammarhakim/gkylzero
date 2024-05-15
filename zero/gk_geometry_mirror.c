#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_comm.h>
#include <gkyl_deflate_geo.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mirror.h>
#include <gkyl_math.h>
#include <gkyl_nodal_ops.h>

#include <gkyl_mirror_geo.h>
#include <gkyl_calc_derived_geo.h>
#include <gkyl_calc_metric.h>
#include <gkyl_calc_bmag.h>


// write out nodal coordinates 
static void
write_nodal_coordinates(const char *nm, struct gkyl_range *nrange,
  struct gkyl_array *nodes)
{
  double lower[3] = { 0.0, 0.0, 0.0 };
  double upper[3] = { 1.0, 1.0, 1.0 };
  int cells[3];
  for (int i=0; i<nrange->ndim; ++i)
    cells[i] = gkyl_range_shape(nrange, i);
  
  struct gkyl_rect_grid grid;
  gkyl_rect_grid_init(&grid, 3, lower, upper, cells);

  gkyl_grid_sub_array_write(&grid, nrange, 0, nodes, nm);
}

struct gk_geometry*
gkyl_gk_geometry_mirror_new(struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry* gk_geom_3d;
  struct gk_geometry* gk_geom;
  gk_geom_3d = gkyl_gk_geometry_mirror_advance(geometry_inp);
  if(geometry_inp->cdim < 3) {
    gkyl_gk_geometry_c2fa_deflate(gk_geom_3d, geometry_inp);
  } else {
    gkyl_gk_geometry_c2fa_acquire(gk_geom_3d, geometry_inp);
  }
  double nonuniform_frac = geometry_inp->nonuniform_map_fraction;
  if (nonuniform_frac > 0.0 & nonuniform_frac <= 1.0) {
    // Copy deflate geometry if necessary
    if(geometry_inp->cdim < 3)
      gk_geom = gkyl_gk_geometry_deflate(gk_geom_3d, geometry_inp);
    else
      gk_geom = gkyl_gk_geometry_acquire(gk_geom_3d);
    struct gkyl_array *bmag_global = gkyl_array_new(GKYL_DOUBLE, geometry_inp->basis.num_basis, geometry_inp->global_ext.volume);
    if (geometry_inp->use_gpu) { // Allgather is only a GPU operation, so we must copy these arrays to GPU, then back to CPU
      struct gkyl_array *bmag_global_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, geometry_inp->basis.num_basis, geometry_inp->global_ext.volume);
      struct gkyl_array *bmag_dev = gkyl_array_cu_dev_new(GKYL_DOUBLE, geometry_inp->basis.num_basis, geometry_inp->local_ext.volume);
      gkyl_array_copy(bmag_dev, gk_geom->bmag);
      gkyl_comm_array_allgather(geometry_inp->comm, &geometry_inp->local, &geometry_inp->global, bmag_dev, bmag_global_dev);
      gkyl_array_copy(bmag_global, bmag_global_dev);
      gkyl_array_release(bmag_global_dev);
      gkyl_array_release(bmag_dev);
    }
    else {
      gkyl_comm_array_allgather(geometry_inp->comm, &geometry_inp->local, &geometry_inp->global, gk_geom->bmag, bmag_global);
    }
    geometry_inp->nonuniform_geom = true;
    geometry_inp->bmag_global = bmag_global;
    struct gkyl_mirror_geo_c2fa_ctx *c2fa_app = geometry_inp->mirror_geo_c2fa_ctx;
    gkyl_gk_geometry_release(gk_geom_3d); // release temporary 3d geometry
    gkyl_gk_geometry_release(gk_geom); // release 3d geometry
    gkyl_array_release(c2fa_app->c2fa);
    gkyl_array_release(c2fa_app->c2fa_deflate);
    gk_geom_3d = gkyl_gk_geometry_mirror_advance(geometry_inp);
    if(geometry_inp->cdim < 3) {
      gkyl_gk_geometry_c2fa_deflate(gk_geom_3d, geometry_inp);
    } else {
      gkyl_gk_geometry_c2fa_acquire(gk_geom_3d, geometry_inp);
    }
    gkyl_array_release(bmag_global);
  }
  else if (nonuniform_frac != 0.0) {
    printf("Invalid non-uniform mapping fraction %f. Must be between 0 and 1", nonuniform_frac);
  }
  return gk_geom_3d;
}


struct gk_geometry*
gkyl_gk_geometry_mirror_advance(struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->basis = geometry_inp->geo_basis;
  up->local = geometry_inp->geo_local;
  up->local_ext = geometry_inp->geo_local_ext;
  up->global = geometry_inp->geo_global;
  up->global_ext = geometry_inp->geo_global_ext;
  up->grid = geometry_inp->geo_grid;
  up->bmag_global = geometry_inp->bmag_global;
  up->decomp_basis = geometry_inp->basis;
  up->decomp_local = geometry_inp->local;
  up->decomp_local_ext = geometry_inp->local_ext;
  up->decomp_global = geometry_inp->global;
  up->decomp_global_ext = geometry_inp->global_ext;
  up->decomp_grid = geometry_inp->grid;

  struct gkyl_range nrange;
  double dzc[3] = {0.0};

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
  int num_fd_nodes = 13;
  struct gkyl_array* mc2p_nodal_fd = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*num_fd_nodes, nrange.volume);
  struct gkyl_array* mc2p_nodal = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim, nrange.volume);
  up->mc2p = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);

  struct gkyl_mirror_geo_c2fa_ctx *c2fa_app = geometry_inp->mirror_geo_c2fa_ctx;
  struct gkyl_array* map_c2fa_nodal_fd = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*num_fd_nodes, nrange.volume);
  struct gkyl_array* map_c2fa_nodal = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim, nrange.volume);
  c2fa_app->c2fa = gkyl_array_new(GKYL_DOUBLE, up->grid.ndim*up->basis.num_basis, up->local_ext.volume);

  // bmag, metrics and derived geo quantities
  up->bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->dxdz = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->dzdx = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->b_i = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
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
  up->bmag_mid = gkyl_array_new(GKYL_DOUBLE, 1, 1);

  const struct gkyl_mirror_geo_efit_inp *inp = geometry_inp->mirror_efit_info;
  struct gkyl_mirror_geo_grid_inp *ginp = geometry_inp->mirror_grid_info;
  bool nonuniform_geom = geometry_inp->nonuniform_geom;
  ginp->nonuniform_mapping_fraction = geometry_inp->nonuniform_map_fraction;
  ginp->cgrid = up->grid;
  ginp->cbasis = up->basis;
  
  struct gkyl_mirror_geo *geo = gkyl_mirror_geo_new(inp);
  // calculate mapc2p
  gkyl_mirror_geo_calc(up, &nrange, dzc, NULL, geo, NULL, ginp, mc2p_nodal_fd, mc2p_nodal, up->mc2p,
    nonuniform_geom, map_c2fa_nodal_fd, map_c2fa_nodal, c2fa_app->c2fa);
  // calculate bmag
  gkyl_calc_bmag *bcalculator = gkyl_calc_bmag_new(&up->basis, &geo->rzbasis, &geo->fbasis, &up->grid, &geo->rzgrid, &geo->fgrid, geo->psisep, false);
  gkyl_calc_bmag_advance(bcalculator, &up->local, &up->local_ext, &up->global, &geo->rzlocal, &geo->rzlocal_ext, &geo->frange, &geo->frange_ext,
    geo->psiRZ, geo->psibyrRZ, geo->psibyr2RZ, up->bmag, geo->fpoldg, up->mc2p, false);
  gkyl_calc_bmag_release(bcalculator);
  // now calculate the metrics
  struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(&up->basis, &up->grid, &up->global, &up->global_ext, &up->local, &up->local_ext, false);
  gkyl_calc_metric_advance(mcalc, &nrange, mc2p_nodal_fd, dzc, up->g_ij, up->dxdz, up->dzdx, &up->local);
  gkyl_calc_metric_release(mcalc);

  // calculate the derived geometric quantities
  struct gkyl_calc_derived_geo *jcalculator = gkyl_calc_derived_geo_new(&up->basis, &up->grid, false);
  gkyl_calc_derived_geo_advance(jcalculator, &up->local, up->g_ij, up->bmag, 
    up->jacobgeo, up->jacobgeo_inv, up->gij, up->b_i, up->cmag, up->jacobtot, up->jacobtot_inv, 
    up->bmag_inv, up->bmag_inv_sq, up->gxxj, up->gxyj, up->gyyj, up->gxzj, up->eps2);
  gkyl_calc_derived_geo_release(jcalculator);

  c2fa_app->grid = geo->rzgrid;
  c2fa_app->cgrid = up->grid;
  c2fa_app->range = geo->rzlocal;
  c2fa_app->crange = up->local;
  c2fa_app->crange_global = up->global;
  c2fa_app->basis = geo->rzbasis;
  c2fa_app->cbasis = up->basis;

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself

  gkyl_mirror_geo_release(geo);
  gkyl_array_release(mc2p_nodal_fd);
  gkyl_array_release(mc2p_nodal);

  gkyl_array_release(map_c2fa_nodal_fd);
  gkyl_array_release(map_c2fa_nodal);

  return up;
}

void
gkyl_gk_geometry_c2fa_deflate(const struct gk_geometry* up_3d, struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  struct gkyl_mirror_geo_c2fa_ctx *c2fa_app = geometry_inp->mirror_geo_c2fa_ctx;
  up->basis = geometry_inp->basis;
  up->local = geometry_inp->local;
  up->local_ext = geometry_inp->local_ext;
  up->grid = geometry_inp->grid;
  c2fa_app->c2fa_deflate = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
  // Now fill the arrays by deflation
  int rem_dirs[3] = {0};
  if (up->grid.ndim==1) {
    rem_dirs[0] = 1;
    rem_dirs[1] = 1;
  }
  else if (up->grid.ndim==2) {
    rem_dirs[1] = 1;
  }
  struct gkyl_deflate_geo* deflator = gkyl_deflate_geo_new(&up_3d->basis, &up->basis, &up_3d->grid, &up->grid, rem_dirs, false);
  c2fa_app->grid_deflate = geometry_inp->grid;
  c2fa_app->basis_deflate = geometry_inp->basis;
  c2fa_app->range_deflate = geometry_inp->local;
  c2fa_app->range_global_deflate = geometry_inp->global;

  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, c2fa_app->c2fa, c2fa_app->c2fa_deflate, 3);
  // Done deflating
  gkyl_deflate_geo_release(deflator);
  free(up);
}

void
gkyl_gk_geometry_c2fa_acquire(const struct gk_geometry* up_3d, struct gkyl_gk_geometry_inp *geometry_inp)
{
  struct gkyl_mirror_geo_c2fa_ctx *c2fa_app = geometry_inp->mirror_geo_c2fa_ctx;
  c2fa_app->c2fa_deflate = gkyl_array_new(GKYL_DOUBLE, 3*geometry_inp->basis.num_basis, geometry_inp->local_ext.volume);
  gkyl_array_copy(c2fa_app->c2fa_deflate, c2fa_app->c2fa);
  c2fa_app->grid_deflate = geometry_inp->grid;
  c2fa_app->basis_deflate = geometry_inp->basis;
  c2fa_app->range_deflate = geometry_inp->local;
  c2fa_app->range_global_deflate = geometry_inp->global;
}
