#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_basis.h>
#include <gkyl_eval_on_nodes.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_gk_geometry_mirror.h>
#include <gkyl_math.h>
#include <gkyl_nodal_ops.h>

#include <gkyl_mirror_geo.h>
#include <gkyl_calc_derived_geo.h>
#include <gkyl_calc_metric.h>
#include <gkyl_calc_bmag.h>

struct gk_geometry*
gkyl_gk_geometry_mirror_new(struct gkyl_gk_geometry_inp *geometry_inp)
{

  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->basis = geometry_inp->geo_basis;
  up->local = geometry_inp->geo_local;
  up->local_ext = geometry_inp->geo_local_ext;
  up->global = geometry_inp->geo_global;
  up->global_ext = geometry_inp->geo_global_ext;
  up->grid = geometry_inp->geo_grid;

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

  // bmag, metrics and derived geo quantities
  up->bmag = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->g_ij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
  up->dxdz = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->dzdx = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->normals = gkyl_array_new(GKYL_DOUBLE, 9*up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->jacobgeo_inv = gkyl_array_new(GKYL_DOUBLE, up->basis.num_basis, up->local_ext.volume);
  up->gij = gkyl_array_new(GKYL_DOUBLE, 6*up->basis.num_basis, up->local_ext.volume);
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
  up->bmag_mid = gkyl_array_new(GKYL_DOUBLE, 1, 1);

  const struct gkyl_mirror_geo_efit_inp inp = geometry_inp->mirror_efit_info;
  struct gkyl_mirror_geo_grid_inp ginp = geometry_inp->mirror_grid_info;
  ginp.cgrid = up->grid;
  ginp.cbasis = up->basis;
  struct gkyl_mirror_geo *geo = gkyl_mirror_geo_new(&inp);
  // calculate mapc2p
  gkyl_mirror_geo_calc(up, &nrange, dzc, geo, &ginp,
    mc2p_nodal_fd, mc2p_nodal, up->mc2p);
  // calculate bmag
  gkyl_calc_bmag *bcalculator = gkyl_calc_bmag_new(&up->basis, &geo->rzbasis, &geo->fbasis, &up->grid, &geo->rzgrid, &geo->fgrid, geo->psisep, false);
  gkyl_calc_bmag_advance(bcalculator, &up->local, &up->local_ext, &up->global, &geo->rzlocal, &geo->rzlocal_ext, &geo->frange, &geo->frange_ext, geo->psiRZ, geo->psibyrRZ, geo->psibyr2RZ, up->bmag, geo->fpoldg, up->mc2p, false);
  gkyl_calc_bmag_release(bcalculator);
  // now calculate the metrics
  struct gkyl_calc_metric* mcalc = gkyl_calc_metric_new(&up->basis, &up->grid, &up->global, &up->global_ext, &up->local, &up->local_ext, false);
  gkyl_calc_metric_advance(mcalc, &nrange, mc2p_nodal_fd, dzc, up->g_ij, up->dxdz, up->dzdx, up->normals, &up->local);

  // calculate the derived geometric quantities
  struct gkyl_calc_derived_geo *jcalculator = gkyl_calc_derived_geo_new(&up->basis, &up->grid, false);
  gkyl_calc_derived_geo_advance(jcalculator, &up->local, up->g_ij, up->bmag, 
    up->jacobgeo, up->jacobgeo_inv, up->gij, up->b_i, up->cmag, up->jacobtot, up->jacobtot_inv, 
    up->bmag_inv, up->bmag_inv_sq, up->gxxj, up->gxyj, up->gyyj, up->gxzj, up->eps2);
  gkyl_calc_derived_geo_release(jcalculator);
  gkyl_calc_metric_advance_bcart(mcalc, &nrange, up->b_i, up->dzdx, up->bcart, &up->local);
  gkyl_calc_metric_release(mcalc);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself

  gkyl_mirror_geo_release(geo);
  gkyl_array_release(mc2p_nodal_fd);
  gkyl_array_release(mc2p_nodal);

  return up;
}

