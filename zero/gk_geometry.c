#include <gkyl_range.h>
#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_rio.h>
#include <gkyl_math.h>
#include <gkyl_basis.h>
#include <gkyl_deflate_geo.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_alloc_flags_priv.h>

struct gk_geometry*
gkyl_gk_geometry_new(struct gk_geometry* geo_host, struct gkyl_gk_geometry_inp *geometry_inp, bool use_gpu)
{

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_gk_geometry_cu_dev_new(geo_host, geometry_inp);
  } 
#endif 

  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->basis = geometry_inp->basis;
  up->local = geometry_inp->local;
  up->local_ext = geometry_inp->local_ext;
  up->global = geometry_inp->global;
  up->global_ext = geometry_inp->global_ext;
  up->grid = geometry_inp->grid;

  // bmag, metrics and derived geo quantities
  up->mc2p = gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
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

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself
                   
  return up;
}

bool
gkyl_gk_geometry_is_cu_dev(const struct gk_geometry* up)
{
  return GKYL_IS_CU_ALLOC(up->flags);
}

struct gkyl_rect_grid gkyl_gk_geometry_augment_grid(struct gkyl_rect_grid grid, struct gkyl_gk_geometry_inp geometry)
{
  struct gkyl_rect_grid augmented_grid;
  int cells[3];
  double lower[3];
  double upper[3];

  if (grid.ndim==1) {
    cells[0] = 1;
    cells[1] = 1;
    cells[2] = grid.cells[0];

    lower[0] = geometry.world[0] - 1e-5;
    lower[1] = geometry.world[1] - 1e-1;
    lower[2] = grid.lower[0];

    upper[0] = geometry.world[0] + 1e-5;
    upper[1] = geometry.world[1] + 1e-1;
    upper[2] = grid.upper[0];
  }
  else if (grid.ndim==2) {
    cells[0] = grid.cells[0];
    cells[1] = 1;
    cells[2] = grid.cells[1];

    lower[0] = grid.lower[0];
    lower[1] = geometry.world[0] - 1e-1;
    lower[2] = grid.lower[1];

    upper[0] = grid.upper[0];
    upper[1] = geometry.world[0] + 1e-1;
    upper[2] = grid.upper[1];
  }

  gkyl_rect_grid_init(&augmented_grid, 3, lower, upper, cells);
  return augmented_grid;
}


void gkyl_gk_geometry_augment_local(const struct gkyl_range *inrange,
  const int *nghost, struct gkyl_range *ext_range, struct gkyl_range *range)
{
  if (inrange->ndim == 2) {
    int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
    
    lower_ext[0] = inrange->lower[0]-nghost[0];
    upper_ext[0] = inrange->upper[0]+nghost[0];
    lower[0] = inrange->lower[0];
    upper[0] = inrange->upper[0];

    lower_ext[1] = 1 - 1;
    upper_ext[1] = 1 + 1;
    lower[1] = 1;
    upper[1] = 1;

    lower_ext[2] = inrange->lower[1]-nghost[1];
    upper_ext[2] = inrange->upper[1]+nghost[1];
    lower[2] = inrange->lower[1];
    upper[2] = inrange->upper[1];


    gkyl_range_init(ext_range, inrange->ndim+1, lower_ext, upper_ext);
    gkyl_sub_range_init(range, ext_range, lower, upper);  
  }
  else if (inrange->ndim == 1) {
    int lower_ext[GKYL_MAX_DIM], upper_ext[GKYL_MAX_DIM];
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
    
    lower_ext[0] = 1 - 1;
    upper_ext[0] = 1 + 1;
    lower[0] = 1;
    upper[0] = 1;

    lower_ext[1] = 1 - 1;
    upper_ext[1] = 1 + 1;
    lower[1] = 1;
    upper[1] = 1;

    lower_ext[2] = inrange->lower[0]-nghost[0];
    upper_ext[2] = inrange->upper[0]+nghost[0];
    lower[2] = inrange->lower[0];
    upper[2] = inrange->upper[0];



    gkyl_range_init(ext_range, inrange->ndim+2, lower_ext, upper_ext);
    gkyl_sub_range_init(range, ext_range, lower, upper);  
  }
}

void gkyl_gk_geometry_bmag_mid(struct gk_geometry* up) {
  int cdim = up->grid.ndim;
  int idx_mid[cdim];
  double xc[cdim];
  for(int i = 0; i <cdim; i++) {
    idx_mid[i] = up->grid.cells[i]/2+1;
    xc[i] = up->grid.cells[i]%2 == 0? -1.0 : 0.0;
  }

  double bmag_mid = 0.0;
  if (gkyl_range_contains_idx(&up->local, idx_mid)) {
    long lidx = gkyl_range_idx(&up->local, idx_mid);
    const double *bcoeffs = gkyl_array_cfetch(up->bmag, lidx);
    double *bmag_mid = gkyl_array_fetch(up->bmag_mid, 0);
    bmag_mid[0] = up->basis.eval_expand(xc, bcoeffs);
  }
}

struct gk_geometry*
gkyl_gk_geometry_deflate(const struct gk_geometry* up_3d, struct gkyl_gk_geometry_inp *geometry_inp)
{

  struct gk_geometry *up = gkyl_malloc(sizeof(struct gk_geometry));
  up->basis = geometry_inp->basis;
  up->local = geometry_inp->local;
  up->local_ext = geometry_inp->local_ext;
  up->grid = geometry_inp->grid;

  // bmag, metrics and derived geo quantities
  up->mc2p= gkyl_array_new(GKYL_DOUBLE, 3*up->basis.num_basis, up->local_ext.volume);
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

  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->mc2p, up->mc2p, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->bmag, up->bmag, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->g_ij, up->g_ij, 6);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->dxdz, up->dxdz, 9);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->dzdx, up->dzdx, 9);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->jacobgeo, up->jacobgeo, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->jacobgeo_inv, up->jacobgeo_inv, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->gij, up->gij, 6);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->b_i, up->b_i, 3);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->cmag, up->cmag, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->jacobtot, up->jacobtot, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->jacobtot_inv, up->jacobtot_inv, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->bmag_inv, up->bmag_inv, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->bmag_inv_sq, up->bmag_inv_sq, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->gxxj, up->gxxj, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->gxyj, up->gxyj, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->gyyj, up->gyyj, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->gxzj, up->gxzj, 1);
  gkyl_deflate_geo_advance(deflator, &up_3d->local, &up->local, up_3d->eps2, up->eps2, 1);
  // Done deflating
  gkyl_deflate_geo_release(deflator);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->ref_count = gkyl_ref_count_init(gkyl_gk_geometry_free);
  up->on_dev = up; // CPU eqn obj points to itself

  return up;
}

void
gkyl_gk_geometry_free(const struct gkyl_ref_count *ref)
{
  struct gk_geometry *up = container_of(ref, struct gk_geometry, ref_count);
  gkyl_array_release(up->mc2p);
  gkyl_array_release(up->bmag);
  gkyl_array_release(up->g_ij);
  gkyl_array_release(up->jacobgeo);
  gkyl_array_release(up->jacobgeo_inv);
  gkyl_array_release(up->dxdz);
  gkyl_array_release(up->dzdx);
  gkyl_array_release(up->gij);
  gkyl_array_release(up->b_i);
  gkyl_array_release(up->cmag);
  gkyl_array_release(up->jacobtot);
  gkyl_array_release(up->jacobtot_inv);
  gkyl_array_release(up->bmag_inv);
  gkyl_array_release(up->bmag_inv_sq);
  gkyl_array_release(up->gxxj);
  gkyl_array_release(up->gxyj);
  gkyl_array_release(up->gyyj);
  gkyl_array_release(up->gxzj);
  gkyl_array_release(up->eps2);
  gkyl_array_release(up->bmag_mid);
  if (gkyl_gk_geometry_is_cu_dev(up)) 
    gkyl_cu_free(up->on_dev); 

  gkyl_free(up);
}

struct gk_geometry*
gkyl_gk_geometry_acquire(const struct gk_geometry* up)
{
  gkyl_ref_count_inc(&up->ref_count);
  return (struct gk_geometry*) up;
}

void
gkyl_gk_geometry_release(const struct gk_geometry *up)
{
  gkyl_ref_count_dec(&up->ref_count);
}



