#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_mom_calc_bcorr.h>
#include <gkyl_mom_calc_bcorr_priv.h>
#include <gkyl_util.h>

void
gkyl_mom_calc_bcorr_advance(gkyl_mom_calc_bcorr *bcorr,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *GKYL_RESTRICT fIn, struct gkyl_array *GKYL_RESTRICT out)
{
  double xc[GKYL_MAX_DIM];
  struct gkyl_range edge_rng;
  struct gkyl_range_iter conf_iter, edge_iter;
  
  int pidx[GKYL_MAX_DIM], eiter_idx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  int elower_idx[GKYL_MAX_DIM], eupper_idx[GKYL_MAX_DIM] = { 0 };
  int conf_idx[GKYL_MAX_CDIM];
  int dim_rng;
  enum gkyl_vel_edge edge;

  for (int dim=0; dim<phase_rng->ndim; ++dim) {
    elower_idx[dim] = phase_rng->lower[dim];
    eupper_idx[dim] = phase_rng->upper[dim];
  }

  // iterate over either velocity space or configuration space
  dim_rng = phase_rng->ndim - conf_rng->ndim;
  if (bcorr->space == 0) {
    dim_rng = conf_rng->ndim;
  }

  for (int d=0; d<dim_rng; ++d) {
    edge = d + GKYL_MAX_CDIM;
    elower_idx[bcorr->space + d] = phase_rng->upper[bcorr->space + d];
    eupper_idx[bcorr->space + d] = phase_rng->upper[bcorr->space + d];
    gkyl_sub_range_init(&edge_rng, phase_rng, elower_idx, eupper_idx);
    gkyl_range_iter_no_split_init(&edge_iter, &edge_rng);

    while (gkyl_range_iter_next(&edge_iter)) {
      gkyl_rect_grid_cell_center(&bcorr->grid, edge_iter.idx, xc);

      for (int i=0; i<conf_rng->ndim; ++i) {
        conf_idx[i] = edge_iter.idx[i];
      }
      long fidx = gkyl_range_idx(&edge_rng, edge_iter.idx);
      long midx = gkyl_range_idx(conf_rng, conf_idx);
      
      gkyl_mom_type_calc(bcorr->momt, xc, bcorr->grid.dx, pidx,
        gkyl_array_cfetch(fIn, fidx), gkyl_array_fetch(out, midx), &edge
      );
    }

    edge = d;
    elower_idx[d] = phase_rng->lower[bcorr->space + d];
    eupper_idx[d] = phase_rng->lower[bcorr->space + d];
    gkyl_sub_range_init(&edge_rng, phase_rng, elower_idx, eupper_idx);
    gkyl_range_iter_no_split_init(&edge_iter, &edge_rng);

    while (gkyl_range_iter_next(&edge_iter)) {
      gkyl_rect_grid_cell_center(&bcorr->grid, edge_iter.idx, xc);

      for (int i=0; i<conf_rng->ndim; ++i) {
        conf_idx[i] = edge_iter.idx[i];
      }
      
      long fidx = gkyl_range_idx(&edge_rng, edge_iter.idx);
      long midx = gkyl_range_idx(conf_rng, conf_idx);
       
      gkyl_mom_type_calc(bcorr->momt, xc, bcorr->grid.dx, pidx,
        gkyl_array_cfetch(fIn, fidx), gkyl_array_fetch(out, midx), &edge
      );
    }
    
    elower_idx[bcorr->space + d] = phase_rng->lower[bcorr->space + d];
    eupper_idx[bcorr->space + d] = phase_rng->upper[bcorr->space + d];
  }
}

gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt, const char *space)
{
  gkyl_mom_calc_bcorr *up = gkyl_malloc(sizeof(gkyl_mom_calc_bcorr));
  up->grid = *grid;
  up->momt = gkyl_mom_type_acquire(momt);
  if (strcmp(space, "conf") == 0) {
    up->space = 0;
  } else {
    up->space = momt->cdim;
  }
  
  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up;
  
  return up;
}

void
gkyl_mom_calc_bcorr_release(gkyl_mom_calc_bcorr* up)
{
  gkyl_mom_type_release(up->momt);
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

void
gkyl_mom_calc_bcorr_advance_cu(gkyl_mom_calc_bcorr *bcorr,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *GKYL_RESTRICT fIn, struct gkyl_array *GKYL_RESTRICT out)
{
  assert(false);
}

gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_mom_type *momt, const char *space)
{
  assert(false);
}

#endif
