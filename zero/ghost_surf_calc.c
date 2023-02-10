#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_ghost_surf_calc.h>
#include <gkyl_ghost_surf_calc_priv.h>
#include <gkyl_util.h>

void
gkyl_ghost_surf_calc_advance(gkyl_ghost_surf_calc *gcalc,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *fIn, struct gkyl_array *rhs)
{
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM], idx_edge[GKYL_MAX_DIM];
  double xcl[GKYL_MAX_DIM], xcc[GKYL_MAX_DIM], xcr[GKYL_MAX_DIM];
  struct gkyl_range edge_rng;
  struct gkyl_range_iter conf_iter, edge_iter;
  
  int pidx[GKYL_MAX_DIM], eiter_idx[GKYL_MAX_DIM], rem_dir[GKYL_MAX_DIM] = { 0 };
  int clower_idx[GKYL_MAX_DIM], cupper_idx[GKYL_MAX_DIM] = { 0 };
  int conf_idx[GKYL_MAX_CDIM];

  for (int d=0; d<phase_rng->ndim; ++d) {
    clower_idx[d] = phase_rng->lower[d];
    cupper_idx[d] = phase_rng->upper[d];
  }

  for (int dir=0; dir<conf_rng->ndim; ++dir) {
    clower_idx[dir] = phase_rng->upper[dir];
    cupper_idx[dir] = phase_rng->upper[dir];
    gkyl_sub_range_init(&edge_rng, phase_rng, clower_idx, cupper_idx);
    gkyl_range_iter_no_split_init(&edge_iter, &edge_rng);

    while (gkyl_range_iter_next(&edge_iter)) {
      gkyl_copy_int_arr(edge_rng.ndim, edge_iter.idx, idxc);
      gkyl_copy_int_arr(edge_rng.ndim, edge_iter.idx, idxl);
      idxl[dir] = idxl[dir]-1;

      gkyl_rect_grid_cell_center(&gcalc->grid, idxc, xcc);
      gkyl_rect_grid_cell_center(&gcalc->grid, idxl, xcl);

      long linc = gkyl_range_idx(&edge_rng, idxc);
      long linl = gkyl_range_idx(&edge_rng, idxl);

      gcalc->equation->boundary_surf_term(gcalc->equation, dir, xcl, xcc,
        gcalc->grid.dx, gcalc->grid.dx, idxc, idxc, 1, 
	gkyl_array_cfetch(fIn, linl), gkyl_array_cfetch(fIn, linc), gkyl_array_fetch(rhs, linc)
      );
    }

    clower_idx[dir] = phase_rng->lower[dir];
    cupper_idx[dir] = phase_rng->lower[dir];
    gkyl_sub_range_init(&edge_rng, phase_rng, clower_idx, cupper_idx);
    gkyl_range_iter_no_split_init(&edge_iter, &edge_rng);

    while (gkyl_range_iter_next(&edge_iter)) {
      gkyl_copy_int_arr(edge_rng.ndim, edge_iter.idx, idxc);
      gkyl_copy_int_arr(edge_rng.ndim, edge_iter.idx, idxr);
      idxr[dir] = idxr[dir]+1;

      gkyl_rect_grid_cell_center(&gcalc->grid, idxc, xcc);
      gkyl_rect_grid_cell_center(&gcalc->grid, idxr, xcr);

      long linc = gkyl_range_idx(&edge_rng, idxc); 
      long linr = gkyl_range_idx(&edge_rng, idxr);

      gcalc->equation->boundary_surf_term(gcalc->equation, dir, xcc, xcr,
        gcalc->grid.dx, gcalc->grid.dx, idxc, idxr, -1,
	gkyl_array_cfetch(fIn, linr), gkyl_array_cfetch(fIn, linc), gkyl_array_fetch(rhs, linc)
      );
    }
    clower_idx[dir] = phase_rng->lower[dir];
    cupper_idx[dir] = phase_rng->upper[dir];
  }
}

gkyl_ghost_surf_calc*
gkyl_ghost_surf_calc_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_dg_eqn *equation)
{
  gkyl_ghost_surf_calc *up = gkyl_malloc(sizeof(gkyl_ghost_surf_calc));
  up->grid = *grid;
  up->equation = gkyl_dg_eqn_acquire(equation);
  
  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  
  up->on_dev = up;
  
  return up;
}

void
gkyl_ghost_surf_calc_release(gkyl_ghost_surf_calc* up)
{
  gkyl_dg_eqn_release(up->equation);
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

void
gkyl_ghost_surf_calc_advance_cu(gkyl_ghost_surf_calc *gcalc,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *fIn, struct gkyl_array *rhs)
{
  assert(false);
}

gkyl_ghost_surf_calc*
gkyl_ghost_surf_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_dg_eqn *equation)
{
  assert(false);
}

#endif
