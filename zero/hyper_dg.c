#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_hyper_dg_priv.h>
#include <gkyl_util.h>

void
gkyl_hyper_dg_set_update_vol(gkyl_hyper_dg *hdg, int update_vol_term)
{
  hdg->update_vol_term = update_vol_term;
}

void
gkyl_hyper_dg_advance(gkyl_hyper_dg *hdg, struct gkyl_range update_range,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  int ndim = hdg->ndim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM], idx_edge[GKYL_MAX_DIM];
  double xcl[GKYL_MAX_DIM], xcc[GKYL_MAX_DIM], xcr[GKYL_MAX_DIM], xc_edge[GKYL_MAX_DIM];
  // integer used for selecting between left-edge zero-flux BCs and right-edge zero-flux BCs
  int edge;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &update_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(ndim, iter.idx, idxc);
    gkyl_rect_grid_cell_center(&hdg->grid, idxc, xcc);

    long linc = gkyl_range_idx(&update_range, idxc);
    if (hdg->update_vol_term) {
      double cflr = hdg->equation->vol_term(
        hdg->equation, xcc, hdg->grid.dx, idxc,
        gkyl_array_cfetch(fIn, linc), gkyl_array_fetch(rhs, linc)
      );
      //printf("linc is %ld\n", linc);
      double *cflrate_d = gkyl_array_fetch(cflrate, linc);
      cflrate_d[0] += cflr; // frequencies are additive
    }
    
    for (int d=0; d<hdg->num_up_dirs; ++d) {
      int dir = hdg->update_dirs[d];
      // TODO: fix for arbitrary subrange
      if (hdg->zero_flux_flags[dir] &&
        (idxc[dir] == update_range.lower[dir] || idxc[dir] == update_range.upper[dir]) ) {
        gkyl_copy_int_arr(ndim, iter.idx, idx_edge);
        edge = (idxc[dir] == update_range.lower[dir]) ? -1 : 1;
        // idx_edge stores interior edge index (first index away from skin cell)
        idx_edge[dir] = idx_edge[dir]-edge;

        gkyl_rect_grid_cell_center(&hdg->grid, idx_edge, xc_edge);
        long lin_edge = gkyl_range_idx(&update_range, idx_edge);

        hdg->equation->boundary_surf_term(hdg->equation,
          dir, xc_edge, xcc, hdg->grid.dx, hdg->grid.dx,
          idx_edge, idxc, edge,
          gkyl_array_cfetch(fIn, lin_edge), gkyl_array_cfetch(fIn, linc),
          gkyl_array_fetch(rhs, linc)
        );
      }
      else {
        gkyl_copy_int_arr(ndim, iter.idx, idxl);
        gkyl_copy_int_arr(ndim, iter.idx, idxr);
        idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;
        
        gkyl_rect_grid_cell_center(&hdg->grid, idxl, xcl);
        gkyl_rect_grid_cell_center(&hdg->grid, idxr, xcr);
        long linl = gkyl_range_idx(&update_range, idxl); 
        long linr = gkyl_range_idx(&update_range, idxr);

        hdg->equation->surf_term(hdg->equation,
          dir, xcl, xcc, xcr, hdg->grid.dx, hdg->grid.dx, hdg->grid.dx,
          idxl, idxc, idxr,
          gkyl_array_cfetch(fIn, linl), gkyl_array_cfetch(fIn, linc), gkyl_array_cfetch(fIn, linr),
          gkyl_array_fetch(rhs, linc)
        );
      }
    }
  }
}

void
gkyl_hyper_dg_advance_no_iter(gkyl_hyper_dg *hdg, struct gkyl_range update_range,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  int ndim = hdg->ndim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  double xcl[GKYL_MAX_DIM], xcc[GKYL_MAX_DIM], xcr[GKYL_MAX_DIM];
  // integer used for selecting between left-edge zero-flux BCs and right-edge zero-flux BCs
  int edge;

  for (long linc1 = 0; linc1 < update_range.volume; linc1++) {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={0,0,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&update_range, linc1, idxc);
    gkyl_rect_grid_cell_center(&hdg->grid, idxc, xcc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&update_range, idxc);

    if (hdg->update_vol_term) {
      double cflr = hdg->equation->vol_term(
        hdg->equation, xcc, hdg->grid.dx, idxc,
        gkyl_array_cfetch(fIn, linc), gkyl_array_fetch(rhs, linc)
      );
      double *cflrate_d = gkyl_array_fetch(cflrate, linc);
      cflrate_d[0] += cflr; // frequencies are additive
    }
    
    for (int d=0; d<hdg->num_up_dirs; ++d) {
      int dir = hdg->update_dirs[d];
      gkyl_copy_int_arr(ndim, idxc, idxl);
      gkyl_copy_int_arr(ndim, idxc, idxr);
      // TODO: fix for arbitrary subrange
      if (hdg->zero_flux_flags[dir] &&
        (idxc[dir] == update_range.lower[dir] || idxc[dir] == update_range.upper[dir])) {
        edge = (idxc[dir] == update_range.lower[dir]) ? -1 : 1;
        // use idxl to store interior edge index (first index away from skin cell)
        idxl[dir] = idxl[dir]-edge;

        gkyl_rect_grid_cell_center(&hdg->grid, idxl, xcl);
        long linl = gkyl_range_idx(&update_range, idxl);

        hdg->equation->boundary_surf_term(hdg->equation,
          dir, xcl, xcc, hdg->grid.dx, hdg->grid.dx,
          idxl, idxc, edge,
          gkyl_array_cfetch(fIn, linl), gkyl_array_cfetch(fIn, linc),
          gkyl_array_fetch(rhs, linc)
        );     
      }
      else {
        idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;
        gkyl_rect_grid_cell_center(&hdg->grid, idxl, xcl);
        gkyl_rect_grid_cell_center(&hdg->grid, idxr, xcr);
        long linl = gkyl_range_idx(&update_range, idxl); 
        long linr = gkyl_range_idx(&update_range, idxr);

        hdg->equation->surf_term(hdg->equation,
          dir, xcl, xcc, xcr, hdg->grid.dx, hdg->grid.dx, hdg->grid.dx,
          idxl, idxc, idxr,
          gkyl_array_cfetch(fIn, linl), gkyl_array_cfetch(fIn, linc), gkyl_array_cfetch(fIn, linr),
          gkyl_array_fetch(rhs, linc)
        );
      }
    }
  }
}

gkyl_hyper_dg*
gkyl_hyper_dg_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_dg_eqn *equation,
  int num_up_dirs, int update_dirs[GKYL_MAX_DIM], int zero_flux_flags[GKYL_MAX_DIM],
  int update_vol_term)
{
  gkyl_hyper_dg *up = gkyl_malloc(sizeof(gkyl_hyper_dg));

  up->grid = *grid;
  up->ndim = basis->ndim;
  up->num_basis = basis->num_basis;
  up->num_up_dirs = num_up_dirs;

  for (int i=0; i<num_up_dirs; ++i)
    up->update_dirs[i] = update_dirs[i];
  for (int i=0; i<GKYL_MAX_DIM; ++i)
    up->zero_flux_flags[i] = zero_flux_flags[i];
    
  up->update_vol_term = update_vol_term;
  up->equation = gkyl_dg_eqn_acquire(equation);

  return up;
}

void gkyl_hyper_dg_release(gkyl_hyper_dg* hdg)
{
  gkyl_dg_eqn_release(hdg->equation);
  gkyl_free(hdg);
}

#ifndef GKYL_HAVE_CUDA

// default functions
gkyl_hyper_dg*
gkyl_hyper_dg_cu_dev_new(const struct gkyl_rect_grid *grid_cu,
  const struct gkyl_basis *basis, const struct gkyl_dg_eqn *equation_cu,
  int num_up_dirs, int update_dirs[], int zero_flux_flags[], int update_vol_term)
{
  assert(false);
  return 0;
}

void gkyl_hyper_dg_advance_cu(gkyl_hyper_dg* hdg, struct gkyl_range update_range,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs)
{
  assert(false);
}

void
gkyl_hyper_dg_set_update_vol_cu(gkyl_hyper_dg *hdg, int update_vol_term)
{
  assert(false);
}

#endif
