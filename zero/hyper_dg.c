#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_util.h>

struct gkyl_hyper_dg {
    struct gkyl_rect_grid grid; // grid object
    int ndim; // number of dimensions
    int numBasis; // number of basis functions
    int num_up_dirs; // number of update directions
    int update_dirs[GKYL_MAX_DIM]; // directions to update
    int zero_flux_flags[GKYL_MAX_DIM]; // directions with zero flux
    int update_vol_term; // should we update volume term?
    const struct gkyl_dg_eqn *equation; // equation object
};

void
gkyl_hyper_dg_advance(const gkyl_hyper_dg *hdg, const struct gkyl_range *update_range,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs, double *maxs)
{
  int ndim = hdg->ndim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  double xcl[GKYL_MAX_DIM], xcc[GKYL_MAX_DIM], xcr[GKYL_MAX_DIM];
  // integer used for selecting between left-edge zero-flux BCs and right-edge zero-flux BCs
  int edge;

  double maxs_old[GKYL_MAX_DIM];
  for (int i=0; i<hdg->ndim; ++i)
    maxs_old[i] = maxs[i];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(ndim, iter.idx, idxc);
    gkyl_rect_grid_cell_center(&hdg->grid, idxc, xcc);

    long linc = gkyl_range_idx(update_range, idxc);
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
      gkyl_copy_int_arr(ndim, iter.idx, idxl);
      gkyl_copy_int_arr(ndim, iter.idx, idxr);
      // TODO: fix for arbitrary subrange
      if (hdg->zero_flux_flags[d] && (idxc[dir] == update_range->lower[dir] || idxc[dir] == update_range->upper[dir])) {
        edge = (idxc[dir] == update_range->lower[dir]) ? -1 : 1;
        // use idxl to store interior edge index (first index away from skin cell)
        idxl[dir] = idxl[dir]-edge;

        gkyl_rect_grid_cell_center(&hdg->grid, idxl, xcl);
        long linl = gkyl_range_idx(update_range, idxl);

        double mdir = hdg->equation->boundary_surf_term(hdg->equation,
          dir, xcl, xcc, hdg->grid.dx, hdg->grid.dx,
          maxs_old[dir], idxl, idxc, edge,
          gkyl_array_cfetch(fIn, linl), gkyl_array_cfetch(fIn, linc),
          gkyl_array_fetch(rhs, linc)
        );
        maxs[dir] = fmax(maxs[dir], mdir);         
      }
      else {
        idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;
        gkyl_rect_grid_cell_center(&hdg->grid, idxl, xcl);
        gkyl_rect_grid_cell_center(&hdg->grid, idxr, xcr);
        long linl = gkyl_range_idx(update_range, idxl); 
        long linr = gkyl_range_idx(update_range, idxr);

        double mdir = hdg->equation->surf_term(hdg->equation,
          dir, xcl, xcc, xcr, hdg->grid.dx, hdg->grid.dx, hdg->grid.dx,
          maxs_old[dir], idxl, idxc, idxr,
          gkyl_array_cfetch(fIn, linl), gkyl_array_cfetch(fIn, linc), gkyl_array_cfetch(fIn, linr),
          gkyl_array_fetch(rhs, linc)
        );
        maxs[dir] = fmax(maxs[dir], mdir);
      }
    }
  }
}

gkyl_hyper_dg*
gkyl_hyper_dg_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_dg_eqn *equation,
  int num_up_dirs, int update_dirs[], int zero_flux_flags[], int update_vol_term)
{
  gkyl_hyper_dg *up = gkyl_malloc(sizeof(gkyl_hyper_dg));

  up->grid = *grid;
  up->ndim = basis->ndim;
  up->numBasis = basis->numBasis;
  up->num_up_dirs = num_up_dirs;

  for (int i=0; i<num_up_dirs; ++i) {
    up->update_dirs[i] = update_dirs[i];
    up->zero_flux_flags[i] = zero_flux_flags[i];
  }
  up->update_vol_term = update_vol_term;
  up->equation = gkyl_dg_eqn_aquire(equation);

  return up;
}

void
gkyl_hyper_dg_set_update_vol(gkyl_hyper_dg *hdg, int update_vol_term)
{
  hdg->update_vol_term = update_vol_term;
}

void
gkyl_hyper_dg_release(gkyl_hyper_dg* up)
{
  gkyl_dg_eqn_release(up->equation);
  free(up);
}
