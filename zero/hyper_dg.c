#include <math.h>

#include <gkyl_alloc.h>
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

    double maxs[GKYL_MAX_DIM];
    double maxs_old[GKYL_MAX_DIM], maxs_local[GKYL_MAX_DIM];
};

static inline void
copy_int_arr(int n, const int *restrict inp, int *restrict out)
{
  for (int i=0; i<n; ++i)  out[i] = inp[i];
}

void
gkyl_hyper_dg_advance(gkyl_hyper_dg *hdg, const struct gkyl_range *update_range,
  const struct gkyl_array *fIn, const struct gkyl_array *rhs)
{
  int ndim = hdg->ndim;
  int firstDir = 1;
  int idxm[GKYL_MAX_DIM], idxp[GKYL_MAX_DIM];
  double xcm[GKYL_MAX_DIM], xcp[GKYL_MAX_DIM];
  
  // we need to use previous step values for relaxation parameter
  for (int i=0; i<hdg->ndim; ++i)
    hdg->maxs_old[i] = hdg->maxs[i];
  
  for (int d=0; d<hdg->num_up_dirs; ++d) {
    int dir = hdg->update_dirs[d];

    int dirLoIdx = update_range->lower[dir];
    int dirUpIdx = update_range->upper[dir]+1; // one more edge than cells

    if (hdg->zero_flux_flags[dir]) {
      dirLoIdx = dirLoIdx+1;
      dirUpIdx = dirUpIdx-1;
    }

    // create range perpendicular to 'dir'
    struct gkyl_range perp_range;
    gkyl_range_shorten(&perp_range, update_range, dir, 1);
    struct gkyl_range_iter iter;
    gkyl_range_iter_init(&iter, &perp_range);

    while (gkyl_range_iter_next(&iter)) {
      copy_int_arr(ndim, iter.idx, idxm);
      copy_int_arr(ndim, iter.idx, idxp);

      for (int i=dirLoIdx; i<=dirUpIdx; ++i) { // note dirUpIdx is inclusive
        idxm[dir] = i-1; idxp[dir] = i;

        gkyl_rect_grid_cell_center(&hdg->grid, idxm, xcm);
        gkyl_rect_grid_cell_center(&hdg->grid, idxp, xcp);

        long linIdxm = gkyl_range_idx(update_range, idxm);
        long linIdxp = gkyl_range_idx(update_range, idxp);

        if (firstDir && hdg->update_vol_term && i<=dirUpIdx-1) {
          double cflRate = hdg->equation->vol_term(
            hdg->equation, xcp, hdg->grid.dx, idxp,
            gkyl_array_fetch(fIn, linIdxp), gkyl_array_fetch(rhs, linIdxp)
          );
        }

        double maxs = hdg->equation->surf_term(hdg->equation,
          dir, xcm, xcp, hdg->grid.dx, hdg->grid.dx,
          hdg->maxs_old[dir], idxm, idxp,
          gkyl_array_fetch(fIn, linIdxm), gkyl_array_fetch(fIn, linIdxp),
          gkyl_array_fetch(rhs, linIdxm), gkyl_array_fetch(rhs, linIdxp)
        );
        hdg->maxs[dir] = fmax(hdg->maxs[dir], maxs);
      }
    }
    firstDir = 0;
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

  for (int i=0; i<up->ndim; ++i)
    up->maxs[i] = up->maxs_old[i] = up->maxs_local[i] = 0.0;

  return up;
}

void
gkyl_hyper_dg_release(gkyl_hyper_dg* up)
{
  gkyl_dg_eqn_release(up->equation);
  free(up);
}
