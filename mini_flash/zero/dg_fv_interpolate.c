#include <gkyl_alloc.h>
#include <gkyl_basis.h>
#include <gkyl_dg_fv_interpolate.h>
#include <gkyl_dg_fv_interpolate_priv.h>

#include <assert.h>

struct gkyl_dg_fv_interpolate*
gkyl_dg_fv_interpolate_new(struct gkyl_basis *basis, 
  struct gkyl_rect_grid *grid_dg, struct gkyl_rect_grid *grid_fv)
{
  // Allocate memory for updater.
  struct gkyl_dg_fv_interpolate *up = gkyl_malloc(sizeof(*up));

  up->num_basis = basis->num_basis;

  // Verify that grid dimensions match.
  assert(grid_dg->ndim == grid_fv->ndim);
  up->ndim = grid_dg->ndim;

  // Verify that grid resolutions are correct.
  for (int i=0; i<up->ndim; ++i) {
    assert(2*grid_dg->cells[i] == grid_fv->cells[i]);
  }
  up->grid_dg = grid_dg;
  up->grid_fv = grid_fv;

  // Set function pointers.
  up->kernels.modal_to_fv = choose_dg_fv_interp_modal_to_fv_kernel(up->ndim);
  up->kernels.fv_to_modal = choose_dg_fv_interp_fv_to_modal_kernel(up->ndim);
  up->kernels.fv_to_modal_recov = choose_dg_fv_interp_fv_to_modal_recov_kernel(up->ndim);

  return up;
}

void
gkyl_dg_fv_interpolate_fv_to_modal_advance(gkyl_dg_fv_interpolate *up, 
  struct gkyl_range *range_dg, struct gkyl_range *range_fv,
  const struct gkyl_array *fv_in, int comp_in,
  struct gkyl_array *dg_out, int comp_out)
{
  // Array of pointers to finite volume data. Up to 8 because we only ever
  // need at most 8 FV cells per DG cell (3D).
  const double *fv_data[8];

  // Sub range to iterate over corresponding cells in FV grid.
  struct gkyl_range sub_range;
  int sub_lower[GKYL_MAX_DIM], sub_upper[GKYL_MAX_DIM];

  // Iterate over coarse DG grid.
  struct gkyl_range_iter iter, sub_iter;
  gkyl_range_iter_init(&iter, range_dg);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(range_dg, iter.idx);
    double *dg_out_d = &(((double *)gkyl_array_fetch(dg_out, linc))[comp_out*up->num_basis]);

    // Compute indices of cells in refined grid at lower and 
    // upper corners of the coarse cell.
    for (int i=0; i<up->ndim; ++i) {
      sub_lower[i] = iter.idx[i]*2 - 1;
      sub_upper[i] = sub_lower[i] + 1;
    }

    // Fetch pointers to corresponding refined grid cells.
    gkyl_sub_range_init(&sub_range, range_fv, sub_lower, sub_upper);
    gkyl_range_iter_init(&sub_iter, &sub_range);
    int count = 0;
    while (gkyl_range_iter_next(&sub_iter)) {
      long sub_linc = gkyl_range_idx(range_fv, sub_iter.idx);

      fv_data[count] = &(((const double*)gkyl_array_cfetch(fv_in, sub_linc))[comp_in]);
      count++;
    }

    // Convert finite volume data to p=1 modal data.
    up->kernels.fv_to_modal(fv_data, dg_out_d);
  }
}


void
gkyl_dg_fv_interpolate_fv_to_modal_recov_advance(gkyl_dg_fv_interpolate *up, 
  struct gkyl_range *range_dg, struct gkyl_range *range_fv,
  const struct gkyl_array *fv_in, int comp_in,
  struct gkyl_array *dg_out, int comp_out)
{
  // Array of pointers to finite volume data. Up to 8 because we only ever
  // need at most 8 FV cells per DG cell (3D).
  const double *fv_data[8];

  // Sub range to iterate over corresponding cells in FV grid.
  struct gkyl_range sub_range;
  int sub_lower[GKYL_MAX_DIM], sub_upper[GKYL_MAX_DIM];

  // Iterate over coarse DG grid.
  struct gkyl_range_iter iter, sub_iter;
  gkyl_range_iter_init(&iter, range_dg);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(range_dg, iter.idx);

    double *dg_out_d = &(((double *)gkyl_array_fetch(dg_out, linc))[comp_out*up->num_basis]);

    // Compute indices of cells in refined grid at lower and 
    // upper corners of the coarse cell.
    for (int i=0; i<up->ndim; ++i) {
      sub_lower[i] = iter.idx[i]*2 - 1;
      sub_upper[i] = sub_lower[i] + 1;
    }

    // Fetch pointers to corresponding refined grid cells.
    gkyl_sub_range_init(&sub_range, range_fv, sub_lower, sub_upper);
    gkyl_range_iter_init(&sub_iter, &sub_range);
    int count = 0;
    while (gkyl_range_iter_next(&sub_iter)) {
      long sub_linc = gkyl_range_idx(range_fv, sub_iter.idx);
      
      fv_data[count] = &(((const double*)gkyl_array_cfetch(fv_in, sub_linc))[comp_in]);
      count++;
    }

    // Convert finite volume data to p=1 modal data.
    up->kernels.fv_to_modal_recov(fv_data, dg_out_d);
  }
}

void
gkyl_dg_fv_interpolate_modal_to_fv_advance(gkyl_dg_fv_interpolate *up, 
  struct gkyl_range *range_dg, struct gkyl_range *range_fv,
  const struct gkyl_array *dg_in, int comp_in,
  struct gkyl_array *fv_out, int comp_out) {
  // Array of pointers to finite volume data. Up to 8 because we only ever
  // need at most 8 cells (3D).
  double *fv_data[8];

  // Sub range to iterate over corresponding cells in FV grid.
  struct gkyl_range sub_range;
  int sub_lower[GKYL_MAX_DIM], sub_upper[GKYL_MAX_DIM];

  // Iterate over coarse DG grid.
  struct gkyl_range_iter iter, sub_iter;
  gkyl_range_iter_init(&iter, range_dg);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(range_dg, iter.idx);
    const double *dg_in_d = &(((double *)gkyl_array_cfetch(dg_in, linc))[comp_in*up->num_basis]);

    // Compute indices of cells in refined grid at lower and 
    // upper corners of the coarse cell.
    for (int i=0; i<up->ndim; ++i) {
      sub_lower[i] = iter.idx[i]*2 - 1;
      sub_upper[i] = sub_lower[i] + 1;
    }

    // Fetch pointers to corresponding refined grid cells.
    gkyl_sub_range_init(&sub_range, range_fv, sub_lower, sub_upper);
    gkyl_range_iter_init(&sub_iter, &sub_range);
    int count = 0;
    while (gkyl_range_iter_next(&sub_iter)) {
      long sub_linc = gkyl_range_idx(range_fv, sub_iter.idx);

      fv_data[count] = &(((double*)gkyl_array_fetch(fv_out, sub_linc))[comp_out]);
      count++;
    }

    // Convert finite volume data to p=1 modal data.
    up->kernels.modal_to_fv(dg_in_d, fv_data);
  }
}

void
gkyl_dg_fv_interpolate_release(gkyl_dg_fv_interpolate *up)
{
  gkyl_free(up);
}