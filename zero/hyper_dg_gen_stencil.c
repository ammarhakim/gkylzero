#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_hyper_dg_gen_stencil.h>
#include <gkyl_hyper_dg_gen_stencil_priv.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

static void
create_offsets(const struct gkyl_range *range, long offsets[])
{
  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, (int[]) { -1, -1, -1 }, (int[]) { 1, 1, 1 });
  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);
  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);
}

void
gkyl_hyper_dg_gen_stencil_advance(gkyl_hyper_dg_gen_stencil *hdg, const struct gkyl_range *update_range,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  int ndim = hdg->ndim;
  long sz[] = { 3, 9, 27 };
  long offsets[sz[ndim-1]];
  create_offsets(update_range, offsets);

  // idx, xc, and dx for volume update
  int idxc[GKYL_MAX_DIM];
  double xcc[GKYL_MAX_DIM];

  // idx, xc, and dx for generic surface update
  int idx[sz[ndim-1]][GKYL_MAX_DIM];
  double xc[sz[ndim-1]][GKYL_MAX_DIM];
  double dx[sz[ndim-1]][GKYL_MAX_DIM];
  const double* fIn_d[sz[ndim-1]];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, update_range);
  while (gkyl_range_iter_next(&iter)) {
    long linc = gkyl_range_idx(update_range, iter.idx);

    // Call volume kernel and get CFL rate
    gkyl_copy_int_arr(ndim, iter.idx, idxc);
    gkyl_rect_grid_cell_center(&hdg->grid, idxc, xcc);
    double cflr = hdg->equation->vol_term(
      hdg->equation, xcc, hdg->grid.dx, idxc,
      gkyl_array_cfetch(fIn, linc), gkyl_array_fetch(rhs, linc)
    );
    double *cflrate_d = gkyl_array_fetch(cflrate, linc);
    cflrate_d[0] += cflr; // frequencies are additive

    // Get pointers to all neighbor values (i.e., 9 cells in 2D, 27 cells in 3D)
    for (int i=0; i<sz[ndim-1]; ++i) {
      gkyl_sub_range_inv_idx(update_range, linc+offsets[i], idx[i]);
      gkyl_rect_grid_cell_center(&hdg->grid, idx[i], xc[i]);
      for (int j=0; j<ndim; ++j)
        dx[i][j] = hdg->grid.dx[j];
      fIn_d[i] = gkyl_array_cfetch(fIn, linc + offsets[i]);
    }

    // Loop over surfaces and update using any/all neighbors needed
    for (int d=0; d<hdg->num_up_dirs; ++d) {
      int dir = hdg->update_dirs[d];
      hdg->equation->gen_term(hdg->equation,
        dir, xc, dx, idx, fIn_d,
        gkyl_array_fetch(rhs, linc)
      );
    }
  }
}

gkyl_hyper_dg_gen_stencil*
gkyl_hyper_dg_gen_stencil_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_dg_eqn *equation,
  int num_up_dirs, int update_dirs[GKYL_MAX_DIM], bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_hyper_dg_gen_stencil_cu_dev_new(grid, basis, equation, num_up_dirs, update_dirs);
  } 
#endif
  gkyl_hyper_dg_gen_stencil *up = gkyl_malloc(sizeof(gkyl_hyper_dg_gen_stencil));

  up->grid = *grid;
  up->ndim = basis->ndim;
  up->num_basis = basis->num_basis;
  up->num_up_dirs = num_up_dirs;

  for (int i=0; i<num_up_dirs; ++i)
    up->update_dirs[i] = update_dirs[i];
  
  up->equation = gkyl_dg_eqn_acquire(equation);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  
  up->on_dev = up; // on host, on_dev points to itself

  return up;
}

void gkyl_hyper_dg_gen_stencil_release(gkyl_hyper_dg_gen_stencil* hdg)
{
  gkyl_dg_eqn_release(hdg->equation);
  if (GKYL_IS_CU_ALLOC(hdg->flags))
    gkyl_cu_free(hdg->on_dev);
  gkyl_free(hdg);
}

#ifndef GKYL_HAVE_CUDA

// default functions
gkyl_hyper_dg_gen_stencil*
gkyl_hyper_dg_gen_stencil_cu_dev_new(const struct gkyl_rect_grid *grid_cu,
  const struct gkyl_basis *basis, const struct gkyl_dg_eqn *equation_cu,
  int num_up_dirs, int update_dirs[])
{
  assert(false);
  return 0;
}

void gkyl_hyper_dg_gen_stencil_advance_cu(gkyl_hyper_dg_gen_stencil* hdg, const struct gkyl_range *update_range,
  const struct gkyl_array* GKYL_RESTRICT fIn, struct gkyl_array* GKYL_RESTRICT cflrate,
  struct gkyl_array* GKYL_RESTRICT rhs)
{
  assert(false);
}

#endif
