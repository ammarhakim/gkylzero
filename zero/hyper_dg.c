#include <assert.h>
#include <math.h>
#include <time.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_hyper_dg.h>
#include <gkyl_hyper_dg_priv.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

static void
create_offsets(gkyl_hyper_dg *hdg, const int num_up_dirs, const int update_dirs[], 
  const struct gkyl_range *range, const int idxc[GKYL_MAX_DIM], long offsets[9])
{
  // Check if we're at an upper or lower edge in each direction
  bool is_edge_upper[2], is_edge_lower[2];
  for (int i=0; i<num_up_dirs; ++i) {
    is_edge_lower[i] = idxc[update_dirs[i]] == range->lower[update_dirs[i]];
    is_edge_upper[i] = idxc[update_dirs[i]] == range->upper[update_dirs[i]];
  }

  // Construct the offsets *only* in the directions being updated.
  // No need to load the neighbors that are not needed for the update.
  int lower_offset[GKYL_MAX_DIM] = {0};
  int upper_offset[GKYL_MAX_DIM] = {0};
  for (int d=0; d<num_up_dirs; ++d) {
    int dir = update_dirs[d];
    lower_offset[dir] = -1 + is_edge_lower[d];
    upper_offset[dir] = 1 - is_edge_upper[d];
  }  

  // box spanning stencil
  struct gkyl_range box3;
  gkyl_range_init(&box3, range->ndim, lower_offset, upper_offset);
  struct gkyl_range_iter iter3;
  gkyl_range_iter_init(&iter3, &box3);
  // construct list of offsets
  int count = 0;
  while (gkyl_range_iter_next(&iter3))
    offsets[count++] = gkyl_range_offset(range, iter3.idx);

}

static int idx_to_inloup_ker(int dim, const int *idx, const int *dirs, const int *num_cells) {
  int iout = 0;

  for (int d=0; d<dim; ++d) {
    if (idx[dirs[d]] == 1) {
      iout = 2*iout+(int)(pow(3,d)+0.5);
    } else if (idx[dirs[d]] == num_cells[dirs[d]]) {
      iout = 2*iout+(int)(pow(3,d)+0.5)+1;
    }
  }
  return iout;
}

void
gkyl_hyper_dg_set_update_vol(gkyl_hyper_dg *hdg, int update_vol_term)
{
  hdg->update_vol_term = update_vol_term;
}

void
gkyl_hyper_dg_advance(struct gkyl_hyper_dg *hdg, const struct gkyl_range *update_range,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  int ndim = hdg->ndim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM], idx_edge[GKYL_MAX_DIM];
  double xcl[GKYL_MAX_DIM], xcc[GKYL_MAX_DIM], xcr[GKYL_MAX_DIM], xc_edge[GKYL_MAX_DIM];
  // integer used for selecting between left-edge zero-flux BCs and right-edge zero-flux BCs
  int edge;

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
      double cfls = 0.0;
      // Assumes update_range owns lower and upper edges of the domain
      if (hdg->zero_flux_flags[dir] &&
        (idxc[dir] == update_range->lower[dir] || idxc[dir] == update_range->upper[dir]) ) {
        gkyl_copy_int_arr(ndim, iter.idx, idx_edge);
        edge = (idxc[dir] == update_range->lower[dir]) ? -1 : 1;
        // idx_edge stores interior edge index (first index away from skin cell)
        idx_edge[dir] = idx_edge[dir]-edge;

        gkyl_rect_grid_cell_center(&hdg->grid, idx_edge, xc_edge);
        long lin_edge = gkyl_range_idx(update_range, idx_edge);

        cfls = hdg->equation->boundary_surf_term(hdg->equation,
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
        long linl = gkyl_range_idx(update_range, idxl); 
        long linr = gkyl_range_idx(update_range, idxr);

        cfls = hdg->equation->surf_term(hdg->equation,
          dir, xcl, xcc, xcr, hdg->grid.dx, hdg->grid.dx, hdg->grid.dx,
          idxl, idxc, idxr,
          gkyl_array_cfetch(fIn, linl), gkyl_array_cfetch(fIn, linc), gkyl_array_cfetch(fIn, linr),
          gkyl_array_fetch(rhs, linc)
        );
      }
      double *cflrate_d = gkyl_array_fetch(cflrate, linc);
      cflrate_d[0] += cfls; // frequencies are additive      
    }
  }
}

void
gkyl_hyper_dg_gen_stencil_advance(gkyl_hyper_dg *hdg, const struct gkyl_range *update_range,
  const struct gkyl_array *fIn, struct gkyl_array *cflrate, struct gkyl_array *rhs)
{
  int ndim = hdg->ndim;

  // idxc, xc, and dx for volume update
  int idxc[GKYL_MAX_DIM] = {0};
  double xcc[GKYL_MAX_DIM] = {0.0};

  // idx, xc, and dx for generic surface update
  int idx[9][GKYL_MAX_DIM] = {0};
  double xc[9][GKYL_MAX_DIM] = {0.0};
  double dx[9][GKYL_MAX_DIM] = {0.0};
  const double* fIn_d[9];

  // bool for checking if index is in the domain
  int in_grid = 1;

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
    cflrate_d[0] += cflr;

    for (int d1=0; d1<hdg->num_up_dirs; ++d1) {
      for (int d2=0; d2<hdg->num_up_dirs; ++d2) {
        int dir1 = hdg->update_dirs[d1];
        int dir2 = hdg->update_dirs[d2];
        int update_dirs[] = {dir1, dir2};

        long offsets[9] = {0};
        int keri = 0;

        // Create offsets for 2D stencil
        if (dir1 != dir2) {
          int num_up_dirs = 2;
          create_offsets(hdg, num_up_dirs, update_dirs, update_range, idxc, offsets);

          // Index into kernel list
          keri = idx_to_inloup_ker(num_up_dirs, idxc, update_dirs, update_range->upper);
        } 
        else {
          int num_up_dirs = 1;
          create_offsets(hdg, num_up_dirs, update_dirs, update_range, idxc, offsets);

          // Index into kernel list
          keri = idx_to_inloup_ker(num_up_dirs, idxc, update_dirs, update_range->upper);
        }

        // Get pointers to all neighbor values
        for (int i=0; i<9; ++i) {
          gkyl_range_inv_idx(update_range, linc+offsets[i], idx[i]);
    
          // Check if index is in the domain
          // Assumes update_range owns lower and upper edges of the domain
          for (int d=0; d<hdg->num_up_dirs; ++d) {
            int dir = hdg->update_dirs[d];
            if (idx[i][dir] < update_range->lower[dir] || idx[i][dir] > update_range->upper[dir]) {
              in_grid = 0;
            }
          }

          // If the index is in the domain, fetch the pointer
          // otherwise, point to the central cell
          if (in_grid) {
            fIn_d[i] = gkyl_array_cfetch(fIn, linc + offsets[i]);
          }
          else {
            fIn_d[i] = gkyl_array_cfetch(fIn, linc);
          }
          // reset in_grid for next neighbor value check
          in_grid = 1;
        }

        // Domain stencil location is handled by the kernel selectors
        // gen_surf_term contains both surf and boundary surf kernels
        hdg->equation->gen_surf_term(hdg->equation,
          dir1, dir2, xcc, hdg->grid.dx, idxc,
          keri, idx, fIn_d,
          gkyl_array_fetch(rhs, linc));
      }
    }
  }
}

gkyl_hyper_dg*
gkyl_hyper_dg_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_dg_eqn *equation,
  int num_up_dirs, int update_dirs[GKYL_MAX_DIM], int zero_flux_flags[GKYL_MAX_DIM],
  int update_vol_term, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_hyper_dg_cu_dev_new(grid, basis, equation, num_up_dirs, update_dirs, zero_flux_flags, update_vol_term);
  } 
#endif
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

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  
  up->on_dev = up; // on host, on_dev points to itself

  return up;
}

void gkyl_hyper_dg_release(struct gkyl_hyper_dg* hdg)
{
  gkyl_dg_eqn_release(hdg->equation);
  if (GKYL_IS_CU_ALLOC(hdg->flags))
    gkyl_cu_free(hdg->on_dev);
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

void gkyl_hyper_dg_advance_cu(gkyl_hyper_dg* hdg, const struct gkyl_range *update_range,
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
