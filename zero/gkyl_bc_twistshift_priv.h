#pragma once

// Private header for bc_twistshift updater, not for direct use in user code.

#include <gkyl_bc_twistshift.h>
#include <gkyl_bc_twistshift_gyrokinetic_kernels.h>
#include <assert.h>
#include <gkyl_mat.h>
#include <string.h> // memcpy

// Function pointer type for twistshift kernels.
typedef void (*twistshift_xlimdg_t)(double sFac, const double *xLimLo,
  const double *xLimUp, double yLimLo, double yLimUp,
  double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);

typedef void (*twistshift_ylimdg_t)(double sFac, double xLimLo,
  double xLimUp, const double *yLimLo, const double *yLimUp,
  double dyDo, double yOff, const double *ySh, struct gkyl_mat *tsmat);

typedef void (*twistshift_fullcell_t)(double dyDo, double yOff,
  const double *ySh, struct gkyl_mat *tsmat);

typedef struct { twistshift_xlimdg_t kernels[3]; }   twistshift_xlimdg_kern_list;  // For use in kernel tables.
typedef struct { twistshift_ylimdg_t kernels[3]; }   twistshift_ylimdg_kern_list;  // For use in kernel tables.
typedef struct { twistshift_fullcell_t kernels[3]; } twistshift_fullcell_kern_list;  // For use in kernel tables.

// Serendipity  kernels.
static const twistshift_xlimdg_kern_list ser_twistshift_xlimdg_list_0v[] = {
  {NULL, twistshift_xlimdg_2x_ser_p1_yshift_p1, NULL,},
  {NULL, twistshift_xlimdg_3x_ser_p1_yshift_p1, NULL,},
};
static const twistshift_ylimdg_kern_list ser_twistshift_ylimdg_list_0v[] = {
  {NULL, twistshift_ylimdg_2x_ser_p1_yshift_p1, NULL,},
  {NULL, twistshift_ylimdg_3x_ser_p1_yshift_p1, NULL,},
};
static const twistshift_fullcell_kern_list ser_twistshift_fullcell_list_0v[] = {
  {NULL, twistshift_fullcell_2x_ser_p1_yshift_p1, NULL,},
  {NULL, twistshift_fullcell_3x_ser_p1_yshift_p1, NULL,},
};

static const twistshift_xlimdg_kern_list ser_twistshift_xlimdg_list_2v[] = {
  {NULL, NULL, NULL,},
  {NULL, twistshift_xlimdg_3x2v_ser_p1_yshift_p1, NULL,},
};
static const twistshift_ylimdg_kern_list ser_twistshift_ylimdg_list_2v[] = {
  {NULL, NULL, NULL,},
  {NULL, twistshift_ylimdg_3x2v_ser_p1_yshift_p1, NULL,},
};
static const twistshift_fullcell_kern_list ser_twistshift_fullcell_list_2v[] = {
  {NULL, NULL, NULL,},
  {NULL, twistshift_fullcell_3x2v_ser_p1_yshift_p1, NULL,},
};

struct gkyl_bc_twistshift_kernels {
  twistshift_xlimdg_t xlimdg;
  twistshift_ylimdg_t ylimdg;
  twistshift_fullcell_t fullcell;
};

// Primary struct in this updater.
struct gkyl_bc_twistshift {
  int dir; // Direction of the BC is applied in.
  int shift_dir; // Direction of the shift.
  int shear_dir; // Direction the shift varies in (shear).
  enum gkyl_edge_loc edge; // Indicates if BC is for lowe/upper edge.
  struct gkyl_basis basis; // Basis the shifted field is defined with.
  struct gkyl_range local_r; // Local range.
  struct gkyl_range local_ext_r; // Local extended range.
  struct gkyl_rect_grid grid; // Grid the shifted field is defined in.
  bool use_gpu; // Whether to apply the BC on the GPU.

  struct gkyl_rect_grid shift_grid; // 1D grid in the direction of the shift.
  struct gkyl_range shift_r; // 1D range in the direction of the shift.

  struct gkyl_rect_grid shear_grid; // 1D grid in the direction of the shear.
  struct gkyl_range shear_r; // 1D range in the direction of the shear.

  struct gkyl_rect_grid ts_grid; // Grid the shift twistshift takes place in.
  struct gkyl_range ts_r; // Range the twistshift takes place in.
  int shift_dir_in_ts_grid; // Dimension the shift is in, in the TS grid.
  int shear_dir_in_ts_grid; // Dimension the shear is in, in the TS grid.

  int shift_poly_order; // Poly order of the DG representation of the shift.
  struct gkyl_basis shift_b; // 1D Basis for the DG shift.
  struct gkyl_array *shift; // DG shift.

  int *num_do; // Number of donors at each cell in shear_dir;
  int *shift_dir_idx_do; // Indices of donor cells, in the direction of the
                         // shift, for each cell in the TS grid.

  // *****    IN PROGRESS   ****
//  struct gkyl_bc_twistshift_kernels *kernels;  // kernels for sub-cell integrals.
//
//  int *ndonors_cum_cu;
//  int *remDir;
//  int *locDir;
//  int *remDir_do;
//  int *locDir_do;
//  long *locs;
//  long *locs_cu;
//  long *tar_locs;
//  long *tar_locs_cu;
//  struct gkyl_range *yrange;
//  struct gkyl_range *xrange;
//  struct gkyl_nmat *matsdo;
//  struct gkyl_nmat *matsdo_ho;
//  struct gkyl_nmat *vecsdo;
//  struct gkyl_nmat *vecstar;
};

void gkyl_bc_twistshift_choose_kernels(struct gkyl_basis basis, int cdim,
  struct gkyl_bc_twistshift_kernels *kers)
{
  int dim = basis.ndim;
  int vdim = dim - cdim;
  enum gkyl_basis_type basis_type = basis.b_type;
  int poly_order = basis.poly_order;
  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->xlimdg   = vdim==0?   ser_twistshift_xlimdg_list_0v[cdim-2].kernels[poly_order] :   ser_twistshift_xlimdg_list_2v[cdim-2].kernels[poly_order];
      kers->ylimdg   = vdim==0?   ser_twistshift_ylimdg_list_0v[cdim-2].kernels[poly_order] :   ser_twistshift_ylimdg_list_2v[cdim-2].kernels[poly_order];
      kers->fullcell = vdim==0? ser_twistshift_fullcell_list_0v[cdim-2].kernels[poly_order] : ser_twistshift_fullcell_list_2v[cdim-2].kernels[poly_order];
      return;
    default:
      assert(false);
      break;
  }
}

double
ts_grid_cell_boundary_in_dir(struct gkyl_rect_grid *grid, const int *idx, enum gkyl_edge_loc edge, int dir)
{
  // Get the coordinate of the cell boundary in specified direction.
  double xc[grid->ndim];
  gkyl_rect_grid_cell_center(grid, idx, xc);
  return edge == GKYL_LOWER_EDGE? xc[dir]-0.5*grid->dx[dir] : xc[dir]+0.5*grid->dx[dir];
}

void
ts_grid_cell_boundaries(struct gkyl_rect_grid *grid, const int *idx, double *cell_bounds)
{
  // Get the cell boundaries in every dimension. The array cell_bounds
  // must be a 2*grid->ndim array.
  for (int d=0; d<grid->ndim; d++) {
    cell_bounds[d*2] = ts_grid_cell_boundary_in_dir(grid, idx, GKYL_LOWER_EDGE, d);
    cell_bounds[d*2+1] = ts_grid_cell_boundary_in_dir(grid, idx, GKYL_UPPER_EDGE, d);
  }
}

static inline double
ts_p2l(double coord, double cell_center, double dx)
{
  // Transform a physical coordinate (coord) to the [-1,1] logical
  // space in a cell centered at cell_center and with length dx. 
  return 2.0*(coord - cell_center)/dx;
}

static inline double
ts_grid_length_in_dir(struct gkyl_rect_grid *grid, int dir)
{
  return grid->upper[dir] - grid->lower[dir];
}

double
ts_wrap_to_range(double val, double lower, double upper, bool pick_upper)
{
  // Wrap a number to range [lower,upper]. If pickUpper=true, output upper when
  // val is a multiple of upper. Otherwise multiples of upper wrap to lower.
  double L        = upper - lower;
  double disp     = fmod(val - lower, L);
  double vwrapped = lower + fmod(L + disp, L);
  double eps      = 1.e-12;
  if ( (lower-eps < vwrapped && vwrapped < lower + eps) ||
       (upper-eps < vwrapped && vwrapped < upper + eps) ) {
    if (pick_upper) 
      return upper;
    else
      return lower;
  }
  else
    return vwrapped;
}

long
ts_shift_dir_idx_do_linidx(const int *num_do,
  int shear_dir_idx, int shift_dir_idx, int shift_dir_num_cells)
{
  // Return the linear index to the first donor for the idx=(i,j) target cell,
  // in the shift_dir_idx_do array. We assume shift_dir_idx_do (whose dimensions
  // are Nx,Ny,num_do(i)) is in row-major order, and that it has num_do donors
  // at each cell in the shear_dir_in_ts_grid direction.
  long linc = 0;
  // Count the number of donors in cells with an idx in the shear dir lower
  // than this one. NOTE: the -1 here is because the idx is often 1-index
  // (since ghost cells are the 0th index) but num_do is only defined on the
  // local range. 
  for (int i=0; i<shear_dir_idx-1; i++)
    linc += num_do[i] * shift_dir_num_cells;

  // Count the number of donors in cells with the same idx in the shear dir.
  // NOTE: the -1 here is because the idx is often 1-index
  // (since ghost cells are the 0th index) but num_do is only defined on the
  // local range.
  linc += num_do[shear_dir_idx-1] * (shift_dir_idx-1);

  return linc;
}

void
ts_check_shifted_test_point(struct gkyl_bc_twistshift *up, const double *test_pt, const double *xc,
  const double *dx, double *shift_c, const int *idx, bool pick_lower, int *num_do_curr, gkyl_mem_buff shift_dir_idx_do_buff)
{
  int shear_idx[] = {idx[up->shear_dir_in_ts_grid]};
  int shift_idx[] = {idx[up->shift_dir_in_ts_grid]};

  int *shift_dir_idx_do_buff_ptr = (int *) gkyl_mem_buff_data(shift_dir_idx_do_buff);

  // Evaluate the shift at this test point.
  double test_pt_in_shear_dir = test_pt[up->shear_dir_in_ts_grid];
  double xc_in_shear_dir = xc[up->shear_dir_in_ts_grid];
  double dx_in_shear_dir = dx[up->shear_dir_in_ts_grid];
  double shift_at_pt = up->shift_b.eval_expand(
    &(double) {ts_p2l(test_pt_in_shear_dir, xc_in_shear_dir, dx_in_shear_dir)}, shift_c);

  // Find the index of the cell that owns the shifted point.
  double shifted_test_pt[] = { ts_wrap_to_range(test_pt[up->shift_dir_in_ts_grid] - shift_at_pt,
    up->ts_grid.lower[up->shift_dir_in_ts_grid], up->ts_grid.upper[up->shift_dir_in_ts_grid],
    false) }; // Shifted test point.
  int shift_dir_idx_test_pt[1];
//  gkyl_rect_grid_coord_idx(&up->shift_grid, shifted_test_pt, shift_dir_idx_test_pt); 
  gkyl_rect_grid_find_cell(&up->shift_grid, shifted_test_pt, pick_lower, (int[]) {-1}, shift_dir_idx_test_pt);

  // Get the linear index to the list of donors for this target.
  long linidx = ts_shift_dir_idx_do_linidx(up->num_do,
    shear_idx[0], shift_idx[0], up->ts_grid.cells[up->shift_dir_in_ts_grid]);
//  printf("idx=%d,%d | test_pt=%g,%g | shift=%g | shifted_test_pt=%g | shift_dir_idx_test_pt=%d | linidx=%ld | num_do=%d\n",idx[0],idx[1],test_pt[0],test_pt[1],shift_at_pt,shifted_test_pt[0],shift_dir_idx_test_pt[0],linidx,up->num_do[shear_idx[0]-1]);
  // If this donor is not in our list of donors, include it.
  bool donor_not_found = true;
  for (int k=0; k<num_do_curr[0]; k++) {
    if (shift_dir_idx_do_buff_ptr[linidx+k] == shift_dir_idx_test_pt[0]) {
      donor_not_found = false;
      break;
    }
  }
  if (donor_not_found) {
    if (num_do_curr[0] > 0) {
      // Insert one more int. Only if num_do_curr>0 because we already
      // allocated space for the first donor.
      size_t new_buff_sz = gkyl_mem_buff_size(shift_dir_idx_do_buff) + sizeof(int);
      shift_dir_idx_do_buff = gkyl_mem_buff_resize(shift_dir_idx_do_buff, new_buff_sz);
    }

    // Get the pointer again in case it changed.
    shift_dir_idx_do_buff_ptr = (int *) gkyl_mem_buff_data(shift_dir_idx_do_buff);
    shift_dir_idx_do_buff_ptr[linidx+num_do_curr[0]] = shift_dir_idx_test_pt[0];

    num_do_curr[0] += 1;
  }
}

void
ts_find_donors(struct gkyl_bc_twistshift *up)
{
  // Find the donor cells for each target cell in the TS grid.

  double delta_frac = 1.e-9; // Distance away from the boundary, as fraction of cell length.
  int num_test_pt[2] = {10, 10}; // Number of test points taken along each side of the cell.

  double step_sz[2] = {0.0}; // Size of the step between test points.
  double delta[2] = {0.0}; // Space between cell boundary and test points.
  for (int d=0; d<2; d++) {
    delta[d] = delta_frac*up->ts_grid.dx[d];
    step_sz[d] = (up->ts_grid.dx[d] - 2.0*delta[d])/(num_test_pt[d]-1);
  }

  // Number of donors at each cell of the shear direction.
  up->num_do = gkyl_malloc(up->shear_r.volume * sizeof(int));
  for (int i=0; i<up->shear_r.volume; i++)
    up->num_do[i] = 0;

  // Temporary buffer to store donors at. One for each cell for now.
  size_t curr_buff_sz = up->ts_r.volume * sizeof(int);
  gkyl_mem_buff shift_dir_idx_do_buff = gkyl_mem_buff_new(curr_buff_sz);
  printf("\n");

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->ts_r);
  while (gkyl_range_iter_next(&iter)) {

    // Get the cell boundaries and cell center.
    double cell_b[4] = {0.0}; // Cell boundaries, x lo and up, y lo and up;
    double xc[2] = {0.0}; // Cell center.
    ts_grid_cell_boundaries(&up->ts_grid, iter.idx, cell_b);
    gkyl_rect_grid_cell_center(&up->ts_grid, iter.idx, xc);

    int shear_idx[] = {iter.idx[up->shear_dir_in_ts_grid]};
    int shift_idx[] = {iter.idx[up->shift_dir_in_ts_grid]};
    long shift_loc = gkyl_range_idx(&up->shift_r, shear_idx);
    double *shift_c = gkyl_array_fetch(up->shift, shift_loc);

    int num_do_curr = 0;

    for (int dC=0; dC<2; dC++) { // dC=0: x=const, dC=1: y=const (boundaries).
      for (int xS=0; xS<2; xS++) { // xS=0: lower, xS=1 upper (boundary).

        double test_pt[2] = {0.0}; // Test point to shift.
        for (int d=0; d<2; d++)
          test_pt[d] = cell_b[2*d]+delta[d];

        // Search first shifted point. Use pick_lower=false in find_cell unless
        // searching for points along a x=const line near the upper y-boundary. 
        bool pick_lower = dC==0 && xS==1;
        test_pt[dC] += xS*(up->ts_grid.dx[dC]-2.0*delta[dC]);

        // Shift the test point, find the cell that contains it, and if we
        // haven't included it yet, add it to our list of donors.
        ts_check_shifted_test_point(up, test_pt, xc, up->ts_grid.dx,
          shift_c, iter.idx, pick_lower, &num_do_curr, shift_dir_idx_do_buff);

        // Search for other shifted points along this line. 
        int step_dim = (dC+1) % 2;
        for (int sI=1; sI<num_test_pt[step_dim]; sI++) {
          test_pt[step_dim] += sI*step_sz[step_dim];

          // Shift the test point, find the cell that contains it, and if we
          // haven't included it yet, add it to our list of donors.
          ts_check_shifted_test_point(up, test_pt, xc, up->ts_grid.dx,
            shift_c, iter.idx, pick_lower, &num_do_curr, shift_dir_idx_do_buff);
        }
      }
    }

    up->num_do[shear_idx[0]-1] = num_do_curr;
  }

  // Copy the donor list to the persistent object and release the buffer.
  long buff_sz = gkyl_mem_buff_size(shift_dir_idx_do_buff);
  up->shift_dir_idx_do = gkyl_malloc(buff_sz);
  int *shift_dir_idx_do_buff_ptr = (int *) gkyl_mem_buff_data(shift_dir_idx_do_buff);
  memcpy(up->shift_dir_idx_do, shift_dir_idx_do_buff_ptr, buff_sz);
  gkyl_mem_buff_release(shift_dir_idx_do_buff);
}

//void gkyl_bc_twistshift_integral_xlimdg(struct gkyl_bc_twistshift *up,
//  double sFac, const double *xLimLo, const double *xLimUp, double yLimLo, double yLimUp,
//  double dyDo, double yOff, const double *ySh, int cellidx, int doidx) {
//  
//  size_t linidx = 0;
//  for(int i = 0; i<cellidx-1; i++)
//    linidx += up->ndonors[i];
//  
//  linidx += doidx;
//  struct gkyl_mat tsmat = gkyl_nmat_get(up->matsdo_ho, linidx);
//  up->kernels->xlimdg(sFac, xLimLo, xLimUp, yLimLo, yLimUp, dyDo, yOff, ySh, &tsmat);
//}
//
//void gkyl_bc_twistshift_integral_ylimdg(struct gkyl_bc_twistshift *up,
//  double sFac, double xLimLo, double xLimUp, const double *yLimLo, const double *yLimUp,
//  double dyDo, double yOff, const double *ySh, int cellidx, int doidx) {
//  size_t linidx = 0;
//  for(int i = 0; i<cellidx-1; i++)
//    linidx += up->ndonors[i];
//  
//  linidx += doidx;
//  struct gkyl_mat tsmat = gkyl_nmat_get(up->matsdo_ho, linidx);
//  up->kernels->ylimdg(sFac, xLimLo, xLimUp, yLimLo, yLimUp, dyDo, yOff, ySh, &tsmat);
//}
//
//void gkyl_bc_twistshift_integral_fullcelllimdg(struct gkyl_bc_twistshift *up,
//  double dyDo, double yOff, const double *ySh, int cellidx, int doidx) {
//  size_t linidx = 0;
//  for(int i = 0; i<cellidx-1; i++)
//    linidx += up->ndonors[i];
//  
//  linidx += doidx;
//  struct gkyl_mat tsmat = gkyl_nmat_get(up->matsdo_ho, linidx);
//  up->kernels->fullcell(dyDo, yOff, ySh, &tsmat);
//}
//
//void gkyl_bc_twistshift_precalc_mats()
//{
//
//  // y-index of the reference target used to precalc matrices. For positive(negative)
//  // yShift idx=1(last) might be better, but ideally it shouldn't matter.
//  int yIdxTar = 1;
//
//  int idxTar[2]; // Target index.
//  idxTar[1] = yIdxTar;
//
//  // Create a range over x.
//  up->remDir_do[0] = 0;
//  for(int d=1; d<up->local_r->ndim; d++)
//    up->remDir_do[d] = 1;
//  for(int d=0; d<up->local_r->ndim; d++){
//    up->locDir_do[d] = up->local_r->lower[d];
//  gkyl_range_deflate(up->xrange, up->local_r, up->remDir_do, up->locDir_do);
//
//  struct gkyl_range_iter iter;
//  gkyl_range_iter_init(&iter, up->xrange);
//  while (gkyl_range_iter_next(&iter)) {
//
//    gkyl_rect_grid_cell_center(up->grid2dconst struct gkyl_rect_grid *grid,
//  const int *idx, double *xc)
//  }
//}
