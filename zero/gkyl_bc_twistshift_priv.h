#pragma once

// Private header for bc_twistshift updater, not for direct use in user code.
//
// Notes:
//   a) Hard-coded parameters:
//     - wrap_to_range: eps.
//     - find_donors: delta_frac, num_test_pt.
//     - find_intersect: tol, max_iter, num_steps.
//     - calc_mats: shift_dir_idx_tar.
//
// Module functions (called in bc_twistshift.c):
//   - ts_find_donors.
//   - ts_calc_mats.
//
// Helper functions (all have ts_ prefix):
//   - grid_cell_boundary_in_dir: cell boundary coordinate in given dir.
//   - grid_cell_boundaries: get all cell boundary coords.
//   - p2l: physical to logical transform.
//   - interval_dx_and_xc: compute length and center of an interval.
//   - grid_length_in_dir: length of the grid in given dir.
//   - wrap_to_range: wrap a number to a range assuming periodicity.
//   - shift_dir_idx_do_linidx: linear index to first donor of a given target
//                              cell in shift_dir_idx_do.
//   - check_shifted_test_point: evaluate a shifted point's cell as a
//                               potential donor cell.
//   - root_find: Finds the root of a given function.
//   - shifted_coord_loss_func: Loss function used to find where yTar-S
//                              intersects yDo.
//   - sign: return the sign of a double.
//   - ts_donor_target_offset: offset between donor and target cells.
//   - find_intersect: finds the intersection of yTar-S and yDo.
//   - comp_to_phys: transform a computational to a physical coord.
//   - nod2mod_proj_1d: evaluate a 1D function at nodes and do a n2m transform
//                      to get the coefficients of the DG representation.
//   - ts_integral_xlimdg: subcell integral with variable x limits.
//   - ts_integral_ylimdg: subcell integral with variable y limits.
//   - ts_integral_fullcelllimdg: integral over the whole cell.
//   - one: return 1 (for projections).
//   - minus_one: return -1 (for projections).
//   - shift_coord_shifted_log: coordinate in shift_dir shifted and transformed
//                              to logical space.
//   

#include <gkyl_bc_twistshift.h>
#include <gkyl_bc_twistshift_gyrokinetic_kernels.h>
#include <assert.h>
#include <gkyl_mat.h>
#include <gkyl_math.h>
#include <gkyl_eval_on_nodes.h>
#include <string.h> // memcpy

// Indices in 4-element cell boundary array.
#define cellb_lo(dir) (2*dir)
#define cellb_up(dir) (2*dir+1)

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
  evalf_t shift_func; // Function defining the shift.
  void *shift_func_ctx; // Context for shift_func.
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

  struct gkyl_bc_twistshift_kernels *kernels;  // kernels for sub-cell integrals.

  // Projection object used in constructing the matrices.
  struct gkyl_eval_on_nodes *ev_on_nod1d;
  // Evaluations of a function at 1D nodes.
  struct gkyl_array *func_nod1d;

  // *****    IN PROGRESS   ****
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

void
ts_interval_dx_and_xc(const double *interval, double *dx, double *xc)
{
  // Compute the lenth (dx) and center (xc) of [interval[0], interval[1]].
  double lo = interval[0], up = interval[1];
  dx[0] = up - lo;
  xc[0] = 0.5*(up + lo);
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

struct gkyl_qr_res
ts_root_find(double (*func)(double,void*), void *ctx, const double *lims, int max_iter, double tol)
{
  // Use a Ridder's root finder to find the root of func in the interval
  // [lims[0],lims[1]] down to a tolerance 'tol'. Return the interval limit
  // if the function is smaller than the tolerance there. Return nil if the
  // function does not change sign in the interval (interval doesn't contain the root).
  double funcLo = func(lims[0], ctx), funcUp = func(lims[1], ctx);
  if (fabs(funcLo) < tol)
    return (struct gkyl_qr_res) {.res=lims[0], .nevals=2, .status=0};
  else if (fabs(funcUp) < tol)
    return (struct gkyl_qr_res) {.res=lims[1], .nevals=2, .status=0};
  else {
    if (funcLo*funcUp < 0)
      return gkyl_ridders(func, ctx, lims[0], lims[1], funcLo, funcUp, max_iter, tol);
    else
      return (struct gkyl_qr_res) {.nevals=2, .status=1};
  }
  return (struct gkyl_qr_res) {.nevals=2, .status=1};
}

struct ts_shifted_coord_loss_func_ctx {
  double shiftCoordTar; // Target coordinate in shift_dir.
  double shiftCoordDo; // Donor coordinate in shift_dir.
  double shiftDirL; // Length of the domain in shift_dir.
  int periodicCopyIdx; // Used to search a periodic copy of the domain (signed).
  evalf_t shift_func; // Shift function.
  void *shift_func_ctx; // Shift function context.
};

double ts_shifted_coord_loss_func(double shearCoord, void *ctx)
{
  // Loss function used to find the shear
  // coord
  struct ts_shifted_coord_loss_func_ctx *tsctx = ctx;
  double shift;
  tsctx->shift_func(0.0, (double[]){shearCoord}, &shift, tsctx->shift_func_ctx);
  return tsctx->shiftCoordTar - shift
    - (tsctx->shiftCoordDo - tsctx->periodicCopyIdx * tsctx->shiftDirL);
}

int static inline
ts_sign(double a)
{
  if (a < 0.0)
    return -1;
  else if (a > 0.0)
    return 1;
  else
    return 0;
}

double
ts_donor_target_offset(struct gkyl_bc_twistshift *up, const double *xc_do, const double *xc_tar) {
  // y-offset between the donor and the target cell (yDo-yTar), in the direction of the shift.
  //   xc_do:  cell center coordinates of donor cell.
  //   xc_tar: cell center coordinates of target cell.
  int shear_dir = up->shear_dir_in_ts_grid;
  int shift_dir = up->shift_dir_in_ts_grid;

  double shift;
  up->shift_func(0.0, (double[]){xc_do[up->shear_dir]}, &shift, up->shift_func_ctx);

  int shift_sign = ts_sign(shift);
  double shift_dir_L = up->ts_grid.upper[up->shift_dir] - up->ts_grid.lower[up->shift_dir];
                             
  // The idea here is that we keep shifting the donor cell center until it is in a
  // periodic copy of our domain which overlaps with the shifted target cell center.
  double xs_shifted_do = xc_do[shift_dir];
  double xs_shifted_tar = xc_tar[shift_dir] - shift;
  bool keep_shifting = true;
  while (keep_shifting) {
    double xs_shifted_dolo = xs_shifted_do - shift_dir_L/2.0;
    double xs_shifted_doup = xs_shifted_do + shift_dir_L/2.0;
    if (xs_shifted_dolo <= xs_shifted_tar && xs_shifted_tar <= xs_shifted_doup) {
      keep_shifting = false;
      break;
    }
    else
      xs_shifted_do = xs_shifted_do - shift_sign*shift_dir_L;
  }
  return xc_tar[shift_dir] - xs_shifted_do;
}

struct gkyl_qr_res
ts_find_intersect(struct gkyl_bc_twistshift *up, double shiftCoordTar, double shiftCoordDo,
  const double *shearDirBounds, const double *shiftDirLimits)
{
  // Given a y-coordinate of the target cell (yTar), and a y-coordinate
  // of the donor cell (yDo), find the x-coordinate of the point where
  // the yTar-yShift(x) and y=yDo lines intersect.
  //  yTar:    target cell y coordinate.
  //  yDo:     donor cell y coordinate.
  //  xBounds: search in the interval [xBounds[1],xBounds[2]].
  //  yLims:   lower and upper limits of the grid.
  // If y-yShift-yDo=0 has no roots, it is possible that y-yShift intersects a periodic
  // copy of this domain. Check for such cases by looking for the roots of
  // yTar-yShift-(yDo-N*Ly)=0 where Ly is the length of the domain along y and N is an integer.
  double tol = 1.e-13;
  int max_iter = 100;

  double shiftDirL = shiftDirLimits[1] - shiftDirLimits[0];

  struct ts_shifted_coord_loss_func_ctx func_ctx = {
    .shiftCoordTar = shiftCoordTar,
    .shiftCoordDo = shiftCoordDo,
    .shiftDirL = shiftDirL,
    .periodicCopyIdx = 0,
    .shift_func = up->shift_func,
    .shift_func_ctx = up->shift_func_ctx,
  };
  struct gkyl_qr_res rootfind_res = ts_root_find(ts_shifted_coord_loss_func, &func_ctx,
    shearDirBounds, max_iter, tol);

  if (rootfind_res.status == 1) {
    // Maybe yTar-ySh intersects y=yDo in a periodic copy of the domain. Find roots of
    // yTar-ySh-(yDo-nP*Ly)=0 where Ly is the y-length of the domain and nP is an integer.

    // Evaluate yTar-ySh at some points in [xBounds.lo,xBounds.up] and obtain potential nP's.
    int num_steps = 10;
    double step_sz = (shearDirBounds[1]-shearDirBounds[0])/num_steps;
    int nP[num_steps+1], num_unique_nP = 0;
    for (int sI=0; sI<num_steps+1; sI++) {
      double xp = shearDirBounds[0] + sI*step_sz;
      double xp_shift;
      up->shift_func(0.0, (double[]){xp}, &xp_shift, up->shift_func_ctx);
      double ex_shift = xp_shift<0.? ((shiftCoordTar-xp_shift)-shiftDirLimits[0])/shiftDirL
                                   : (shiftDirLimits[1]-(shiftCoordTar-xp_shift))/shiftDirL;
      double nP_new = ts_sign(xp_shift)*floor(fabs(ex_shift));
      // If we haven't accounted for this nP, add it to our list.
      bool nP_not_found = true;
      for (int n=0; n<num_unique_nP; n++) {
        if (nP[n] == nP_new) {
          nP_not_found = false;
          break;
        }
      }
      if (nP_not_found) {
        nP[num_unique_nP] = nP_new;
        num_unique_nP++;
      }
    }

    for (int n=0; n<num_unique_nP; n++) {
      func_ctx.periodicCopyIdx = nP[n];
      rootfind_res = ts_root_find(ts_shifted_coord_loss_func, &func_ctx,
        shearDirBounds, max_iter, tol);
      if (rootfind_res.status == 0) break;
    }
  }

  return rootfind_res;
}

struct ts_val_found {
  bool status; // =0 not found, =1 found.
  double value; // value found.
};

static inline void
ts_comp_to_phys(int ndim, const double *eta,
  const double * GKYL_RESTRICT dx, const double * GKYL_RESTRICT xc,
  double* GKYL_RESTRICT xout)
{
  for (int d=0; d<ndim; ++d) xout[d] = 0.5*dx[d]*eta[d]+xc[d];
}

void
ts_nod2mod_proj_1d(struct gkyl_bc_twistshift *up, evalf_t func, void *func_ctx, const double *interval, double *out)
{
  // Project 'func' onto 1D DG basis in [interval[0], interval[1]]
  // evaluating at nodes and transforming to modal representation.
  //   func: 1D scalar function to be projected.
  //   interval: limits of the interval in which to project the function.
  //   out: DG field output.
  double dx[] = {0.0}, xc[] = {0.0}, xmu[] = {0.0};
  ts_interval_dx_and_xc(interval, dx, xc);

  for (int i=0; i<up->shift_b.num_basis; ++i) {
    ts_comp_to_phys(1, gkyl_eval_on_nodes_fetch_node(up->ev_on_nod1d, i),
      dx, xc, xmu);
    func(0.0, xmu, (double *)gkyl_array_fetch(up->func_nod1d,i), func_ctx);
  }
  gkyl_eval_on_nodes_nod2mod(up->ev_on_nod1d, up->func_nod1d, out); 
}

void
ts_integral_xlimdg(struct gkyl_bc_twistshift *up, double sFac, const double *xLimLo,
  const double *xLimUp, double yLimLo, double yLimUp, double dyDo, double yOff,
  const double *ySh, struct gkyl_mat *mat_do) {
  // Populate a matrix (mat_do) with a sub-cell integral that has variably x limits
  // represented by a DG polynomial, and a y-integral that goes from yLimLo to yLimUp.
  //   up: BC updater.
  //   sFac: +/-1 factor to add or subtract this subcell integral.
  //   xLimLo: DG representation of the lower x-limit.
  //   xLimUp: DG representation of the upper x-limit.
  //   yLimLo: lower y-limit.
  //   yLimUp: upper y-limit.
  //   dyDo: Cell length along y.
  //   yOff: Offset along y.
  //   ySh: DG representation of the y shift.
  //   mat_do: donor matrix.
  up->kernels->xlimdg(sFac, xLimLo, xLimUp, yLimLo, yLimUp, dyDo, yOff, ySh, mat_do);
}

void
ts_integral_ylimdg(struct gkyl_bc_twistshift *up, double sFac, double xLimLo, double xLimUp,
  const double *yLimLo, const double *yLimUp, double dyDo, double yOff, const double *ySh,
  struct gkyl_mat *mat_do) {
  // Populate a matrix (mat_do) with a sub-cell integral that has variable y limits
  // represented by a DG polynomial, and a x-integral that goes from xLimLo to xLimUp.
  //   up: BC updater.
  //   sFac: +/-1 factor to add or subtract this subcell integral.
  //   xLimLo: lower x-limit.
  //   xLimUp: upper x-limit.
  //   yLimLo: DG representation of the lower y-limit.
  //   yLimUp: DG representation of the upper y-limit.
  //   dyDo: Cell length along y.
  //   yOff: Offset along y.
  //   ySh: DG representation of the y shift.
  //   mat_do: donor matrix.
  up->kernels->ylimdg(sFac, xLimLo, xLimUp, yLimLo, yLimUp, dyDo, yOff, ySh, mat_do);
}

void
ts_integral_fullcelllimdg(struct gkyl_bc_twistshift *up, double dyDo, double yOff,
  const double *ySh, struct gkyl_mat *mat_do) {
  // Populate a matrix (mat_do) with the full-cell integral.
  //   up: BC updater.
  //   sFac: +/-1 factor to add or subtract this subcell integral.
  //   dyDo: Cell length along y.
  //   yOff: Offset along y.
  //   ySh: DG representation of the y shift.
  //   mat_do: donor matrix.
  up->kernels->fullcell(dyDo, yOff, ySh, mat_do);
}

void
ts_one(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = 1.0;
}

void
ts_minus_one(double t, const double *xn, double *fout, void *ctx)
{
  fout[0] = -1.0;
}

struct ts_shift_coord_shifted_log_ctx {
  double shift_coord_tar; // Target coordinate in shift_dir.
  int shift_sign_fac;
  const double *xc_do, *xc_tar; // Cell centers (donor and target).
  double *dx; // Cell lengths.
  bool pick_upper;
  int shear_dir, shift_dir; // Shear and shift directions.
  double shift_dir_bounds[2]; // Domain boundaries in shift_dir.
  evalf_t shift_func; // Shift function.
  void *shift_func_ctx; // Shift function context.
};

void
ts_shift_coord_shifted_log(double t, const double *xn, double *fout, void *ctx)
{
  // Given a logical space x coordinate (xi) and a (physical) y-coordinate in the target cell,
  // compute the shifted y-coordinate in the logical space of the donor cell (eta \in [-1,1]).
  //   xi:    logical space x coordinate.
  //   yTar:  physical y-coordinate in target cell.
  //   pmSh:  factor multiplying the y-shift (+/- 1).
  //   xcDo:  cell center coordinates of donor cell.
  //   xcTar: cell center coordinates of target cell.
  //   dx:    cell lengths.
  //   pickUpper: boolean indicating if wrapping function should return upper/lower boundary.

  double xi = xn[0];

  struct ts_shift_coord_shifted_log_ctx *tsctx = ctx;
  double shift_coord_tar = tsctx->shift_coord_tar;
  int shift_sign_fac = tsctx->shift_sign_fac;
  const double *xc_do = tsctx->xc_do, *xc_tar = tsctx->xc_tar;
  double *dx = tsctx->dx;
  bool pick_upper = tsctx->pick_upper;
  int shear_dir = tsctx->shear_dir, shift_dir = tsctx->shift_dir;
  double *shift_dir_bounds = tsctx->shift_dir_bounds;

  double shear_coord_phys = xc_tar[shear_dir] + 0.5*dx[shear_dir]*xi;

  double shift;
  tsctx->shift_func(0.0, (double[]){shear_coord_phys}, &shift, tsctx->shift_func_ctx);

  double shift_coord_shifted = shift_coord_tar - shift_sign_fac * shift;
  shift_coord_shifted = ts_wrap_to_range(shift_coord_shifted, shift_dir_bounds[0], shift_dir_bounds[1], pick_upper);

  fout[0] = ts_p2l(shift_coord_shifted, xc_do[shift_dir], dx[shift_dir]);
}

void
ts_subcellint_sxv_sxvi(struct gkyl_bc_twistshift *up, struct ts_val_found *inter_pts,
  const double *xc_do, const double *xc_tar, const double *cellb_do, const double *cellb_tar,
  bool is_upper_shift_dir_cell, const double *shift_c, struct gkyl_mat *mat_do)
{
  // Perform subcell integral sxv or sxvi, using fixed x-limits and variable y limits.
  //   inter_pts: intersections y_{j_tar-/+1/2}-yShift and y_{j_do-/+1/2} (lower/upper y-boundaries of donor cell).
  //   xc_do:     cell center of donor cell.
  //   xc_tar:    cell center of target cell.
  //   cellb_do:  boundaries of target cell.
  //   cellb_tar: boundaries of target cell.
  //   shift_c:   DG coefficients of the shift.
  //   mat_do:    current donor matrix.

  double xi_b[] = {-1.0, 1.0};  // Limits of xi integral.

  double shift_dir_bounds[] = {up->ts_grid.lower[up->shift_dir_in_ts_grid],
                               up->ts_grid.upper[up->shift_dir_in_ts_grid]};
  double shift;
  up->shift_func(0.0, (double[]){xc_do[up->shear_dir_in_ts_grid]}, &shift, up->shift_func_ctx);

  double shifted_coord = cellb_tar[cellb_lo(up->shift_dir_in_ts_grid)] - shift;
  shifted_coord = ts_wrap_to_range(shifted_coord, shift_dir_bounds[0], shift_dir_bounds[1],
    is_upper_shift_dir_cell);

  struct ts_shift_coord_shifted_log_ctx eta_lims_ctx = {
    .shift_sign_fac = 1,
    .xc_do = xc_do,
    .xc_tar = xc_tar,
    .dx = up->ts_grid.dx,
    .pick_upper = is_upper_shift_dir_cell,
    .shear_dir = up->shear_dir_in_ts_grid,
    .shift_dir = up->shift_dir_in_ts_grid,
    .shift_dir_bounds = {up->ts_grid.lower[up->shift_dir_in_ts_grid], 
                         up->ts_grid.upper[up->shift_dir_in_ts_grid]},
    .shift_func     = up->shift_func,
    .shift_func_ctx = up->shift_func_ctx,
  };

  evalf_t eta_lims[2]; // Table of functions definting the limits of eta integral.
  if ( cellb_do[cellb_lo(up->shift_dir_in_ts_grid)] <= shifted_coord &&
       shifted_coord <= cellb_do[cellb_up(up->shift_dir_in_ts_grid)] ) {
    // sxv integral.
    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_lo(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_shift_coord_shifted_log;
    eta_lims[1] = ts_one;
  }
  else {
    // sxvi integral.
    eta_lims_ctx.shift_coord_tar = cellb_tar[cellb_up(up->shift_dir_in_ts_grid)];
    eta_lims[0] = ts_minus_one;
    eta_lims[1] = ts_shift_coord_shifted_log;
  }

  double etalo_xi[up->shift_b.num_basis], etaup_xi[up->shift_b.num_basis];
  ts_nod2mod_proj_1d(up, eta_lims[0], &eta_lims_ctx, xi_b, etalo_xi);
  ts_nod2mod_proj_1d(up, eta_lims[1], &eta_lims_ctx, xi_b, etaup_xi);
  // Offset between cell centers in direction of the shift.
  double xs_off = ts_donor_target_offset(up, xc_do, xc_tar);

  ts_integral_ylimdg(up, 1.0, xi_b[0], xi_b[1], etalo_xi, etaup_xi,
    up->ts_grid.dx[up->shift_dir_in_ts_grid], xs_off, shift_c, mat_do);
}

void ts_calc_mats(struct gkyl_bc_twistshift *up, struct gkyl_nmat *matsdo)
{

  // y-index of the reference target used to precalc matrices. For positive(negative)
  // yShift idx=1(last) might be better, but ideally it shouldn't matter.
  int shift_dir_idx_tar = 1;

  double shift_dir_lims[] = {up->ts_grid.lower[up->shift_dir_in_ts_grid],
                             up->ts_grid.upper[up->shift_dir_in_ts_grid]};

  // Create an eval_on_nodes updater to use its nodes and functions (but not
  // the whole advance method).
  up->ev_on_nod1d = gkyl_eval_on_nodes_new(&up->shear_grid, &up->shift_b, 1, ts_one, NULL);
  // Create an array to store evaluations of a function at 1D nodes.
  up->func_nod1d = gkyl_array_new(GKYL_DOUBLE, 1, up->shift_b.num_basis);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->shear_r);
  while (gkyl_range_iter_next(&iter)) {

    // Get the cell boundaries and cell center.
    int idx_tar[2]; // Target index.
    idx_tar[up->shift_dir_in_ts_grid] = shift_dir_idx_tar;
    idx_tar[up->shear_dir_in_ts_grid] = iter.idx[0];
    double cellb_tar[4] = {0.0}; // Cell boundaries, x lo and up, y lo and up;
    double xc_tar[2] = {0.0}; // Cell center.
    ts_grid_cell_boundaries(&up->ts_grid, idx_tar, cellb_tar);
    gkyl_rect_grid_cell_center(&up->ts_grid, idx_tar, xc_tar);

    long shift_loc = gkyl_range_idx(&up->shift_r, iter.idx);
    double *shift_c = gkyl_array_fetch(up->shift, shift_loc);

    long linidx_do = ts_shift_dir_idx_do_linidx(up->num_do,
      iter.idx[0], shift_dir_idx_tar, up->ts_grid.cells[up->shift_dir_in_ts_grid]);
    int *shift_dir_idx_do_ptr = &up->shift_dir_idx_do[linidx_do];

    for (int iC=0; iC<up->num_do[iter.idx[0]-1]; iC++){
      int idx_do[2]; // Target index.
      idx_do[up->shift_dir_in_ts_grid] = shift_dir_idx_do_ptr[iC];
      idx_do[up->shear_dir_in_ts_grid] = iter.idx[0];

      double cellb_do[4] = {0.0}; // Cell boundaries, x lo and up, y lo and up;
      double xc_do[2] = {0.0}; // Cell center.
      ts_grid_cell_boundaries(&up->ts_grid, idx_do, cellb_do);
      gkyl_rect_grid_cell_center(&up->ts_grid, idx_do, xc_do);

      // Get the matrix we are presently assigning.
      struct gkyl_mat mat_do = gkyl_nmat_get(matsdo, linidx_do+iC);

      // Find the points where y_{j_tar-/+1/2}-yShift intersect the y=y_{j_do-/+1/2} lines.
      // Also record the number and indices of points found/not found.
      struct ts_val_found inter_pts[4] = {};
      int num_inter_pts_found = 0, num_inter_pts_not_found = 4;
      int inter_pts_found_idxs[4];
      int inter_pts_not_found_idxs[4];
      for (int i=0; i<2; i++) { // Loop over j_tar-/+1/2
        for (int j=0; j<2; j++) { // Loop over j_do-/+1/2
          double shift_dir_coord_tar = cellb_tar[2*up->shift_dir_in_ts_grid+i];
          double shift_dir_coord_do = cellb_do[2*up->shift_dir_in_ts_grid+j];
          struct gkyl_qr_res inter_res = ts_find_intersect(up, shift_dir_coord_tar, shift_dir_coord_do,
            (double[]) {cellb_tar[cellb_lo(up->shear_dir_in_ts_grid)],cellb_tar[cellb_up(up->shear_dir_in_ts_grid)]},
            shift_dir_lims);
          int ip_linc = i*2+j;
          inter_pts[ip_linc].status = inter_res.status == 0;
          if (inter_res.status == 0) {
            inter_pts[ip_linc].value = inter_res.res;
            inter_pts_found_idxs[num_inter_pts_found] = ip_linc;
            num_inter_pts_found++;
          }
          else {
            inter_pts_not_found_idxs[num_inter_pts_not_found] = ip_linc;
            num_inter_pts_not_found--;
          }
        }
      }

      bool is_upper_shift_dir_cell = idx_do[up->shift_dir_in_ts_grid] == up->ts_grid.cells[up->shift_dir_in_ts_grid];

      if (num_inter_pts_found == 4) {
        // sN: all intersections are found at this cell.
//        ts_subcellint_sN(up, inter_pts, xc_do, xc_tar, cellb_do, cellb_tar, is_upper_shift_dir_cell, shift_c, &mat_do);
      }
      else {
        if (num_inter_pts_found == 0) {
          // sxv:  y_{j_tar-1/2}-yShift crosses x_{i-/+1/2}.
          // sxvi: y_{j_tar+1/2}-yShift crosses x_{i-/+1/2}.
          ts_subcellint_sxv_sxvi(up, inter_pts, xc_do, xc_tar, cellb_do, cellb_tar, is_upper_shift_dir_cell, shift_c, &mat_do);
        }
      }
    }

  }

  gkyl_array_release(up->func_nod1d);
  gkyl_eval_on_nodes_release(up->ev_on_nod1d);
}
