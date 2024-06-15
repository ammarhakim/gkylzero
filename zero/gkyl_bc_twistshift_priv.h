#pragma once

// Private header for bc_twistshift updater, not for direct use in user code.

#include <gkyl_bc_twistshift.h>
#include <gkyl_bc_twistshift_gyrokinetic_kernels.h>
#include <assert.h>
#include <gkyl_mat.h>

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
  int dir;
  int do_dir;
  int shift_dir;
  enum gkyl_edge_loc edge;
  const struct gkyl_basis *basis;
  bool use_gpu;
  struct gkyl_bc_twistshift_kernels *kernels;  // kernels.
  struct gkyl_bc_twistshift_kernels *kernels_cu;  // device copy.
  const struct gkyl_rect_grid *grid;
  const int *ndonors;
  int *ndonors_cum_cu;
  const int *cells_do; // y indices of donor cells for each x and y
  int *remDir;
  int *locDir;
  int *remDir_do;
  int *locDir_do;
  long *locs;
  long *locs_cu;
  long *tar_locs;
  long *tar_locs_cu;
  const struct gkyl_range *local_range_ext;
  const struct gkyl_range *local_range_update;
  struct gkyl_range *yrange;
  struct gkyl_range *xrange;
  struct gkyl_nmat *matsdo;
  struct gkyl_nmat *matsdo_ho;
  struct gkyl_nmat *vecsdo;
  struct gkyl_nmat *vecstar;
};

void gkyl_bc_twistshift_choose_kernels_cu(const struct gkyl_basis *basis, int cdim,
  struct gkyl_bc_twistshift_kernels *kers)
{
  int dim = basis->ndim;
  int vdim = dim - cdim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;
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

void gkyl_bc_twistshift_choose_kernels(const struct gkyl_basis *basis, int cdim,
  struct gkyl_bc_twistshift_kernels *kers)
{
  int dim = basis->ndim;
  int vdim = dim - cdim;
  enum gkyl_basis_type basis_type = basis->b_type;
  int poly_order = basis->poly_order;
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

// PROTOTYPING code for reconstructing matrix from Lua file
//
// Module functions (to start)
//
// - getDonors
// - set_yShiftF
//
// Local helper functions (to start).
// - p2l
// - wrapNum
// - yShiftedLog
// - subCellInt_sxvORsxvi
//

// Important to add these things first!
// yshift --> eval on nodes (g2 updater)
// ndonors : see tsFun.getDonors() number of donor cells
// cells_do : index of each donor

// Returning double depends on how double is used.
// If array is needed, use void.

// stFun.getDonors returns doCells
int* get_donors(const struct gkyl_rect_grid *grid, const struct gkyl_range *local_range,
		const struct gkyl_array *ySh_arr, const struct gkyl_range *ySh_range,
		const struct gkyl_basis *ySh_basis) {

  // create buffer here
  // copy contents of buffer into regular int array

  int ndim = 2; // dimension of xy grid
  double deltaFrac = 1.e-9;
  int numSteps = 10; // consider 10 points on each boundary

  int ySh_nc = ySh_arr.num_comp;
  double yShBasisEv[ySh_nc];
  double p2l_eval, yShEv;

  double cellLim_lo[] = {0., 0.};
  double cellLim_up[] = {0., 0.};
  double stepSz[] = {0., 0.};
  double delta[] = {0., 0.};
  double evPoint[] = {0., 0.};
  double newP[] = {0., 0.};
  int newIdx[] = {0, 0};
  
  int idxShifted[ndim], idxP[ndim];
  double searchPoint[ndim];

  for (int d=0; d<ndim; ++d) {
    idxShifted[d] = 1;
    idxP[d] = 1;
    searchPoint[d] = grid.lower[d]*deltaFrac*grid.dx[d];
  }

  // create xyrange
  int loc_dim = local_range.ndim; // dimension of grid
  struct gkyl_range xyrange;
  int remDir[loc_dim] = {0};
  int locDir[loc_dim];
  for (int d=ndim; d<loc_dim; ++d) remDir[d] = 1;
  for(int d=0;i<loc_dim;d++){
    up->locDir[d] = local_range->lower[i];
  }
  deflate(xy_range, local_range, remDir, locDir);

  double xc[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, xy_range);
  
  while (gkyl_range_iter_next(&iter)) {

    // local doCellsC = doCells[idx[1]][idx[2]]

    /* for d = 1,2 do */
    /*    cellLim[d].lo, cellLim[d].up = grid:cellLowerInDir(d), grid:cellUpperInDir(d) */
    /*    stepSz[d] = grid:dx(d)/numSteps[d] */
    /*    delta[d]  = deltaFrac*grid:dx(d) */
    /* end */
    // set cellLim_lo, cellLim_up based on current cell
    gkyl_rect_grid_cell_center(grid, iter.idx, xc);
    for (int d=0; d<ndim; ++d) {
      cellLim_lo[d] = xc[d] - 0.5*grid.dx[d];
      cellLim_up[d] = xc[d] + 0.5*grid.dx[d];
      stepSz[d] = grid.dx[d]/numSteps;
      delta[d] = deltaFrac*grid.dx[d];
    }

    // need separate xrange for indexing into yshift
    // yShift:fill(yShIndexer(idx), yShItr)
    long ySh_idx = gkyl_range_idx(ySh_range, iter.idx);
    const double *ySh_arr_d = gkyl_array_cfetch(ySh_arr, ySh_idx);

    for (int dC=0; dC<2; dC++) {   // dC=0 : y=const, dC=2: x=const (boundaries).
      for (int xS=0; dC<2; dC++) { // xS=1: lower, xS=2 upper (boundary).

	/* -- Search first shifted point. Use pickLower=false in findCell unless */
	/* -- searching for points along a x=const line near the upper y-boundary. */
	/* for d = 1,2 do evPoint[d] = cellLim[d].lo+delta[d] end */
	/* evPoint[dC] = evPoint[dC]+(xS-1)*(grid:dx(dC)-2.*delta[dC]) */
	for (int d=0; d<ndim; ++d) evPoint[d] = cellLim_lo[d] + delta[d];
	evPoint[dC] = evPoint[dC] + xS*grid.dx[dC] - 2.*delta[dC];

	/* -- Evaluate yShift at this point. */
	/* local p2l_eval= Lin.Vec(1) */
	/* p2l_eval[1] = p2l(evPoint[1],grid:cellCenterInDir(1),grid:dx(1)) */
	/* yShBasis:evalBasis(p2l_eval:data(), yShBasisEv:data()) */
	/* local yShEv = 0. */
	/* for k = 1,yShNumB do yShEv = yShEv + yShItr[k]*yShBasisEv[k] end */
	double p2l_eval = p2l(evPoint[0], xc[0], grid.dx[0]);
	ySh_basis->eval(p2l_eval, yShBasisEv);
	double yShEv = 0.;
	for (int k=0; k<ySh_nc; k++) yShEv = yShEv + ySh_arr_d[k]*yShBasisEv[k];

	// Remaining for this function:
	// 1. Find index of cell that owns shifted point. use gkyl_rect_grid_find_cell()
	// 2. Search for other shifted points along the line
	// 3. Find number of donor cells for each target and store in array of length NX.
	// 4. 

      }
    }
	

  }
  
  
}

// Transform a physical coordinate (valIn) to the [-1,1] logical
//space in a cell with cell center xcIn and length dxIn.
static inline
double p2l(double valIn, double xcIn, double dxIn) return 2.*(valIn-xcIn)/dxIn;

// Wrap a number to range [lower,upper]. If pickUpper=true, output upper when
// val is a multiple of upper. Otherwise multiples of upper wrap to lower.
double wrapNum(double *val, double *lim_lo, double *lim_up, bool pickUpper) {
  double lower = lim_lo;
  double upper = lim_up;
  double L = upper - lower;
  double disp = (val - lower) % L;
  double newCoord = lower + (L + disp) % L;
  double eps = 1.e-12;
  if ( ((lower-eps < newCoord) && (newCoord < lower + eps)) ||
       ((upper-eps < newCoord) && (newCoord < upper + eps)) ) {
    if pickUpper return upper;
    else return lower;
  }
  else return newCoord;
}

/* // Given a logical space x coordinate (xi) and a (physical) y-coordinate in the target cell,  */
/* // compute the shifted y-coordinate in the logical space of the donor cell (eta \in [-1,1]).  */
/* //   xi:    logical space x coordinate. */
/* //   yTar:  physical y-coordinate in target cell.  */
/* //   pmSh:  factor multiplying the y-shift (+/- 1). */
/* //   xcDo:  cell center coordinates of donor cell. */
/* //   xcTar: cell center coordinates of target cell.  */
/* //   dx:    cell lengths. */
/* //   pickUpper: boolean indicating if wrapping function should return upper/lower boundary. */
/* double yShiftedLog(double xi, double yTar, double pmSh, double xcDo, double xcTar, bool pickUpper) { */
/*   double xPhys = xcTar[0] + 0.5*dx[0]*xi; */
/*   double yS = yTar-pmSh*yShiftF(xPhys); */
/*   yS = wrapNum(yS, domLim[1], pickUpper); */
/*   return p2l(yS, xcDo[1], dx[1]); */
//}
