#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>
#include <gkyl_mat.h>
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <assert.h>

// Object type
typedef struct gkyl_bc_twistshift gkyl_bc_twistshift;

/**
 * Create a new updater to apply twist-shift BCs.
 *
 * @param dir Direction in which to apply BC.
 * @param edge Lower or upper edge at which to apply BC (see gkyl_edge_loc).
 * @param local_range_ext Local extended range.
 * @param num_ghosts Number of ghosts in each dimension.
 * @param basis Basis on which coefficients in array are expanded.
 * @param grid Grid dynamic field is defined on.
 * @param cdim Number of configuration space dimensions.
 * @param yshift Discrete y-shift defined on a 1D x-grid.
 * @param ndonors Number of donors for each x-cell.
 * @param use_gpu Boolean to indicate whether to use the GPU.
 * @return New updater pointer.
 */
struct gkyl_bc_twistshift* gkyl_bc_twistshift_new(int dir, int do_dir, int shift_dir, enum gkyl_edge_loc edge,
  const struct gkyl_range *local_range_ext, const int *num_ghosts, const struct gkyl_basis *basis,
  const struct gkyl_rect_grid *grid, int cdim,
  const struct gkyl_array *yshift, const int *ndonors, const int *cells_do, bool use_gpu);

/**
 * Populate a matrix in mats corresponding to the x-cell with cellidx and
 * the donor cell corresponding to doidx, for the sub-cell integral that
 * has variably x limits represented by a DG polynomial, and a y-integral that
 * goes from yLimLo to yLimUp.
 *
 * @param up BC updater.
 * @param sFac +/-1 factor to add or subtract this subcell integral.
 * @param xLimLo DG representation of the lower x-limit.
 * @param xLimUp DG representation of the upper x-limit.
 * @param yLimLo lower y-limit.
 * @param yLimUp upper y-limit.
 * @param dyDo Cell length along y.
 * @param yOff Offset along y.
 * @param ySh DG representation of the y shift.
 * @param mats Matrices.
 * @param cellidx Cell index along x.
 * @param doidx Donor index.
 */
void gkyl_bc_twistshift_integral_xlimdg(struct gkyl_bc_twistshift *up,
  double sFac, const double *xLimLo, const double *xLimUp, double yLimLo, double yLimUp,
  double dyDo, double yOff, const double *ySh, struct gkyl_nmat *mats, int cellidx, int doidx);

/**
 * Populate a matrix in mats corresponding to the x-cell with cellidx and
 * the donor cell corresponding to doidx, for the sub-cell integral that
 * has variable y limits represented by a DG polynomial, and a x-integral that
 * goes from xLimLo to xLimUp.
 *
 * @param up BC updater.
 * @param sFac +/-1 factor to add or subtract this subcell integral.
 * @param xLimLo lower x-limit.
 * @param xLimUp upper x-limit.
 * @param yLimLo DG representation of the lower y-limit.
 * @param yLimUp DG representation of the upper y-limit.
 * @param dyDo Cell length along y.
 * @param yOff Offset along y.
 * @param ySh DG representation of the y shift.
 * @param mats Matrices.
 * @param cellidx Cell index along x.
 * @param doidx Donor index.
 */
void gkyl_bc_twistshift_integral_ylimdg(struct gkyl_bc_twistshift *up,
  double sFac, double xLimLo, double xLimUp, const double *yLimLo, const double *yLimUp,
  double dyDo, double yOff, const double *ySh, struct gkyl_nmat *mats, int cellidx, int doidx);

/**
 * Populate a matrix in mats corresponding to the x-cell with cellidx and
 * the donor cell corresponding to doidx, for the full-cell integral.
 *
 * @param up BC updater.
 * @param sFac +/-1 factor to add or subtract this subcell integral.
 * @param dyDo Cell length along y.
 * @param yOff Offset along y.
 * @param ySh DG representation of the y shift.
 * @param mats Matrices.
 * @param cellidx Cell index along x.
 * @param doidx Donor index.
 */
void gkyl_bc_twistshift_integral_fullcelllimdg(struct gkyl_bc_twistshift *up,
  double dyDo, double yOff, const double *ySh, struct gkyl_nmat *mats, int cellidx, int doidx);

/**
 * Multiply the donor matrices by the donor cell dg coefficients
 *
 * @param matsdo nmat of donor matrices for all x grid cells and all donor cells
 * @param vecsdo nmat of donor vectors for all x grid cells and all donor cells
 * @param vecstar nmat of target vectors for all x grid cells
 */
void gkyl_bc_twistshift_mv(struct gkyl_bc_twistshift *up, struct gkyl_nmat *matsdo, struct gkyl_array *fdo, struct gkyl_array *ftar);

/**
 * Free memory associated with bc_twistshift updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_twistshift_release(struct gkyl_bc_twistshift *up);
