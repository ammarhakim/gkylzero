#pragma once

#include <gkyl_range.h>
#include <gkyl_basis.h>
#include <gkyl_array.h>
#include <gkyl_rect_grid.h>
#include <gkyl_evalf_def.h>
#include <assert.h>

// Object type
typedef struct gkyl_bc_twistshift gkyl_bc_twistshift;

struct gkyl_bc_twistshift_inp {
  int bc_dir; // Direction in which to apply this BC.
  int shift_dir; // Direction of the shift.
  int shear_dir; // Direction in which the shift varies (shear).
  enum gkyl_edge_loc edge; // Edge of to apply this BC at (lower/upper).
  int cdim; // Configuration space dimensions.
  struct gkyl_range local_ext_r; // Local extended range.
  const int *num_ghost; // Number of ghost cells in each direction.
  struct gkyl_basis basis; // Basis of the field shifted.
  struct gkyl_rect_grid grid; // Grid the field shifted is defined on.
  evalf_t shift_func; // Function defining the shift.
  void *shift_func_ctx; // Context for shift_func.
  bool use_gpu; // Whether to apply the BC using the GPU.
  // Optional inputs:
  int shift_poly_order; // Basis order for the DG representation of the shift.
};

/**
 * Create a new updater to apply twist-shift BCs.
 *
 * @param inp bc_twistshift_inp struct containing the inputs to the updater.
 * @return New updater pointer.
 */
struct gkyl_bc_twistshift* gkyl_bc_twistshift_new(const struct gkyl_bc_twistshift_inp *inp);
 
/**
 * Apply the twist-shift. It assumes that periodicity along bc_dir has been
 * applied to the donor field. Can be used in-place.
 *
 * @param up Twist-shift BC updater object.
 * @param fdo Donor field.
 * @param ftar Target field.
 */
void gkyl_bc_twistshift_advance(struct gkyl_bc_twistshift *up, struct gkyl_array *fdo, struct gkyl_array *ftar);

/**
 * Free memory associated with bc_twistshift updater.
 *
 * @param up BC updater.
 */
void gkyl_bc_twistshift_release(struct gkyl_bc_twistshift *up);

// /**
//  * Populate a matrix in mats corresponding to the x-cell with cellidx and
//  * the donor cell corresponding to doidx, for the sub-cell integral that
//  * has variably x limits represented by a DG polynomial, and a y-integral that
//  * goes from yLimLo to yLimUp.
//  *
//  * @param up BC updater.
//  * @param sFac +/-1 factor to add or subtract this subcell integral.
//  * @param xLimLo DG representation of the lower x-limit.
//  * @param xLimUp DG representation of the upper x-limit.
//  * @param yLimLo lower y-limit.
//  * @param yLimUp upper y-limit.
//  * @param dyDo Cell length along y.
//  * @param yOff Offset along y.
//  * @param ySh DG representation of the y shift.
//  * @param mats Matrices.
//  * @param cellidx Cell index along x.
//  * @param doidx Donor index.
//  */
// void gkyl_bc_twistshift_integral_xlimdg(struct gkyl_bc_twistshift *up,
//   double sFac, const double *xLimLo, const double *xLimUp, double yLimLo, double yLimUp,
//   double dyDo, double yOff, const double *ySh, int cellidx, int doidx);
// 
// /**
//  * Populate a matrix in mats corresponding to the x-cell with cellidx and
//  * the donor cell corresponding to doidx, for the sub-cell integral that
//  * has variable y limits represented by a DG polynomial, and a x-integral that
//  * goes from xLimLo to xLimUp.
//  *
//  * @param up BC updater.
//  * @param sFac +/-1 factor to add or subtract this subcell integral.
//  * @param xLimLo lower x-limit.
//  * @param xLimUp upper x-limit.
//  * @param yLimLo DG representation of the lower y-limit.
//  * @param yLimUp DG representation of the upper y-limit.
//  * @param dyDo Cell length along y.
//  * @param yOff Offset along y.
//  * @param ySh DG representation of the y shift.
//  * @param mats Matrices.
//  * @param cellidx Cell index along x.
//  * @param doidx Donor index.
//  */
// void gkyl_bc_twistshift_integral_ylimdg(struct gkyl_bc_twistshift *up,
//   double sFac, double xLimLo, double xLimUp, const double *yLimLo, const double *yLimUp,
//   double dyDo, double yOff, const double *ySh, int cellidx, int doidx);
// 
// /**
//  * Populate a matrix in mats corresponding to the x-cell with cellidx and
//  * the donor cell corresponding to doidx, for the full-cell integral.
//  *
//  * @param up BC updater.
//  * @param sFac +/-1 factor to add or subtract this subcell integral.
//  * @param dyDo Cell length along y.
//  * @param yOff Offset along y.
//  * @param ySh DG representation of the y shift.
//  * @param mats Matrices.
//  * @param cellidx Cell index along x.
//  * @param doidx Donor index.
//  */
// void gkyl_bc_twistshift_integral_fullcelllimdg(struct gkyl_bc_twistshift *up,
//   double dyDo, double yOff, const double *ySh, int cellidx, int doidx);
// 
// /**
//  * Copy donor matrices to device if necessary
//  */
// void gkyl_bc_twistshift_copy_matsdo(struct gkyl_bc_twistshift *up);
// 
// 
// 
// /**
//  * Increment target field by target vectors. Parallelized over x and dg coeff
//  * @param ftar target field at specific location
//  * @param tar_locs target locations
//  * @param num_tar_locs number of target locations
//  * @param vecstar vectors to be accumulated into ftar
//  * @param ndonors_cum cumulative list of number of donors
//  * @param
//  */
// 
// void gkyl_bc_twistshift_inc_cu(const struct gkyl_array* ftar, long* tar_locs, int num_tar_locs, struct gkyl_nmat* vecstar, int* ndonors_cum);
// 
// /**
//  * Zero out target field
//  * @param ftar target array
//  * @param loc locations coefficient to zero
//  * @param num_locs number of locations coefficient to zero
//  */
// void gkyl_bc_twistshift_clear_cu(struct gkyl_array* ftar, long* locs, int num_locs);
// 
// /**
//  * Fill donor vecs from donor field
//  * @param fdo donor field
//  * @param locs donor locations
//  * @param vecsdo donor vectors
//  */
// void gkyl_bc_twistshift_set_vecsdo_cu(struct gkyl_array* fdo, long* locs, struct gkyl_nmat* vecsdo);
