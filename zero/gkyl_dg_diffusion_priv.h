#pragma once

#include <gkyl_dg_diffusion_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_diffusion_priv_kerdecl.h>
#include <gkyl_dg_diffusion_priv_kerdecl_vlasov.h>
#include <gkyl_dg_diffusion_priv_kerdecl_gyrokinetic.h>

static inline int numeq_from_diffid(enum gkyl_diffusion_id diffusion_id) {
  // Number of equations based on diffusion id.
  int num_equations = 1;
  if ((diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_EULER) || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_VAR_EULER))
    num_equations = 4;
  else if ((diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_EULER_ISO) || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_VAR_EULER_ISO))
    num_equations = 5;
  if ((diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_PKPM) || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_VAR_PKPM))
    num_equations = 3;
  return num_equations;
}

static inline bool constcoeff_from_diffid(enum gkyl_diffusion_id diffusion_id) {
  // Determine of diffusion coefficient is constant.
  bool const_coeff = (   (diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST)
                      || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_EULER)
                      || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_EULER_ISO)
                      || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_VLASOV)
                      || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_GYROKINETIC)
                      || (diffusion_id == GKYL_DIFFUSION_DIAGONAL_CONST_PKPM) );
  return const_coeff;
}

static inline int diffdirs_linidx(const bool *isdirdiff, int cdim) {
  // Compute the linear index into the array of volume kernels (one
  // kernel for each combination of diffusive directions).
  bool diff_in_dir[GKYL_MAX_CDIM];
  if (isdirdiff)
    for (size_t d=0; d<cdim; d++) diff_in_dir[d] = isdirdiff[d];
  else
    for (size_t d=0; d<cdim; d++) diff_in_dir[d] = true;

  // Linear index into list of volume kernels.
  int dirs_bin_key[] = {1,2,4,8,16,32}; // Binary: 000001, 000010, 000100, 001000, 010000, 100000.
  int dirs_linidx = 0; // Binary 000000.
  for (int d=0; d<cdim; d++) {
     if (diff_in_dir[d]) dirs_linidx = dirs_linidx | dirs_bin_key[d];
  }
  dirs_linidx -= 1;
  return dirs_linidx;
}

// private header for use in diffusion DG equation object creation
// functions

#define SURFKERIDX(cdim,vdim) cdim-1+vdim-1+GKYL_MIN(1,cdim-1)
#define SURFKERIDXGK(cdim,vdim) (cdim-1+vdim-1)*2-(vdim-1)

// Macro for choosing surface kernel for fluid diffusion.
#define CKSURFCONF(lst,diff_order,cdim,poly_order) lst[diff_order/2-1].list[SURFKERIDX(cdim,cdim)].kernels[poly_order-1]

// Macro for choosing surface kernel for vlasov diffusion.
#define CKSURFVLASOV(lst,diff_order,cdim,vdim,poly_order) lst[diff_order/2-1].list[SURFKERIDX(cdim,vdim)].kernels[poly_order-1]

// Macro for choosing surface kernel for gyrokinetic diffusion.
#define CKSURFGYROKINETIC(lst,diff_order,cdim,vdim,poly_order) lst[diff_order/2-1].list[SURFKERIDXGK(cdim,vdim)].kernels[poly_order-1]

// Macro for choosing volume and surface kernels.
#define CKVOL(lst,cdim,diff_order,poly_order,diffdir_linidx) lst[cdim-1].list[diff_order/2-1].list[poly_order-1].kernels[diffdir_linidx]
#define CKSURF(lst,diff_order,cdim,vdim,poly_order,diffid) ((diffid==GKYL_DIFFUSION_DIAGONAL_CONST_VLASOV) || (diffid==GKYL_DIFFUSION_DIAGONAL_VAR_VLASOV))? CKSURFVLASOV(lst,diff_order,cdim,vdim,poly_order) : (((diffid==GKYL_DIFFUSION_DIAGONAL_CONST_GYROKINETIC) || (diffid==GKYL_DIFFUSION_DIAGONAL_VAR_GYROKINETIC))? CKSURFGYROKINETIC(lst,diff_order,cdim,vdim,poly_order) : CKSURFCONF(lst,diff_order,cdim,poly_order))

/**
 * Free diffusion equation object
 *
 * @param ref Reference counter for constant diffusion equation
 */
void gkyl_diffusion_free(const struct gkyl_ref_count* ref);

GKYL_CU_D static double surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR,
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  if (diffusion->diff_in_dir[dir]) {
    for (size_t c=0; c<diffusion->num_equations; c++) {
      int off = c*diffusion->num_basis;
      diffusion->surf[dir](xcC, dxC, _cfD(idxC), qInL+off, qInC+off, qInR+off, qRhsOut+off);
    }
  }
  return 0.;  // CFL frequency computed in volume term.
}

GKYL_CU_D static double boundary_surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcEdge, const double* xcSkin, const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{ 
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  if (diffusion->diff_in_dir[dir]) {
    for (size_t c=0; c<diffusion->num_equations; c++) {
      int off = c*diffusion->num_basis;
      diffusion->boundary_surf[dir](xcSkin, dxSkin, _cfD(idxSkin), edge, qInSkin+off, qInEdge+off, qRhsOut+off);
    }
  }
  return 0.;  // CFL frequency computed in volume term.
}

#undef _cfD

#ifdef GKYL_HAVE_CUDA
/**
 * Create a new diffusion equation object on the device.
 *
 * @param basis Basis functions of the equation system.
 * @param cbasis Configuration space basis.
 * @param diffusion_id Diffusion type (constant/varying & fluid/vlasov/etc).
 * @param diff_in_dir Whether to apply diffusion in each direction.
 * @param diff_order Diffusion order.
 * @param conf_range Configuration space range (to index diff coefficient).
 * @param use_gpu Whether to run on host or device.
 * @return Pointer to diffusion equation object
 */
struct gkyl_dg_eqn*
gkyl_dg_diffusion_cu_dev_new(const struct gkyl_basis *basis, const struct gkyl_basis *cbasis,
  enum gkyl_diffusion_id diffusion_id, const bool *diff_in_dir, int diff_order,
  const struct gkyl_range *conf_range);
#endif
