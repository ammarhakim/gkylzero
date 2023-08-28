#pragma once

#include <gkyl_dg_diffusion_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_diffusion_priv_kerdecl.h>
#include <gkyl_dg_diffusion_priv_kerdecl_vlasov.h>

// private header for use in diffusion DG equation object creation
// functions

#define SURFKERIDX(cdim,vdim) cdim-1+vdim-1+GKYL_MIN(1,cdim-1)

// Macro for choosing surface kernel for fluid diffusion.
#define CKSURFFLUID(lst,diff_order,cdim,poly_order) lst[diff_order/2-1].list[SURFKERIDX(cdim,cdim)].kernels[poly_order-1]

// Macro for choosing surface kernel for vlasov diffusion.
#define CKSURFVLASOV(lst,diff_order,cdim,vdim,poly_order) lst[diff_order/2-1].list[SURFKERIDX(cdim,vdim)].kernels[poly_order-1]

// Macro for choosing volume and surface kernels.
#define CKVOL(lst,cdim,diff_order,poly_order,diffdir_linidx) lst[cdim-1].list[diff_order/2-1].list[poly_order-1].kernels[diffdir_linidx]
#define CKSURF(lst,diff_order,cdim,vdim,poly_order,diffid) ((diffid==GKYL_DIFFUSION_DIAGONAL_CONST_VLASOV) || (diffid==GKYL_DIFFUSION_DIAGONAL_VAR_VLASOV))? CKSURFVLASOV(lst,diff_order,cdim,vdim,poly_order) : CKSURFFLUID(lst,diff_order,cdim,poly_order)

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
    long cidx = gkyl_range_idx(&diffusion->conf_range, idxC);
    for (size_t c=0; c<diffusion->num_equations; c++) {
      int off = c*diffusion->num_basis;
      diffusion->surf[dir](xcC, dxC, _cfD(cidx), qInL+off, qInC+off, qInR+off, qRhsOut+off);
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
    long cidx = gkyl_range_idx(&diffusion->conf_range, idxSkin);
    for (size_t c=0; c<diffusion->num_equations; c++) {
      int off = c*diffusion->num_basis;
      diffusion->boundary_surf[dir](xcSkin, dxSkin, _cfD(cidx), edge, qInSkin+off, qInEdge+off, qRhsOut+off);
    }
  }
  return 0.;  // CFL frequency computed in volume term.
}

#undef _cfD
