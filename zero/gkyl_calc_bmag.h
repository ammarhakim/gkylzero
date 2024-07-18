#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_tok_geo.h>

// Object type
typedef struct gkyl_calc_bmag gkyl_calc_bmag;
typedef struct bmag_ctx bmag_ctx;

/**
 * Create new updater to compute the bmag on the compuational grid 
 *
 * @param cbasis Basis object ( computational configuration space).
 * @ param pbasis (physical RZ basis)
 * @param flux basis (Poloidal flux basis)
 * @param cgrid computational grid
 * @param pgrid physical RZ grid
 * @param fgrid poloidal flux grid from psi_min to psi_max. Usually from EFIT.
 * @param psisep poloidal flux at separatrix
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_calc_bmag* 
gkyl_calc_bmag_new(const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, const struct gkyl_basis *fbasis,
  const struct gkyl_rect_grid *cgrid, const struct gkyl_rect_grid *pgrid, const struct gkyl_rect_grid *fgrid, double psisep, bool use_gpu);


/**
 * Advance calc_bmag (compute bmag given dg fields Psi Psi/R and Psi/R^2 on the RZ grid 
 * as well as Fpol = RB_phi on the poloidal flux grid).
 *
 * @param up calc_bmag updater object.
 * @param crange, _ext computational local Config-space range and extended range.
 * @param crange_global, _ext computational global Config-space range and extended range.
 * @param prange, _ext physical RZ range and extended range.
 * @param frange, _ext poloidal flux range and extended range.
 * @param psidg, psibyrdg, psibyr2dg: DG Psi(R,Z), Psi(R,Z)/R, Psi(R,Z)/R^2 on the RZ grid
 * @param bmag_compdg output field where DG bmag on the computational basis/grid will be plaed
 * @param fpol_dg DG rep of fpol = RB_phi on the poloidal flux grid
 * @param mapc2p DG rep of mapc2p on the computational grid
 * @param calc_bphi whether or not to calculate bphi. set false for mirrors
 * @param mapc2 field containing DG rep of cylindrical coordinates
 * @param gFld output field where metric coefficients will be placed
 */

void gkyl_calc_bmag_advance(const gkyl_calc_bmag *up, const struct gkyl_range *crange, const struct gkyl_range *crange_ext,  const struct gkyl_range *crange_global, const struct gkyl_range *prange, const struct gkyl_range *prange_ext, const struct gkyl_range *frange, const struct gkyl_range* frange_ext, const struct gkyl_array *psidg, const struct gkyl_array *psibyrdg, const struct gkyl_array *psibyr2dg, struct gkyl_array* bmag_compdg, const struct gkyl_array* fpoldg, struct gkyl_array* mapc2p, bool calc_bphi);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_calc_bmag_release(gkyl_calc_bmag* up);
