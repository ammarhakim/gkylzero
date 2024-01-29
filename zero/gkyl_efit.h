#pragma once

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>

// Object type
typedef struct gkyl_efit gkyl_efit;

struct gkyl_efit{
  const char* filepath;
  int nr, nz;
  double rdim, zdim, rcentr, rleft, zmid, rmaxis, zmaxis, simag, sibry, bcentr, current, xdum;
  double rmin, rmax, zmin, zmax;

  struct gkyl_basis *rzbasis;
  struct gkyl_rect_grid *rzgrid;
  struct gkyl_range *rzlocal;
  struct gkyl_range *rzlocal_ext;
  struct gkyl_array *psizr;
  struct gkyl_array *psibyrzr;
  struct gkyl_array *psibyr2zr;

  struct gkyl_basis *fluxbasis;
  struct gkyl_rect_grid *fluxgrid;
  struct gkyl_range *fluxlocal;
  struct gkyl_range *fluxlocal_ext;
  struct gkyl_array* fpolflux;
  struct gkyl_array* qflux;
  bool use_gpu;
};

/**
 * Create new updater to project psi, psi/R, psi/R^2
 * new method fills info to contruct a grid
 *
 * @param rzbasis Basis object 
 * @param rz grid to be filled from efit
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_efit* gkyl_efit_new(const char *filepath, int rz_poly_order, enum gkyl_basis_type rz_basis_type, int flux_poly_order, bool use_gpu);



/**
 * Project psi, psi/R, psi/R^2
 *
 * @param rzbasis Basis object 
 * @param rz grid to be filled from efit
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */

void gkyl_efit_advance(gkyl_efit* up, struct gkyl_rect_grid* rzgrid, struct gkyl_rect_grid* fluxgrid, struct gkyl_range* rzlocal, struct gkyl_range* rzlocal_ext, struct gkyl_array* psizr, struct gkyl_array* psibyrzr,struct gkyl_array* psibyr2zr, struct gkyl_range* fluxlocal, struct gkyl_range* fluxlocal_ext, struct gkyl_array* fpolflux, struct gkyl_array* qflux);

void gkyl_efit_release(gkyl_efit* up);
