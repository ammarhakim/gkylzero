#pragma once

#include <gkyl_alloc.h>
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

  double psisep; // Separatrix psi for our DG representation
                 // Can differ from sibry, but we need to keep sibry
                 // because fpol, q, etc. are defined based on it

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

  int num_xpts; // Number of X-points
  double *Rxpt; // R coordinates of X points
  double *Zxpt; // Z coordinates of X-points

  bool reflect;
  bool use_gpu;
};

/**
 * Create new updater which reads in magnetic equilibrium 
 * parameters from a geqdsk file,
 * projects the poloidal flux psi(R,Z), psi/R, psi/R^2 on 
 * the RZ grid,
 * and projects F(psi)= R*B_phi on a poloidal flux grid
 *
 * @param filepath path to eqdsk file
 * @param rz_basis_type RZ basis to use for DG rep of psi, psi/R, and psi/R^2
 * @param rz_poly_order poly order for DG rep of psi, psi/R, and psi/R^2
 * @param flux_poly_order poly order to use for DG rep of F(psi)
 * @param reflect boolean indicating whether to reflect psi across R-axis
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_efit* gkyl_efit_new(const char *filepath, int rz_poly_order, enum gkyl_basis_type rz_basis_type, int flux_poly_order, bool reflect, bool use_gpu);


void gkyl_efit_release(gkyl_efit* up);
