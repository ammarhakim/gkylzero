#pragma once

// Private header for bc_see updater, not for direct use in user code.

#include <gkyl_bc_see.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_updater_moment.h>
#include <assert.h>

// Function pointer type for sheath reflection kernels.
typedef void (*see_func_t)(double *denominator, double *numerator, double *tot_flux, double *in_flux, int dir, double *gain, double *elastic, struct gkyl_basis *cbasis, enum gkyl_edge_loc edge, const double *inp, double *out, void *ctx, int idx[GKYL_MAX_DIM], double xc[GKYL_MAX_DIM]);

// Primary struct in this updater.
struct gkyl_bc_see {
  int dir, cdim;
  enum gkyl_edge_loc edge;
  see_func_t func;
  struct gkyl_rect_grid *grid; // grid object
  void *ctx;
  void *ctx_on_dev;
  double mass, charge, phi;
  struct gkyl_range skin_r, ghost_r, pos_r, neg_r, cskin_r, cghost_r;
  struct gkyl_dg_updater_moment *mcalc;
  struct gkyl_array *flux;
  struct gkyl_array *f_store;
  const struct gkyl_basis *cbasis;
  double *gain;
  double *elastic;
  gkyl_proj_on_basis *boundary_proj;
  bool use_gpu;
};

// context for use in BCs
struct bc_see_ctx {
  int dir; // direction for BCs.
  int cdim; // config-space dimensions.
  int ncomp; // number of components within a cell.
  double *gain;
  double *elastic;
  double mass, charge, phi, k;
  const struct gkyl_basis *basis; // basis function.
  const struct gkyl_basis *cbasis; // basis function.
};

GKYL_CU_D
static void
bc_see(double *denominator, double *numerator, double *tot_flux, double *in_flux, int dir, double *gain, double *elastic, struct gkyl_basis *cbasis, enum gkyl_edge_loc edge, const double *inp, double *out, void *ctx, int idx[GKYL_MAX_DIM], double xc[GKYL_MAX_DIM])
{
  struct bc_see_ctx *mc = (struct bc_see_ctx*) ctx;
  int cdim = mc->cdim;
  int nbasis = mc->basis->num_basis;

  mc->basis->flip_odd_sign(dir, inp, out);
  for (int c=0; c<nbasis; ++c) out[c] = elastic[idx[1] - 1]*out[c];
  mc->basis->flip_odd_sign(dir+cdim, out, out);

  if ((edge == GKYL_LOWER_EDGE && xc[cdim+dir] < 0) || (edge == GKYL_UPPER_EDGE && xc[cdim+dir] > 0)) { 
    denominator[0] += fabs(xc[cdim+dir])*inp[0]*gain[idx[1] - 1];
    numerator[0] += fabs(xc[cdim+dir])*inp[0];
  }
}
