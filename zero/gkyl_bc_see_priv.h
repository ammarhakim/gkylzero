#pragma once

// Private header for bc_see updater, not for direct use in user code.

#include <gkyl_bc_see.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_updater_moment.h>
#include <math.h>
#include <assert.h>

// CHANGE NEEDED - too many arguments, fix issues with accessing ctx instead
typedef void (*see_func_t)(double *out, const double *inp, void *ctx, int idx[GKYL_MAX_DIM], double xc[GKYL_MAX_DIM]);

// CHANGE NEEDED - too many elements, fix issues with accessing ctx and store them there instead
struct gkyl_bc_see {
  int dir, cdim;
  enum gkyl_edge_loc edge;
  see_func_t func;
  struct gkyl_rect_grid *grid; // grid object
  void *ctx;
  void *ctx_on_dev;
  struct gkyl_range skin_r, ghost_r, pos_r, neg_r, cskin_r, cghost_r;
  struct gkyl_dg_updater_moment *mcalc;
  struct gkyl_array *flux;
  struct gkyl_array *f_proj;
  double tot_flux;
  const struct gkyl_basis *cbasis;
  gkyl_proj_on_basis *boundary_proj;
  bool use_gpu;
};

struct bc_see_ctx {
  int dir; // direction for BCs.
  int cdim, vdim; // config-space dimensions.
  int ncomp; // number of components within a cell.
  double *gain;
  double *elastic;
  double mass, charge, phi, k;
  double numerator, denominator;
  enum gkyl_edge_loc edge;
  const struct gkyl_basis *basis; // basis function.
};


GKYL_CU_D
static void
bc_see(double *out, const double *inp, void *ctx, int idx[GKYL_MAX_DIM], double xc[GKYL_MAX_DIM])
{
  struct bc_see_ctx *mc = (struct bc_see_ctx*) ctx;
  int cdim = mc->cdim;
  int dir = mc->dir;
  int nbasis = mc->basis->num_basis;

  // CHANGE NEEDED - elastic backscattering should be optional/separated from the inelastic
  mc->basis->flip_odd_sign(dir, inp, out);
  for (int c=0; c<nbasis; ++c) out[c] = mc->elastic[idx[1] - 1]*out[c];
  mc->basis->flip_odd_sign(dir+cdim, out, out);

  // CHANGE NEEDED - verify this is in fact correct, might be flipped in the positive/negative plane
  if ((mc->edge == GKYL_LOWER_EDGE && xc[cdim+dir] < 0) || (mc->edge == GKYL_UPPER_EDGE && xc[cdim+dir] > 0)) { 
    mc->numerator += fabs(xc[cdim+dir])*inp[0]*mc->gain[idx[1] - 1]; // CHANGE NEEDED - store coefficients in gkyl_array and index more rigorously
    mc->denominator += fabs(xc[cdim+dir])*inp[0];
  }
}
