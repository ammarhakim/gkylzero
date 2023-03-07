#pragma once

// Private header for bc_emission updater, not for direct use in user code.

#include <gkyl_bc_emission.h>
#include <gkyl_array_ops.h>
#include <assert.h>

// Function pointer type for sheath reflection kernels.
typedef void (*emission_func_t)(double *out, const double *inp, void *ctx, int in_idx[GKYL_MAX_DIM], int out_idx[GKYL_MAX_DIM]);

// Primary struct in this updater.
struct gkyl_bc_emission {
  int dir, cdim;
  enum gkyl_edge_loc edge;
  emission_func_t func;
  void *ctx;
  void *ctx_on_dev;
  enum gkyl_bc_emission_type bctype;
  struct gkyl_range skin_r, ghost_r;
  int in_idx[GKYL_MAX_DIM], out_idx[GKYL_MAX_DIM];
  bool use_gpu;
};

// context for use in BCs
struct bc_gain_ctx {
  int dir; // direction for BCs.
  int cdim; // config-space dimensions.
  int ncomp; // number of components within a cell.
  double gain;
  const struct gkyl_basis *basis; // basis function.
};

GKYL_CU_D
static void
bc_emission_gain(double *out, const double *inp, void *ctx, int in_idx[GKYL_MAX_DIM], int out_idx[GKYL_MAX_DIM])
{
  struct bc_gain_ctx *mc = (struct bc_gain_ctx*) ctx;
  int dir = mc->dir, cdim = mc->cdim;
  int nbasis = mc->basis->num_basis;

  mc->basis->flip_odd_sign(dir, inp, out);
  for (int c=0; c<nbasis; ++c) out[c] = mc->gain*out[c];
  mc->basis->flip_odd_sign(dir+cdim, out, out);
}
