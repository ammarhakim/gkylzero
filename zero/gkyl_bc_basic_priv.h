#pragma once

// Private header for bc_basic updater, not for direct use in user code.

#include <gkyl_bc_basic.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <assert.h>

// Primary struct in this updater.
struct gkyl_bc_basic {
  int dir;
  enum gkyl_edge_loc edge;
  struct gkyl_range skin_r, ghost_r;
  struct gkyl_array_copy_func *array_copy_func;
};

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set up function to apply boundary conditions.

 * @param dir Direction in which to apply BC .
 * @param cdim Number of configuration space dimensions.
 * @param bctype Type of BC .
 * @param basis Basis in which to expand coefficients in array we apply BC to.
 * @return Pointer to array_copy_func which can be passed to array_copy_fn methods.
 */
struct gkyl_array_copy_func* gkyl_bc_basic_create_arr_copy_func_cu(int dir, int cdim,
  enum gkyl_bc_basic_type bctype, const struct gkyl_basis *basis);

#endif

// context for use in BCs
struct dg_bc_ctx {
  int dir; // direction for BCs
  int cdim; // config-space dimensions
  const struct gkyl_basis *basis; // basis function
};

GKYL_CU_D
static void
species_absorb_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int nbasis = mc->basis->num_basis;
  for (int c=0; c<nbasis; ++c) out[c] = 0.0;
}

GKYL_CU_D
static void
species_reflect_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir, cdim = mc->cdim;

  mc->basis->flip_odd_sign(dir, inp, out);
  mc->basis->flip_odd_sign(dir+cdim, out, out);
}

