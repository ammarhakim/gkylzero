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

struct gkyl_array_copy_func*
gkyl_bc_basic_create_arr_copy_func(int dir, int cdim, enum gkyl_bc_basic_type bctype,
  const struct gkyl_basis *basis)
{
#ifdef GKYL_HAVE_CUDA
//  This needs to be updated. This is from old vlasov code.
//  if (gkyl_dg_eqn_is_cu_dev(eqn)) {
//    return gkyl_vlasov_absorb_bc_create_cu(dir, basis);
//  }
#endif

  struct dg_bc_ctx *ctx = (struct dg_bc_ctx*) gkyl_malloc(sizeof(struct dg_bc_ctx));
  ctx->basis = basis;
  ctx->dir = dir;
  ctx->cdim = cdim;

  struct gkyl_array_copy_func *fout = (struct gkyl_array_copy_func*) gkyl_malloc(sizeof(struct gkyl_array_copy_func));
  switch (bctype) {
    case BC_ABSORB:
      fout->func = species_absorb_bc;
      break;

    case BC_REFLECT:
      fout->func = species_reflect_bc;
      break;

    default:
      assert(false);
      break;
  }
  fout->ctx = ctx;
  fout->ctx_on_dev = fout->ctx;

  fout->flags = 0;
  GKYL_CLEAR_CU_ALLOC(fout->flags);
  fout->on_dev = fout; // CPU function obj points to itself.
  return fout;
}

