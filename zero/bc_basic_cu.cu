/* -*- c++ -*- */

extern "C" {
#include <gkyl_bc_basic.h>
#include <gkyl_bc_basic_priv.h>
#include <gkyl_alloc.h>
// #include <gkyl_alloc_flags_priv.h>
}

__global__ static void
gkyl_bc_basic_create_set_cu_dev_ptrs(int dir, int cdim, enum gkyl_bc_basic_type bctype,
  const struct gkyl_basis* basis, struct dg_bc_ctx *ctx, struct gkyl_array_copy_func *fout)
{
  ctx->dir = dir;
  ctx->cdim = cdim;
  ctx->basis = basis;

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
}

struct gkyl_array_copy_func*
gkyl_bc_basic_create_arr_copy_func_cu(int dir, int cdim, enum gkyl_bc_basic_type bctype,
  const struct gkyl_basis *basis)
{
  // create host context and bc func structs
  struct dg_bc_ctx *ctx = (struct dg_bc_ctx*) gkyl_malloc(sizeof(struct dg_bc_ctx));
  struct gkyl_array_copy_func *fout = (struct gkyl_array_copy_func*) gkyl_malloc(sizeof(struct gkyl_array_copy_func));
  fout->ctx = ctx;

  fout->flags = 0;
  GKYL_SET_CU_ALLOC(fout->flags);

  // create device context and bc func structs
  struct dg_bc_ctx *ctx_cu = (struct dg_bc_ctx*) gkyl_cu_malloc(sizeof(struct dg_bc_ctx));
  struct gkyl_array_copy_func *fout_cu = (struct gkyl_array_copy_func*) gkyl_cu_malloc(sizeof(struct gkyl_array_copy_func));

  gkyl_cu_memcpy(ctx_cu, ctx, sizeof(struct dg_bc_ctx), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(fout_cu, fout, sizeof(struct gkyl_array_copy_func), GKYL_CU_MEMCPY_H2D);

  fout->ctx_on_dev = ctx_cu;

  gkyl_bc_basic_create_set_cu_dev_ptrs<<<1,1>>>(dir, cdim, bctype, basis, ctx_cu, fout_cu);

  // set parent on_dev pointer
  fout->on_dev = fout_cu;
  return fout;
}
