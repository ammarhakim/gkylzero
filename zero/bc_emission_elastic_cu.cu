/* -*- c++ -*- */
extern "C" {
#include <gkyl_bc_emission_elastic.h>
#include <gkyl_bc_emission_elastic_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
}

__global__ static void
gkyl_bc_emission_elastic_create_set_cu_dev_ptrs(int dir, int cdim, const struct gkyl_basis* basis, int ncomp, struct bc_elastic_ctx *ctx, struct gkyl_array_copy_func *fout)
{
  ctx->dir = dir;
  ctx->cdim = cdim;
  ctx->basis = basis;
  ctx->ncomp = ncomp;

  fout->func = reflection;
  fout->ctx = ctx;
}

// CUDA kernel to set device pointers to function for SEY calculation.
__global__ static void
gkyl_bc_emission_elastic_set_cu_yield_func_ptrs(enum gkyl_bc_emission_elastic_type yield_type,
  struct gkyl_bc_emission_elastic_funcs *funcs, void *elastic_param_cu)
{
  switch (yield_type) {
    case GKYL_BS_FURMAN_PIVI:
      funcs->yield = furman_pivi_yield;
      funcs->elastic_param = (struct gkyl_bc_emission_elastic_furman_pivi *) elastic_param_cu;
      break;
    case GKYL_BS_CAZAUX:
      funcs->yield = cazaux_yield;
      funcs->elastic_param = (struct gkyl_bc_emission_elastic_cazaux *) elastic_param_cu;
      break;
    case GKYL_BS_CONSTANT:
      funcs->yield = constant_yield;
      funcs->elastic_param = (struct gkyl_bc_emission_elastic_constant *) elastic_param_cu;
      break; 
    default:
      assert(false);
      break;
  }
};

void
gkyl_bc_emission_elastic_choose_elastic_cu(enum gkyl_bc_emission_elastic_type elastic_type,
  struct gkyl_bc_emission_elastic_funcs *funcs, void *elastic_param)
{
  void *elastic_param_cu;
  switch (elastic_type) {
    case GKYL_BS_FURMAN_PIVI:
      elastic_param_cu = (struct gkyl_bc_emission_elastic_furman_pivi *)
        gkyl_cu_malloc(sizeof(struct gkyl_bc_emission_elastic_furman_pivi));
      gkyl_cu_memcpy(elastic_param_cu,
        (struct gkyl_bc_emission_elastic_furman_pivi *) elastic_param,
        sizeof(struct gkyl_bc_emission_elastic_furman_pivi), GKYL_CU_MEMCPY_H2D);
      break;
    case GKYL_BS_CAZAUX:
      elastic_param_cu = (struct gkyl_bc_emission_elastic_cazaux *)
        gkyl_cu_malloc(sizeof(struct gkyl_bc_emission_elastic_cazaux));
      gkyl_cu_memcpy(elastic_param_cu,
        (struct gkyl_bc_emission_elastic_cazaux *) elastic_param,
        sizeof(struct gkyl_bc_emission_elastic_cazaux), GKYL_CU_MEMCPY_H2D);
      break;
    case GKYL_BS_CONSTANT:
      elastic_param_cu = (struct gkyl_bc_emission_elastic_constant *)
        gkyl_cu_malloc(sizeof(struct gkyl_bc_emission_elastic_constant));
      gkyl_cu_memcpy(elastic_param_cu,
        (struct gkyl_bc_emission_elastic_constant *) elastic_param,
        sizeof(struct gkyl_bc_emission_elastic_constant), GKYL_CU_MEMCPY_H2D);
      break;
    default:
      assert(false);
      break;
  }
  gkyl_bc_emission_elastic_set_cu_yield_func_ptrs<<<1,1>>>(elastic_type, funcs, elastic_param_cu);
};

struct gkyl_array_copy_func*
gkyl_bc_emission_elastic_create_arr_copy_func_cu(int dir, int cdim, const struct gkyl_basis *basis,
  int ncomp)
{
  struct bc_elastic_ctx *ctx = (struct bc_elastic_ctx*) gkyl_malloc(sizeof(struct bc_elastic_ctx));
  struct gkyl_array_copy_func *fout = (struct gkyl_array_copy_func*)
    gkyl_malloc(sizeof(struct gkyl_array_copy_func));
  fout->ctx = ctx;
  
  fout->flags = 0;
  GKYL_SET_CU_ALLOC(fout->flags);
  
  struct bc_elastic_ctx *ctx_cu = (struct bc_elastic_ctx*)
    gkyl_cu_malloc(sizeof(struct bc_elastic_ctx));
  struct gkyl_array_copy_func *fout_cu = (struct gkyl_array_copy_func*)
    gkyl_cu_malloc(sizeof(struct gkyl_array_copy_func));

  gkyl_cu_memcpy(ctx_cu, ctx, sizeof(struct bc_elastic_ctx), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(fout_cu, fout, sizeof(struct gkyl_array_copy_func), GKYL_CU_MEMCPY_H2D);

  fout->ctx_on_dev = ctx_cu;

  gkyl_bc_emission_elastic_create_set_cu_dev_ptrs<<<1,1>>>(dir, cdim, basis, ncomp, ctx_cu,
    fout_cu);

  fout->on_dev = fout_cu;
  return fout;
}
