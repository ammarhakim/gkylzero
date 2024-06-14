/* -*- c++ -*- */
extern "C" {
#include <gkyl_bc_emission_elastic.h>
#include <gkyl_bc_emission_elastic_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
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
