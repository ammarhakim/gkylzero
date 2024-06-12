/* -*- c++ -*- */
extern "C" {
#include <gkyl_bc_emission_elastic.h>
#include <gkyl_bc_emission_elastic_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
}

// start ID for use in various loops
#define START_ID (threadIdx.x + blockIdx.x*blockDim.x)

// CUDA kernel to set device pointers to function for SEY calculation.
__global__ static void
gkyl_bc_emission_elastic_set_cu_yield_func_ptrs(enum gkyl_bc_emission_elastic_type yield_type,
  struct gkyl_bc_emission_elastic_funcs *funcs)
{
  switch (yield_type) {
    case GKYL_BS_FURMAN_PIVI:
      funcs->yield = furman_pivi_yield;
      break;
    case GKYL_BS_CAZAUX:
      funcs->yield = cazaux_yield;
      break;
    case GKYL_BS_CONSTANT:
      funcs->yield = constant_yield;
      break; 
    default:
      assert(false);
      break;
  }
};

void
gkyl_bc_emission_elastic_choose_func_cu(enum gkyl_bc_emission_elastic_type yield_type,
  struct gkyl_bc_emission_elastic_funcs *funcs)
{
  gkyl_bc_emission_elastic_set_cu_yield_func_ptrs<<<1,1>>>(yield_type, funcs);
}

void
gkyl_bc_emission_elastic_advance_cu(const struct gkyl_bc_emission_elastic *up,
  struct gkyl_range *emit_skin_r, struct gkyl_array *buff_arr, struct gkyl_array *f_skin,
  struct gkyl_array *f_emit, struct gkyl_array *elastic_yield, struct gkyl_basis *basis)
{
  int nblocks = emit_skin_r->nblocks, nthreads = emit_skin_r->nthreads;

  // gkyl_array_clear_range(weight, 0.0, conf_r);

  // // Calculate weighted mean numerator and denominator
  // gkyl_bc_emission_spectrum_advance_cu_weight_ker<<<nblocks, nthreads>>>(up->cdim, up->dir, up->edge, f_skin->on_dev, weight->on_dev, *grid, gamma->on_dev, *skin_r, *ghost_r, *conf_r, up->funcs_cu, up->bc_param_cu);

  // nblocks = buff_r->nblocks;
  // nthreads = buff_r->nthreads;

  // // Finish weighted mean calculation and accumulate to buffer
  // gkyl_bc_emission_spectrum_advance_cu_accumulate_ker<<<nblocks, nthreads>>>(f_proj->on_dev, f_buff->on_dev, weight->on_dev, k->on_dev, flux->on_dev, *buff_r, *conf_r, up->funcs_cu, up->bc_param_cu);
}
