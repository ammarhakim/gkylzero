/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_sr_euler.h>    
#include <gkyl_wv_sr_euler_priv.h>
}

#include <cassert>

// CUDA kernel to set device pointers to special relativistic euler kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
wv_sr_euler_set_cu_dev_ptrs(struct wv_sr_euler *sr_euler)
{
  sr_euler->eqn.waves_func = wave_roe;
  sr_euler->eqn.qfluct_func = qfluct_roe;
  sr_euler->eqn.max_speed_func = max_speed;
  sr_euler->eqn.rotate_to_local_func = rot_to_local;
  sr_euler->eqn.rotate_to_global_func = rot_to_global;
}

struct gkyl_wv_eqn*
gkyl_wv_sr_euler_cu_dev_new(double gas_gamma)
{
  struct wv_sr_euler *sr_euler = (struct wv_sr_euler*) gkyl_malloc(sizeof(struct wv_sr_euler));

  sr_euler->eqn.type = GKYL_EQN_SR_EULER;
  sr_euler->eqn.num_equations = 5;
  sr_euler->eqn.num_waves = 3;
  sr_euler->gas_gamma = gas_gamma;

  sr_euler->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(sr_euler->eqn.flags);
  sr_euler->eqn.ref_count = gkyl_ref_count_init(gkyl_sr_euler_free);

  // copy the host struct to device struct
  struct wv_sr_euler *sr_euler_cu = (struct wv_sr_euler*) gkyl_cu_malloc(sizeof(struct wv_sr_euler));
  gkyl_cu_memcpy(sr_euler_cu, sr_euler, sizeof(struct wv_sr_euler), GKYL_CU_MEMCPY_H2D);

  wv_sr_euler_set_cu_dev_ptrs<<<1,1>>>(sr_euler_cu);

  sr_euler->eqn.on_dev = &sr_euler_cu->eqn; // CPU eqn obj points to itself
  return &sr_euler->eqn;
}
