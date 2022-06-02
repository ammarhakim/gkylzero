/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_euler.h>    
#include <gkyl_wv_euler_priv.h>
}

#include <cassert>

// CUDA kernel to set device pointers to isothermal euler kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
wv_euler_set_cu_dev_ptrs(struct wv_euler *euler)
{
  euler->eqn.waves_func = wave_roe;
  euler->eqn.qfluct_func = qfluct_roe;
  euler->eqn.max_speed_func = max_speed;
  euler->eqn.rotate_to_local_func = rot_to_local;
  euler->eqn.rotate_to_global_func = rot_to_global;

  euler->eqn.wall_bc_func = euler_wall;
}

struct gkyl_wv_eqn*
gkyl_wv_euler_cu_dev_new(double gas_gamma)
{
  struct wv_euler *euler = (struct wv_euler*) gkyl_malloc(sizeof(struct wv_euler));

  euler->eqn.type = GKYL_EQN_EULER;
  euler->eqn.num_equations = 5;
  euler->eqn.num_waves = 3;
  euler->gas_gamma = gas_gamma;

  euler->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(euler->eqn.flags);
  euler->eqn.ref_count = gkyl_ref_count_init(gkyl_euler_free);

  // copy the host struct to device struct
  struct wv_euler *euler_cu = (struct wv_euler*) gkyl_cu_malloc(sizeof(struct wv_euler));
  gkyl_cu_memcpy(euler_cu, euler, sizeof(struct wv_euler), GKYL_CU_MEMCPY_H2D);

  wv_euler_set_cu_dev_ptrs<<<1,1>>>(euler_cu);

  euler->eqn.on_dev = &euler_cu->eqn; // CPU eqn obj points to itself
  return &euler->eqn;
}
