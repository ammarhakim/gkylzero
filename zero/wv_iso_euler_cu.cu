/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_iso_euler.h>    
#include <gkyl_wv_iso_euler_priv.h>
}

#include <cassert>

// CUDA kernel to set device pointers to isothermal euler kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
wv_iso_euler_set_cu_dev_ptrs(struct wv_iso_euler *iso_euler)
{
  iso_euler->eqn.waves_func = wave_roe;
  iso_euler->eqn.qfluct_func = qfluct_roe;
  iso_euler->eqn.max_speed_func = max_speed;
  iso_euler->eqn.rotate_to_local_func = rot_to_local;
  iso_euler->eqn.rotate_to_global_func = rot_to_global;
}

struct gkyl_wv_eqn*
gkyl_wv_iso_euler_cu_dev_new(double vt)
{
  struct wv_iso_euler *iso_euler = (struct wv_iso_euler*) gkyl_malloc(sizeof(struct wv_iso_euler));

  iso_euler->eqn.type = GKYL_EQN_ISO_EULER;
  iso_euler->eqn.num_equations = 4;
  iso_euler->eqn.num_waves = 3;
  iso_euler->vt = vt;

  iso_euler->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(iso_euler->eqn.flags);
  iso_euler->eqn.ref_count = gkyl_ref_count_init(gkyl_iso_euler_free);

  // copy the host struct to device struct
  struct wv_iso_euler *iso_euler_cu = (struct wv_iso_euler*) gkyl_cu_malloc(sizeof(struct wv_iso_euler));
  gkyl_cu_memcpy(iso_euler_cu, iso_euler, sizeof(struct wv_iso_euler), GKYL_CU_MEMCPY_H2D);

  wv_iso_euler_set_cu_dev_ptrs<<<1,1>>>(iso_euler_cu);

  iso_euler->eqn.on_dev = &iso_euler_cu->eqn; // CPU eqn obj points to itself
  return &iso_euler->eqn;
}
