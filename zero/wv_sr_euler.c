#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_sr_euler.h>
#include <gkyl_wv_sr_euler_priv.h>

void
gkyl_sr_euler_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  
  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct wv_sr_euler *sr_euler = container_of(base->on_dev, struct wv_sr_euler, eqn);
    gkyl_cu_free(sr_euler);
  }
  
  struct wv_sr_euler *sr_euler = container_of(base, struct wv_sr_euler, eqn);
  gkyl_free(sr_euler);
}

struct gkyl_wv_eqn*
gkyl_wv_sr_euler_new(double gas_gamma)
{
  struct wv_sr_euler *sr_euler = gkyl_malloc(sizeof(struct wv_sr_euler));

  sr_euler->eqn.type = GKYL_EQN_SR_EULER;
  sr_euler->eqn.num_equations = 5;
  sr_euler->eqn.num_waves = 3;
  sr_euler->gas_gamma = gas_gamma;
  sr_euler->eqn.waves_func = wave_roe;
  sr_euler->eqn.qfluct_func = qfluct_roe;
  sr_euler->eqn.max_speed_func = max_speed;
  sr_euler->eqn.rotate_to_local_func = rot_to_local;
  sr_euler->eqn.rotate_to_global_func = rot_to_global;

  sr_euler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(sr_euler->eqn.flags);
  sr_euler->eqn.ref_count = gkyl_ref_count_init(gkyl_sr_euler_free);

  sr_euler->eqn.on_dev = &sr_euler->eqn; // CPU eqn obj points to itself
  return &sr_euler->eqn;
}

double
gkyl_wv_sr_euler_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_sr_euler *sr_euler = container_of(eqn, struct wv_sr_euler, eqn);
  return sr_euler->gas_gamma;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_wv_eqn*
gkyl_wv_sr_euler_cu_dev_new(double gas_gamma)
{
  assert(false);
  return 0;
}

#endif
