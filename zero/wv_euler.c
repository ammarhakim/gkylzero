#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_euler.h>
#include <gkyl_wv_euler_priv.h>

void
gkyl_euler_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  
  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct wv_euler *euler = container_of(base->on_dev, struct wv_euler, eqn);
    gkyl_cu_free(euler);
  }
  
  struct wv_euler *euler = container_of(base, struct wv_euler, eqn);
  gkyl_free(euler);
}

struct gkyl_wv_eqn*
gkyl_wv_euler_new(double gas_gamma)
{
  struct wv_euler *euler = gkyl_malloc(sizeof(struct wv_euler));

  euler->eqn.type = GKYL_EQN_EULER;
  euler->eqn.num_equations = 5;
  euler->eqn.num_waves = 3;
  euler->gas_gamma = gas_gamma;
  euler->eqn.waves_func = wave_roe;
  euler->eqn.qfluct_func = qfluct_roe;
  euler->eqn.max_speed_func = max_speed;

  euler->eqn.rotate_to_local_func = rot_to_local;
  euler->eqn.rotate_to_global_func = rot_to_global;

  euler->eqn.wall_bc_func = euler_wall;

  euler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(euler->eqn.flags);
  euler->eqn.ref_count = gkyl_ref_count_init(gkyl_euler_free);

  return &euler->eqn;
}

double
gkyl_wv_euler_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  return euler->gas_gamma;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_wv_eqn*
gkyl_wv_euler_cu_dev_new(double gas_gamma)
{
  assert(false);
  return 0;
}

#endif
