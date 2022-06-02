#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_ten_moment.h>
#include <gkyl_wv_ten_moment_priv.h>

void
gkyl_ten_moment_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  
  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct wv_ten_moment *ten_moment = container_of(base->on_dev, struct wv_ten_moment, eqn);
    gkyl_cu_free(ten_moment);
  }
  
  struct wv_ten_moment *ten_moment = container_of(base, struct wv_ten_moment, eqn);
  gkyl_free(ten_moment);
}

struct gkyl_wv_eqn*
gkyl_wv_ten_moment_new(double k0)
{
  struct wv_ten_moment *ten_moment = gkyl_malloc(sizeof(struct wv_ten_moment));

  ten_moment->eqn.type = GKYL_EQN_TEN_MOMENT;
  ten_moment->eqn.num_equations = 10;
  ten_moment->eqn.num_waves = 5;
  ten_moment->k0 = k0;
  ten_moment->eqn.waves_func = wave_roe;
  ten_moment->eqn.qfluct_func = qfluct_roe;
  ten_moment->eqn.max_speed_func = max_speed;
  
  ten_moment->eqn.rotate_to_local_func = rot_to_local;
  ten_moment->eqn.rotate_to_global_func = rot_to_global;

  ten_moment->eqn.wall_bc_func = ten_moment_wall;

  ten_moment->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(ten_moment->eqn.flags);
  ten_moment->eqn.ref_count = gkyl_ref_count_init(gkyl_ten_moment_free);

  return &ten_moment->eqn;
}

double
gkyl_wv_ten_moment_k0(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_ten_moment *tm = container_of(eqn, struct wv_ten_moment, eqn);
  return tm->k0;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_wv_eqn*
gkyl_wv_ten_moment_cu_dev_new(double gas_gamma)
{
  assert(false);
  return 0;
}

#endif