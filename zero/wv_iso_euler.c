#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_iso_euler.h>
#include <gkyl_wv_iso_euler_priv.h>

void
gkyl_iso_euler_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  
  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct wv_iso_euler *iso_euler = container_of(base->on_dev, struct wv_iso_euler, eqn);
    gkyl_cu_free(iso_euler);
  }
  
  struct wv_iso_euler *iso_euler = container_of(base, struct wv_iso_euler, eqn);
  gkyl_free(iso_euler);
}

struct gkyl_wv_eqn*
gkyl_wv_iso_euler_new(double vt)
{
  struct wv_iso_euler *iso_euler = gkyl_malloc(sizeof(struct wv_iso_euler));

  iso_euler->eqn.type = GKYL_EQN_ISO_EULER;
  iso_euler->eqn.num_equations = 4;
  iso_euler->eqn.num_waves = 3;
  iso_euler->vt = vt;
  iso_euler->eqn.waves_func = wave_roe;
  iso_euler->eqn.qfluct_func = qfluct_roe;
  iso_euler->eqn.max_speed_func = max_speed;
  iso_euler->eqn.rotate_to_local_func = rot_to_local;
  iso_euler->eqn.rotate_to_global_func = rot_to_global;

  iso_euler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(iso_euler->eqn.flags);
  iso_euler->eqn.ref_count = gkyl_ref_count_init(gkyl_iso_euler_free);

  iso_euler->eqn.on_dev = &iso_euler->eqn; // CPU eqn obj points to itself
  return &iso_euler->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_wv_eqn*
gkyl_wv_iso_euler_cu_dev_new(double vt)
{
  assert(false);
  return 0;
}

#endif
