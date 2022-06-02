#include <math.h>
#include <float.h>
#include <string.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_wv_mhd.h>
#include <gkyl_wv_mhd_priv.h>

struct gkyl_wv_eqn*
gkyl_wv_mhd_new(double gas_gamma, const char *divergence_constraint)
{
  struct wv_mhd *mhd = gkyl_malloc(sizeof(struct wv_mhd));

  mhd->eqn.type = GKYL_EQN_MHD;
  mhd->gas_gamma = gas_gamma;
  mhd->eqn.waves_func = wave_roe;
  mhd->eqn.qfluct_func = qfluct_roe;
  mhd->eqn.max_speed_func = max_speed;
  mhd->eqn.rotate_to_local_func = rot_to_local_rect;
  mhd->eqn.rotate_to_global_func = rot_to_global_rect;

  if (strcmp(divergence_constraint, "none") == 0)
  {
    mhd->divergence_constraint = DIVB_NONE;
    mhd->eqn.num_equations = 8;
    mhd->eqn.num_waves = 7;
  }
  else if (strcmp(divergence_constraint, "eight_waves") == 0)
  {
    mhd->divergence_constraint = DIVB_EIGHT_WAVES;
    mhd->eqn.num_equations = 8;
    mhd->eqn.num_waves = 7;  // will merge the entropy wave and divB wave
  }
  else if (strcmp(divergence_constraint, "glm") == 0)
  {
    mhd->divergence_constraint = DIVB_GLM;
    mhd->eqn.num_equations = 9;
    mhd->eqn.num_waves = 9;
    mhd->eqn.rotate_to_local_func = rot_to_local_rect_glm;
    mhd->eqn.rotate_to_global_func = rot_to_global_rect_glm;
  }
  else {
    // Do not constrain divergence by default. TODO: Warn or throw an error
    mhd->divergence_constraint = DIVB_NONE;
    mhd->eqn.num_equations = 8;
    mhd->eqn.num_waves = 7;
  }

  mhd->eqn.ref_count = gkyl_ref_count_init(gkyl_wv_mhd_free);

  return &mhd->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_wv_eqn*
gkyl_wv_mhd_cu_dev_new(double gas_gamma, const char *divergence_constraint)
{
  assert(false);
  return 0;
}

#endif

double
gkyl_wv_mhd_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_mhd *mhd = container_of(eqn, struct wv_mhd, eqn);
  return mhd->gas_gamma;
}
