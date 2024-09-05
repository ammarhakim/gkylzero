#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_ultra_rel_euler.h>
#include <gkyl_wv_gr_ultra_rel_euler_priv.h>

void
gkyl_gr_ultra_rel_euler_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = container_of(base->on_dev, struct wv_gr_ultra_rel_euler, eqn);
    gkyl_cu_free(gr_ultra_rel_euler);
  }

  struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = container_of(base, struct wv_gr_ultra_rel_euler, eqn);
  gkyl_free(gr_ultra_rel_euler);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_ultra_rel_euler_new(double gas_gamma, struct gkyl_gr_spacetime* spacetime, bool use_gpu)
{
  return gkyl_wv_gr_ultra_rel_euler_inew(&(struct gkyl_wv_gr_ultra_rel_euler_inp) {
      .gas_gamma = gas_gamma,
      .spacetime = spacetime,
      .rp_type = WV_GR_ULTRA_REL_EULER_RP_ROE,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_ultra_rel_euler_inew(const struct gkyl_wv_gr_ultra_rel_euler_inp* inp)
{
  struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = gkyl_malloc(sizeof(struct wv_gr_ultra_rel_euler));

  gr_ultra_rel_euler->eqn.type = GKYL_EQN_GR_ULTRA_REL_EULER;
  gr_ultra_rel_euler->eqn.num_equations = 27;
  gr_ultra_rel_euler->eqn.num_diag = 5;

  gr_ultra_rel_euler->gas_gamma = inp->gas_gamma;
  gr_ultra_rel_euler->spacetime = inp->spacetime;

  gr_ultra_rel_euler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_ultra_rel_euler->eqn.flags);
  gr_ultra_rel_euler->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_ultra_rel_euler_free);
  gr_ultra_rel_euler->eqn.on_dev = &gr_ultra_rel_euler->eqn; // On the CPU, the equation object points to itself.

  return &gr_ultra_rel_euler->eqn;
}

double
gkyl_wv_gr_ultra_rel_euler_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = container_of(eqn, struct wv_gr_ultra_rel_euler, eqn);
  double gas_gamma = gr_ultra_rel_euler->gas_gamma;

  return gas_gamma;
}

struct gkyl_gr_spacetime*
gkyl_wv_gr_ultra_rel_euler_spacetime(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_ultra_rel_euler *gr_ultra_rel_euler = container_of(eqn, struct wv_gr_ultra_rel_euler, eqn);
  struct gkyl_gr_spacetime *spacetime = gr_ultra_rel_euler->spacetime;

  return spacetime;
}