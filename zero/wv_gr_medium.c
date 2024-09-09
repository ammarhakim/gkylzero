#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_medium.h>
#include <gkyl_wv_gr_medium_priv.h>

void
gkyl_gr_medium_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_medium *gr_medium = container_of(base->on_dev, struct wv_gr_medium, eqn);
    gkyl_cu_free(gr_medium);
  }

  struct wv_gr_medium *gr_medium = container_of(base, struct wv_gr_medium, eqn);
  gkyl_free(gr_medium);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_medium_new(double gas_gamma, double kappa, bool use_gpu)
{
  return gkyl_wv_gr_medium_inew(&(struct gkyl_wv_gr_medium_inp) {
      .gas_gamma = gas_gamma,
      .kappa = kappa,
      .rp_type = WV_GR_MEDIUM_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_medium_inew(const struct gkyl_wv_gr_medium_inp* inp)
{
  struct wv_gr_medium *gr_medium = gkyl_malloc(sizeof(struct wv_gr_medium));

  gr_medium->eqn.type = GKYL_EQN_GR_MEDIUM;
  gr_medium->eqn.num_equations = 15;
  gr_medium->eqn.num_diag = 15;

  gr_medium->gas_gamma = inp->gas_gamma;
  gr_medium->kappa = inp->kappa;

  if (inp->rp_type == WV_GR_MEDIUM_RP_LAX) {
    gr_medium->eqn.num_waves = 2;
  }

  gr_medium->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_medium->eqn.flags);
  gr_medium->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_medium_free);
  gr_medium->eqn.on_dev = &gr_medium->eqn; // On the CPU, the equation object points to itself.

  return &gr_medium->eqn;
}

double
gkyl_wv_gr_medium_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_medium *gr_medium = container_of(eqn, struct wv_gr_medium, eqn);
  double gas_gamma = gr_medium->gas_gamma;

  return gas_gamma;
}

double
gkyl_wv_gr_medium_kappa(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_medium *gr_medium = container_of(eqn, struct wv_gr_medium, eqn);
  double kappa = gr_medium->kappa;

  return kappa;
}