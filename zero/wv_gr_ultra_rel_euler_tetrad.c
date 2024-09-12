#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_ultra_rel_euler_tetrad.h>
#include <gkyl_wv_gr_ultra_rel_euler_tetrad_priv.h>

static inline void
gr_ultra_rel_euler_tetrad_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 27; i++) {
    diag[i] = qin[i];
  }
}

static inline void
gr_ultra_rel_euler_tetrad_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  // TODO: Rewrite source solver(s).
  for (int i = 0; i < 27; i++) {
    sout[i] = 0.0;
  }
}

void
gkyl_gr_ultra_rel_euler_tetrad_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_ultra_rel_euler_tetrad *gr_ultra_rel_euler_tetrad = container_of(base->on_dev, struct wv_gr_ultra_rel_euler_tetrad, eqn);
    gkyl_cu_free(gr_ultra_rel_euler_tetrad);
  }

  struct wv_gr_ultra_rel_euler_tetrad *gr_ultra_rel_euler_tetrad = container_of(base, struct wv_gr_ultra_rel_euler_tetrad, eqn);
  gkyl_free(gr_ultra_rel_euler_tetrad);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_ultra_rel_euler_tetrad_new(double gas_gamma, struct gkyl_gr_spacetime* spacetime, bool use_gpu)
{
  return gkyl_wv_gr_ultra_rel_euler_tetrad_inew(&(struct gkyl_wv_gr_ultra_rel_euler_tetrad_inp) {
      .gas_gamma = gas_gamma,
      .spacetime = spacetime,
      .rp_type = WV_GR_ULTRA_REL_EULER_TETRAD_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_ultra_rel_euler_tetrad_inew(const struct gkyl_wv_gr_ultra_rel_euler_tetrad_inp* inp)
{
  struct wv_gr_ultra_rel_euler_tetrad *gr_ultra_rel_euler_tetrad = gkyl_malloc(sizeof(struct wv_gr_ultra_rel_euler_tetrad));

  gr_ultra_rel_euler_tetrad->eqn.type = GKYL_EQN_GR_ULTRA_REL_EULER_TETRAD;
  gr_ultra_rel_euler_tetrad->eqn.num_equations = 27;
  gr_ultra_rel_euler_tetrad->eqn.num_diag = 5;

  gr_ultra_rel_euler_tetrad->gas_gamma = inp->gas_gamma;
  gr_ultra_rel_euler_tetrad->spacetime = inp->spacetime;

  gr_ultra_rel_euler_tetrad->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_ultra_rel_euler_tetrad->eqn.flags);
  gr_ultra_rel_euler_tetrad->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_ultra_rel_euler_tetrad_free);
  gr_ultra_rel_euler_tetrad->eqn.on_dev = &gr_ultra_rel_euler_tetrad->eqn; // On the CPU, the equation object points to itself.

  return &gr_ultra_rel_euler_tetrad->eqn;
}

double
gkyl_wv_gr_ultra_rel_euler_tetrad_gas_gamma(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_ultra_rel_euler_tetrad *gr_ultra_rel_euler_tetrad = container_of(eqn, struct wv_gr_ultra_rel_euler_tetrad, eqn);
  double gas_gamma = gr_ultra_rel_euler_tetrad->gas_gamma;

  return gas_gamma;
}

struct gkyl_gr_spacetime*
gkyl_wv_gr_ultra_rel_euler_tetrad_spacetime(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_ultra_rel_euler_tetrad *gr_ultra_rel_euler_tetrad = container_of(eqn, struct wv_gr_ultra_rel_euler_tetrad, eqn);
  struct gkyl_gr_spacetime *spacetime = gr_ultra_rel_euler_tetrad->spacetime;

  return spacetime;
}