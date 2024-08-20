#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_gr_maxwell.h>
#include <gkyl_wv_gr_maxwell_priv.h>

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 19; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 19; i++) {
    qout[i] = win[i];
  }
}

static inline void
gr_maxwell_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  for (int i = 0; i < 19; i++) {
    diag[i] = qin[i];
  }
}

static inline void
gr_maxwell_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  for (int i = 0; i < 19; i++) {
    sout[i] = 0.0;
  }
}

void
gkyl_gr_maxwell_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_gr_maxwell *gr_maxwell = container_of(base->on_dev, struct wv_gr_maxwell, eqn);
    gkyl_cu_free(gr_maxwell);
  }

  struct wv_gr_maxwell *gr_maxwell = container_of(base, struct wv_gr_maxwell, eqn);
  gkyl_free(gr_maxwell);
}

struct gkyl_wv_eqn*
gkyl_wv_gr_maxwell_new(double light_speed, double e_fact, double b_fact, struct gkyl_gr_spacetime* spacetime, bool use_gpu)
{
  return gkyl_wv_gr_maxwell_inew(&(struct gkyl_wv_gr_maxwell_inp) {
      .light_speed = light_speed,
      .e_fact = e_fact,
      .b_fact = b_fact,
      .spacetime = spacetime,
      .rp_type = WV_GR_MAXWELL_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_gr_maxwell_inew(const struct gkyl_wv_gr_maxwell_inp* inp)
{
  struct wv_gr_maxwell *gr_maxwell = gkyl_malloc(sizeof(struct wv_gr_maxwell));

  gr_maxwell->eqn.type = GKYL_EQN_GR_MAXWELL;
  gr_maxwell->eqn.num_equations = 19;
  gr_maxwell->eqn.num_diag = 6;

  gr_maxwell->light_speed = inp->light_speed;
  gr_maxwell->e_fact = inp->e_fact;
  gr_maxwell->b_fact = inp->b_fact;
  gr_maxwell->spacetime = inp->spacetime;

  if (inp->rp_type == WV_GR_MAXWELL_RP_LAX) {
    gr_maxwell->eqn.num_waves = 2;
  }

  gr_maxwell->eqn.cons_to_riem = cons_to_riem;
  gr_maxwell->eqn.riem_to_cons = riem_to_cons;

  gr_maxwell->eqn.cons_to_diag = gr_maxwell_cons_to_diag;

  gr_maxwell->eqn.source_func = gr_maxwell_source;

  gr_maxwell->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gr_maxwell->eqn.flags);
  gr_maxwell->eqn.ref_count = gkyl_ref_count_init(gkyl_gr_maxwell_free);
  gr_maxwell->eqn.on_dev = &gr_maxwell->eqn; // On the CPU, the equation object points to itself.

  return &gr_maxwell->eqn;
}

double
gkyl_wv_gr_maxwell_light_speed(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_maxwell *gr_maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);
  double light_speed = gr_maxwell->light_speed;

  return light_speed;
}

double
gkyl_wv_gr_maxwell_e_fact(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_maxwell *gr_maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);
  double e_fact = gr_maxwell->e_fact;

  return e_fact;
}

double
gkyl_wv_gr_maxwell_b_fact(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_gr_maxwell *gr_maxwell = container_of(eqn, struct wv_gr_maxwell, eqn);
  double b_fact = gr_maxwell->b_fact;

  return b_fact;
}