#include <gkyl_alloc.h>
#include <gkyl_wv_euler.h>

struct wv_euler {
    struct gkyl_wv_eqn eqn; // base object
    double gas_gamma; // gas adiabatic constant
};

static void
euler_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_wv_eqn *base = container_of(ref, struct gkyl_wv_eqn, ref_count);
  struct wv_euler *euler = container_of(base, struct wv_euler, eqn);
  gkyl_free(euler);
}

// Waves and speeds using Roe averaging
static double
wave_roe(const struct gkyl_wv_eqn *eqn, 
  int dir, const double *ql, const double *qr, double *waves, double *speeds)
{
  const struct wv_euler *euler = container_of(eqn, struct wv_euler, eqn);
  double gas_gamma = euler->gas_gamma;
  

  return 0;
}

struct gkyl_wv_eqn*
wv_euler_new(double gas_gamma)
{
  struct wv_euler *euler = gkyl_malloc(sizeof(struct wv_euler));

  euler->eqn.num_equations = 5;
  euler->eqn.num_waves = 3;
  euler->gas_gamma = gas_gamma;
  euler->eqn.wave_func = wave_roe;

  euler->eqn.ref_count = (struct gkyl_ref_count) { euler_free, 1 };

  return &euler->eqn;
}
