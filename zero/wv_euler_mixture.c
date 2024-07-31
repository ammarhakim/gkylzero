#include <assert.h>
#include <math.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_wv_euler_mixture.h>
#include <gkyl_wv_euler_mixture_priv.h>

void
gkyl_euler_mixture_prim_vars(int num_species, double* gas_gamma_s, const double* q, double* v)
{
  double *rho_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_s[i] = q[i];
  }
  double momx = q[num_species];
  double momy = q[num_species + 1];
  double momz = q[num_species + 2];
  double Etot = q[num_species + 3];

  double rho_total = 0.0;
  for (int i = 0; i < num_species; i++) {
    rho_total += rho_s[i];
  }

  double vx = momx / rho_total;
  double vy = momy / rho_total;
  double vz = momz / rho_total;

  double *p_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    p_s[i] = (gas_gamma_s[i] - 1.0) * ((rho_s[i] * (Etot / rho_total)) - (0.5 * rho_s[i] * ((vx * vx) + (vy * vy) + (vz * vz))));
  }

  double p_total = 0.0;
  for (int i = 0; i < num_species; i++) {
    p_total += p_s[i];
  }

  for (int i = 0; i < num_species; i++) {
    v[i] = rho_s[i];
  }
  v[num_species] = vx;
  v[num_species + 1] = vy;
  v[num_species + 2] = vz;
  v[num_species + 3] = p_total;
}

static inline double
gkyl_euler_mixture_max_abs_speed(int num_species, double* gas_gamma_s, const double* q)
{
  double *v = gkyl_malloc(sizeof(double[4 + num_species]));
  gkyl_euler_mixture_prim_vars(num_species, gas_gamma_s, q, v);

  double rho_total = 0.0;
  for (int i = 0; i < num_species; i++) {
    rho_total += v[i];
  }

  double vx = v[num_species];
  double vy = v[num_species + 1];
  double vz = v[num_species + 2];
  double p_total = v[num_species + 3];

  double v_mag = sqrt((vx * vx) + (vy * vy) + (vz * vz));

  double max_abs_speed = 0.0;
  for (int i = 0; i < num_species; i++) {
    if (fabs(v_mag) + sqrt(gas_gamma_s[i] * (p_total / rho_total)) > max_abs_speed) {
      max_abs_speed = fabs(v_mag) + sqrt(gas_gamma_s[i] * (p_total / rho_total));
    }
  }

  return max_abs_speed;
}

void
gkyl_euler_mixture_flux(int num_species, double* gas_gamma_s, const double* q, double* flux)
{
  double *v = gkyl_malloc(sizeof(double[4 + num_species]));
  gkyl_euler_mixture_prim_vars(num_species, gas_gamma_s, q, v);

  double *rho_s = gkyl_malloc(sizeof(double[num_species]));
  for (int i = 0; i < num_species; i++) {
    rho_s[i] = v[i];
  }

  double rho_total = 0.0;
  for (int i = 0; i < num_species; i++) {
    rho_total += rho_s[i];
  }

  double vx = v[num_species];
  double vy = v[num_species + 1];
  double vz = v[num_species + 2];
  double p_total = v[num_species + 3];

  double Etot = q[num_species + 3];

  for (int i = 0; i < num_species; i++) {
    flux[i] = rho_s[i] * vx;
  }
  flux[num_species] = (rho_total * (vx * vx)) + p_total;
  flux[num_species + 1] = rho_total * (vx * vy);
  flux[num_species + 2] = rho_total * (vx * vz);
  flux[num_species + 3] = (Etot + p_total) * vx;
}

static inline void
cons_to_riem(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* qin, double* wout)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 4 + num_species; i++) {
    wout[i] = qin[i];
  }
}

static inline void
riem_to_cons(const struct gkyl_wv_eqn* eqn, const double* qstate, const double* win, double* qout)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  // TODO: This should use a proper L matrix.
  for (int i = 0; i < 4 + num_species; i++) {
    qout[i] = win[i];
  }
}

static void
euler_mixture_wall(double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  long num_species = 2;

  for (int i = 0; i < 4 + num_species; i++) {
    ghost[i] = skin[i];
  }

  ghost[num_species] = -ghost[num_species];
}

static void
euler_mixture_no_slip(double t, int nc, const double* skin, double* GKYL_RESTRICT ghost, void* ctx)
{
  long num_species = 2;

  for (int i = 0; i < num_species; i++) {
    ghost[i] = skin[i];
  }

  for (int i = num_species; i < 3 + num_species; i++) {
    ghost[i] = -skin[i];
  }

  ghost[3 + num_species] = skin[3 + num_species];
}

static inline void
rot_to_local(const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qglobal, double* GKYL_RESTRICT qlocal)
{
  long num_species = 2;

  for (int i = 0; i < num_species; i++) {
    qlocal[i] = qglobal[i];
  }

  qlocal[num_species] = (qglobal[num_species] * norm[0]) + (qglobal[num_species + 1] * norm[1]) + (qglobal[num_species + 2] * norm[2]);
  qlocal[num_species + 1] = (qglobal[num_species] * tau1[0]) + (qglobal[num_species + 1] * tau1[1]) + (qglobal[num_species + 2] * tau1[2]);
  qlocal[num_species + 2] = (qglobal[num_species] * tau2[0]) + (qglobal[num_species + 1] * tau2[1]) + (qglobal[num_species + 2] * tau2[2]);
  qlocal[num_species + 3] = qglobal[num_species + 3];
}

static inline void
rot_to_global(const double* tau1, const double* tau2, const double* norm, const double* GKYL_RESTRICT qlocal, double* GKYL_RESTRICT qglobal)
{
  long num_species = 2;

  for (int i = 0; i < num_species; i++) {
    qglobal[i] = qlocal[i];
  }

  qglobal[num_species] = (qlocal[num_species] * norm[0]) + (qlocal[num_species + 1] * tau1[0]) + (qlocal[num_species + 2] * tau2[0]);
  qglobal[num_species + 1] = (qlocal[num_species] * norm[1]) + (qlocal[num_species + 1] * tau1[1]) + (qlocal[num_species + 2] * tau2[1]);
  qglobal[num_species + 2] = (qlocal[num_species] * norm[2]) + (qlocal[num_species + 1] * tau1[2]) + (qlocal[num_species + 2] * tau2[2]);
  qglobal[num_species + 3] = qlocal[num_species + 3];
}

static double
wave_lax(const struct gkyl_wv_eqn* eqn, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;
  double* gas_gamma_s = euler_mixture->gas_gamma_s;

  double sl = gkyl_euler_mixture_max_abs_speed(num_species, gas_gamma_s, ql);
  double sr = gkyl_euler_mixture_max_abs_speed(num_species, gas_gamma_s, qr);
  double amax = fmax(sl, sr);

  double *fl = gkyl_malloc(sizeof(double[4 + num_species]));
  double *fr = gkyl_malloc(sizeof(double[4 + num_species]));
  gkyl_euler_mixture_flux(num_species, gas_gamma_s, ql, fl);
  gkyl_euler_mixture_flux(num_species, gas_gamma_s, qr, fr);

  double *w0 = &waves[0], *w1 = &waves[4 + num_species];
  for (int i = 0; i < 4 + num_species; i++) {
    w0[i] = 0.5 * ((qr[i] - ql[i]) - (fr[i] - fl[i]) / amax);
    w1[i] = 0.5 * ((qr[i] - ql[i]) + (fr[i] - fl[i]) / amax);
  }

  s[0] = -amax;
  s[1] = amax;

  return s[1];
}

static void
qfluct_lax(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, const double* waves, const double* s, double* amdq, double* apdq)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  const double *w0 = &waves[0], *w1 = &waves[4 + num_species];
  double s0m = fmin(0.0, s[0]), s1m = fmin(0.0, s[1]);
  double s0p = fmax(0.0, s[0]), s1p = fmax(0.0, s[1]);

  for (int i = 0; i < 4 + num_species; i++) {
    amdq[i] = (s0m * w0[i]) + (s1m * w1[i]);
    apdq[i] = (s0p * w0[i]) + (s1p * w1[i]);
  }
}

static double
wave_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* delta, const double* ql, const double* qr, double* waves, double* s)
{
  return wave_lax(eqn, delta, ql, qr, waves, s);
}

static void
qfluct_lax_l(const struct gkyl_wv_eqn* eqn, enum gkyl_wv_flux_type type, const double* ql, const double* qr, const double* waves, const double* s,
  double* amdq, double* apdq)
{
  return qfluct_lax(eqn, ql, qr, waves, s, amdq, apdq);
}

static double
flux_jump(const struct gkyl_wv_eqn* eqn, const double* ql, const double* qr, double* flux_jump)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;
  double *gas_gamma_s = euler_mixture->gas_gamma_s;

  double *fr = gkyl_malloc(sizeof(double[4 + num_species]));
  double *fl = gkyl_malloc(sizeof(double[4 + num_species]));
  gkyl_euler_mixture_flux(num_species, gas_gamma_s, ql, fl);
  gkyl_euler_mixture_flux(num_species, gas_gamma_s, qr, fr);

  for (int m = 0; m < 4 + num_species; m++) {
    flux_jump[m] = fr[m] - fl[m];
  }

  double amaxl = gkyl_euler_mixture_max_abs_speed(num_species, gas_gamma_s, ql);
  double amaxr = gkyl_euler_mixture_max_abs_speed(num_species, gas_gamma_s, qr);

  return fmax(amaxl, amaxr);
}

static bool
check_inv(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;
  double *gas_gamma_s = euler_mixture->gas_gamma_s;

  double *v = gkyl_malloc(sizeof(double[4 + num_species]));
  gkyl_euler_mixture_prim_vars(num_species, gas_gamma_s, q, v);

  for (int i = 0; i < num_species; i++) {
    if (v[i] < 0.0) {
      return false;
    }
  }

  if (v[num_species + 3] < 0.0) {
    return false;
  }
  else {
    return true;
  }
}

static double
max_speed(const struct gkyl_wv_eqn* eqn, const double* q)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;
  double* gas_gamma_s = euler_mixture->gas_gamma_s;

  return gkyl_euler_mixture_max_abs_speed(num_species, gas_gamma_s, q);
}

static inline void
euler_mixture_cons_to_diag(const struct gkyl_wv_eqn* eqn, const double* qin, double* diag)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  for (int i = 0; i < 4 + num_species; i++) {
    diag[i] = qin[i];
  }
}

static inline void
euler_mixture_source(const struct gkyl_wv_eqn* eqn, const double* qin, double* sout)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  for (int i = 0; i < 4 + num_species; i++) {
    sout[i] = 0.0;
  }
}

void
gkyl_euler_mixture_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_wv_eqn* base = container_of(ref, struct gkyl_wv_eqn, ref_count);

  if (gkyl_wv_eqn_is_cu_dev(base)) {
    // Free inner on_dev object.
    struct wv_euler_mixture *euler_mixture = container_of(base->on_dev, struct wv_euler_mixture, eqn);
    gkyl_cu_free(euler_mixture);
  }

  struct wv_euler_mixture *euler_mixture = container_of(base, struct wv_euler_mixture, eqn);
  gkyl_free(euler_mixture);
}

struct gkyl_wv_eqn*
gkyl_wv_euler_mixture_new(int num_species, double* gas_gamma_s, bool use_gpu)
{
  return gkyl_wv_euler_mixture_inew(&(struct gkyl_wv_euler_mixture_inp) {
      .num_species = num_species,
      .gas_gamma_s = gas_gamma_s,
      .rp_type = WV_EULER_MIXTURE_RP_LAX,
      .use_gpu = use_gpu,
    }
  );
}

struct gkyl_wv_eqn*
gkyl_wv_euler_mixture_inew(const struct gkyl_wv_euler_mixture_inp* inp)
{
  struct wv_euler_mixture *euler_mixture = gkyl_malloc(sizeof(struct wv_euler_mixture));

  euler_mixture->eqn.type = GKYL_EQN_EULER_MIXTURE;
  euler_mixture->eqn.num_equations = 4 + inp->num_species;
  euler_mixture->eqn.num_diag = 4 + inp->num_species;

  euler_mixture->num_species = inp->num_species;
  euler_mixture->gas_gamma_s = inp->gas_gamma_s;

  if (inp->rp_type == WV_EULER_MIXTURE_RP_LAX) {
    euler_mixture->eqn.num_waves = 2;
    euler_mixture->eqn.waves_func = wave_lax_l;
    euler_mixture->eqn.qfluct_func = qfluct_lax_l;
  }

  euler_mixture->eqn.flux_jump = flux_jump;
  euler_mixture->eqn.check_inv_func = check_inv;
  euler_mixture->eqn.max_speed_func = max_speed;
  euler_mixture->eqn.rotate_to_local_func = rot_to_local;
  euler_mixture->eqn.rotate_to_global_func = rot_to_global;
  
  euler_mixture->eqn.wall_bc_func = euler_mixture_wall;
  euler_mixture->eqn.no_slip_bc_func = euler_mixture_no_slip;

  euler_mixture->eqn.cons_to_riem = cons_to_riem;
  euler_mixture->eqn.riem_to_cons = riem_to_cons;

  euler_mixture->eqn.cons_to_diag = euler_mixture_cons_to_diag;

  euler_mixture->eqn.source_func = euler_mixture_source;

  euler_mixture->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(euler_mixture->eqn.flags);
  euler_mixture->eqn.ref_count = gkyl_ref_count_init(gkyl_euler_mixture_free);
  euler_mixture->eqn.on_dev = &euler_mixture->eqn; // On the CPU, the equation object points ot itself.

  return &euler_mixture->eqn;
}

int
gkyl_wv_euler_mixture_num_species(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  int num_species = euler_mixture->num_species;

  return num_species;
}

double*
gkyl_wv_euler_mixture_gas_gamma_s(const struct gkyl_wv_eqn* eqn)
{
  const struct wv_euler_mixture *euler_mixture = container_of(eqn, struct wv_euler_mixture, eqn);
  double *gas_gamma_s = euler_mixture->gas_gamma_s;

  return gas_gamma_s;
}