#include <gkyl_yield_model.h>
#include <math.h>
#include <gkyl_alloc.h>

static void
furman_pivi_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_yield_model *yield = container_of(ref, struct gkyl_yield_model, ref_count);
  struct gkyl_yield_furman_pivi *model = container_of(yield, struct gkyl_yield_furman_pivi, yield);
  gkyl_free(model);
}

static void
schou_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_yield_model *yield = container_of(ref, struct gkyl_yield_model, ref_count);
  struct gkyl_yield_schou *model = container_of(yield, struct gkyl_yield_schou, yield);
  gkyl_free(model);
}

static void
constant_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_yield_model *yield = container_of(ref, struct gkyl_yield_model, ref_count);
  struct gkyl_yield_constant *model = container_of(yield, struct gkyl_yield_constant, yield);
  gkyl_free(model);
}

// Furman-Pivi SEY calculation */
GKYL_CU_D
static void
furman_pivi_yield(double *out, struct gkyl_yield_model *yield, double xc[GKYL_MAX_DIM])
// Electron impact model adapted from https://link.aps.org/doi/10.1103/PhysRevSTAB.5.124404
{
  const struct gkyl_yield_furman_pivi *model = container_of(yield,
    struct gkyl_yield_furman_pivi, yield);
  
  int cdim = yield->cdim;
  int vdim = yield->vdim;
  double mass = yield->mass;
  double charge = yield->charge;

  double deltahat_ts = model->deltahat_ts;
  double Ehat_ts = model->Ehat_ts;
  double t1 = model->t1;
  double t2 = model->t2;
  double t3 = model->t3;
  double t4 = model->t4;
  double s = model->s;
  
  double E = 0.0;
  double mu = 1.0; // currently hardcoded to normal, will add angular dependence later
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xc[cdim+d]*xc[cdim+d]/fabs(charge); // Calculate energy in eV
  }
  double deltahat = deltahat_ts*(1 + t1*(1 - pow(mu, t2)));
  double Ehat = Ehat_ts*(1 + t3*(1 - pow(mu, t4)));
  double x = E/Ehat;
  out[0] = deltahat*s*x/(s - 1 + pow(x, s));
}

// Schou SEY calculation
GKYL_CU_D
static void
schou_yield(double *out, struct gkyl_yield_model *yield, double xc[GKYL_MAX_DIM])
// Ion impact model adapted from https://doi.org/10.1103/PhysRevB.22.2141
{ // No angular dependence atm. Will have to add later
  const struct gkyl_yield_schou *model = container_of(yield,
    struct gkyl_yield_schou, yield);
  int cdim = yield->cdim;
  int vdim = yield->vdim;
  double mass = yield->mass;
  double charge = yield->charge;
  double int_wall = model->int_wall;
  double A2 = model->a2;  // Note: Starts at 2 to match notation from source: https://doi.org/10.1093/jicru_os25.2.18
  double A3 = model->a3;
  double A4 = model->a4;
  double A5 = model->a5;
  double nw = model->nw;   // Number density of wall material in m^-3

  double E = 0.0;
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xc[cdim+d]*xc[cdim+d]/fabs(charge)/1000;  // Calculate energy in keV
  }

  double Es = E*1.66053906660e-27/1.67262192369e-27;  // Scale energy by ratio of atomic mass unit to proton mass
  double eps_low = A2*pow(Es,0.45);
  double eps_high = 0.0;
  if (Es != 0) { // Divide by zero error catching. If Es is not 0, then do the normal calculation.
    eps_high = (A3/Es)*log(1+A4/Es+A5*Es);
  }
  double eps = eps_low*eps_high/(eps_low+eps_high);

  out[0] =  eps*nw*fabs(charge)*int_wall/1.0e19;
}

// Fixed constant SEY
GKYL_CU_D
static void
constant_yield(double *out, struct gkyl_yield_model *yield, double xc[GKYL_MAX_DIM])
{
  const struct gkyl_yield_constant *model = container_of(yield,
    struct gkyl_yield_constant, yield);
  double delta = model->delta;

  out[0] = delta;
}

struct gkyl_yield_model*
gkyl_yield_furman_pivi_new(double deltahat_ts, double Ehat_ts, double t1, double t2, double t3,
  double t4, double s)
{
  struct gkyl_yield_furman_pivi *model = gkyl_malloc(sizeof(struct gkyl_yield_furman_pivi));
  
  model->deltahat_ts = deltahat_ts;
  model->Ehat_ts = Ehat_ts;
  model->t1 = t1;
  model->t2 = t2;
  model->t3 = t3;
  model->t4 = t4;
  model->s = s;
  model->yield.function = furman_pivi_yield;

  model->yield.ref_count = gkyl_ref_count_init(furman_pivi_free);

  return &model->yield;
}

struct gkyl_yield_model*
gkyl_yield_schou_new(double int_wall, double a2, double a3, double a4, double a5, double nw)
{
  struct gkyl_yield_schou *model = gkyl_malloc(sizeof(struct gkyl_yield_schou));
  
  model->int_wall = int_wall;
  model->a2 = a2;
  model->a3 = a3;
  model->a4 = a4;
  model->a5 = a5;
  model->nw = nw;
  model->yield.function = schou_yield;

  model->yield.ref_count = gkyl_ref_count_init(schou_free);

  return &model->yield;
}

struct gkyl_yield_model*
gkyl_yield_constant_new(double delta)
{
  struct gkyl_yield_constant *model = gkyl_malloc(sizeof(struct gkyl_yield_constant));
  
  model->delta = delta;
  model->yield.function = constant_yield;

  model->yield.ref_count = gkyl_ref_count_init(constant_free);

  return &model->yield;
}

struct gkyl_yield_model*
gkyl_yield_model_acquire(const struct gkyl_yield_model* model)
{
  gkyl_ref_count_inc(&model->ref_count);
  return (struct gkyl_yield_model*) model;
}

void
gkyl_yield_model_release(const struct gkyl_yield_model* model)
{
  gkyl_ref_count_dec(&model->ref_count);
}
