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
