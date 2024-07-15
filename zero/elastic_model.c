#include <gkyl_elastic_model.h>
#include <math.h>
#include <gkyl_alloc.h>

static void
furman_pivi_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_elastic_model *elastic = container_of(ref, struct gkyl_elastic_model, ref_count);
  struct gkyl_elastic_furman_pivi *model = container_of(elastic, struct gkyl_elastic_furman_pivi, elastic);
  gkyl_free(model);
}

// Furman-Pivi SEY calculation
GKYL_CU_D
static void
furman_pivi_yield(double t, const double *xn, double *fout, void *ctx)
// Electron impact model adapted from https://link.aps.org/doi/10.1103/PhysRevSTAB.5.124404
{
  struct gkyl_elastic_model *elastic = (struct gkyl_elastic_model *) ctx;
  const struct gkyl_elastic_furman_pivi *model = container_of(elastic,
    struct gkyl_elastic_furman_pivi, elastic);

  int cdim = elastic->cdim;
  int vdim = elastic->vdim;
  double mass = elastic->mass;
  double charge = elastic->charge;

  double P1_inf = model->P1_inf;
  double P1_hat = model->P1_hat;
  double E_hat = model->E_hat;
  double W = model->W;
  double p = model->p;

  double E = 0.0;
  double mu = 1.0; // currently hardcoded to normal, will add angular dependence later
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xn[cdim+d]*xn[cdim+d]/fabs(charge); // Calculate energy in eV
  }
  fout[0] = P1_inf + (P1_hat - P1_inf)*exp(pow(-fabs(E - E_hat)/W, p)/p);
}

struct gkyl_elastic_model*
gkyl_elastic_furman_pivi_new(double P1_inf, double P1_hat, double E_hat, double W, double p)
{
  struct gkyl_elastic_furman_pivi *model = gkyl_malloc(sizeof(struct gkyl_elastic_furman_pivi));
  
  model->P1_inf = P1_inf;
  model->P1_hat = P1_hat;
  model->E_hat = E_hat;
  model->W = W;
  model->p = p;
  model->elastic.function = furman_pivi_yield;

  model->elastic.ref_count = gkyl_ref_count_init(furman_pivi_free);

  return &model->elastic;
}

struct gkyl_elastic_model*
gkyl_elastic_model_acquire(const struct gkyl_elastic_model* model)
{
  gkyl_ref_count_inc(&model->ref_count);
  return (struct gkyl_elastic_model*) model;
}

void
gkyl_elastic_model_release(const struct gkyl_elastic_model* model)
{
  gkyl_ref_count_dec(&model->ref_count);
}
