#include <gkyl_elastic_model.h>
#include <math.h>
#include <gkyl_alloc.h>

struct gkyl_elastic_model*
gkyl_elastic_furman_pivi_new(double charge, double P1_inf, double P1_hat, double E_hat, double W, double p, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_elastic_furman_pivi_cu_dev_new(charge, P1_inf, P1_hat, E_hat, W, p);
  }
#endif
  struct gkyl_elastic_furman_pivi *model = gkyl_malloc(sizeof(struct gkyl_elastic_furman_pivi));

  model->P1_inf = P1_inf;
  model->P1_hat = P1_hat;
  model->E_hat = E_hat;
  model->W = W;
  model->p = p;
  model->elastic.charge = charge;
  model->elastic.function = gkyl_elastic_furman_pivi_yield;

  model->elastic.ref_count = gkyl_ref_count_init(gkyl_elastic_furman_pivi_free);

  return &model->elastic;
}

struct gkyl_elastic_model*
gkyl_elastic_cazaux_new(double charge, double E_f, double phi, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_elastic_cazaux_cu_dev_new(charge, E_f, phi);
  }
#endif
  struct gkyl_elastic_cazaux *model = gkyl_malloc(sizeof(struct gkyl_elastic_cazaux));
  
  model->E_f = E_f;
  model->phi = phi;
  model->elastic.charge = charge;
  model->elastic.function = gkyl_elastic_cazaux_yield;

  model->elastic.ref_count = gkyl_ref_count_init(gkyl_elastic_cazaux_free);

  return &model->elastic;
}

struct gkyl_elastic_model*
gkyl_elastic_constant_new(double charge, double delta, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_elastic_constant_cu_dev_new(charge, delta);
  }
#endif
  struct gkyl_elastic_constant *model = gkyl_malloc(sizeof(struct gkyl_elastic_constant));
  
  model->delta = delta;
  model->elastic.charge = charge;
  model->elastic.function = gkyl_elastic_constant_yield;

  model->elastic.ref_count = gkyl_ref_count_init(gkyl_elastic_constant_free);

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
