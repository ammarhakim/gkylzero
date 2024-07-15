#include <gkyl_yield_model.h>
#include <math.h>
#include <gkyl_alloc.h>

struct gkyl_yield_model*
gkyl_yield_furman_pivi_new(double charge, double deltahat_ts, double Ehat_ts, double t1, double t2, double t3,
  double t4, double s, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_yield_furman_pivi_cu_dev_new(charge, deltahat_ts, Ehat_ts, t1, t2, t3, t4, s);
  }
#endif
  struct gkyl_yield_furman_pivi *model = gkyl_malloc(sizeof(struct gkyl_yield_furman_pivi));
  
  model->deltahat_ts = deltahat_ts;
  model->Ehat_ts = Ehat_ts;
  model->t1 = t1;
  model->t2 = t2;
  model->t3 = t3;
  model->t4 = t4;
  model->s = s;
  model->yield.charge = charge;
  model->yield.function = furman_pivi_yield;

  model->yield.ref_count = gkyl_ref_count_init(furman_pivi_free);

  return &model->yield;
}

struct gkyl_yield_model*
gkyl_yield_schou_new(double charge, double int_wall, double a2, double a3, double a4, double a5,
  double nw, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_yield_schou_cu_dev_new(charge, int_wall, a2, a3, a4, a5, nw);
  }
#endif
  struct gkyl_yield_schou *model = gkyl_malloc(sizeof(struct gkyl_yield_schou));
  
  model->int_wall = int_wall;
  model->a2 = a2;
  model->a3 = a3;
  model->a4 = a4;
  model->a5 = a5;
  model->nw = nw;
  model->yield.charge = charge;
  model->yield.function = schou_yield;

  model->yield.ref_count = gkyl_ref_count_init(schou_free);

  return &model->yield;
}

struct gkyl_yield_model*
gkyl_yield_constant_new(double charge, double delta, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_yield_constant_cu_dev_new(charge, delta);
  }
#endif
  struct gkyl_yield_constant *model = gkyl_malloc(sizeof(struct gkyl_yield_constant));
  
  model->delta = delta;
  model->yield.charge = charge;
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
