/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_elastic_model.h> 
}

#include <cassert>  

__global__ static void
furman_pivi_set_cu_dev_ptrs(struct gkyl_elastic_furman_pivi *model)
{
  model->elastic.function = gkyl_elastic_furman_pivi_yield;
}

__global__ static void
cazaux_set_cu_dev_ptrs(struct gkyl_elastic_cazaux *model)
{
  model->elastic.function = gkyl_elastic_cazaux_yield;
}

__global__ static void
constant_set_cu_dev_ptrs(struct gkyl_elastic_constant *model)
{
  model->elastic.function = gkyl_elastic_constant_yield;
}

struct gkyl_elastic_model*
gkyl_elastic_furman_pivi_cu_dev_new(double charge, double P1_inf, double P1_hat, double E_hat, double W, double p)
{
  struct gkyl_elastic_furman_pivi *model = (struct gkyl_elastic_furman_pivi*) gkyl_malloc(sizeof(struct gkyl_elastic_furman_pivi));

  model->P1_inf = P1_inf;
  model->P1_hat = P1_hat;
  model->E_hat = E_hat;
  model->W = W;
  model->p = p;
  model->elastic.charge = charge;

  model->elastic.ref_count = gkyl_ref_count_init(gkyl_elastic_furman_pivi_free);

  struct gkyl_elastic_furman_pivi *model_cu = (struct gkyl_elastic_furman_pivi*)
    gkyl_cu_malloc(sizeof(struct gkyl_elastic_furman_pivi));
  gkyl_cu_memcpy(model_cu, model, sizeof(struct gkyl_elastic_furman_pivi), GKYL_CU_MEMCPY_H2D);

  furman_pivi_set_cu_dev_ptrs<<<1,1>>>(model_cu);

  model->elastic.on_dev = &model_cu->elastic;

  return &model->elastic;
}

struct gkyl_elastic_model*
gkyl_elastic_cazaux_cu_dev_new(double charge, double E_f, double phi)
{
  struct gkyl_elastic_cazaux *model = (struct gkyl_elastic_cazaux*) gkyl_malloc(sizeof(struct gkyl_elastic_cazaux));

  model->E_f = E_f;
  model->phi = phi;
  model->elastic.charge = charge;

  model->elastic.ref_count = gkyl_ref_count_init(gkyl_elastic_cazaux_free);

  struct gkyl_elastic_cazaux *model_cu = (struct gkyl_elastic_cazaux*)
    gkyl_cu_malloc(sizeof(struct gkyl_elastic_cazaux));
  gkyl_cu_memcpy(model_cu, model, sizeof(struct gkyl_elastic_cazaux), GKYL_CU_MEMCPY_H2D);

  cazaux_set_cu_dev_ptrs<<<1,1>>>(model_cu);

  model->elastic.on_dev = &model_cu->elastic;

  return &model->elastic;
}

struct gkyl_elastic_model*
gkyl_elastic_constant_cu_dev_new(double charge, double delta)
{
  struct gkyl_elastic_constant *model = (struct gkyl_elastic_constant*) gkyl_malloc(sizeof(struct gkyl_elastic_constant));

  model->delta = delta;
  model->elastic.charge = charge;

  model->elastic.ref_count = gkyl_ref_count_init(gkyl_elastic_constant_free);

  struct gkyl_elastic_constant *model_cu = (struct gkyl_elastic_constant*)
    gkyl_cu_malloc(sizeof(struct gkyl_elastic_constant));
  gkyl_cu_memcpy(model_cu, model, sizeof(struct gkyl_elastic_constant), GKYL_CU_MEMCPY_H2D);

  constant_set_cu_dev_ptrs<<<1,1>>>(model_cu);

  model->elastic.on_dev = &model_cu->elastic;

  return &model->elastic;
}
