/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_spectrum_model.h> 
}

#include <cassert>  

__global__ static void
chung_everhart_set_cu_dev_ptrs(struct gkyl_spectrum_chung_everhart *model)
{
  model->spectrum.distribution = chung_everhart_dist;
  model->spectrum.normalization = chung_everhart_norm;
}

__global__ static void
gaussian_set_cu_dev_ptrs(struct gkyl_spectrum_gaussian *model)
{
  model->spectrum.distribution = gaussian_dist;
  model->spectrum.normalization = gaussian_norm;
}

__global__ static void
maxwellian_set_cu_dev_ptrs(struct gkyl_spectrum_maxwellian *model)
{
  model->spectrum.distribution = maxwellian_dist;
  model->spectrum.normalization = maxwellian_norm;
}

struct gkyl_spectrum_model*
gkyl_spectrum_chung_everhart_cu_dev_new(double charge, double phi)
{
  struct gkyl_spectrum_chung_everhart *model = (struct gkyl_spectrum_chung_everhart*) gkyl_malloc(sizeof(struct gkyl_spectrum_chung_everhart));

  model->phi = phi;
  model->spectrum.charge = charge;

  model->spectrum.ref_count = gkyl_ref_count_init(chung_everhart_free);

  struct gkyl_spectrum_chung_everhart *model_cu = (struct gkyl_spectrum_chung_everhart*) gkyl_cu_malloc(sizeof(struct gkyl_spectrum_chung_everhart));
  gkyl_cu_memcpy(model_cu, model, sizeof(struct gkyl_spectrum_chung_everhart), GKYL_CU_MEMCPY_H2D);

  chung_everhart_set_cu_dev_ptrs<<<1,1>>>(model_cu);

  model->spectrum.on_dev = &model_cu->spectrum;

  return &model->spectrum;
}

struct gkyl_spectrum_model*
gkyl_spectrum_gaussian_cu_dev_new(double charge, double E_0, double tau)
{
  struct gkyl_spectrum_gaussian *model = (struct gkyl_spectrum_gaussian*) gkyl_malloc(sizeof(struct gkyl_spectrum_gaussian));

  model->E_0 = E_0;
  model->tau = tau;
  model->spectrum.charge = charge;

  model->spectrum.ref_count = gkyl_ref_count_init(gaussian_free);

  struct gkyl_spectrum_gaussian *model_cu = (struct gkyl_spectrum_gaussian*) gkyl_cu_malloc(sizeof(struct gkyl_spectrum_gaussian));
  gkyl_cu_memcpy(model_cu, model, sizeof(struct gkyl_spectrum_gaussian), GKYL_CU_MEMCPY_H2D);

  gaussian_set_cu_dev_ptrs<<<1,1>>>(model_cu);

  model->spectrum.on_dev = &model_cu->spectrum;

  return &model->spectrum;
}

struct gkyl_spectrum_model*
gkyl_spectrum_maxwellian_cu_dev_new(double charge, double vt)
{
  struct gkyl_spectrum_maxwellian *model = (struct gkyl_spectrum_maxwellian*) gkyl_malloc(sizeof(struct gkyl_spectrum_maxwellian));

  model->vt = vt;
  model->spectrum.charge = charge;

  model->spectrum.ref_count = gkyl_ref_count_init(maxwellian_free);

  struct gkyl_spectrum_maxwellian *model_cu = (struct gkyl_spectrum_maxwellian*) gkyl_cu_malloc(sizeof(struct gkyl_spectrum_maxwellian));
  gkyl_cu_memcpy(model_cu, model, sizeof(struct gkyl_spectrum_maxwellian), GKYL_CU_MEMCPY_H2D);

  maxwellian_set_cu_dev_ptrs<<<1,1>>>(model_cu);

  model->spectrum.on_dev = &model_cu->spectrum;

  return &model->spectrum;
}
