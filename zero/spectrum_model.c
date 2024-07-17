#include <math.h>
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_spectrum_model.h>

struct gkyl_spectrum_model*
gkyl_spectrum_chung_everhart_new(double charge, double phi, bool use_gpu)
{
  struct gkyl_spectrum_chung_everhart *model = gkyl_malloc(sizeof(struct gkyl_spectrum_chung_everhart));
  
  model->phi = phi;
  model->spectrum.charge = charge;
  model->spectrum.distribution = chung_everhart_dist;
  model->spectrum.normalization = chung_everhart_norm;

  model->spectrum.flags = 0;
  GKYL_CLEAR_CU_ALLOC(model->spectrum.flags);
  model->spectrum.ref_count = gkyl_ref_count_init(chung_everhart_free);

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    model->spectrum.on_dev = gkyl_spectrum_chung_everhart_cu_dev_new(model, charge, phi);
  }
#endif

  return &model->spectrum;
}

struct gkyl_spectrum_model*
gkyl_spectrum_gaussian_new(double charge, double E_0, double tau, bool use_gpu)
{
  struct gkyl_spectrum_gaussian *model = gkyl_malloc(sizeof(struct gkyl_spectrum_gaussian));

  model->E_0 = E_0;
  model->tau = tau;
  model->spectrum.charge = charge;
  model->spectrum.distribution = gaussian_dist;
  model->spectrum.normalization = gaussian_norm;

  model->spectrum.flags = 0;
  GKYL_CLEAR_CU_ALLOC(model->spectrum.flags);
  model->spectrum.ref_count = gkyl_ref_count_init(gaussian_free);

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    model->spectrum.on_dev = gkyl_spectrum_gaussian_cu_dev_new(model, charge, E_0, tau);
  }
#endif

  return &model->spectrum;
}

struct gkyl_spectrum_model*
gkyl_spectrum_maxwellian_new(double charge, double vt, bool use_gpu)
{
  struct gkyl_spectrum_maxwellian *model = gkyl_malloc(sizeof(struct gkyl_spectrum_maxwellian));

  model->vt = vt;
  model->spectrum.charge = charge;
  model->spectrum.distribution = maxwellian_dist;
  model->spectrum.normalization = maxwellian_norm;

  model->spectrum.flags = 0;
  GKYL_CLEAR_CU_ALLOC(model->spectrum.flags);
  model->spectrum.ref_count = gkyl_ref_count_init(maxwellian_free);

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    model->spectrum.on_dev = gkyl_spectrum_maxwellian_cu_dev_new(model, charge, vt);
  }
#endif

  return &model->spectrum;
}

bool
gkyl_spectrum_model_is_cu_dev(const struct gkyl_spectrum_model *model)
{
  return GKYL_IS_CU_ALLOC(model->flags);
}

struct gkyl_spectrum_model*
gkyl_spectrum_model_acquire(const struct gkyl_spectrum_model* spectrum)
{
  gkyl_ref_count_inc(&spectrum->ref_count);
  return (struct gkyl_spectrum_model*) spectrum;
}

void
gkyl_spectrum_model_release(const struct gkyl_spectrum_model* spectrum)
{
  gkyl_ref_count_dec(&spectrum->ref_count);
}
