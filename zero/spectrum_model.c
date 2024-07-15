#include <gkyl_spectrum_model.h>
#include <math.h>
#include <gkyl_alloc.h>

struct gkyl_spectrum_model*
gkyl_spectrum_chung_everhart_new(double charge, double phi, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_spectrum_chung_everhart_cu_dev_new(charge, phi);
  }
#endif
  struct gkyl_spectrum_chung_everhart *model = gkyl_malloc(sizeof(struct gkyl_spectrum_chung_everhart));
  
  model->phi = phi;
  model->spectrum.charge = charge;
  model->spectrum.distribution = chung_everhart_dist;
  model->spectrum.normalization = chung_everhart_norm;

  model->spectrum.ref_count = gkyl_ref_count_init(chung_everhart_free);

  return &model->spectrum;
}

struct gkyl_spectrum_model*
gkyl_spectrum_gaussian_new(double charge, double E_0, double tau, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_spectrum_gaussian_cu_dev_new(charge, E_0, tau);
  }
#endif
  struct gkyl_spectrum_gaussian *model = gkyl_malloc(sizeof(struct gkyl_spectrum_gaussian));

  model->E_0 = E_0;
  model->tau = tau;
  model->spectrum.charge = charge;
  model->spectrum.distribution = gaussian_dist;
  model->spectrum.normalization = gaussian_norm;

  model->spectrum.ref_count = gkyl_ref_count_init(gaussian_free);

  return &model->spectrum;
}

struct gkyl_spectrum_model*
gkyl_spectrum_maxwellian_new(double charge, double vt, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_spectrum_maxwellian_cu_dev_new(charge, vt);
  }
#endif
  struct gkyl_spectrum_maxwellian *model = gkyl_malloc(sizeof(struct gkyl_spectrum_maxwellian));

  model->vt = vt;
  model->spectrum.charge = charge;
  model->spectrum.distribution = maxwellian_dist;
  model->spectrum.normalization = maxwellian_norm;

  model->spectrum.ref_count = gkyl_ref_count_init(maxwellian_free);

  return &model->spectrum;
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
