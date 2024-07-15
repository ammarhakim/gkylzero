#include <gkyl_spectrum_model.h>
#include <math.h>
#include <gkyl_alloc.h>

static void
chung_everhart_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_spectrum_model *spectrum = container_of(ref, struct gkyl_spectrum_model, ref_count);
  struct gkyl_spectrum_chung_everhart *model = container_of(spectrum, struct gkyl_spectrum_chung_everhart, spectrum);
  gkyl_free(model);
}

GKYL_CU_D
static void
chung_everhart_dist(double t, const double *xn, double *fout, void *ctx)
{
  struct gkyl_spectrum_model *spectrum = (struct gkyl_spectrum_model *) ctx;
  const struct gkyl_spectrum_chung_everhart *model = container_of(spectrum,
    struct gkyl_spectrum_chung_everhart, spectrum);
  
  int cdim = spectrum->cdim;
  int vdim = spectrum->vdim;
  double mass = spectrum->mass;
  double charge = spectrum->charge;
  double phi = model->phi;

  double E = 0.0;
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xn[cdim+d]*xn[cdim+d]/fabs(charge);
  }

  fout[0] = E/pow(E + phi, 4);
}

// Chung-Everhart normalization factor
GKYL_CU_D
static void
chung_everhart_norm(double *out, struct gkyl_spectrum_model *spectrum, const double *flux,
  double effective_delta)
{
  const struct gkyl_spectrum_chung_everhart *model = container_of(spectrum,
    struct gkyl_spectrum_chung_everhart, spectrum);
  double mass = spectrum->mass;
  double charge = spectrum->charge;
  double phi = model->phi;
  
  out[0] = 6.0*effective_delta*flux[0]*phi*phi*mass/fabs(charge);
}

struct gkyl_spectrum_model*
gkyl_spectrum_chung_everhart_new(double phi)
{
  struct gkyl_spectrum_chung_everhart *model = gkyl_malloc(sizeof(struct gkyl_spectrum_chung_everhart));
  
  model->phi = phi;
  model->spectrum.distribution = chung_everhart_dist;
  model->spectrum.normalization = chung_everhart_norm;

  model->spectrum.ref_count = gkyl_ref_count_init(chung_everhart_free);

  return &model->spectrum;
}

/* struct gkyl_spectrum_model* */
/* gkyl_spectrum_gaussian_new(double E_0, double tau); */
/* { */
/*   struct gkyl_spectrum_gaussian *spectrum = gkyl_malloc(sizeof(struct gkyl_spectrum_gaussian)); */
/* } */

/* struct gkyl_spectrum_model* */
/* gkyl_spectrum_maxwellian_new(double vt); */
/* { */
/*   struct gkyl_spectrum_maxwellian *spectrum = gkyl_malloc(sizeof(struct gkyl_spectrum_maxwellian)); */
/* } */

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
