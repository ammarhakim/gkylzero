
#include <math.h>
#include <stddef.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_dist_type.h>
#include <gkyl_dist_type_priv.h>


struct gkyl_vlasov_kappa {
  struct gkyl_dist_proj_type dist;

  enum gkyl_model_id model_id;
  int vdim;
};

/**
 * Free the Kappa distribution projection object.
 */
static void
vlasov_kappa_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_vlasov_kappa *kappa =
    (struct gkyl_vlasov_kappa *)
    ((char *) ref - offsetof(struct gkyl_vlasov_kappa, dist.ref_count));

  gkyl_free(kappa);
}


/**
 * Non-relativistic isotropic Kappa distribution.
 *
 * moms layout:
 *
 *   moms[0]              = n
 *   moms[1 ... vdim]     = V_drift
 *   moms[vdim+1]         = T/m
 *
 * param layout:
 *
 *   param[0]             = kappa
 */
static void
vlasov_kappa_default(
  const struct gkyl_dist_proj_type *dist_proj,
  const double *xc,
  const double *dx,
  const int *idx,
  const double *moms,
  const double *param,
  double *out)
{
  const struct gkyl_vlasov_kappa *kappa = (const struct gkyl_vlasov_kappa *) dist_proj;

  int cdim = dist_proj->cdim;
  int vdim = kappa->vdim;

  double n = moms[0];
  const double *V_drift = &moms[1];
  double T_over_m = moms[vdim+1];
  double kappa_param = param[0];

  if (n <= 0.0 || T_over_m <= 0.0 || kappa_param <= 0.5*vdim) {
    out[0] = 0.0;
    return;
  }

  double efact = 0.0;

  for (int d = 0; d < vdim; ++d) {
    double v = xc[cdim+d] - V_drift[d];

    efact += v*v;
  }

  double alpha = (2.0*kappa_param-vdim)*T_over_m;

  double amplitude = n * exp(
    lgamma(kappa_param+1.0)
    - lgamma(kappa_param+1.0-0.5*vdim)
    - 0.5*vdim*log(GKYL_PI*alpha)
  );

  out[0] = amplitude * pow(1.0 + efact/alpha, -(kappa_param+1.0));
}


struct gkyl_dist_proj_type*
gkyl_vlasov_kappa_new(
  int cdim,
  int pdim,
  int poly_order,
  int num_config,
  int num_phase,
  enum gkyl_model_id model_id)
{
  assert(pdim > cdim);

  struct gkyl_vlasov_kappa *kappa =
    gkyl_malloc(sizeof(*kappa));

  kappa->model_id = model_id;
  kappa->vdim = pdim-cdim;

  kappa->dist.dist_id = GKYL_DIST_TYPE_KAPPA;

  kappa->dist.cdim = cdim;

  kappa->dist.pdim = pdim;

  kappa->dist.poly_order = poly_order;

  kappa->dist.num_config = num_config;
  
  kappa->dist.num_phase = num_phase;

  // n + V_drift[vdim] + T/m
  kappa->dist.num_mom = kappa->vdim + 2;

  switch (model_id) {
    default:
      kappa->dist.proj_kernel = vlasov_kappa_default;
      break;
  }

  kappa->dist.ref_count =
    gkyl_ref_count_init(vlasov_kappa_free);

  kappa->dist.flags = 0;
  kappa->dist.on_dev = &kappa->dist;

  return &kappa->dist;
}
