#include <math.h>
#include <stddef.h>
#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_const.h>
#include <gkyl_dist_type.h>
#include <gkyl_dist_type_priv.h>


struct gkyl_vlasov_multi_temp_max {
  struct gkyl_dist_proj_type dist;

  enum gkyl_model_id model_id;
  int vdim;
};

/**
 * Free the multi-temperature Maxwellian distribution projection object.
 */
static void
vlasov_multi_temp_max_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_vlasov_multi_temp_max *multi_temp_max =
    (struct gkyl_vlasov_multi_temp_max *)
    ((char *) ref- offsetof(struct gkyl_vlasov_multi_temp_max, dist.ref_count));

  gkyl_free(multi_temp_max);
}


/**
 * Non-relativistic multi-temperature Maxwellian distribution.
 *
 * moms layout:
 *
 *   moms[0]                   = n
 *   moms[1 ... vdim]          = V_drift
 *   moms[vdim+1 ... 2*vdim]   = T_i/m
 */
static void
vlasov_multi_temp_max_default(
  const struct gkyl_dist_proj_type *dist_proj,
  const double *xc,
  const double *dx,
  const int *idx,
  const double *moms,
  const double *param,
  double *out)
{
  const struct gkyl_vlasov_multi_temp_max *multi_temp_max = (const struct gkyl_vlasov_multi_temp_max *) dist_proj;

  int cdim = dist_proj->cdim;
  int vdim = multi_temp_max->vdim;

  double n = moms[0];
  const double *V_drift = &moms[1];
  const double *T_over_m = &moms[vdim+1];

  if (n <= 0.0) {
    out[0] = 0.0;
    return;
  }

  double efact = 0.0;
  double amplitude = n;

  for (int d = 0; d < vdim; ++d) {
    if (T_over_m[d] <= 0.0) {
      out[0] = 0.0;
      return;
    }

    double v = xc[cdim+d] - V_drift[d];

    efact += v*v/(2.0*T_over_m[d]);

    amplitude /= sqrt(2.0*GKYL_PI*T_over_m[d]);
  }

  out[0] = amplitude*exp(-efact);
}

struct gkyl_dist_proj_type*
gkyl_vlasov_multi_temp_max_new(
  int cdim,
  int pdim,
  int poly_order,
  int num_config,
  int num_phase,
  enum gkyl_model_id model_id)
{
  struct gkyl_vlasov_multi_temp_max *multi_temp_max = gkyl_malloc(sizeof(*multi_temp_max));

  multi_temp_max->model_id = model_id;
  multi_temp_max->vdim = pdim-cdim;

  multi_temp_max->dist.dist_id = GKYL_DIST_TYPE_MULTI_TEMP_MAX;

  multi_temp_max->dist.cdim = cdim;

  multi_temp_max->dist.pdim = pdim;

  multi_temp_max->dist.poly_order = poly_order;

  multi_temp_max->dist.num_config =  num_config;

  multi_temp_max->dist.num_phase = num_phase;

  // n + V_drift[vdim] + T_i/m[vdim]
  multi_temp_max->dist.num_mom = 2*multi_temp_max->vdim + 1;

  switch (model_id) {
    default:
      multi_temp_max->dist.proj_kernel = vlasov_multi_temp_max_default;
      break;
  }

  multi_temp_max->dist.ref_count = gkyl_ref_count_init(vlasov_multi_temp_max_free);

  multi_temp_max->dist.flags = 0;
  multi_temp_max->dist.on_dev = &multi_temp_max->dist;

  return &multi_temp_max->dist;
}