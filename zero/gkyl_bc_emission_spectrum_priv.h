#pragma once

// Private header for bc_emission_spectrum updater, not for direct use in user code.

#include <gkyl_bc_emission_spectrum.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_updater_moment.h>
#include <math.h>
#include <assert.h>

typedef void (*emission_spectrum_func_t)(const double *inp, void *ctx, int idx[GKYL_MAX_DIM], double xc[GKYL_MAX_DIM], double *gain, double *weight);
typedef double (*emission_spectrum_norm_func_t)(double flux, double *param, double effective_gamma);

struct gkyl_bc_emission_spectrum {
  int dir, cdim;
  enum gkyl_edge_loc edge;
  emission_spectrum_func_t func;
  emission_spectrum_norm_func_t norm;
  void *ctx;
  void *ctx_on_dev;
  struct gkyl_range ghost_r;
  const struct gkyl_basis *cbasis;
  bool use_gpu;
};

struct bc_emission_spectrum_ctx {
  int dir; // direction for BCs.
  int cdim, vdim; // config-space dimensions.
  enum gkyl_edge_loc edge;
  const struct gkyl_basis *basis; // basis function.
};

GKYL_CU_D
void
bc_weighted_gamma(const double *inp, void *ctx, int idx[GKYL_MAX_DIM], double xc[GKYL_MAX_DIM], double *gain, double *weight)
{
  struct bc_emission_spectrum_ctx *mc = (struct bc_emission_spectrum_ctx*) ctx;
  int cdim = mc->cdim;
  int dir = mc->dir;
  int nbasis = mc->basis->num_basis;

  // CHANGE NEEDED - verify this is in fact correct, might be flipped in the positive/negative plane
  if ((mc->edge == GKYL_LOWER_EDGE && xc[cdim+dir] < 0) || (mc->edge == GKYL_UPPER_EDGE && xc[cdim+dir] > 0)) { 
    weight[0] += fabs(xc[cdim+dir])*inp[0]*gain[idx[1] - 1]; // CHANGE NEEDED - store coefficients in gkyl_array and index more rigorously
    weight[1] += fabs(xc[cdim+dir])*inp[0];
  }
}

GKYL_CU_D
double
chung_everhart_norm(double flux, double *param, double effective_gamma)
{
  double phi = param[0];
  double mass = param[1];
  double charge = param[2];
  
  return 6.0*effective_gamma*flux*phi*phi*mass/fabs(charge);
}

GKYL_CU_D
double
gaussian_norm(double flux, double *param, double effective_gamma)
{
  double E_0 = param[0];
  double tau = param[1];
  double mass = param[2];
  double charge = param[3];
  
  return effective_gamma*flux*mass/(sqrt(2.0*M_PI)*E_0*tau*exp(tau*tau/2.0)*fabs(charge));
}
