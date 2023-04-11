#pragma once

// Private header for bc_emission_spectrum updater, not for direct use in user code.

#include <gkyl_bc_emission_spectrum.h>
#include <gkyl_array_ops.h>
#include <gkyl_proj_on_basis.h>
#include <gkyl_dg_updater_moment.h>
#include <math.h>
#include <assert.h>

typedef void (*emission_spectrum_func_t)(double *out, const double *inp, void *ctx, int idx[GKYL_MAX_DIM], double xc[GKYL_MAX_DIM]);
typedef void (*emission_spectrum_proj_func_t)(double t, const double *xn, double *out, void *ctx);
typedef double (*emission_spectrum_norm_func_t)(double *tot_flux, void *ctx);

struct gkyl_bc_emission_spectrum {
  int dir, cdim;
  enum gkyl_edge_loc edge;
  emission_spectrum_func_t func;
  emission_spectrum_proj_func_t proj;
  emission_spectrum_norm_func_t norm;
  struct gkyl_rect_grid *grid; // grid object
  void *ctx;
  void *ctx_on_dev;
  struct gkyl_range skin_r, ghost_r, pos_r, neg_r, cskin_r, cghost_r;
  struct gkyl_dg_updater_moment *mcalc;
  struct gkyl_array *flux;
  struct gkyl_array *f_proj;
  double tot_flux;
  const struct gkyl_basis *cbasis;
  gkyl_proj_on_basis *boundary_proj;
  bool use_gpu;
};

struct bc_emission_spectrum_ctx {
  int dir; // direction for BCs.
  int cdim, vdim; // config-space dimensions.
  double *gain;
  double *elastic;
  double *param;
  double numerator, denominator;
  enum gkyl_edge_loc edge;
  const struct gkyl_basis *basis; // basis function.
};

GKYL_CU_D
static void
bc_weighted_gamma(double *out, const double *inp, void *ctx, int idx[GKYL_MAX_DIM], double xc[GKYL_MAX_DIM])
{
  struct bc_emission_spectrum_ctx *mc = (struct bc_emission_spectrum_ctx*) ctx;
  int cdim = mc->cdim;
  int dir = mc->dir;
  int nbasis = mc->basis->num_basis;

  if (mc->elastic) {
    mc->basis->flip_odd_sign(dir, inp, out);
    for (int c=0; c<nbasis; ++c) out[c] = mc->elastic[idx[1] - 1]*out[c];
    mc->basis->flip_odd_sign(dir+cdim, out, out);
  }

  // CHANGE NEEDED - verify this is in fact correct, might be flipped in the positive/negative plane
  if ((mc->edge == GKYL_LOWER_EDGE && xc[cdim+dir] < 0) || (mc->edge == GKYL_UPPER_EDGE && xc[cdim+dir] > 0)) { 
    mc->numerator += fabs(xc[cdim+dir])*inp[0]*mc->gain[idx[1] - 1]; // CHANGE NEEDED - store coefficients in gkyl_array and index more rigorously
    mc->denominator += fabs(xc[cdim+dir])*inp[0];
  }
}

GKYL_CU_D
double
chung_everhart_norm(double *tot_flux, void *ctx)
{
  struct bc_emission_spectrum_ctx *mc = (struct bc_emission_spectrum_ctx*) ctx;

  double phi = mc->param[0];
  double mass = mc->param[1];
  double charge = mc->param[2];
  
  double effective_gamma = mc->numerator/mc->denominator;
  
  return 6.0*effective_gamma*tot_flux[0]*phi*phi*mass/fabs(charge);
}

GKYL_CU_D
static void 
chung_everhart_proj(double t, const double *xn, double *out, void *ctx)
{
  struct bc_emission_spectrum_ctx *mc = (struct bc_emission_spectrum_ctx*) ctx;

  int cdim = mc->cdim;
  int vdim = mc->vdim;
  double phi = mc->param[0];
  double mass = mc->param[1];
  double charge = mc->param[2];
  
  double E = 0.0;
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xn[cdim+d]*xn[cdim+d]/fabs(charge);
  }
  out[0] = E/((E + phi)*(E + phi)*(E + phi)*(E + phi));
}

GKYL_CU_D
double
gaussian_norm(double *tot_flux, void *ctx)
{
  struct bc_emission_spectrum_ctx *mc = (struct bc_emission_spectrum_ctx*) ctx;

  double E_0 = mc->param[0];
  double tau = mc->param[1];
  double mass = mc->param[2];
  double charge = mc->param[3];
  
  double effective_gamma = mc->numerator/mc->denominator;
  
  return effective_gamma*tot_flux[0]*mass/(sqrt(2.0*M_PI)*E_0*tau*exp(tau*tau/2.0)*fabs(charge));
}

GKYL_CU_D
static void 
gaussian_proj(double t, const double *xn, double *out, void *ctx)
{
  struct bc_emission_spectrum_ctx *mc = (struct bc_emission_spectrum_ctx*) ctx;
  
  int cdim = mc->cdim;
  int vdim = mc->vdim;
  double E_0 = mc->param[0];
  double tau = mc->param[1];
  double mass = mc->param[2];
  double charge = mc->param[3];
  
  double E = 0.0;
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xn[cdim+d]*xn[cdim+d]/fabs(charge);
  }
  out[0] = exp(-pow(log(E/E_0), 2)/(2.0*tau*tau));
}
