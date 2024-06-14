#pragma once

// Private header for bc_emission_spectrum updater, not for direct use in user code.

#include <gkyl_bc_emission_spectrum.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_updater_moment.h>
#include <math.h>
#include <assert.h>
#include <gkyl_alloc.h>

typedef void (*emission_spectrum_func_t)(const double *inp, int cdim, int dir, enum gkyl_edge_loc edge, double xc[GKYL_MAX_DIM], const double *gain, double *weight);
typedef void (*emission_spectrum_spec_func_t)(double t, const double *xn, double *fout, void *ctx);
typedef void (*emission_spectrum_norm_func_t)(double *out, const double *flux, void *norm_param, double effective_delta);
typedef void (*emission_spectrum_yield_func_t)(double *out, int cdim, int vdim, double xc[GKYL_MAX_DIM], void *yield_param);

struct gkyl_bc_emission_spectrum_funcs {
  emission_spectrum_func_t func;
  emission_spectrum_spec_func_t spec;
  emission_spectrum_norm_func_t norm;
  emission_spectrum_yield_func_t yield;
  void *norm_param;
  void *yield_param;
};

// Primary struct for the updater
struct gkyl_bc_emission_spectrum {
  int dir, cdim, vdim;
  enum gkyl_edge_loc edge;
  double charge;
  double mass;
  struct gkyl_rect_grid *grid;
  struct gkyl_bc_emission_spectrum_funcs *funcs;
  struct gkyl_bc_emission_spectrum_funcs *funcs_cu;
  struct gkyl_array *flux;
  struct gkyl_range *pos_r;
  struct gkyl_range *neg_r;
  bool use_gpu;
};

GKYL_CU_D
static void
chung_everhart_spec(double t, const double *xn, double *fout, void *ctx)
{
  struct gkyl_bc_emission_spectrum *bc_ctx =
    (struct gkyl_bc_emission_spectrum *) ctx;
  struct gkyl_bc_emission_spectrum_norm_chung_everhart *param =
    (struct gkyl_bc_emission_spectrum_norm_chung_everhart *) bc_ctx->funcs->norm_param;
  int cdim = bc_ctx->cdim;
  int vdim = bc_ctx->vdim;
  double mass = param->mass;
  double charge = param->charge;
  double phi = param->phi;

  double E = 0.0;
  double mu = 1.0; // currently hardcoded to normal, will add angular dependence later
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xn[cdim+d]*xn[cdim+d]/fabs(charge);
  }

  fout[0] = E/pow(E + phi, 4);
}

GKYL_CU_D
static void
gaussian_spec(double t, const double *xn, double *fout, void *ctx)
{
  struct gkyl_bc_emission_spectrum *bc_ctx =
    (struct gkyl_bc_emission_spectrum *) ctx;
  struct gkyl_bc_emission_spectrum_norm_gaussian *param =
    (struct gkyl_bc_emission_spectrum_norm_gaussian *) bc_ctx->funcs->norm_param;
  int cdim = bc_ctx->cdim;
  int vdim = bc_ctx->vdim;
  double mass = param->mass;
  double charge = param->charge;
  double E_0 = param->E_0;
  double tau = param->tau;

  double E = 0.0;
  double mu = 1.0; // currently hardcoded to normal, will add angular dependence later
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xn[cdim+d]*xn[cdim+d]/fabs(charge);
  }

  fout[0] = exp(-pow(log(E/E_0), 2)/(2.0*pow(tau, 2)));
}

GKYL_CU_D
static void
maxwellian_spec(double t, const double *xn, double *fout, void *ctx)
{
  struct gkyl_bc_emission_spectrum *bc_ctx = 
    (struct gkyl_bc_emission_spectrum *) ctx;
  struct gkyl_bc_emission_spectrum_norm_maxwellian *param = 
    (struct gkyl_bc_emission_spectrum_norm_maxwellian *) bc_ctx->funcs->norm_param;
  int cdim = bc_ctx->cdim;
  int vdim = bc_ctx->vdim;
  double mass = param->mass;
  double charge = param->charge;
  double vt = param->vt;

  double v_sq = 0.0;
  for (int d=0; d<vdim; d++) {
    v_sq += xn[cdim+d]*xn[cdim+d];
  }

  fout[0] = exp(-v_sq/(2.0*pow(vt, 2)));
}

// Function to calculate the weighted mean of the SE yield
GKYL_CU_D
static void
bc_weighted_delta(const double *inp, int cdim, int dir, enum gkyl_edge_loc edge, double xc[GKYL_MAX_DIM], const double *gain, double *weight)
{
  if ((edge == GKYL_LOWER_EDGE && xc[cdim+dir] < 0) || (edge == GKYL_UPPER_EDGE && xc[cdim+dir] > 0)) { 
    weight[0] += inp[0]*gain[0];
    weight[1] += inp[0];
  }
}

// Chung-Everhart normalization factor
GKYL_CU_D
static void
chung_everhart_norm(double *out, const double *flux, void *norm_param, double effective_delta)
{
  struct gkyl_bc_emission_spectrum_norm_chung_everhart *param = 
    (struct gkyl_bc_emission_spectrum_norm_chung_everhart *) norm_param;
  double mass = param->mass;
  double charge = param->charge;
  double phi = param->phi;
  
  out[0] = 6.0*effective_delta*flux[0]*phi*phi*mass/fabs(charge);
}

// Gaussian normalization factor
GKYL_CU_D
static void
gaussian_norm(double *out, const double *flux, void *norm_param, double effective_delta)
{
  struct gkyl_bc_emission_spectrum_norm_gaussian *param = 
    (struct gkyl_bc_emission_spectrum_norm_gaussian *) norm_param;
  double mass = param->mass;
  double charge = param->charge;
  double E_0 = param->E_0;
  double tau = param->tau;

  out[0] = effective_delta*flux[0]*mass/(sqrt(2.0*M_PI)*E_0*tau*exp(tau*tau/2.0)*fabs(charge));
}

// Maxwellian normalization factor
GKYL_CU_D
static void
maxwellian_norm(double *out, const double *flux, void *norm_param, double effective_delta)
{
  struct gkyl_bc_emission_spectrum_norm_maxwellian *param = 
    (struct gkyl_bc_emission_spectrum_norm_maxwellian *) norm_param;
  double vt = param->vt;
  int vdim = param->vdim;
  
  out[0] = effective_delta*flux[0]/(pow(2.0*M_PI, (vdim - 1)/2.0)*pow(vt, vdim + 1));
}

// Furman-Pivi SEY calculation
GKYL_CU_D
static void
furman_pivi_yield(double *out, int cdim, int vdim, double xc[GKYL_MAX_DIM],
  void *yield_param)
// Electron impact model adapted from https://link.aps.org/doi/10.1103/PhysRevSTAB.5.124404
{
  struct gkyl_bc_emission_spectrum_yield_furman_pivi *param = 
    (struct gkyl_bc_emission_spectrum_yield_furman_pivi *) yield_param;
  double mass = param->mass;
  double charge = param->charge;
  double deltahat_ts = param->deltahat_ts;
  double Ehat_ts = param->Ehat_ts;
  double t1 = param->t1;
  double t2 = param->t2;
  double t3 = param->t3;
  double t4 = param->t4;
  double s = param->s;
  
  double E = 0.0;
  double mu = 1.0; // currently hardcoded to normal, will add angular dependence later
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xc[cdim+d]*xc[cdim+d]/fabs(charge);
  }
  double deltahat = deltahat_ts*(1 + t1*(1 - pow(mu, t2)));
  double Ehat = Ehat_ts*(1 + t3*(1 - pow(mu, t4)));
  double x = E/Ehat;
  out[0] = deltahat*s*x/(s - 1 + pow(x, s));
}

// Schou SEY calculation
GKYL_CU_D
static void
schou_yield(double *out, int cdim, int vdim, double xc[GKYL_MAX_DIM],
  void *yield_param)
// Ion impact model adapted from https://doi.org/10.1103/PhysRevB.22.2141
{ // No angular dependence atm. Will have to add later
  struct gkyl_bc_emission_spectrum_yield_schou *param =
    (struct gkyl_bc_emission_spectrum_yield_schou *) yield_param;
  double mass = param->mass;   
  double charge = param->charge;
  double int_wall = param->int_wall;
  double A2 = param->a2;  // Note: Starts at 2 to match notation from source: https://doi.org/10.1093/jicru_os25.2.18
  double A3 = param->a3;
  double A4 = param->a4;
  double A5 = param->a5;
  double nw = param->nw;   // Number density of wall material in m^-3

  double E = 0.0;   
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xc[cdim+d]*xc[cdim+d]/fabs(charge)/1000;  // Calculate energy in keV
  }

  double Es = E*1.66053906660e-27/1.67262192369e-27;  // Scale energy by ratio of atomic mass unit to proton mass
  double eps_low = A2*pow(Es,0.45);
  double eps_high = 0.0;
  if (Es != 0) { // Divide by zero error catching. If Es is not 0, then do the normal calculation.
    eps_high = (A3/Es)*log(1+A4/Es+A5*Es);
  }
  double eps = eps_low*eps_high/(eps_low+eps_high);

  out[0] =  eps*nw*fabs(charge)*int_wall/1.0e19;
}

// Fixed constant SEY
GKYL_CU_D
static void
constant_yield(double *out, int cdim, int vdim, double xc[GKYL_MAX_DIM], void *yield_param)
{
  struct gkyl_bc_emission_spectrum_yield_constant *param =
    (struct gkyl_bc_emission_spectrum_yield_constant *) yield_param;
  double delta = param->delta;

  out[0] = delta;
}

void
gkyl_bc_emission_spectrum_choose_func_cu(enum gkyl_bc_emission_spectrum_norm_type norm_type,
  enum gkyl_bc_emission_spectrum_yield_type yield_type, struct gkyl_bc_emission_spectrum_funcs *funcs);

void
gkyl_bc_emission_spectrum_choose_norm_cu(enum gkyl_bc_emission_spectrum_norm_type norm_type,
  struct gkyl_bc_emission_spectrum_funcs *funcs, void *norm_param);

void
gkyl_bc_emission_spectrum_choose_yield_cu(enum gkyl_bc_emission_spectrum_yield_type yield_type, struct gkyl_bc_emission_spectrum_funcs *funcs, void *yield_param);

GKYL_CU_D
static emission_spectrum_spec_func_t
bc_emission_spectrum_choose_spec_func(enum gkyl_bc_emission_spectrum_norm_type norm_type)
{
  switch (norm_type) {
    case GKYL_SEE_CHUNG_EVERHART:
      return chung_everhart_spec;
    case GKYL_SEE_GAUSSIAN:
      return gaussian_spec;
    case GKYL_SEE_MAXWELLIAN:
      return maxwellian_spec;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static emission_spectrum_norm_func_t
bc_emission_spectrum_choose_norm_func(enum gkyl_bc_emission_spectrum_norm_type norm_type)
{
  switch (norm_type) {
    case GKYL_SEE_CHUNG_EVERHART:
      return chung_everhart_norm;
    case GKYL_SEE_GAUSSIAN:
      return gaussian_norm;
    case GKYL_SEE_MAXWELLIAN:
      return maxwellian_norm;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static emission_spectrum_yield_func_t
bc_emission_spectrum_choose_yield_func(enum gkyl_bc_emission_spectrum_yield_type yield_type)
{
  switch (yield_type) {
    case GKYL_SEE_FURMAN_PIVI:
      return furman_pivi_yield;
    case GKYL_SEE_SCHOU:
      return schou_yield;
    case GKYL_SEE_CONSTANT:
      return constant_yield;
    default:
      assert(false);
      break;
  }
}

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set up function to apply boundary conditions.
 *
 * @param up BC updater
 * @param f_skin Skin cell distribution
 * @param f_proj Projected spectrum distribution
 * @param f_buff Distribution buffer array
 * @param weight Weighting coefficients
 * @param k Normalization factor
 * @param flux Flux into boundary
 * @param grid Domain grid
 * @param gamma SE yield values on incoming ghost space
 * @param skin_r Incoming skin space range
 * @param ghost_r Incoming ghost space range
 * @param conf_r Configuration space range
 * @param buff_r Buffer array range
 */
void gkyl_bc_emission_spectrum_advance_cu(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_range *impact_buff_r, struct gkyl_range *impact_cbuff_r,
  struct gkyl_range *emit_buff_r, struct gkyl_array *bflux, struct gkyl_array *f_emit,
  struct gkyl_array *yield, struct gkyl_array *spectrum, struct gkyl_array *weight,
  struct gkyl_array *flux, struct gkyl_array *k);

/**
 * CUDA device function to set up function to calculate SEY
 *
 * @param up BC updater
 * @param grid Domain grid
 * @param gamma SE yield values on incoming ghost space
 * @param ghost_r Incoming ghost space range
 */
void gkyl_bc_emission_spectrum_sey_calc_cu(const struct gkyl_bc_emission_spectrum *up, struct gkyl_array *yield, struct gkyl_rect_grid *grid, const struct gkyl_range *gamma_r);

#endif
