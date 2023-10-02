#pragma once

// Private header for bc_emission_spectrum updater, not for direct use in user code.

#include <gkyl_bc_emission_spectrum.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_updater_moment.h>
#include <math.h>
#include <assert.h>

typedef void (*emission_spectrum_func_t)(const double *inp, int cdim, int dir, enum gkyl_edge_loc edge, double xc[GKYL_MAX_DIM], const double *gain, double *weight);
typedef double (*emission_spectrum_norm_func_t)(double *out, const double *flux, double *param, double effective_gamma);
typedef void (*emission_spectrum_gamma_func_t)(double *out, int cdim, int vdim, double xc[GKYL_MAX_DIM], double *param);

struct gkyl_bc_emission_spectrum_funcs {
  emission_spectrum_func_t func;
  emission_spectrum_norm_func_t norm;
  emission_spectrum_gamma_func_t gamma;
};

// Primary struct for the updater
struct gkyl_bc_emission_spectrum {
  int dir, cdim, vdim;
  enum gkyl_edge_loc edge;
  double *bc_param;
  double *bc_param_cu;
  double *sey_param;
  double *sey_param_cu;
  struct gkyl_bc_emission_spectrum_funcs *funcs;
  struct gkyl_bc_emission_spectrum_funcs *funcs_cu;
  struct gkyl_array *gamma;
  bool use_gpu;
};

// Function to calculate the weighted mean of the SE yield
GKYL_CU_D
static void
bc_weighted_gamma(const double *inp, int cdim, int dir, enum gkyl_edge_loc edge, double xc[GKYL_MAX_DIM], const double *gain, double *weight)
{
  if ((edge == GKYL_LOWER_EDGE && xc[cdim+dir] < 0) || (edge == GKYL_UPPER_EDGE && xc[cdim+dir] > 0)) { 
    weight[0] += fabs(xc[cdim+dir])*inp[0]*gain[0];
    weight[1] += fabs(xc[cdim+dir])*inp[0];
  }
}

// Chung-Everhart normalization factor
GKYL_CU_D
static double
chung_everhart_norm(double *out, const double *flux, double *param, double effective_gamma)
  
{
  double phi = param[2];
  double mass = param[0];
  double charge = param[1];
  
  out[0] = 6.0*effective_gamma*flux[0]*phi*phi*mass/fabs(charge);
}

// Gaussian normalization factor
GKYL_CU_D
static double
gaussian_norm(double *out, const double *flux, double *param, double effective_gamma)
{
  double E_0 = param[2];
  double tau = param[3];
  double mass = param[0];
  double charge = param[1];
  
  out[0] = effective_gamma*flux[0]*mass/(sqrt(2.0*M_PI)*E_0*tau*exp(tau*tau/2.0)*fabs(charge));
}

// Maxwellian normalization factor
GKYL_CU_D
static double
maxwellian_norm(double *out, const double *flux, double *param, double effective_gamma)
{
  double vt = param[2];
  
  out[0] = effective_gamma*flux[0]/(vt*vt);
}

// Furman-Pivi SEY calculation
GKYL_CU_D
static void
furman_pivi_gamma(double *out, int cdim, int vdim, double xc[GKYL_MAX_DIM], double *param) // Electron impact model adapted from https://link.aps.org/doi/10.1103/PhysRevSTAB.5.124404
{
  double mass = param[0];
  double charge = param[1];
  double gammahat_ts = param[2];
  double Ehat_ts = param[3];
  double t1 = param[4];
  double t2 = param[5];
  double t3 = param[6];
  double t4 = param[7];
  double s = param[8];

  double E = 0.0;
  double mu = 1.0; // currently hardcoded to normal, will add angular dependence later
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xc[cdim+d]*xc[cdim+d]/fabs(charge);
  }
  double gammahat = gammahat_ts*(1 + t1*(1 - pow(mu, t2)));
  double Ehat = Ehat_ts*(1 + t3*(1 - pow(mu, t4)));
  double x = E/Ehat;

  out[0] = gammahat*s*x/(s - 1 + pow(x, s));
}

// Schou SEY calculation
GKYL_CU_D
static void
schou_gamma(double *out, int cdim, int vdim, double xc[GKYL_MAX_DIM], double *param) // Ion impact model adapted from https://doi.org/10.1103/PhysRevB.22.2141
{ // No angular dependence atm. Will have to add later

  double mass = param[0];   
  double charge = param[1];
  double intWall = param[2];
  double A2 = param[3];  // Note: Starts at 2 to match notation from source: https://doi.org/10.1093/jicru_os25.2.18
  double A3 = param[4];
  double A4 = param[5];
  double A5 = param[6];
  double nw = param[7];   // Number density of wall material in m^-3

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

  out[0] =  eps*nw*fabs(charge)*intWall/1.0e19;
}

// Fixed constant SEY
GKYL_CU_D
static void
constant_gamma(double *out, int cdim, int vdim, double xc[GKYL_MAX_DIM], double *param)
{
  double gain = param[2];

  out[0] = gain;
}

void
gkyl_bc_emission_spectrum_choose_func_cu(enum gkyl_bc_emission_spectrum_type bctype,
  enum gkyl_bc_emission_spectrum_gamma_type gammatype, struct gkyl_bc_emission_spectrum_funcs *funcs);

GKYL_CU_D
static emission_spectrum_norm_func_t
bc_emission_spectrum_choose_norm_func(enum gkyl_bc_emission_spectrum_type bctype)
{
  switch (bctype) {
    case GKYL_BC_CHUNG_EVERHART:
      return chung_everhart_norm;
    case GKYL_BC_GAUSSIAN:
      return gaussian_norm;
    case GKYL_BC_MAXWELLIAN:
      return maxwellian_norm;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static emission_spectrum_gamma_func_t
bc_emission_spectrum_choose_gamma_func(enum gkyl_bc_emission_spectrum_gamma_type gammatype)
{
  switch (gammatype) {
    case GKYL_BC_FURMAN_PIVI:
      return furman_pivi_gamma;
    case GKYL_BC_SCHOU:
      return schou_gamma;
    case GKYL_BC_CONSTANT:
      return constant_gamma;
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
  const struct gkyl_array *f_skin, const struct gkyl_array *f_proj, struct gkyl_array *f_buff,
  struct gkyl_array *weight, struct gkyl_array *k,
  const struct gkyl_array *flux, struct gkyl_rect_grid *grid, struct gkyl_array *gamma,
  const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r, const struct gkyl_range *conf_r,
  const struct gkyl_range *buff_r);

/**
 * CUDA device function to set up function to calculate SEY
 *
 * @param up BC updater
 * @param grid Domain grid
 * @param gamma SE yield values on incoming ghost space
 * @param ghost_r Incoming ghost space range
 */
void gkyl_bc_emission_spectrum_sey_calc_cu(const struct gkyl_bc_emission_spectrum *up,
  struct gkyl_array *gamma, struct gkyl_rect_grid *grid, const struct gkyl_range *ghost_r);

#endif
