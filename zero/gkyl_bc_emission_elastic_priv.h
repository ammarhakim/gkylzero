#pragma once

// Private header for bc_emission_spectrum updater, not for direct use in user code.

#include <gkyl_bc_emission_elastic.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_updater_moment.h>
#include <math.h>
#include <assert.h>

typedef void (*emission_elastic_yield_func_t)(double t, const double *xn, double *fout, void *ctx);

struct gkyl_bc_emission_elastic_funcs {
  emission_elastic_yield_func_t yield;
  void *elastic_param;
};

// Primary struct for the updater
struct gkyl_bc_emission_elastic {
  int dir, cdim, vdim;
  enum gkyl_edge_loc edge;
  struct gkyl_array_copy_func *reflect_func;
  struct gkyl_bc_emission_elastic_funcs *funcs;
  struct gkyl_bc_emission_elastic_funcs *funcs_cu;
  bool use_gpu;
};

struct bc_elastic_ctx {
  int dir; // direction for BCs.
  int cdim; // config-space dimensions.
  int ncomp; // number of components within a cell.
  const struct gkyl_basis *basis; // basis function.
};

GKYL_CU_D
static void
reflection(size_t nc, double *out, const double *inp, void *ctx)
{
  struct bc_elastic_ctx *bc_ctx = (struct bc_elastic_ctx *) ctx;
  int dir = bc_ctx->dir, cdim = bc_ctx->cdim;
  
  bc_ctx->basis->flip_odd_sign(dir, inp, out);
  bc_ctx->basis->flip_odd_sign(dir+cdim, out, out);
}

// Furman-Pivi SEY calculation
GKYL_CU_D
static void
furman_pivi_yield(double t, const double *xn, double *fout, void *ctx)
// Electron impact model adapted from https://link.aps.org/doi/10.1103/PhysRevSTAB.5.124404
{
  struct gkyl_bc_emission_elastic *bc_ctx = 
    (struct gkyl_bc_emission_elastic *) ctx;
  struct gkyl_bc_emission_elastic_furman_pivi *param =
    (struct gkyl_bc_emission_elastic_furman_pivi *) bc_ctx->funcs->elastic_param;
  int cdim = bc_ctx->cdim;
  int vdim = bc_ctx->vdim;
  double mass = param->mass;
  double charge = param->charge;
  double P1_inf = param->P1_inf;
  double P1_hat = param->P1_hat;
  double E_hat = param->E_hat;
  double W = param->W;
  double p = param->p;

  double E = 0.0;
  double mu = 1.0; // currently hardcoded to normal, will add angular dependence later
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xn[cdim+d]*xn[cdim+d]/fabs(charge); // Calculate energy in eV
  }
  fout[0] = P1_inf + (P1_hat - P1_inf)*exp(pow(-fabs(E - E_hat)/W, p)/p);
}

// Cazaux backscattering
GKYL_CU_D
static void
cazaux_yield(double t, const double *xn, double *fout, void *ctx)
// Low-energy backscattering model adapted from https://doi.org/10.1063/1.3691956
{ 
  struct gkyl_bc_emission_elastic *bc_ctx = 
    (struct gkyl_bc_emission_elastic *) ctx;
  struct gkyl_bc_emission_elastic_cazaux *param = 
    (struct gkyl_bc_emission_elastic_cazaux *)bc_ctx->funcs->elastic_param;
  int cdim = bc_ctx->cdim;
  int vdim = bc_ctx->vdim;
  double mass = param->mass;   
  double charge = param->charge;
  double E_f = param->E_f;
  double phi = param->phi;

  double E = 0.0;   
  for (int d=0; d<vdim; d++) {
    E += 0.5*mass*xn[cdim+d]*xn[cdim+d]/fabs(charge);  // Calculate energy in eV
  }  
  double E_s = E + E_f + phi;
  double G = 1 + (E_s - E)/E;

  fout[0] = pow(1 - sqrt(G), 2)/pow(1 + sqrt(G), 2);
}

// Fixed constant reflection
GKYL_CU_D
static void
constant_yield(double t, const double *xn, double *fout, void *ctx)
{
  struct gkyl_bc_emission_elastic *bc_ctx = 
    (struct gkyl_bc_emission_elastic *) ctx;
  struct gkyl_bc_emission_elastic_constant *param = 
    (struct gkyl_bc_emission_elastic_constant *) bc_ctx->funcs->elastic_param;
  double delta = param->delta;

  fout[0] = delta;
}

struct gkyl_array_copy_func*
gkyl_bc_emission_elastic_create_arr_copy_func_cu(int dir, int cdim, const struct gkyl_basis *basis,
  int ncomp);

void
gkyl_bc_emission_elastic_choose_elastic_cu(enum gkyl_bc_emission_elastic_type elastic_type,
  struct gkyl_bc_emission_elastic_funcs *funcs, void *elastic_param);

GKYL_CU_D
static emission_elastic_yield_func_t
bc_emission_elastic_choose_yield_func(enum gkyl_bc_emission_elastic_type yield_type)
{
  switch (yield_type) {
    case GKYL_BS_FURMAN_PIVI:
      return furman_pivi_yield;
    case GKYL_BS_CAZAUX:
      return cazaux_yield;
    case GKYL_BS_CONSTANT:
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
void gkyl_bc_emission_elastic_advance_cu(const struct gkyl_bc_emission_elastic *up,
  struct gkyl_range *emit_skin_r, struct gkyl_array *buff_arr, struct gkyl_array *f_skin,
  struct gkyl_array *f_emit, struct gkyl_array *elastic_yield, struct gkyl_basis *basis);
#endif
