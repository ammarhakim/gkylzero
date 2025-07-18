// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gk_neut_fluid_prim_vars_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_wv_eqn.h>
#include <assert.h>

typedef void (*gk_nf_udrift_set_prob_t)(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs,
  const double *moms);

typedef void (*gk_nf_udrift_get_sol_t)(int count, struct gkyl_nmat *xsol,
  double* GKYL_RESTRICT out);

typedef void (*gk_nf_pressure_t)(double gas_gamma, const double *moms, const double *u,
  double* GKYL_RESTRICT out);

typedef void (*gk_nf_temp_set_prob_t)(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs,
  const double *moms, double gas_gamma, double mass);

typedef void (*gk_nf_temp_get_sol_t)(int count, struct gkyl_nmat *xsol,
  double* GKYL_RESTRICT out);

typedef void (*gk_nf_udrift_temp_set_prob_t)(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs,
  const double *moms, double gas_gamma, double mass);

typedef void (*gk_nf_udrift_temp_get_sol_t)(int count, struct gkyl_nmat *xsol,
  double* GKYL_RESTRICT out);

// For use in kernel tables.
typedef struct { gk_nf_udrift_set_prob_t kernels[4]; } gkyl_gk_nf_prim_vars_udrift_set_prob_kern_list;
typedef struct { gk_nf_udrift_get_sol_t kernels[4]; } gkyl_gk_nf_prim_vars_udrift_get_sol_kern_list;
typedef struct { gk_nf_pressure_t kernels[4]; } gkyl_gk_nf_prim_vars_pressure_kern_list;
typedef struct { gk_nf_temp_set_prob_t kernels[4]; } gkyl_gk_nf_prim_vars_temp_set_prob_kern_list;
typedef struct { gk_nf_temp_get_sol_t kernels[4]; } gkyl_gk_nf_prim_vars_temp_get_sol_kern_list;
typedef struct { gk_nf_udrift_temp_set_prob_t kernels[4]; } gkyl_gk_nf_prim_vars_udrift_temp_set_prob_kern_list;
typedef struct { gk_nf_udrift_temp_get_sol_t kernels[4]; } gkyl_gk_nf_prim_vars_udrift_temp_get_sol_kern_list;

struct gkyl_gk_neut_fluid_prim_vars {
  double gas_gamma; // Adiabatic index.
  double mass; // Species mass.

  int cdim; // Configuration space dimensionality.
  int poly_order; // Polynomial order.
  int num_basis; // Number of basis functions.
  struct gkyl_range mem_range; // Configuration space range for linear solve.

  struct gkyl_nmat *As, *xs; // matrices for LHS and RHS
  gkyl_nmat_mem *mem; // memory for use in batched linear solve
  int udrift_ncomp; // number of components in the linear solve (6 variables being solved for)

  gk_nf_udrift_set_prob_t udrift_set_prob_ker; // Kernel for setting matrices for linear solve.
  gk_nf_udrift_get_sol_t udrift_get_sol_ker; // Kernel for copying solution to output.
  gk_nf_pressure_t pressure_ker; // Kernel for computing pressure.
  gk_nf_temp_set_prob_t temp_set_prob_ker; // Kernel for setting matrices for linear solve.
  gk_nf_temp_get_sol_t temp_get_sol_ker; // Kernel for copying solution to output.
  gk_nf_udrift_temp_set_prob_t udrift_temp_set_prob_ker; // Kernel for setting matrices for linear solve.
  gk_nf_udrift_temp_get_sol_t udrift_temp_get_sol_ker; // Kernel for copying solution to output.

  uint32_t flags;
  struct gkyl_gk_neut_fluid_prim_vars *on_dev; // pointer to itself or device data

  // Method chosen at runtime.
  void (*advance_func)(struct gkyl_gk_neut_fluid_prim_vars *up,
    const struct gkyl_array* moms, struct gkyl_array *out, int out_coff);
};

// Set matrices for computing fluid flow velocity (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_nf_prim_vars_udrift_set_prob_kern_list ser_gk_nf_prim_vars_udrift_set_prob_kernels[] = {
  { gk_neut_fluid_prim_vars_udrift_set_prob_1x_ser_p1, gk_neut_fluid_prim_vars_udrift_set_prob_1x_ser_p2, NULL },
  { gk_neut_fluid_prim_vars_udrift_set_prob_2x_ser_p1, NULL, NULL },
  { gk_neut_fluid_prim_vars_udrift_set_prob_3x_ser_p1, NULL, NULL },
};

// Copy solution for fluid flow velocity (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_nf_prim_vars_udrift_get_sol_kern_list ser_gk_nf_prim_vars_udrift_get_sol_kernels[] = {
  { gk_neut_fluid_prim_vars_udrift_get_sol_1x_ser_p1, gk_neut_fluid_prim_vars_udrift_get_sol_1x_ser_p2, NULL },
  { gk_neut_fluid_prim_vars_udrift_get_sol_2x_ser_p1, NULL, NULL },
  { gk_neut_fluid_prim_vars_udrift_get_sol_3x_ser_p1, NULL, NULL },
};

// Scalar pressure p = (gas_gamma - 1)*(E - 1/2 rho u^2) (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_nf_prim_vars_pressure_kern_list ser_gk_nf_prim_vars_pressure_kernels[] = {
  { gk_neut_fluid_prim_vars_pressure_1x_ser_p1, gk_neut_fluid_prim_vars_pressure_1x_ser_p2, NULL },
  { gk_neut_fluid_prim_vars_pressure_2x_ser_p1, NULL, NULL },
  { gk_neut_fluid_prim_vars_pressure_3x_ser_p1, NULL, NULL },
};

// Set matrices for computing temperature (Serendipity kernels).
GKYL_CU_D
static const gkyl_gk_nf_prim_vars_temp_set_prob_kern_list ser_gk_nf_prim_vars_temp_set_prob_kernels[] = {
  { gk_neut_fluid_prim_vars_temp_set_prob_1x_ser_p1, gk_neut_fluid_prim_vars_temp_set_prob_1x_ser_p2, NULL },
  { gk_neut_fluid_prim_vars_temp_set_prob_2x_ser_p1, NULL, NULL },
  { gk_neut_fluid_prim_vars_temp_set_prob_3x_ser_p1, NULL, NULL },
};

// Copy solution for temperature. (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_nf_prim_vars_temp_get_sol_kern_list ser_gk_nf_prim_vars_temp_get_sol_kernels[] = {
  { gk_neut_fluid_prim_vars_temp_get_sol_1x_ser_p1, gk_neut_fluid_prim_vars_temp_get_sol_1x_ser_p2, NULL },
  { gk_neut_fluid_prim_vars_temp_get_sol_2x_ser_p1, NULL, NULL },
  { gk_neut_fluid_prim_vars_temp_get_sol_3x_ser_p1, NULL, NULL },
};

// Set matrices for computing temperature (Serendipity kernels).
GKYL_CU_D
static const gkyl_gk_nf_prim_vars_udrift_temp_set_prob_kern_list ser_gk_nf_prim_vars_udrift_temp_set_prob_kernels[] = {
  { gk_neut_fluid_prim_vars_udrift_temp_set_prob_1x_ser_p1, gk_neut_fluid_prim_vars_udrift_temp_set_prob_1x_ser_p2, NULL },
  { gk_neut_fluid_prim_vars_udrift_temp_set_prob_2x_ser_p1, NULL, NULL },
  { gk_neut_fluid_prim_vars_udrift_temp_set_prob_3x_ser_p1, NULL, NULL },
};

// Copy solution for temperature. (Serendipity kernels)
GKYL_CU_D
static const gkyl_gk_nf_prim_vars_udrift_temp_get_sol_kern_list ser_gk_nf_prim_vars_udrift_temp_get_sol_kernels[] = {
  { gk_neut_fluid_prim_vars_udrift_temp_get_sol_1x_ser_p1, gk_neut_fluid_prim_vars_udrift_temp_get_sol_1x_ser_p2, NULL },
  { gk_neut_fluid_prim_vars_udrift_temp_get_sol_2x_ser_p1, NULL, NULL },
  { gk_neut_fluid_prim_vars_udrift_temp_get_sol_3x_ser_p1, NULL, NULL },
};

GKYL_CU_D
static gk_nf_udrift_set_prob_t
choose_udrift_set_prob_ker(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_gk_nf_prim_vars_udrift_set_prob_kernels[cdim-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static gk_nf_udrift_get_sol_t
choose_udrift_get_sol_ker(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_gk_nf_prim_vars_udrift_get_sol_kernels[cdim-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static gk_nf_pressure_t
choose_pressure_ker(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_gk_nf_prim_vars_pressure_kernels[cdim-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static gk_nf_temp_set_prob_t
choose_temp_set_prob_ker(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_gk_nf_prim_vars_temp_set_prob_kernels[cdim-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static gk_nf_temp_get_sol_t
choose_temp_get_sol_ker(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_gk_nf_prim_vars_temp_get_sol_kernels[cdim-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static gk_nf_udrift_temp_set_prob_t
choose_udrift_temp_set_prob_ker(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_gk_nf_prim_vars_udrift_temp_set_prob_kernels[cdim-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}

GKYL_CU_D
static gk_nf_udrift_temp_get_sol_t
choose_udrift_temp_get_sol_ker(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_gk_nf_prim_vars_udrift_temp_get_sol_kernels[cdim-1].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}
