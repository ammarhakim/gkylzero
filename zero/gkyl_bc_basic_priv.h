#pragma once

// Private header for bc_basic updater, not for direct use in user code.

#include <gkyl_bc_basic.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <assert.h>

// Primary struct in this updater.
struct gkyl_bc_basic {
  int dir, cdim;
  enum gkyl_edge_loc edge;
  enum gkyl_bc_basic_type bctype;
  const struct gkyl_range *skin_r, *ghost_r;
  struct gkyl_array_copy_func *array_copy_func;
  bool use_gpu;
};

#ifdef GKYL_HAVE_CUDA

/**
 * CUDA device function to set up function to apply boundary conditions.

 * @param dir Direction in which to apply BC .
 * @param cdim Number of configuration space dimensions.
 * @param bctype Type of BC .
 * @param basis Basis in which to expand coefficients in array we apply BC to.
 * @param num_comp Number of components (DOFs) within a cell.
 * @return Pointer to array_copy_func which can be passed to array_copy_fn methods.
 */
struct gkyl_array_copy_func* gkyl_bc_basic_create_arr_copy_func_cu(int dir, int cdim,
  enum gkyl_bc_basic_type bctype, const struct gkyl_basis *basis, int num_comp);

#endif

// context for use in BCs
struct dg_bc_ctx {
  int dir; // direction for BCs.
  int cdim; // config-space dimensions.
  int ncomp; // number of components within a cell.
  const struct gkyl_basis *basis; // basis function.
};

GKYL_CU_D
static void
copy_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int num_comp = mc->ncomp;
  for (int c=0; c<num_comp; ++c) out[c] = inp[c];
}

GKYL_CU_D
static void
species_absorb_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int num_comp = mc->ncomp;
  for (int c=0; c<num_comp; ++c) out[c] = 0.0;
}

GKYL_CU_D
static void
species_reflect_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir, cdim = mc->cdim;

  mc->basis->flip_odd_sign(dir, inp, out);
  mc->basis->flip_odd_sign(dir+cdim, out, out);
}

enum { M_EX, M_EY, M_EZ, M_BX, M_BY, M_BZ }; // components of EM field
GKYL_CU_D static const int m_flip_even[3][3] = { // zero tangent E and zero normal B
  {M_BX, M_EY, M_EZ},
  {M_BY, M_EX, M_EZ},
  {M_BZ, M_EX, M_EY},
};
GKYL_CU_D static const int m_flip_odd[3][3] = { // zero gradient
  { M_EX, M_BY, M_BZ },
  { M_EY, M_BX, M_BZ },
  { M_EZ, M_BX, M_BY },
};

// Maxwell's perfect electrical conductor (zero normal B and zero tangent E)
GKYL_CU_D
static void
maxwell_pec_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir;
  int nbasis = mc->basis->num_basis;

  const int *feven = m_flip_even[dir];
  const int *fodd = m_flip_odd[dir];

  for (int i=0; i<3; ++i) {
    int eloc = nbasis*feven[i], oloc = nbasis*fodd[i];
    mc->basis->flip_even_sign(dir, &inp[eloc], &out[eloc]);
    mc->basis->flip_odd_sign(dir, &inp[oloc], &out[oloc]);
  }
  // correction potentials
  int eloc = nbasis*6, oloc = nbasis*7;
  mc->basis->flip_even_sign(dir, &inp[eloc], &out[eloc]);
  mc->basis->flip_odd_sign(dir, &inp[oloc], &out[oloc]);
}

GKYL_CU_D static const int m_sym_flip_even[3][3] = { // zero tangent B and zero normal E
  { M_EX, M_BY, M_BZ },
  { M_EY, M_BX, M_BZ },
  { M_EZ, M_BX, M_BY },
};
GKYL_CU_D static const int m_sym_flip_odd[3][3] = { // zero gradient
  { M_BX, M_EY, M_EZ },
  { M_BY, M_EX, M_EZ },
  { M_BZ, M_EX, M_EY },
};

// Maxwell's symmetry BC (zero normal E and zero tangent B)
GKYL_CU_D
static void
maxwell_sym_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir;
  int nbasis = mc->basis->num_basis;

  const int *feven = m_sym_flip_even[dir];
  const int *fodd = m_sym_flip_odd[dir];

  for (int i=0; i<3; ++i) {
    int eloc = nbasis*feven[i], oloc = nbasis*fodd[i];
    mc->basis->flip_even_sign(dir, &inp[eloc], &out[eloc]);
    mc->basis->flip_odd_sign(dir, &inp[oloc], &out[oloc]);
  }
  // correction potentials
  int eloc = nbasis*6, oloc = nbasis*7;
  mc->basis->flip_even_sign(dir, &inp[eloc], &out[eloc]);
  mc->basis->flip_odd_sign(dir, &inp[oloc], &out[oloc]);
}

// Reservoir Maxwell's BCs for heat flux problem
// Based on Roberg-Clark et al. PRL 2018
// NOTE: ONLY WORKS WITH X BOUNDARY 
GKYL_CU_D
static void
maxwell_reservoir_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir;
  int nbasis = mc->basis->num_basis;

  // Zero gradient for Ex, Ez, Bx, Bz
  mc->basis->flip_odd_sign(dir, &inp[nbasis*0], &out[nbasis*0]);
  mc->basis->flip_odd_sign(dir, &inp[nbasis*2], &out[nbasis*2]);
  mc->basis->flip_odd_sign(dir, &inp[nbasis*3], &out[nbasis*3]);
  mc->basis->flip_odd_sign(dir, &inp[nbasis*5], &out[nbasis*5]);

  // Zero Ey and By
  mc->basis->flip_even_sign(dir, &inp[nbasis*1], &out[nbasis*1]);
  mc->basis->flip_even_sign(dir, &inp[nbasis*4], &out[nbasis*4]);

  // correction potentials
  int eloc = nbasis*6, oloc = nbasis*7;
  mc->basis->flip_even_sign(dir, &inp[eloc], &out[eloc]);
  mc->basis->flip_odd_sign(dir, &inp[oloc], &out[oloc]);
}

// Reflecting wall BCs for PKPM momentum
GKYL_CU_D
static void
pkpm_mom_reflect_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir;
  int nbasis = mc->basis->num_basis;

  // reflect normal component (zero normal) and zero gradient in other components
  for (int i=0; i<3; ++i) {
    int loc = nbasis*i;
    if (i == dir) {
      mc->basis->flip_even_sign(dir, &inp[loc], &out[loc]);
    }
    else {
      mc->basis->flip_odd_sign(dir, &inp[loc], &out[loc]);
    }
  }
}

// Line-tied wall BCs for PKPM momentum
GKYL_CU_D
static void
pkpm_mom_line_tied_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir;
  int nbasis = mc->basis->num_basis;

  // reflect tangential components (zero tangent) and zero gradient in normal
  for (int i=0; i<3; ++i) {
    int loc = nbasis*i;
    if (i == dir) {
      mc->basis->flip_odd_sign(dir, &inp[loc], &out[loc]);
    }
    else {
      mc->basis->flip_even_sign(dir, &inp[loc], &out[loc]);
    }
  }
}

// No-slip wall BCs for PKPM momentum
GKYL_CU_D
static void
pkpm_mom_no_slip_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir;
  int nbasis = mc->basis->num_basis;

  // zero normal and zero tangent
  for (int i=0; i<3; ++i) {
    int loc = nbasis*i;
    mc->basis->flip_even_sign(dir, &inp[loc], &out[loc]);
  }  
}

GKYL_CU_D
static void
pkpm_species_reflect_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir, cdim = mc->cdim;
  int nbasis = mc->basis->num_basis;

  int f_loc = 0*nbasis;
  int g_loc = 1*nbasis;

  // reflect F_0
  mc->basis->flip_odd_sign(dir, &inp[f_loc], &out[f_loc]);
  mc->basis->flip_odd_sign(dir+cdim, &out[f_loc], &out[f_loc]);
  // reflect G
  mc->basis->flip_odd_sign(dir, &inp[g_loc], &out[g_loc]);
  mc->basis->flip_odd_sign(dir+cdim, &out[g_loc], &out[g_loc]);
}

// Reflecting wall BCs for Euler equations
GKYL_CU_D
static void
euler_reflect_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir;
  int nbasis = mc->basis->num_basis;

  // Copy BCs for density and energy
  for (int c=0; c<nbasis; ++c) {
    out[c] = inp[c];
  }
  for (int c=4*nbasis; c<5*nbasis; ++c) {
    out[c] = inp[c];
  }

  // reflect normal component (zero normal) and zero gradient in other components
  for (int i=1; i<4; ++i) {
    int loc = nbasis*i;
    if (i == dir) {
      mc->basis->flip_even_sign(dir, &inp[loc], &out[loc]);
    }
    else {
      mc->basis->flip_odd_sign(dir, &inp[loc], &out[loc]);
    }
  }
}

// Line-tied wall BCs for Euler equations
GKYL_CU_D
static void
euler_line_tied_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir;
  int nbasis = mc->basis->num_basis;

  // Copy BCs for density and energy
  for (int c=0; c<nbasis; ++c) {
    out[c] = inp[c];
  }
  for (int c=4*nbasis; c<5*nbasis; ++c) {
    out[c] = inp[c];
  }

  // reflect tangential components (zero tangent) and zero gradient in normal
  for (int i=1; i<4; ++i) {
    int loc = nbasis*i;
    if (i == dir) {
      mc->basis->flip_odd_sign(dir, &inp[loc], &out[loc]);
    }
    else {
      mc->basis->flip_even_sign(dir, &inp[loc], &out[loc]);
    }
  }
}

// No-slip wall BCs for Euler equations
GKYL_CU_D
static void
euler_no_slip_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir;
  int nbasis = mc->basis->num_basis;

  // Copy BCs for density and energy
  for (int c=0; c<nbasis; ++c) {
    out[c] = inp[c];
  }
  for (int c=4*nbasis; c<5*nbasis; ++c) {
    out[c] = inp[c];
  }

  // zero normal and zero tangent
  for (int i=1; i<4; ++i) {
    int loc = nbasis*i;
    mc->basis->flip_even_sign(dir, &inp[loc], &out[loc]);
  }  
}
