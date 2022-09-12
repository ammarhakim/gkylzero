#pragma once

// Private header for bc_basic updater, not for direct use in user code.

#include <gkyl_bc_basic.h>
#include <gkyl_bc_external.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <assert.h>

// Primary struct in this updater.
struct gkyl_bc_basic {
  int dir, cdim;
  enum gkyl_edge_loc edge;
  enum gkyl_bc_basic_type bctype;
  struct gkyl_range skin_r, ghost_r;
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
  double *external_field; // pointer to external field for some BCs.
  long out_idx;
  long in_idx;
};

GKYL_CU_D
static void
copy_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int nbasis = mc->basis->num_basis;
  for (int c=0; c<nbasis; ++c) out[c] = inp[c];
}

GKYL_CU_D
static void
species_absorb_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int nbasis = mc->basis->num_basis;
  for (int c=0; c<nbasis; ++c) out[c] = 0.0;
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

GKYL_CU_D
static void
species_gain_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir, cdim = mc->cdim;
  double *gain = mc->external_field;

  mc->basis->flip_odd_sign(dir, inp, out);
  out[0] = gain[0]*out[0];
  mc->basis->flip_odd_sign(dir+cdim, out, out);
}

GKYL_CU_D
static void
species_furman_pivi_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir, cdim = mc->cdim;
  double *reflect = mc->external_field;
  long out_idx = mc->out_idx;
  long in_idx = mc->in_idx;
  double coeffs[8];
  for(int a=0; a<8; a++){
    coeffs[a] = 0.;
    for(int b=0; b<8; b++){
      coeffs[a] += inp[b]*reflect[b + 8*(a + 8*((out_idx - 1) + 32*(in_idx - 1)))];
    }
  }
  // mc->basis->flip_odd_sign(dir, coeffs, coeffs);
  // mc->basis->flip_odd_sign(dir+cdim, coeffs, coeffs);
  for(int a=0; a<8; a++){
    out[a] += coeffs[a];
  }
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

// Perfect electrical conductor
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
