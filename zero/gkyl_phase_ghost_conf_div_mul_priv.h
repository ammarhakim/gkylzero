#pragma once

// Private header for phase_ghost_conf_div_mul updater, not for direct use in user code.

#include <gkyl_phase_ghost_conf_div_mul.h>
#include <gkyl_mom_gyrokinetic_kernels.h>
#include <gkyl_mom_gyrokinetic_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_util.h>
#include <assert.h>

struct gkyl_phase_ghost_conf_div_mul_kernels {
  inv_op_t conf_inv_op; // Conf-space weak inversion (1/A) kernel (p=1 only).
  mul_op_t conf_mul_op; // Conf-space weak multiplication kernel.
  mul_op_t conf_phase_mul_op; // Conf-phase weak multiplication kernel.
};

// Primary struct in this updater.
struct gkyl_phase_ghost_conf_div_mul {
  int dir; // Direction perpendicular to the boundary.
  enum gkyl_edge_loc edge; // Boundary edge (lower/upper).
  const struct gkyl_basis *conf_basis; // Configuration space basis object.
  const struct gkyl_range *conf_skin_r; // Conf-space skin range.
  const struct gkyl_range *conf_ghost_r; // Conf-space ghost range.
  const struct gkyl_range *phase_ghost_r; // Phase-space ghost range.
  bool use_gpu; // Whether to run on the GPU or not.
  struct gkyl_phase_ghost_conf_div_mul_kernels *kernels;
};

#ifdef GKYL_HAVE_CUDA
// Declaration of cuda device functions.
void
phase_ghost_conf_div_mul_kernel_cu(struct gkyl_phase_ghost_conf_div_mul_kernels *kernels,
  struct gkyl_basis *cbasis, struct gkyl_basis *pbasis);

void
gkyl_phase_ghost_conf_div_mul_advance_cu(const struct gkyl_phase_ghost_conf_div_mul *up,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_range *conf_skin_r, const struct gkyl_range *conf_ghost_r,
  const struct gkyl_range *phase_ghost_r, struct gkyl_array *GKYL_RESTRICT jac, struct gkyl_array *GKYL_RESTRICT jf);
#endif

GKYL_CU_D
static void phase_ghost_conf_div_mul_choose_kernel(struct gkyl_phase_ghost_conf_div_mul_kernels *kernels,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    pos_shift_gk_choose_shift_kernel_cu(kernels, cbasis, pbasis);
    return;
  }
#endif

  enum gkyl_basis_type cbasis_type = cbasis->b_type, pbasis_type = pbasis->b_type;
  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int poly_order = pbasis->poly_order;

  switch (pbasis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->conf_phase_mul_op = choose_mul_conf_phase_kern(pbasis_type, cdim, pdim-cdim, poly_order);
      break;
    default:
      assert(false);
      break;
  }

  switch (cbasis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->conf_inv_op = choose_ser_inv_kern(cdim, poly_order);
      kernels->conf_mul_op = choose_ser_mul_kern(cdim, poly_order);
      break;
    default:
      assert(false);
      break;
  }

}
