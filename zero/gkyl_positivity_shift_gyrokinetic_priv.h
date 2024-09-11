#pragma once

// Private header for translate_dim_gyrokinetic updater, not for direct use in user code.

#include <gkyl_translate_dim_gyrokinetic.h>
#include <gkyl_translate_dim_gyrokinetic_kernels.h>
#include <gkyl_util.h>
#include <assert.h>

// Function pointer type for sheath reflection kernels.
typedef bool (*translate_dim_t)(const double *fdo, double *ftar);

typedef struct { translate_dim_t kernels[3]; } trans_dim_gk_kern_list_shift;  // For use in kernel tables.

// Serendipity  kernels.
GKYL_CU_D
static const trans_dim_gk_kern_list_shift trans_dim_gk_kern_list_shift_ser[] = {
  { translate_dim_gyrokinetic_2x2v_ser_p1_from_1x2v_p1, NULL, NULL },
  { translate_dim_gyrokinetic_3x2v_ser_p1_from_1x2v_p1, NULL, NULL },
  { translate_dim_gyrokinetic_3x2v_ser_p1_from_2x2v_p1, NULL, NULL },
};

struct gkyl_translate_dim_gyrokinetic_kernels {
  translate_dim_t translate;  // Kernel that translate the DG coefficients.
};

// Primary struct in this updater.
struct gkyl_translate_dim_gyrokinetic {
  int cdim_do; // Configuration space dimension of donor field.
  int vdim_do; // Velocity space dimension of donor field.
  int cdim_tar; // Configuration space dimension of target field.
  int vdim_tar; // Velocity space dimension of target field.
  bool use_gpu;
  struct gkyl_translate_dim_gyrokinetic_kernels *kernels;
};

#ifdef GKYL_HAVE_CUDA
// Declaration of cuda device functions.
void trans_dim_gk_choose_shift_kernel_cu(struct gkyl_translate_dim_gyrokinetic_kernels *kernels,
  struct gkyl_basis pbasis);

void gkyl_translate_dim_gyrokinetic_advance_cu(gkyl_translate_dim_gyrokinetic* up,
  const struct gkyl_range *phase_rng_do, const struct gkyl_range *phase_rng_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar);
#endif

GKYL_CU_D
static void trans_dim_gk_choose_shift_kernel(struct gkyl_translate_dim_gyrokinetic_kernels *kernels,
  int cdim_do, struct gkyl_basis pbasis_do, int cdim_tar, struct gkyl_basis pbasis_tar, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    trans_dim_gk_choose_shift_kernel_cu(kernels, cdim_do, pbasis_do, cdim_tar, pbasis_tar);
    return;
  }
#endif

  enum gkyl_basis_type basis_type = pbasis_tar.b_type;
  int poly_order = pbasis_tar.poly_order;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->translate = trans_dim_gk_kern_list_shift_ser[cdim_tar+cdim_do-3].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}
