#pragma once

// Private header for dg_interpolate updater, not for direct use in user code.

#include <gkyl_dg_interpolate.h>
#include <gkyl_dg_interpolate_kernels.h>
#include <gkyl_util.h>
#include <assert.h>

struct dg_interp_stencils {
  int sz[GKYL_MAX_DIM]; // Size of stencil in each direction.
  struct gkyl_range offsets; // Index offset relative to the lower-left
                             // corner cell of the stencil.
};

// Function that translates a grid index into a stencil index.
typedef int (*dg_interp_grid2stencilIdx_t)(int dir, int idx,
  int num_cells, double dx_rat, int stencilIn);

// Function pointer type for sheath reflection kernels.
typedef void (*dg_interp_t)(const double *wDo, const double *wTar,
  const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

typedef struct {dg_interp_t kernels[2];} dg_interp_kern_list_gk;  // For use in kernel tables.

// Primary struct in this updater.
struct gkyl_dg_interpolate {
  int ndim; // Phase space dimensionality of the fields.
  bool use_gpu; // Whether to use the GPU.
  struct gkyl_rect_grid grid_do; // Donor grid.
  struct gkyl_rect_grid grid_tar; // Target grid.
  double dxRat[GKYL_MAX_DIM]; // Ratio of donor to target cell lengths in each direction.

  struct dg_interp_stencils *stencils; // Size of stencil in region.

  dg_interp_t interp;  // Kernel that performs the interpolation.
  dg_interp_grid2stencilIdx_t grid2stencil[GKYL_MAX_DIM]; // Translate grid to stencil index.
};

// Serendipity  kernels.
GKYL_CU_D
static const dg_interp_kern_list_gk dg_interp_kern_list_gk_ser[] = {
  { dg_interpolate_gyrokinetic_1x1v_ser_p1, NULL },
  { NULL, NULL },
  { NULL, NULL },
  { NULL, NULL },
};

#ifdef GKYL_HAVE_CUDA
// Declaration of cuda device functions.
void dg_interp_gk_choose_shift_kernel_cu(struct gkyl_dg_interpolate_kernels *kernels,
  int cdim, struct gkyl_basis pbasis);

void gkyl_dg_interpolate_advance_cu(gkyl_dg_interpolate* up,
  const struct gkyl_range *phase_rng_do, const struct gkyl_range *phase_rng_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar);
#endif

GKYL_CU_D
static dg_interp_t
dg_interp_choose_gk_interp_kernel(const struct gkyl_basis *pbasis)
{
  enum gkyl_basis_type basis_type = pbasis->b_type;
  int pdim = pbasis->ndim;
  int poly_order = pbasis->poly_order;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return dg_interp_kern_list_gk_ser[pdim-2].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }
}
