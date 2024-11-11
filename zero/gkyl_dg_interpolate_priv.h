#pragma once

// Private header for dg_interpolate updater, not for direct use in user code.

#include <gkyl_dg_interpolate.h>
#include <gkyl_dg_interpolate_kernels.h>
#include <gkyl_util.h>
#include <assert.h>

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
  int *offset_upper; // Upper index of the offset in each direction for each stencil.

  dg_interp_t interp;  // Kernel that performs the interpolation.
  dg_interp_grid2stencilIdx_t grid2stencil[GKYL_MAX_DIM]; // Translate grid to stencil index.

  uint32_t flags;
  struct gkyl_dg_interpolate *on_dev; // pointer to itself or device data
};

// Serendipity  kernels.
GKYL_CU_D
static const dg_interp_kern_list_gk dg_interp_kern_list_gk_ser[] = {
  { dg_interpolate_gyrokinetic_1x1v_ser_p1, NULL },
  { dg_interpolate_gyrokinetic_1x2v_ser_p1, NULL },
  { NULL, NULL },
  { NULL, NULL },
};

#ifdef GKYL_HAVE_CUDA
// Declaration of cuda device functions.
void gkyl_dg_interpolate_new_cu(struct gkyl_dg_interpolate* up, struct gkyl_basis pbasis);

void gkyl_dg_interpolate_advance_cu(gkyl_dg_interpolate* up,
  const struct gkyl_range *phase_rng_do, const struct gkyl_range *phase_rng_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar);
#endif

GKYL_CU_D
static dg_interp_t
dg_interp_choose_gk_interp_kernel(struct gkyl_basis pbasis)
{
  enum gkyl_basis_type basis_type = pbasis.b_type;
  int pdim = pbasis.ndim;
  int poly_order = pbasis.poly_order;

  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return dg_interp_kern_list_gk_ser[pdim-2].kernels[poly_order-1];
      break;
    default:
      assert(false);
      break;
  }

  return 0;
}

static int dg_interp_prime_factors(int n, int *pfs, int pfs_size)
{
  // Find the prime factos of number `n`, and put them into `pfs`. We assume
  // there are fewer than `pfs_size` prime factors.
  int pf_count = 0;
  int c = 2;
  while (n > 1) {
    if (n % c == 0) {
      pfs[pf_count] = c;
      pf_count++;
      assert(pf_count < pfs_size);
      n /= c;
    }
    else
      c++;
  }
  return pf_count;
}

GKYL_CU_DH 
static int dg_interp_index_stencil_map_refine(int dir, int idx,
  int num_cells, double dx_rat, int stencilIn)
{
  // Given an index 'idx' to a cell in the coarse grid with 'num_cells'
  // cells, return the index of the refinement stencil needed, within
  // the table that holds stencils. Here 'dx_rat' is the ratio of
  // donor to target cell length in the 'dir' direction.
  double remDecL = (idx-1)*dx_rat-floor((idx-1)*dx_rat);
  double remDecU = ceil(idx*dx_rat)-idx*dx_rat;
  int stencilOut = stencilIn;
  if ((idx == 1) ||   // First cell.
      (remDecL == 0) ||    // Interior cell with a left-boundary-like stencil.
      ((remDecL <= 0.5) && (remDecU <= 0.5))) {
    stencilOut = 2*stencilIn + pow(3,dir-1) - 1;
  }
  else if ((idx == num_cells) || // Last cell.
          (remDecU == 0)) { // Interior cell with a right-boundary-like stencil.
    stencilOut = 2*stencilIn + pow(3,dir-1);
  }
  return stencilOut;
}

GKYL_CU_DH 
static int dg_interp_index_stencil_map_coarsen(int dir, int idx,
  int num_cells, double dx_rat, int stencilIn)
{
  // Given an index 'idx' to a cell in the fine grid with 'num_cells'
  // cells, return the index of the coarsening stencil needed, within
  // the table that holds stencils. Here 'dx_rat' is the ratio of
  // donor to target cell length in the 'dir' direction.
  double remDecL = (idx-1)*dx_rat-floor(idx*dx_rat);
  double remDecU = ceil(idx*dx_rat)-idx*dx_rat;
  int stencilOut = stencilIn;
  if ((idx == 1) || // First cell.
      (remDecL == 0) || // Interior cell with a left-boundary-like stencil.
      ((remDecL > 0) && (remDecU > 0))) {
     stencilOut = 2*stencilIn + pow(3,dir-1) - 1;
  }
  else if ((idx == num_cells) || // Last cell.
          (remDecU == 0)) { // Interior cell with a right-boundary-like stencil.
     stencilOut = 2*stencilIn + pow(3,dir-1);
  }
  return stencilOut;
}

