#pragma once

// Private header for dg_interpolate updater, not for direct use in user code.

#include <gkyl_dg_interpolate.h>
#include <gkyl_dg_interpolate_kernels.h>
#include <gkyl_util.h>
#include <assert.h>

// Function that translates a grid index into a stencil index.
typedef int (*dg_interp_grid2stencilIdx_t)(int idx,
  int num_cells, double dx_rat);

// Function pointer type for sheath reflection kernels.
typedef void (*dg_interp_t)(const double *wDo, const double *wTar,
  const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);

// For use in kernel tables.
typedef struct {dg_interp_t dirs[6];} dg_interp_kern_dir_list;
typedef struct {dg_interp_kern_dir_list list[2];} dg_interp_kern_p_list;
typedef struct {dg_interp_kern_p_list vdim[3];} dg_interp_kern_p_list_vlasov;

struct gkyl_dg_interpolate_kernels {
  dg_interp_t interp;  // Kernel that performs the interpolation.
  dg_interp_grid2stencilIdx_t grid2stencil; // Translate grid to stencil index.
};

// Primary struct in this updater.
struct gkyl_dg_interpolate {
  int ndim; // Phase space dimensionality of the fields.
  bool use_gpu; // Whether to use the GPU.
  int num_interp_dirs; // Number of directions to interpolate in.
  int interp_dirs[GKYL_MAX_DIM]; // Directions to interpolate in.
  struct gkyl_rect_grid *grids; // Target grids, one for each interp_dir.
  struct gkyl_dg_interpolate **interp_ops; // Interp updater for each dimension.
  struct gkyl_range *ranges; // Ranges to interpolate in for each dimension.
  struct gkyl_array **fields; // Field at each grid.
  struct gkyl_rect_grid grid_do; // Donor grid.
  struct gkyl_rect_grid grid_tar; // Target grid.
  int dir; // Interpolation direction.
  double dxRat; // Ratio of donor to target cell length in interpolation dir.
  int *offset_upper; // Upper index of the offset in each direction for each stencil.
  struct gkyl_dg_interpolate_kernels *kernels;
};

// Serendipity  kernels.
GKYL_CU_D
static const dg_interp_kern_p_list dg_interp_kern_list_ser[] = {
  // 1x
  { .list = {
      { dg_interpolate_1x_ser_p1_x, NULL, NULL, NULL, NULL, NULL, },
      { dg_interpolate_1x_ser_p2_x, NULL, NULL, NULL, NULL, NULL, },
    },
  },
  // 2x
  { .list = {
      { dg_interpolate_2x_ser_p1_x, dg_interpolate_2x_ser_p1_y, NULL, NULL, NULL, NULL, },
      { dg_interpolate_2x_ser_p2_x, dg_interpolate_2x_ser_p2_y, NULL, NULL, NULL, NULL, },
    },
  },
  // 3x
  { .list = {
      { dg_interpolate_3x_ser_p1_x, dg_interpolate_3x_ser_p1_y, dg_interpolate_3x_ser_p1_z, NULL, NULL, NULL, },
      { dg_interpolate_3x_ser_p2_x, dg_interpolate_3x_ser_p2_y, dg_interpolate_3x_ser_p2_z, NULL, NULL, NULL, },
    },
  },
};

GKYL_CU_D
static const dg_interp_kern_p_list dg_interp_kern_list_gk_ser[] = {
  // 1x1v
  { .list = {
      { dg_interpolate_gyrokinetic_1x1v_ser_p1_x, dg_interpolate_gyrokinetic_1x1v_ser_p1_vpar, NULL, NULL, NULL, NULL, },
      { dg_interpolate_gyrokinetic_1x1v_ser_p2_x, dg_interpolate_gyrokinetic_1x1v_ser_p2_vpar, NULL, NULL, NULL, NULL, },
    },
  },
  // 1x2v
  { .list = {
      { dg_interpolate_gyrokinetic_1x2v_ser_p1_x, dg_interpolate_gyrokinetic_1x2v_ser_p1_vpar, dg_interpolate_gyrokinetic_1x2v_ser_p1_mu, NULL, NULL, NULL, },
      { dg_interpolate_gyrokinetic_1x2v_ser_p2_x, dg_interpolate_gyrokinetic_1x2v_ser_p2_vpar, dg_interpolate_gyrokinetic_1x2v_ser_p2_mu, NULL, NULL, NULL, },
    },
  },
  // 2x2v
  { .list = {
      { dg_interpolate_gyrokinetic_2x2v_ser_p1_x, dg_interpolate_gyrokinetic_2x2v_ser_p1_z, dg_interpolate_gyrokinetic_2x2v_ser_p1_vpar, dg_interpolate_gyrokinetic_2x2v_ser_p1_mu, NULL, NULL, },
      { dg_interpolate_gyrokinetic_2x2v_ser_p2_x, dg_interpolate_gyrokinetic_2x2v_ser_p2_z, dg_interpolate_gyrokinetic_2x2v_ser_p2_vpar, dg_interpolate_gyrokinetic_2x2v_ser_p2_mu, NULL, NULL, },
    },
  },
  // 3x2v
  { .list = {
      { dg_interpolate_gyrokinetic_3x2v_ser_p1_x, dg_interpolate_gyrokinetic_3x2v_ser_p1_y, dg_interpolate_gyrokinetic_3x2v_ser_p1_z, dg_interpolate_gyrokinetic_3x2v_ser_p1_vpar, dg_interpolate_gyrokinetic_3x2v_ser_p1_mu, NULL, },
      { NULL, NULL, NULL, NULL, NULL, NULL, },
    },
  },
};

GKYL_CU_D
static const dg_interp_kern_p_list_vlasov dg_interp_kern_list_vlasov_ser[] = {
  // 1x
  { .vdim = {
      { .list = {
          { dg_interpolate_vlasov_1x1v_ser_p1_x, dg_interpolate_vlasov_1x1v_ser_p1_vx, NULL, NULL, NULL, NULL, },
          { dg_interpolate_vlasov_1x1v_ser_p2_x, dg_interpolate_vlasov_1x1v_ser_p2_vx, NULL, NULL, NULL, NULL, },
        },
      },
      { .list = {
          { dg_interpolate_vlasov_1x2v_ser_p1_x, dg_interpolate_vlasov_1x2v_ser_p1_vx, dg_interpolate_vlasov_1x2v_ser_p1_vy, NULL, NULL, NULL, },
          { dg_interpolate_vlasov_1x2v_ser_p2_x, dg_interpolate_vlasov_1x2v_ser_p2_vx, dg_interpolate_vlasov_1x2v_ser_p2_vy, NULL, NULL, NULL, },
        },
      },
      { .list = {
          { dg_interpolate_vlasov_1x3v_ser_p1_x, dg_interpolate_vlasov_1x3v_ser_p1_vx, dg_interpolate_vlasov_1x2v_ser_p1_vy, dg_interpolate_vlasov_1x3v_ser_p1_vz, NULL, NULL, },
          { dg_interpolate_vlasov_1x3v_ser_p2_x, dg_interpolate_vlasov_1x3v_ser_p2_vx, dg_interpolate_vlasov_1x2v_ser_p2_vy, dg_interpolate_vlasov_1x3v_ser_p2_vz, NULL, NULL, },
        },
      },
    },
  },
  // 2x
  { .vdim = {
      { .list = {
          { NULL, NULL, NULL, NULL, NULL, NULL, },
          { NULL, NULL, NULL, NULL, NULL, NULL, },
        },
      },
      { .list = {
          { dg_interpolate_vlasov_2x2v_ser_p1_x, dg_interpolate_vlasov_2x2v_ser_p1_y, dg_interpolate_vlasov_2x2v_ser_p1_vx, dg_interpolate_vlasov_2x2v_ser_p1_vy, NULL, NULL, },
          { dg_interpolate_vlasov_2x2v_ser_p2_x, dg_interpolate_vlasov_2x2v_ser_p2_y, dg_interpolate_vlasov_2x2v_ser_p2_vx, dg_interpolate_vlasov_2x2v_ser_p2_vy, NULL, NULL, },
        },
      },
      { .list = {
          { dg_interpolate_vlasov_2x3v_ser_p1_x, dg_interpolate_vlasov_2x3v_ser_p1_y, dg_interpolate_vlasov_2x3v_ser_p1_vx, dg_interpolate_vlasov_2x2v_ser_p1_vy, dg_interpolate_vlasov_2x3v_ser_p1_vz, NULL, },
          { dg_interpolate_vlasov_2x3v_ser_p2_x, dg_interpolate_vlasov_2x3v_ser_p2_y, dg_interpolate_vlasov_2x3v_ser_p2_vx, dg_interpolate_vlasov_2x2v_ser_p2_vy, dg_interpolate_vlasov_2x3v_ser_p2_vz, NULL, },
        },
      },
    },
  },
  // 3x
  { .vdim = {
      { .list = {
          { NULL, NULL, NULL, NULL, NULL, NULL, },
          { NULL, NULL, NULL, NULL, NULL, NULL, },
        },
      },
      { .list = {
          { NULL, NULL, NULL, NULL, NULL, NULL, },
          { NULL, NULL, NULL, NULL, NULL, NULL, },
        },
      },
      { .list = {
          { dg_interpolate_vlasov_3x3v_ser_p1_x, dg_interpolate_vlasov_3x3v_ser_p1_y, dg_interpolate_vlasov_3x3v_ser_p1_z, dg_interpolate_vlasov_3x3v_ser_p1_vx, dg_interpolate_vlasov_3x3v_ser_p1_vy, dg_interpolate_vlasov_3x3v_ser_p1_vz, },
          { NULL, NULL, NULL, NULL, NULL, NULL, },
        },
      },
    },
  },
};

#ifdef GKYL_HAVE_CUDA
// Declaration of cuda device functions.
void dg_interp_choose_kernel_cu(struct gkyl_dg_interpolate_kernels *kernels,
  int cdim, struct gkyl_basis basis, int dir, double dxRat);

void gkyl_dg_interpolate_advance_1x_cu(gkyl_dg_interpolate* up,
  const struct gkyl_range *phase_rng_do, const struct gkyl_range *phase_rng_tar,
  const struct gkyl_array *GKYL_RESTRICT fdo, struct gkyl_array *GKYL_RESTRICT ftar);
#endif

GKYL_CU_D
static dg_interp_t
dg_interp_choose_gk_interp_kernel(int cdim, struct gkyl_basis basis, int dir)
{
  enum gkyl_basis_type basis_type = basis.b_type;
  int ndim = basis.ndim;
  int vdim = ndim - cdim;
  int poly_order = basis.poly_order;

  if (vdim == 0) {
    switch (basis_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
        return dg_interp_kern_list_ser[ndim-1].list[poly_order-1].dirs[dir];
        break;
      default:
        assert(false);
        break;
    }
  }
  else {
    switch (basis_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
        return dg_interp_kern_list_vlasov_ser[cdim-1].vdim[vdim-1].list[poly_order-1].dirs[dir];
        break;
      case GKYL_BASIS_MODAL_HYBRID:
        return dg_interp_kern_list_vlasov_ser[cdim-1].vdim[vdim-1].list[poly_order-1].dirs[dir];
        break;
      case GKYL_BASIS_MODAL_GKHYBRID:
        return dg_interp_kern_list_gk_ser[ndim-2].list[poly_order-1].dirs[dir];
        break;
      default:
        assert(false);
        break;
    }
  }

  return 0;
}

static int dg_interp_prime_factors(int n, int *pfs, int pfs_size)
{
  // Find the prime factors of number `n`, and put them into `pfs`. We assume
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
static int dg_interp_index_stencil_map_refine(int idx,
  int num_cells, double dx_rat)
{
  // Given an index 'idx' to a cell in the coarse grid with 'num_cells'
  // cells, return the index of the refinement stencil needed, within
  // the table that holds stencils. Here 'dx_rat' is the ratio of
  // donor to target cell length in the interpolating direction.
  double remDecL = (idx-1)*dx_rat-floor((idx-1)*dx_rat);
  double remDecU = ceil(idx*dx_rat)-idx*dx_rat;
  int stencilOut = 1;
  if ((idx == 1) ||   // First cell.
      (remDecL == 0) ||    // Interior cell with a left-boundary-like stencil.
      ((remDecL <= 0.5) && (remDecU <= 0.5))) {
    stencilOut = 2*stencilOut - 1;
  }
  else if ((idx == num_cells) || // Last cell.
          (remDecU == 0)) { // Interior cell with a right-boundary-like stencil.
    stencilOut = 2*stencilOut;
  }
  return stencilOut;
}

GKYL_CU_DH 
static int dg_interp_index_stencil_map_coarsen(int idx,
  int num_cells, double dx_rat)
{
  // Given an index 'idx' to a cell in the fine grid with 'num_cells'
  // cells, return the index of the coarsening stencil needed, within
  // the table that holds stencils. Here 'dx_rat' is the ratio of
  // donor to target cell length in the interpolating direction.
  double remDecL = (idx-1)*dx_rat-floor(idx*dx_rat);
  double remDecU = ceil(idx*dx_rat)-idx*dx_rat;
  int stencilOut = 1;
  if ((idx == 1) || // First cell.
      (remDecL == 0) || // Interior cell with a left-boundary-like stencil.
      ((remDecL > 0) && (remDecU > 0))) {
     stencilOut = 2*stencilOut - 1;
  }
  else if ((idx == num_cells) || // Last cell.
          (remDecU == 0)) { // Interior cell with a right-boundary-like stencil.
     stencilOut = 2*stencilOut;
  }
  return stencilOut;
}

