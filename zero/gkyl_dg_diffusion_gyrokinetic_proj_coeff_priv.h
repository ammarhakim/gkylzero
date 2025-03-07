#pragma once

#include <gkyl_dg_diffusion_gyrokinetic_proj_coeff.h>
#include <gkyl_dg_diffusion_gyrokinetic_proj_coeff_kernels.h>
#include <assert.h>

// Types for various kernels
typedef void (*diffusion_proj_coeff_t)(const double *xc, const double *dx, const double *nu,
  const double *xi, double mass, double vtsq_min, const double *gijJ, const double *bmag,
  const double *vtsq, const double *vmap, const double *vmapSq, double *out);

struct gkyl_dg_diffusion_gyrokinetic_proj_coeff_kernels {
  diffusion_proj_coeff_t projC; // Kernel projecting the diff coeff onto basis.
};

struct gkyl_dg_diffusion_gyrokinetic_proj_coeff {
  bool use_gpu; // Whether to use the GPU.
  int cdim; // Conf space dimensionality.
  int pdim; // Phase space dimensionality.
  struct gkyl_rect_grid *grid; // Grid object.
  double mass; // Species mass.
  double vtsq_min; // Minimum thermal speed squared supported by the grid.
  double nu[3], xi[3]; // Coefficients in the model.
  struct gkyl_dg_diffusion_gyrokinetic_proj_coeff_kernels *kernels;
};

// For use in kernel tables.
typedef struct {diffusion_proj_coeff_t kernels[7];} diff_proj_coeff_kern_list;
typedef struct {diff_proj_coeff_kern_list list[2];} diff_proj_coeff_kern_list_pOrder;

// Serendipity kernels.
GKYL_CU_D
static const diff_proj_coeff_kern_list_pOrder ser_diff_proj_coeff_kernels[] = {
  // 1x1v
  {.list = 
    {
      {dg_diffusion_gyrokinetic_proj_coeff_1x1v_ser_p1_diffdirsx, NULL, NULL, NULL, NULL, NULL, NULL},
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
    },
  },
  // 1x2v
  {.list = 
    {
      {dg_diffusion_gyrokinetic_proj_coeff_1x2v_ser_p1_diffdirsx, NULL, NULL, NULL, NULL, NULL, NULL},
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
    },
  },
  // 2x2v
  {.list = 
    {
      {dg_diffusion_gyrokinetic_proj_coeff_2x2v_ser_p1_diffdirsx, dg_diffusion_gyrokinetic_proj_coeff_2x2v_ser_p1_diffdirsz, dg_diffusion_gyrokinetic_proj_coeff_2x2v_ser_p1_diffdirsxz, NULL, NULL, NULL, NULL},
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
    },
  },
  // 3x2v
  {.list = 
    {
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
      {NULL, NULL, NULL, NULL, NULL, NULL, NULL},
    },
  },
};

// Macro for choosing kernels.
#define CKVOL(lst,cdim,vdim,poly_order,diffdir_linidx) lst[cdim+vdim-2].list[poly_order-1].kernels[diffdir_linidx]

static inline int dg_diff_gk_projC_diffdirs_linidx(const bool *isdirdiff, int cdim) {
  // Compute the linear index into the array of volume kernels (one
  // kernel for each combination of diffusive directions).
  bool diff_in_dir[GKYL_MAX_CDIM];
  if (isdirdiff)
    for (int d=0; d<cdim; d++) diff_in_dir[d] = isdirdiff[d];
  else
    for (int d=0; d<cdim; d++) diff_in_dir[d] = true;

  // Linear index into list of volume kernels.
  int dirs_bin_key[] = {1,2,4,8,16,32}; // Binary: 000001, 000010, 000100, 001000, 010000, 100000.
  int dirs_linidx = 0; // Binary 000000.
  for (int d=0; d<cdim; d++) {
     if (diff_in_dir[d]) dirs_linidx = dirs_linidx | dirs_bin_key[d];
  }
  dirs_linidx -= 1;
  return dirs_linidx;
}

GKYL_CU_D
static void dg_diff_gk_projC_choose_kernel(struct gkyl_dg_diffusion_gyrokinetic_proj_coeff_kernels *kernels,
  int cdim, struct gkyl_basis pbasis, const bool *diff_in_dir, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    dg_diff_gk_projC_choose_kernel_cu(kernels, cdim, pbasis, diff_in_dir);
    return;
  }
#endif

  int pdim = pbasis.ndim;
  enum gkyl_basis_type basis_type = pbasis.b_type;
  int poly_order = pbasis.poly_order;

  int dirs_linidx = dg_diff_gk_projC_diffdirs_linidx(diff_in_dir, cdim);

  switch (basis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->projC = ser_diff_proj_coeff_kernels[pdim-2].list[poly_order-1].kernels[dirs_linidx];
      break;
    default:
      assert(false);
      break;
  }
}
