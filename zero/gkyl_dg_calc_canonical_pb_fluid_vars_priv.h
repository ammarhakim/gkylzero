// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_canonical_pb_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef int (*canonical_pb_fluid_alpha_surf_t)(const double *w, const double *dxv, const double *phi,
  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 
typedef void (*canonical_pb_fluid_source_t)(const double *dxv, double alpha, double kappa, 
  const double *background_n_gradient, const double *phi, const double *fluid, double* GKYL_RESTRICT rhs); 

// for use in kernel tables
typedef struct { canonical_pb_fluid_alpha_surf_t kernels[3]; } gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list;
typedef struct { canonical_pb_fluid_source_t kernels[3]; } gkyl_dg_canonical_pb_fluid_source_kern_list;

struct gkyl_dg_calc_canonical_pb_fluid_vars {
  struct gkyl_rect_grid conf_grid; // Phase space grid for cell spacing and cell center
  int cdim; // Configuration space dimensionality
  double alpha; // Adiabaticity parameter for adiabatic coupling of vorticity and density.
  double kappa; // Constant density gradient scale length (for turbulence drive). 
  bool is_modified; // Boolean parameter for if we are doing the modified Hasegawa-Wakatani. 
  canonical_pb_fluid_alpha_surf_t alpha_surf[3]; // kernel for computing surface expansion of configuration space flux alpha
  canonical_pb_fluid_alpha_surf_t alpha_edge_surf[3]; // kernel for computing surface expansion of configuration space flux alpha
                                                // at upper configuration space edge
  canonical_pb_fluid_source_t canonical_pb_fluid_source; // Canonical pb fluid variables source function
  uint32_t flags;
  struct gkyl_dg_calc_canonical_pb_fluid_vars *on_dev; // pointer to itself or device data
};

// Canonical PB fluid alpha surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list ser_canonical_pb_fluid_alpha_surfx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_surfx_2x_ser_p1, canonical_pb_alpha_surfx_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB fluid alpha upper edge of configuration space surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list ser_canonical_pb_fluid_alpha_edge_surfx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_edge_surfx_2x_ser_p1, canonical_pb_alpha_edge_surfx_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB fluid alpha surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list ser_canonical_pb_fluid_alpha_surfy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_surfy_2x_ser_p1, canonical_pb_alpha_surfy_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB fluid alpha upper edge of configuration space surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list ser_canonical_pb_fluid_alpha_edge_surfy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_edge_surfy_2x_ser_p1, canonical_pb_alpha_edge_surfy_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB fluid alpha surface expansions in x (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list tensor_canonical_pb_fluid_alpha_surfx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_surfx_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB fluid alpha upper edge of configuration space surface expansions in x (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list tensor_canonical_pb_fluid_alpha_edge_surfx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_edge_surfx_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB fluid alpha surface expansions in y (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list tensor_canonical_pb_fluid_alpha_surfy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_surfy_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB fluid alpha upper edge of configuration space surface expansions in y (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_alpha_surf_kern_list tensor_canonical_pb_fluid_alpha_edge_surfy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_alpha_edge_surfy_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB Hasegawa-Mima source update (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_source_kern_list ser_canonical_pb_fluid_hasegawa_mima_source_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_fluid_hasegawa_mima_source_2x_ser_p1, canonical_pb_fluid_hasegawa_mima_source_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB Hasegawa-Mima source update (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_source_kern_list tensor_canonical_pb_fluid_hasegawa_mima_source_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_fluid_hasegawa_mima_source_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB Hasegawa-Wakatani source update (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_source_kern_list ser_canonical_pb_fluid_hasegawa_wakatani_source_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_fluid_hasegawa_wakatani_source_2x_ser_p1, canonical_pb_fluid_hasegawa_wakatani_source_2x_ser_p2 }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

// Canonical PB Hasegawa-Wakatani source update (Tensor kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_fluid_source_kern_list tensor_canonical_pb_fluid_hasegawa_wakatani_source_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_fluid_hasegawa_wakatani_source_2x_ser_p1, NULL }, // 1
  // 3x kernels
  { NULL, NULL, NULL }, // 2
};

GKYL_CU_D
static canonical_pb_fluid_alpha_surf_t
choose_canonical_pb_fluid_alpha_surf_kern(enum gkyl_basis_type b_type, int dir, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (dir == 0)
        return ser_canonical_pb_fluid_alpha_surfx_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return ser_canonical_pb_fluid_alpha_surfy_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      if (dir == 0)
        return tensor_canonical_pb_fluid_alpha_surfx_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return tensor_canonical_pb_fluid_alpha_surfy_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static canonical_pb_fluid_alpha_surf_t
choose_canonical_pb_fluid_alpha_edge_surf_kern(enum gkyl_basis_type b_type, int dir, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      if (dir == 0)
        return ser_canonical_pb_fluid_alpha_edge_surfx_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return ser_canonical_pb_fluid_alpha_edge_surfy_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    case GKYL_BASIS_MODAL_TENSOR:
      if (dir == 0)
        return tensor_canonical_pb_fluid_alpha_edge_surfx_kernels[cdim-1].kernels[poly_order];
      else if (dir == 1)
        return tensor_canonical_pb_fluid_alpha_edge_surfy_kernels[cdim-1].kernels[poly_order];
      else
        return NULL;
      break;
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static canonical_pb_fluid_source_t
choose_canonical_pb_fluid_hasegawa_mima_source_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_canonical_pb_fluid_hasegawa_mima_source_kernels[cdim-1].kernels[poly_order];
      break; 
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_canonical_pb_fluid_hasegawa_mima_source_kernels[cdim-1].kernels[poly_order];
      break; 
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static canonical_pb_fluid_source_t
choose_canonical_pb_fluid_hasegawa_wakatani_source_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_canonical_pb_fluid_hasegawa_wakatani_source_kernels[cdim-1].kernels[poly_order];
      break; 
    case GKYL_BASIS_MODAL_TENSOR:
      return tensor_canonical_pb_fluid_hasegawa_wakatani_source_kernels[cdim-1].kernels[poly_order];
      break; 
    default:
      assert(false);
      break;  
  }
}

GKYL_CU_D
static canonical_pb_fluid_source_t
choose_canonical_pb_fluid_default_source_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  return canonical_pb_fluid_default_source;
}
