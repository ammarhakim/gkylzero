// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_dg_vlasov_gen_geo_alpha_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef int (*vlasov_gen_geo_alpha_surf_t)(const double *w, const double *dxv, 
  const double *tvComp, const double *gij, 
  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 

typedef void (*vlasov_gen_geo_cot_vec_t)(const double *tvComp, const double *gij, 
  double* GKYL_RESTRICT cot_vec);

// for use in kernel tables
typedef struct { vlasov_gen_geo_alpha_surf_t kernels[3]; } gkyl_dg_vlasov_gen_geo_alpha_surf_kern_list;
typedef struct { vlasov_gen_geo_cot_vec_t kernels[3]; } gkyl_dg_vlasov_gen_geo_cot_vec_kern_list;
typedef struct { vlasov_bmag_cart_t kernels[3]; } gkyl_dg_vlasov_bmag_cart_kern_list;

struct gkyl_dg_calc_vlasov_gen_geo_vars {
  struct gkyl_rect_grid phase_grid; // Phase space grid for cell spacing and cell center
  int cdim; // Configuration space dimensionality
  int pdim; // Phase space dimensionality
  vlasov_gen_geo_alpha_surf_t alpha_surf[3]; // kernel for computing surface expansion of phase space flux alpha
  vlasov_gen_geo_alpha_surf_t alpha_edge_surf[3]; // kernel for computing surface expansion of phase space flux alpha
                                               // at upper configuration space edge
  vlasov_gen_geo_cot_vec_t calc_cot_vec; // kernel for computing volume expansion of cotangent vectors e^i
  
  const struct gk_geometry *gk_geom; // Pointer to geometry struct

  uint32_t flags;
  struct gkyl_dg_calc_vlasov_gen_geo_vars *on_dev; // pointer to itself or device data
};

//
// Serendipity surface kernels general geometry
//
// Vlasov general geometry phase space flux alpha surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_vlasov_gen_geo_alpha_surf_kern_list ser_vlasov_gen_geo_alpha_surfx_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_gen_geo_alpha_surfx_3x3v_ser_p1, NULL }, // 2
};

// Vlasov general geometry phase space flux alpha edge surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_vlasov_gen_geo_alpha_surf_kern_list ser_vlasov_gen_geo_alpha_edge_surfx_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_gen_geo_alpha_edge_surfx_3x3v_ser_p1, NULL }, // 2
};

// Vlasov general geometry phase space flux alpha surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_vlasov_gen_geo_alpha_surf_kern_list ser_vlasov_gen_geo_alpha_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_gen_geo_alpha_surfy_3x3v_ser_p1, NULL }, // 2
};

// Vlasov general geometry phase space flux alpha edge surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_vlasov_gen_geo_alpha_surf_kern_list ser_vlasov_gen_geo_alpha_edge_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_gen_geo_alpha_edge_surfy_3x3v_ser_p1, NULL }, // 2
};

// Vlasov general geometry phase space flux alpha surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_vlasov_gen_geo_alpha_surf_kern_list ser_vlasov_gen_geo_alpha_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_gen_geo_alpha_surfz_3x3v_ser_p1, NULL }, // 2
};

// Vlasov general geometry phase space flux alpha edge surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_vlasov_gen_geo_alpha_surf_kern_list ser_vlasov_gen_geo_alpha_edge_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_gen_geo_alpha_edge_surfz_3x3v_ser_p1, NULL }, // 2
};

// Vlasov general geometry cotangent vector (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_vlasov_gen_geo_cot_vec_kern_list ser_vlasov_gen_geo_cot_vec_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_gen_geo_cot_vec_3x_ser_p1, NULL }, // 2
};

GKYL_CU_D
static vlasov_gen_geo_alpha_surf_t
choose_vlasov_gen_geo_alpha_surf_kern(int dir, int cdim, int poly_order)
{
  if (dir == 0)
    return ser_vlasov_gen_geo_alpha_surfx_kernels[cdim-1].kernels[poly_order];
  else if (dir == 1)
    return ser_vlasov_gen_geo_alpha_surfy_kernels[cdim-1].kernels[poly_order];
  else if (dir == 2)
    return ser_vlasov_gen_geo_alpha_surfz_kernels[cdim-1].kernels[poly_order];
  else
    return NULL;
}

GKYL_CU_D
static vlasov_gen_geo_alpha_surf_t
choose_vlasov_gen_geo_alpha_edge_surf_kern(int dir, int cdim, int poly_order)
{
  if (dir == 0)
    return ser_vlasov_gen_geo_alpha_edge_surfx_kernels[cdim-1].kernels[poly_order];
  else if (dir == 1)
    return ser_vlasov_gen_geo_alpha_edge_surfy_kernels[cdim-1].kernels[poly_order];
  else if (dir == 2)
    return ser_vlasov_gen_geo_alpha_edge_surfz_kernels[cdim-1].kernels[poly_order];
  else
    return NULL;
}

GKYL_CU_D
static vlasov_gen_geo_cot_vec_t
choose_vlasov_gen_geo_cot_vec_kern(int cdim, int poly_order)
{
  return ser_vlasov_gen_geo_cot_vec_kernels[cdim-1].kernels[poly_order];
}
