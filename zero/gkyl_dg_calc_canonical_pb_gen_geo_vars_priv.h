// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_canonical_pb_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef int (*canonical_pb_gen_geo_alpha_surf_t)(const double *w, const double *dxv, 
  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf); 

// for use in kernel tables
typedef struct { canonical_pb_gen_geo_alpha_surf_t kernels[3]; } gkyl_dg_canonical_pb_gen_geo_alpha_surf_kern_list;

struct gkyl_dg_calc_canonical_pb_gen_geo_vars {
  struct gkyl_rect_grid phase_grid; // Phase space grid for cell spacing and cell center
  int cdim; // Configuration space dimensionality
  int pdim; // Phase space dimensionality
  canonical_pb_gen_geo_alpha_surf_t alpha_surf[3]; // kernel for computing surface expansion of phase space flux alpha
  canonical_pb_gen_geo_alpha_surf_t alpha_edge_surf[3]; // kernel for computing surface expansion of phase space flux alpha
                                               // at upper configuration space edge

  uint32_t flags;
  struct gkyl_dg_calc_canonical_pb_gen_geo_vars *on_dev; // pointer to itself or device data
};

//
// Serendipity surface kernels general geometry
//
// canonical_pb general geometry phase space flux alpha surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_gen_geo_alpha_surf_kern_list ser_canonical_pb_gen_geo_alpha_surfx_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, canonical_pb_alpha_surfx_3x3v_ser_p1, NULL }, // 2
};

// canonical_pb general geometry phase space flux alpha edge surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_gen_geo_alpha_surf_kern_list ser_canonical_pb_gen_geo_alpha_edge_surfx_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, canonical_pb_alpha_edge_surfx_3x3v_ser_p1, NULL }, // 2
};

// canonical_pb general geometry phase space flux alpha surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_gen_geo_alpha_surf_kern_list ser_canonical_pb_gen_geo_alpha_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, canonical_pb_alpha_surfy_3x3v_ser_p1, NULL }, // 2
};

// canonical_pb general geometry phase space flux alpha edge surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_gen_geo_alpha_surf_kern_list ser_canonical_pb_gen_geo_alpha_edge_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, canonical_pb_alpha_edge_surfy_3x3v_ser_p1, NULL }, // 2
};

// canonical_pb general geometry phase space flux alpha surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_gen_geo_alpha_surf_kern_list ser_canonical_pb_gen_geo_alpha_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, canonical_pb_alpha_surfz_3x3v_ser_p1, NULL }, // 2
};

// canonical_pb general geometry phase space flux alpha edge surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_canonical_pb_gen_geo_alpha_surf_kern_list ser_canonical_pb_gen_geo_alpha_edge_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, canonical_pb_alpha_edge_surfz_3x3v_ser_p1, NULL }, // 2
};

GKYL_CU_D
static canonical_pb_gen_geo_alpha_surf_t
choose_canonical_pb_gen_geo_alpha_surf_kern(int dir, int cdim, int poly_order)
{
  if (dir == 0)
    return ser_canonical_pb_gen_geo_alpha_surfx_kernels[cdim].kernels[poly_order];
  else if (dir == 1)
    return ser_canonical_pb_gen_geo_alpha_surfy_kernels[cdim].kernels[poly_order];
  else if (dir == 2)
    return ser_canonical_pb_gen_geo_alpha_surfz_kernels[cdim].kernels[poly_order];
  else
    return NULL;
}

GKYL_CU_D
static canonical_pb_gen_geo_alpha_surf_t
choose_canonical_pb_gen_geo_alpha_edge_surf_kern(int dir, int cdim, int poly_order)
{
  if (dir == 0)
    return ser_canonical_pb_gen_geo_alpha_edge_surfx_kernels[cdim].kernels[poly_order];
  else if (dir == 1)
    return ser_canonical_pb_gen_geo_alpha_edge_surfy_kernels[cdim].kernels[poly_order];
  else if (dir == 2)
    return ser_canonical_pb_gen_geo_alpha_edge_surfz_kernels[cdim].kernels[poly_order];
  else
    return NULL;
}