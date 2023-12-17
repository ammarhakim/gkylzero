// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

typedef void (*gyrokinetic_Bstar_Bmag_t)(const double *w, const double *dxv, 
  const double q_, const double m_, 
  const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
  double* GKYL_RESTRICT Bstar_Bmag); 

typedef void (*gyrokinetic_alpha_surf_t)(const double *w, const double *dxv, 
  const double q_, const double m_, 
  const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
  const double *phi, const double *Bstar_Bmag, double* GKYL_RESTRICT alpha_surf); 

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below.
GKYL_CU_D
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense.
  {-1,  0,  1}, // 1x kernel indices.
  {-1, -1,  2}, // 2x kernel indices.
  {-1, -1,  3}, // 3x kernel indices.
};

// for use in kernel tables
typedef struct { gyrokinetic_Bstar_Bmag_t kernels[3]; } gkyl_dg_gyrokinetic_Bstar_Bmag_kern_list;
typedef struct { gyrokinetic_alpha_surf_t kernels[3]; } gkyl_dg_gyrokinetic_alpha_surf_kern_list;

struct gkyl_dg_calc_gyrokinetic_vars {
  struct gkyl_rect_grid phase_grid; // Phase space grid for cell spacing and cell center
  int cdim; // Configuration space dimensionality
  int pdim; // Phase space dimensionality
  gyrokinetic_Bstar_Bmag_t Bstar_Bmag; // kernel for computing surface expansion of phase space flux alpha
  gyrokinetic_alpha_surf_t alpha_surf[4]; // kernel for computing surface expansion of phase space flux alpha
  gyrokinetic_alpha_surf_t alpha_edge_surf[3]; // kernel for computing surface expansion of phase space flux alpha
                                               // at upper configuration space edge
  double charge, mass;
  const struct gk_geometry *gk_geom; // Pointer to geometry struct

  uint32_t flags;
  struct gkyl_dg_calc_gyrokinetic_vars *on_dev; // pointer to itself or device data
};

// Gyrokinetic Bstar/Bmag volume expansion (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_Bstar_Bmag_kern_list ser_gyrokinetic_Bstar_Bmag_kernels[] = {
  { NULL, gyrokinetic_Bstar_Bmag_1x_ser_p1, gyrokinetic_Bstar_Bmag_1x_ser_p2 }, // 0
  { NULL, gyrokinetic_Bstar_Bmag_2x_ser_p1, gyrokinetic_Bstar_Bmag_2x_ser_p2 }, // 1
  { NULL, gyrokinetic_Bstar_Bmag_3x_ser_p1, gyrokinetic_Bstar_Bmag_3x_ser_p2 }, // 2
};

// Gyrokinetic phase space flux alpha surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_alpha_surf_kern_list ser_gyrokinetic_alpha_surfx_kernels[] = {
  { NULL, gyrokinetic_alpha_surfx_1x1v_ser_p1, gyrokinetic_alpha_surfx_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_alpha_surfx_1x2v_ser_p1, gyrokinetic_alpha_surfx_1x2v_ser_p2 }, // 1
  { NULL, gyrokinetic_alpha_surfx_2x2v_ser_p1, gyrokinetic_alpha_surfx_2x2v_ser_p2 }, // 2
  { NULL, gyrokinetic_alpha_surfx_3x2v_ser_p1, gyrokinetic_alpha_surfx_3x2v_ser_p2 }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in x (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_alpha_surf_kern_list ser_gyrokinetic_alpha_edge_surfx_kernels[] = {
  { NULL, gyrokinetic_alpha_edge_surfx_1x1v_ser_p1, gyrokinetic_alpha_edge_surfx_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_alpha_edge_surfx_1x2v_ser_p1, gyrokinetic_alpha_edge_surfx_1x2v_ser_p2 }, // 1
  { NULL, gyrokinetic_alpha_edge_surfx_2x2v_ser_p1, gyrokinetic_alpha_edge_surfx_2x2v_ser_p2 }, // 2
  { NULL, gyrokinetic_alpha_edge_surfx_3x2v_ser_p1, gyrokinetic_alpha_edge_surfx_3x2v_ser_p2 }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_alpha_surf_kern_list ser_gyrokinetic_alpha_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gyrokinetic_alpha_surfy_2x2v_ser_p1, gyrokinetic_alpha_surfy_2x2v_ser_p2 }, // 2
  { NULL, gyrokinetic_alpha_surfy_3x2v_ser_p1, gyrokinetic_alpha_surfy_3x2v_ser_p2 }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in y (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_alpha_surf_kern_list ser_gyrokinetic_alpha_edge_surfy_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, gyrokinetic_alpha_edge_surfy_2x2v_ser_p1, gyrokinetic_alpha_edge_surfy_2x2v_ser_p2 }, // 2
  { NULL, gyrokinetic_alpha_edge_surfy_3x2v_ser_p1, gyrokinetic_alpha_edge_surfy_3x2v_ser_p2 }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_alpha_surf_kern_list ser_gyrokinetic_alpha_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gyrokinetic_alpha_surfz_3x2v_ser_p1, gyrokinetic_alpha_surfz_3x2v_ser_p2 }, // 3
};

// Gyrokinetic phase space flux alpha edge surface expansions in z (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_alpha_surf_kern_list ser_gyrokinetic_alpha_edge_surfz_kernels[] = {
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  { NULL, gyrokinetic_alpha_edge_surfz_3x2v_ser_p1, gyrokinetic_alpha_edge_surfz_3x2v_ser_p2 }, // 3
};

// Gyrokinetic phase space flux alpha surface expansions in vpar (Serendipity kernels)
GKYL_CU_D
static const gkyl_dg_gyrokinetic_alpha_surf_kern_list ser_gyrokinetic_alpha_surfvpar_kernels[] = {
  { NULL, gyrokinetic_alpha_surfvpar_1x1v_ser_p1, gyrokinetic_alpha_surfvpar_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_alpha_surfvpar_1x2v_ser_p1, gyrokinetic_alpha_surfvpar_1x2v_ser_p2 }, // 1
  { NULL, gyrokinetic_alpha_surfvpar_2x2v_ser_p1, gyrokinetic_alpha_surfvpar_2x2v_ser_p2 }, // 2
  { NULL, gyrokinetic_alpha_surfvpar_3x2v_ser_p1, gyrokinetic_alpha_surfvpar_3x2v_ser_p2 }, // 3
};

GKYL_CU_D
static gyrokinetic_Bstar_Bmag_t
choose_gyrokinetic_Bstar_Bmag_kern(int cdim, int poly_order)
{
  return ser_gyrokinetic_Bstar_Bmag_kernels[cdim-1].kernels[poly_order];
}

GKYL_CU_D
static gyrokinetic_alpha_surf_t
choose_gyrokinetic_alpha_surf_conf_kern(int dir, int cdim, int vdim, int poly_order)
{
  if (dir == 0)
    return ser_gyrokinetic_alpha_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 1)
    return ser_gyrokinetic_alpha_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 2)
    return ser_gyrokinetic_alpha_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else
    return NULL;
}

GKYL_CU_D
static gyrokinetic_alpha_surf_t
choose_gyrokinetic_alpha_edge_surf_conf_kern(int dir, int cdim, int vdim, int poly_order)
{
  if (dir == 0)
    return ser_gyrokinetic_alpha_edge_surfx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 1)
    return ser_gyrokinetic_alpha_edge_surfy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else if (dir == 2)
    return ser_gyrokinetic_alpha_edge_surfz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  else
    return NULL;
}

GKYL_CU_D
static gyrokinetic_alpha_surf_t
choose_gyrokinetic_alpha_surf_vpar_kern(int cdim, int vdim, int poly_order)
{
  return ser_gyrokinetic_alpha_surfvpar_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}
