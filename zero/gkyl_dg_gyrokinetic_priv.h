#pragma once

// Private header, not for direct use in user code.

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_gk_geometry.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_dg_gyrokinetic.h>

// Types for various kernels.
typedef double (*gyrokinetic_step2_vol_t)(const double *w, const double *dxv, const double q_, const double m_,
  const double *apardot, const double *f, double* GKYL_RESTRICT out);

typedef double (*gyrokinetic_surf_t)(const double *w, const double *dxv, 
  const double *alpha_surf_l, const double *alpha_surf_r, 
  const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
  const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
  const double *fL, const double *fC, const double *fR, double* GKYL_RESTRICT out);

typedef double (*gyrokinetic_boundary_surf_t)(const double *w, const double *dxv, 
  const double *alpha_surf_edge, const double *alpha_surf_skin, 
  const double *sgn_alpha_surf_edge, const double *sgn_alpha_surf_skin, 
  const int *const_sgn_alpha_edge, const int *const_sgn_alpha_skin, 
  const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below.
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense.
  {-1,  0,  1}, // 1x kernel indices.
  {-1, -1,  2}, // 2x kernel indices.
  {-1, -1,  3}, // 3x kernel indices.
};

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_gyrokinetic_vol_kern_list;
typedef struct { gyrokinetic_step2_vol_t kernels[3]; } gkyl_dg_gyrokinetic_step2_vol_kern_list;
typedef struct { gyrokinetic_surf_t kernels[3]; } gkyl_dg_gyrokinetic_surf_kern_list;
typedef struct { gyrokinetic_boundary_surf_t kernels[3]; } gkyl_dg_gyrokinetic_boundary_surf_kern_list;

struct dg_gyrokinetic {
  struct gkyl_dg_eqn eqn; // Base object.
  int cdim; // Config-space dimensions.
  int pdim; // Phase-space dimensions.
  gyrokinetic_step2_vol_t step2_vol; // Volume kernel.
  gyrokinetic_surf_t surf[4]; // Surface terms.
  gyrokinetic_boundary_surf_t boundary_surf[4]; // Surface terms for velocity boundary.
  struct gkyl_range conf_range; // Configuration space range.
  struct gkyl_range phase_range; // Phase space range.
  double charge, mass;
  const struct gk_geometry *gk_geom; // Pointer to geometry struct
  struct gkyl_dg_gyrokinetic_auxfields auxfields; // Auxiliary fields.
};

//
// Serendipity volume kernels general geometry
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_gyrokinetic_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  long cidx = gkyl_range_idx(&gyrokinetic->conf_range, idx);
  long pidx = gkyl_range_idx(&gyrokinetic->phase_range, idx);
  return gyrokinetic_vol_1x1v_ser_p1(xc, dx,
    gyrokinetic->charge, gyrokinetic->mass,
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->bmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->jacobtot_inv, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->cmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->b_i, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.phi, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apar, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apardot, cidx), 
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_gyrokinetic_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  long cidx = gkyl_range_idx(&gyrokinetic->conf_range, idx);
  long pidx = gkyl_range_idx(&gyrokinetic->phase_range, idx);
  return gyrokinetic_vol_1x1v_ser_p2(xc, dx,
    gyrokinetic->charge, gyrokinetic->mass,
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->bmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->jacobtot_inv, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->cmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->b_i, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.phi, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apar, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apardot, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_gyrokinetic_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  long cidx = gkyl_range_idx(&gyrokinetic->conf_range, idx);
  long pidx = gkyl_range_idx(&gyrokinetic->phase_range, idx);
  return gyrokinetic_vol_1x2v_ser_p1(xc, dx,
    gyrokinetic->charge, gyrokinetic->mass,
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->bmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->jacobtot_inv, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->cmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->b_i, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.phi, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apar, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apardot, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_gyrokinetic_vol_1x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  long cidx = gkyl_range_idx(&gyrokinetic->conf_range, idx);
  long pidx = gkyl_range_idx(&gyrokinetic->phase_range, idx);
  return gyrokinetic_vol_1x2v_ser_p2(xc, dx,
    gyrokinetic->charge, gyrokinetic->mass,
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->bmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->jacobtot_inv, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->cmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->b_i, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.phi, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apar, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apardot, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_gyrokinetic_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  long cidx = gkyl_range_idx(&gyrokinetic->conf_range, idx);
  long pidx = gkyl_range_idx(&gyrokinetic->phase_range, idx);
  return gyrokinetic_vol_2x2v_ser_p1(xc, dx,
    gyrokinetic->charge, gyrokinetic->mass,
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->bmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->jacobtot_inv, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->cmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->b_i, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.phi, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apar, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apardot, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_gyrokinetic_vol_2x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  long cidx = gkyl_range_idx(&gyrokinetic->conf_range, idx);
  long pidx = gkyl_range_idx(&gyrokinetic->phase_range, idx);
  return gyrokinetic_vol_2x2v_ser_p2(xc, dx,
    gyrokinetic->charge, gyrokinetic->mass,
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->bmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->jacobtot_inv, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->cmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->b_i, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.phi, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apar, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apardot, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_gyrokinetic_vol_3x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  long cidx = gkyl_range_idx(&gyrokinetic->conf_range, idx);
  long pidx = gkyl_range_idx(&gyrokinetic->phase_range, idx);
  return gyrokinetic_vol_3x2v_ser_p1(xc, dx,
    gyrokinetic->charge, gyrokinetic->mass,
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->bmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->jacobtot_inv, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->cmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->b_i, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.phi, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apar, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apardot, cidx),
    qIn, qRhsOut);
}

// Volume kernel list.
GKYL_CU_D
static const gkyl_dg_gyrokinetic_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_gyrokinetic_vol_1x1v_ser_p1, kernel_gyrokinetic_vol_1x1v_ser_p2 }, // 0
  { NULL, kernel_gyrokinetic_vol_1x2v_ser_p1, kernel_gyrokinetic_vol_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_gyrokinetic_vol_2x2v_ser_p1, kernel_gyrokinetic_vol_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_gyrokinetic_vol_3x2v_ser_p1, NULL }, // 3
};

// Step 2 (for electromagnetics) volume kernel list.
GKYL_CU_D
static const gkyl_dg_gyrokinetic_step2_vol_kern_list ser_step2_vol_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_step2_vol_1x1v_ser_p1, gyrokinetic_step2_vol_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_step2_vol_1x2v_ser_p1, gyrokinetic_step2_vol_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, NULL, NULL }, // 2
  // 3x kernels
  { NULL, NULL, NULL }, // 3
};

//
// Serendipity volume kernels general geometry, no toroidal field (by=0)
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_gyrokinetic_no_by_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  long cidx = gkyl_range_idx(&gyrokinetic->conf_range, idx);
  long pidx = gkyl_range_idx(&gyrokinetic->phase_range, idx);
  return gyrokinetic_no_by_vol_2x2v_ser_p1(xc, dx,
    gyrokinetic->charge, gyrokinetic->mass,
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->bmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->jacobtot_inv, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->cmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->b_i, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.phi, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apar, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apardot, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_gyrokinetic_no_by_vol_2x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  long cidx = gkyl_range_idx(&gyrokinetic->conf_range, idx);
  long pidx = gkyl_range_idx(&gyrokinetic->phase_range, idx);
  return gyrokinetic_no_by_vol_2x2v_ser_p2(xc, dx,
    gyrokinetic->charge, gyrokinetic->mass,
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->bmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->jacobtot_inv, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->cmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->b_i, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.phi, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apar, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apardot, cidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_gyrokinetic_no_by_vol_3x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  long cidx = gkyl_range_idx(&gyrokinetic->conf_range, idx);
  long pidx = gkyl_range_idx(&gyrokinetic->phase_range, idx);
  return gyrokinetic_no_by_vol_3x2v_ser_p1(xc, dx,
    gyrokinetic->charge, gyrokinetic->mass,
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->bmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->jacobtot_inv, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->cmag, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->gk_geom->b_i, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.phi, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apar, cidx),
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apardot, cidx),
    qIn, qRhsOut);
}

// Volume kernel list.
GKYL_CU_D
static const gkyl_dg_gyrokinetic_vol_kern_list ser_no_by_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_gyrokinetic_vol_1x1v_ser_p1, kernel_gyrokinetic_vol_1x1v_ser_p2 }, // 0
  { NULL, kernel_gyrokinetic_vol_1x2v_ser_p1, kernel_gyrokinetic_vol_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, kernel_gyrokinetic_no_by_vol_2x2v_ser_p1, kernel_gyrokinetic_no_by_vol_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, kernel_gyrokinetic_no_by_vol_3x2v_ser_p1, NULL }, // 3
};

//
// Serendipity surface kernels general geometry
//
// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_surf_kern_list ser_surf_x_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_surfx_1x1v_ser_p1, gyrokinetic_surfx_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_surfx_1x2v_ser_p1, gyrokinetic_surfx_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_surfx_2x2v_ser_p1, gyrokinetic_surfx_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_surfx_3x2v_ser_p1, NULL }, // 3
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_surf_kern_list ser_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  // 2x kernels
  { NULL, gyrokinetic_surfy_2x2v_ser_p1, gyrokinetic_surfy_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_surfy_3x2v_ser_p1, NULL }, // 3
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_surf_kern_list ser_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  // 2x kernels
  { NULL, NULL, NULL }, // 2
  // 3x kernels
  { NULL, gyrokinetic_surfz_3x2v_ser_p1, NULL }, // 3
};

// Acceleration surface kernel list: vpar-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_surf_kern_list ser_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_surfvpar_1x1v_ser_p1, gyrokinetic_surfvpar_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_surfvpar_1x2v_ser_p1, gyrokinetic_surfvpar_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_surfvpar_2x2v_ser_p1, gyrokinetic_surfvpar_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_surfvpar_3x2v_ser_p1, NULL }, // 3
};

// Conf-space advection boundary surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_boundary_surf_kern_list ser_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_boundary_surfx_1x1v_ser_p1, gyrokinetic_boundary_surfx_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_boundary_surfx_1x2v_ser_p1, gyrokinetic_boundary_surfx_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_boundary_surfx_2x2v_ser_p1, gyrokinetic_boundary_surfx_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_boundary_surfx_3x2v_ser_p1, NULL }, // 3
};

// Conf-space advection boundary surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_boundary_surf_kern_list ser_boundary_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  // 2x kernels
  { NULL, gyrokinetic_boundary_surfy_2x2v_ser_p1, gyrokinetic_boundary_surfy_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_boundary_surfy_3x2v_ser_p1, NULL }, // 3
};

// Conf-space advection boundary surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_boundary_surf_kern_list ser_boundary_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  // 2x kernels
  { NULL, NULL, NULL }, // 2
  // 3x kernels
  { NULL, gyrokinetic_boundary_surfz_3x2v_ser_p1, NULL }, // 3
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vpar-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_boundary_surf_kern_list ser_boundary_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_boundary_surfvpar_1x1v_ser_p1, gyrokinetic_boundary_surfvpar_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_boundary_surfvpar_1x2v_ser_p1, gyrokinetic_boundary_surfvpar_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_boundary_surfvpar_2x2v_ser_p1, gyrokinetic_boundary_surfvpar_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_boundary_surfvpar_3x2v_ser_p1, NULL }, // 3
};

//
// Serendipity surface kernels general geometry, no toroidal field (by=0)
//
// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_surf_kern_list ser_no_by_surf_x_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_surfx_1x1v_ser_p1, gyrokinetic_surfx_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_surfx_1x2v_ser_p1, gyrokinetic_surfx_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_no_by_surfx_2x2v_ser_p1, gyrokinetic_no_by_surfx_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_no_by_surfx_3x2v_ser_p1, NULL }, // 3
};

// Surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_surf_kern_list ser_no_by_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  // 2x kernels
  { NULL, gyrokinetic_no_by_surfy_2x2v_ser_p1, gyrokinetic_no_by_surfy_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_no_by_surfy_3x2v_ser_p1, NULL }, // 3
};

// Surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_surf_kern_list ser_no_by_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  // 2x kernels
  { NULL, NULL, NULL }, // 2
  // 3x kernels
  { NULL, gyrokinetic_no_by_surfz_3x2v_ser_p1, NULL }, // 3
};

// Acceleration surface kernel list: vpar-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_surf_kern_list ser_no_by_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_surfvpar_1x1v_ser_p1, gyrokinetic_surfvpar_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_surfvpar_1x2v_ser_p1, gyrokinetic_surfvpar_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_no_by_surfvpar_2x2v_ser_p1, gyrokinetic_no_by_surfvpar_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_no_by_surfvpar_3x2v_ser_p1, NULL }, // 3
};

// Conf-space advection boundary surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_boundary_surf_kern_list ser_no_by_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_boundary_surfx_1x1v_ser_p1, gyrokinetic_boundary_surfx_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_boundary_surfx_1x2v_ser_p1, gyrokinetic_boundary_surfx_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_no_by_boundary_surfx_2x2v_ser_p1, gyrokinetic_no_by_boundary_surfx_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_no_by_boundary_surfx_3x2v_ser_p1, NULL }, // 3
};

// Conf-space advection boundary surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_boundary_surf_kern_list ser_no_by_boundary_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  // 2x kernels
  { NULL, gyrokinetic_no_by_boundary_surfy_2x2v_ser_p1, gyrokinetic_no_by_boundary_surfy_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_no_by_boundary_surfy_3x2v_ser_p1, NULL }, // 3
};

// Conf-space advection boundary surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_boundary_surf_kern_list ser_no_by_boundary_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  // 2x kernels
  { NULL, NULL, NULL }, // 2
  // 3x kernels
  { NULL, gyrokinetic_no_by_boundary_surfz_3x2v_ser_p1, NULL }, // 3
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vpar-direction
GKYL_CU_D
static const gkyl_dg_gyrokinetic_boundary_surf_kern_list ser_no_by_boundary_surf_vpar_kernels[] = {
  // 1x kernels
  { NULL, gyrokinetic_boundary_surfvpar_1x1v_ser_p1, gyrokinetic_boundary_surfvpar_1x1v_ser_p2 }, // 0
  { NULL, gyrokinetic_boundary_surfvpar_1x2v_ser_p1, gyrokinetic_boundary_surfvpar_1x2v_ser_p2 }, // 1
  // 2x kernels
  { NULL, gyrokinetic_no_by_boundary_surfvpar_2x2v_ser_p1, gyrokinetic_no_by_boundary_surfvpar_2x2v_ser_p2 }, // 2
  // 3x kernels
  { NULL, gyrokinetic_no_by_boundary_surfvpar_3x2v_ser_p1, NULL }, // 3
};

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,cdim,vd,poly_order) lst[cv_index[cdim].vdim[vd]].kernels[poly_order]

/**
 * Free gyrokinetic eqn object.
 *
 * @param ref Reference counter for gyrokinetic eqn
 */
void gkyl_gyrokinetic_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
vol_step2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx,
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  long cidx = gkyl_range_idx(&gyrokinetic->conf_range, idx);
  return gyrokinetic->step2_vol(xc, dx,
    gyrokinetic->charge, gyrokinetic->mass,
    (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.apardot, cidx),
    qIn, qRhsOut);
}


GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double* xcL, const double* xcC, const double* xcR,
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  // Only in x,y,z,vpar directions.
  if (dir <= gyrokinetic->cdim) {
    // Each cell owns the *lower* edge surface alpha
    // Since alpha is continuous, fetch alpha_surf in center cell for lower edge
    // and fetch alpha_surf in right cell for upper edge
    long pidxC = gkyl_range_idx(&gyrokinetic->phase_range, idxC);
    long pidxR = gkyl_range_idx(&gyrokinetic->phase_range, idxR);
    return gyrokinetic->surf[dir](xcC, dxC, 
      (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.alpha_surf, pidxC), 
      (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.alpha_surf, pidxR), 
      (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.sgn_alpha_surf, pidxC), 
      (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.sgn_alpha_surf, pidxR), 
      (const int*) gkyl_array_cfetch(gyrokinetic->auxfields.const_sgn_alpha, pidxC), 
      (const int*) gkyl_array_cfetch(gyrokinetic->auxfields.const_sgn_alpha, pidxR), 
      qInL, qInC, qInR, qRhsOut);
  }
  return 0.;
}

GKYL_CU_D
static double
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double* xcEdge, const double* xcSkin,
  const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);

  // Only in x,y,z,vpar directions.
  if (dir <= gyrokinetic->cdim) {
    // Each cell owns the *lower* edge surface alpha
    long pidxEdge = gkyl_range_idx(&gyrokinetic->phase_range, idxEdge);
    long pidxSkin = gkyl_range_idx(&gyrokinetic->phase_range, idxSkin);
    return gyrokinetic->boundary_surf[dir](xcSkin, dxSkin, 
      (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.alpha_surf, pidxEdge), 
      (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.alpha_surf, pidxSkin), 
      (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.sgn_alpha_surf, pidxEdge), 
      (const double*) gkyl_array_cfetch(gyrokinetic->auxfields.sgn_alpha_surf, pidxSkin), 
      (const int*) gkyl_array_cfetch(gyrokinetic->auxfields.const_sgn_alpha, pidxEdge), 
      (const int*) gkyl_array_cfetch(gyrokinetic->auxfields.const_sgn_alpha, pidxSkin), 
      edge, qInEdge, qInSkin, qRhsOut);
  }
  return 0.;
}
