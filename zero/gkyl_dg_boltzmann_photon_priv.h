#pragma once

// Private header, not for direct use in user code.

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_boltzmann_photon_kernels.h>
#include <gkyl_dg_boltzmann_photon.h>

// Types for various kernels.
typedef double (*boltzmann_photon_surf_t)(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*boltzmann_photon_boundary_surf_t)(const double *w, const double *dxv,  double light_speed, double rho_curv, 
  const double *jacob_vel_inv, const double *kpar_abs, 
  const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_boltzmann_photon_vol_kern_list;
typedef struct { boltzmann_photon_surf_t kernels[3]; } gkyl_dg_boltzmann_photon_surf_kern_list;
typedef struct { boltzmann_photon_boundary_surf_t kernels[3]; } gkyl_dg_boltzmann_photon_boundary_surf_kern_list;

struct dg_boltzmann_photon {
  struct gkyl_dg_eqn eqn; // Base object.
  int cdim; // Config-space dimensions.
  int pdim; // Phase-space dimensions.
  boltzmann_photon_surf_t surf[3]; // Surface terms.
  boltzmann_photon_boundary_surf_t boundary_surf[3]; // Boundary surface terms.
  double light_speed; // Speed of light. 
  double rho_curv; // Radius of curvature for magnetic field photons are propagating along. 
  struct gkyl_range vel_range; // Velocity-space range.
  struct gkyl_dg_boltzmann_photon_auxfields auxfields; // Auxiliary fields.
};

//
// Serendipity volume kernels with potentially mapped momentum-space grid
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_boltzmann_photon_vol_1x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_boltzmann_photon *boltzmann_photon = container_of(eqn, struct dg_boltzmann_photon, eqn);

  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<boltzmann_photon->pdim-boltzmann_photon->cdim; ++i)
    idx_vel[i] = idx[boltzmann_photon->cdim+i];

  long vidx = gkyl_range_idx(&boltzmann_photon->vel_range, idx_vel);

  return boltzmann_photon_vol_1x2v_ser_p1(xc, dx, 
    boltzmann_photon->light_speed, boltzmann_photon->rho_curv, 
    (const double*) gkyl_array_cfetch( boltzmann_photon->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch( boltzmann_photon->auxfields.kpar_abs, vidx),
    qIn, qRhsOut);
}

GKYL_CU_DH
static double
kernel_boltzmann_photon_vol_1x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_boltzmann_photon *boltzmann_photon = container_of(eqn, struct dg_boltzmann_photon, eqn);

  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<boltzmann_photon->pdim-boltzmann_photon->cdim; ++i)
    idx_vel[i] = idx[boltzmann_photon->cdim+i];

  long vidx = gkyl_range_idx(&boltzmann_photon->vel_range, idx_vel);

  return boltzmann_photon_vol_1x2v_ser_p2(xc, dx, 
    boltzmann_photon->light_speed, boltzmann_photon->rho_curv, 
    (const double*) gkyl_array_cfetch( boltzmann_photon->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch( boltzmann_photon->auxfields.kpar_abs, vidx),
    qIn, qRhsOut);
}

// Volume kernel list with potentially mapped momentum-space grid (Serendipity basis)
GKYL_CU_D
static const gkyl_dg_boltzmann_photon_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, kernel_boltzmann_photon_vol_1x2v_ser_p1, kernel_boltzmann_photon_vol_1x2v_ser_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

GKYL_CU_DH
static double
kernel_boltzmann_photon_vol_1x2v_tensor_p2(const struct gkyl_dg_eqn *eqn, const double* xc, const double* dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_boltzmann_photon *boltzmann_photon = container_of(eqn, struct dg_boltzmann_photon, eqn);

  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<boltzmann_photon->pdim-boltzmann_photon->cdim; ++i)
    idx_vel[i] = idx[boltzmann_photon->cdim+i];

  long vidx = gkyl_range_idx(&boltzmann_photon->vel_range, idx_vel);

  return boltzmann_photon_vol_1x2v_tensor_p2(xc, dx, 
    boltzmann_photon->light_speed, boltzmann_photon->rho_curv, 
    (const double*) gkyl_array_cfetch( boltzmann_photon->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch( boltzmann_photon->auxfields.kpar_abs, vidx),
    qIn, qRhsOut);
}

// Volume kernel list with potentially mapped momentum-space grid (Tensor basis)
GKYL_CU_D
static const gkyl_dg_boltzmann_photon_vol_kern_list tensor_vol_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, kernel_boltzmann_photon_vol_1x2v_tensor_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

//
// Serendipity surface kernels with potentially mapped momentum-space grids
//

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_boltzmann_photon_surf_kern_list ser_surf_x_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, boltzmann_photon_surfx_1x2v_ser_p1, boltzmann_photon_surfx_1x2v_ser_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// Surface kernel list: vx-direction (kpar direction)
GKYL_CU_D
static const gkyl_dg_boltzmann_photon_surf_kern_list ser_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, boltzmann_photon_surfvx_1x2v_ser_p1, boltzmann_photon_surfvx_1x2v_ser_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// Surface kernel list: vy-direction (kperp direction)
GKYL_CU_D
static const gkyl_dg_boltzmann_photon_surf_kern_list ser_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, boltzmann_photon_surfvy_1x2v_ser_p1, boltzmann_photon_surfvy_1x2v_ser_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

//
// Tensor surface kernels with potentially mapped momentum-space grids
//

// Surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_boltzmann_photon_surf_kern_list tensor_surf_x_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, boltzmann_photon_surfx_1x2v_tensor_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// Surface kernel list: vx-direction (kpar direction)
GKYL_CU_D
static const gkyl_dg_boltzmann_photon_surf_kern_list tensor_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, boltzmann_photon_surfvx_1x2v_tensor_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// Surface kernel list: vy-direction (kperp direction)
GKYL_CU_D
static const gkyl_dg_boltzmann_photon_surf_kern_list tensor_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, boltzmann_photon_surfvy_1x2v_tensor_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

//
// Serendipity boundary surface kernels with potentially mapped momentum-space grids
//

// Boundary surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_boltzmann_photon_boundary_surf_kern_list ser_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, boltzmann_photon_boundary_surfx_1x2v_ser_p1, boltzmann_photon_boundary_surfx_1x2v_ser_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// Boundary surface kernel list: vx-direction (kpar direction)
GKYL_CU_D
static const gkyl_dg_boltzmann_photon_boundary_surf_kern_list ser_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, boltzmann_photon_boundary_surfvx_1x2v_ser_p1, boltzmann_photon_boundary_surfvx_1x2v_ser_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// Boundary surface kernel list: vy-direction (kperp direction)
GKYL_CU_D
static const gkyl_dg_boltzmann_photon_boundary_surf_kern_list ser_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, boltzmann_photon_boundary_surfvy_1x2v_ser_p1, boltzmann_photon_boundary_surfvy_1x2v_ser_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

//
// Tensor boundary surface kernels with potentially mapped momentum-space grids
//

// Boundary surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_boltzmann_photon_boundary_surf_kern_list tensor_boundary_surf_x_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, boltzmann_photon_boundary_surfx_1x2v_tensor_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// Boundary surface kernel list: vx-direction (kpar direction)
GKYL_CU_D
static const gkyl_dg_boltzmann_photon_boundary_surf_kern_list tensor_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, boltzmann_photon_boundary_surfvx_1x2v_tensor_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// Boundary surface kernel list: vy-direction (kperp direction)
GKYL_CU_D
static const gkyl_dg_boltzmann_photon_boundary_surf_kern_list tensor_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, boltzmann_photon_boundary_surfvy_1x2v_tensor_p2 }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,cdim,vd,poly_order) lst[cv_index[cdim].vdim[vd]].kernels[poly_order]

/**
 * Free Boltzmann photon eqn object.
 *
 * @param ref Reference counter for Boltzmann photon eqn
 */
void gkyl_boltzmann_photon_free(const struct gkyl_ref_count *ref);


GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double* xcL, const double* xcC, const double* xcR,
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_boltzmann_photon *boltzmann_photon = container_of(eqn, struct dg_boltzmann_photon, eqn);

  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<boltzmann_photon->pdim-boltzmann_photon->cdim; ++i)
    idx_vel[i] = idxC[boltzmann_photon->cdim+i];

  long vidx = gkyl_range_idx(&boltzmann_photon->vel_range, idx_vel);

  return boltzmann_photon->surf[dir](xcC, dxC, 
    boltzmann_photon->light_speed, boltzmann_photon->rho_curv, 
    (const double*) gkyl_array_cfetch(boltzmann_photon->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch( boltzmann_photon->auxfields.kpar_abs, vidx),
    qInL, qInC, qInR, qRhsOut);
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
  struct dg_boltzmann_photon *boltzmann_photon = container_of(eqn, struct dg_boltzmann_photon, eqn);

  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<boltzmann_photon->pdim-boltzmann_photon->cdim; ++i)
    idx_vel[i] = idxSkin[boltzmann_photon->cdim+i];

  long vidx = gkyl_range_idx(&boltzmann_photon->vel_range, idx_vel);

  return boltzmann_photon->boundary_surf[dir](xcSkin, dxSkin, 
    boltzmann_photon->light_speed, boltzmann_photon->rho_curv, 
    (const double*) gkyl_array_cfetch(boltzmann_photon->auxfields.jacob_vel_inv, vidx),
    (const double*) gkyl_array_cfetch( boltzmann_photon->auxfields.kpar_abs, vidx),
    edge, qInEdge, qInSkin, qRhsOut);
}
