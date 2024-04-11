#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_canonical_pb_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

// Types for various kernels
typedef double (*canonical_pb_stream_surf_t)(const double *w, const double *dxv,
  const double *hamil, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*canonical_pb_accel_surf_t)(const double *w, const double *dxv,
  const double *hamil, const double *qmem, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef double (*canonical_pb_accel_boundary_surf_t)(const double *w, const double *dxv,
  const double *hamil, const double *qmem, 
  const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_canonical_pb_stream_vol_kern_list;
typedef struct { vol_termf_t kernels[3]; } gkyl_dg_canonical_pb_vol_kern_list;
typedef struct { canonical_pb_stream_surf_t kernels[3]; } gkyl_dg_canonical_pb_stream_surf_kern_list;
typedef struct { canonical_pb_accel_surf_t kernels[3]; } gkyl_dg_canonical_pb_accel_surf_kern_list;
typedef struct { canonical_pb_accel_boundary_surf_t kernels[3]; } gkyl_dg_canonical_pb_accel_boundary_surf_kern_list;

struct dg_canonical_pb {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  canonical_pb_stream_surf_t stream_surf[3]; // Surface terms for streaming
  canonical_pb_accel_surf_t accel_surf[3]; // Surface terms for acceleration
  canonical_pb_accel_boundary_surf_t accel_boundary_surf[3]; // Surface terms for acceleration
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_range vel_range; // velocity space range
  struct gkyl_range phase_range; // velocity space range
  struct gkyl_dg_canonical_pb_auxfields auxfields; // Auxiliary fields.
};

//
// Serendipity volume kernels
// Need to be separated like this for GPU build
//

GKYL_CU_DH
static double
kernel_canonical_pb_vol_1x1v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_1x1v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_1x1v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_1x1v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_2x2v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_2x2v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_2x2v_ser_p2(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_2x2v_ser_p2(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

GKYL_CU_DH
static double
kernel_canonical_pb_vol_3x3v_ser_p1(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  return canonical_pb_vol_3x3v_ser_p1(xc, dx,
    (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
    qIn, qRhsOut);  
}

// Volume kernel list for relativistic streaming + EM
GKYL_CU_D
static const gkyl_dg_canonical_pb_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, kernel_canonical_pb_vol_1x1v_ser_p1, kernel_canonical_pb_vol_1x1v_ser_p2 }, // 0 (cdim)
  // 2x kernels
  { NULL, kernel_canonical_pb_vol_2x2v_ser_p1, kernel_canonical_pb_vol_2x2v_ser_p2 }, // 1
  // 3x kernels
  { NULL, kernel_canonical_pb_vol_3x3v_ser_p1, NULL               }, // 2
};

//
// Serendipity surface kernels
//

// Streaming surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_surf_kern_list ser_stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_surfx_1x1v_ser_p1, canonical_pb_surfx_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, canonical_pb_surfx_2x2v_ser_p1, canonical_pb_surfx_2x2v_ser_p2 }, // 1
  // 3x kernels
  { NULL, canonical_pb_surfx_3x3v_ser_p1, NULL                  }, // 2
};

// Streaming surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_surf_kern_list ser_stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_surfy_2x2v_ser_p1, canonical_pb_surfy_2x2v_ser_p2 }, // 1
  // 3x kernels
  { NULL, canonical_pb_surfy_3x3v_ser_p1, NULL                  }, // 2
};

// Streaming surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_stream_surf_kern_list ser_stream_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, NULL, NULL }, // 1
  // 3x kernels
  { NULL, canonical_pb_surfz_3x3v_ser_p1, NULL }, // 2
};

// Acceleration surface kernel list: vx-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_surf_kern_list ser_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_surfvx_1x1v_ser_p1, canonical_pb_surfvx_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, canonical_pb_surfvx_2x2v_ser_p1, canonical_pb_surfvx_2x2v_ser_p2 }, // 1
  // 3x kernels
  { NULL, canonical_pb_surfvx_3x3v_ser_p1, NULL                   }, // 2
};

// Acceleration surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_surf_kern_list ser_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_surfvy_2x2v_ser_p1, canonical_pb_surfvy_2x2v_ser_p2 }, // 1
  // 3x kernels
  { NULL, canonical_pb_surfvy_3x3v_ser_p1, NULL                   }, // 2
};

// Acceleration surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_surf_kern_list ser_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, NULL, NULL }, // 1
  // 3x kernels
  { NULL, canonical_pb_surfvz_3x3v_ser_p1, NULL }, // 2
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vx-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_boundary_surf_kern_list ser_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, canonical_pb_boundary_surfvx_1x1v_ser_p1, canonical_pb_boundary_surfvx_1x1v_ser_p2 }, // 0
  // 2x kernels
  { NULL, canonical_pb_boundary_surfvx_2x2v_ser_p1, canonical_pb_boundary_surfvx_2x2v_ser_p2 }, // 1
  // 3x kernels
  { NULL, canonical_pb_boundary_surfvx_3x3v_ser_p1, NULL                   }, // 2
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vy-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_boundary_surf_kern_list ser_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, canonical_pb_boundary_surfvy_2x2v_ser_p1, canonical_pb_boundary_surfvy_2x2v_ser_p2 }, // 1
  // 3x kernels
  { NULL, canonical_pb_boundary_surfvy_3x3v_ser_p1, NULL                   }, // 2
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vz-direction
GKYL_CU_D
static const gkyl_dg_canonical_pb_accel_boundary_surf_kern_list ser_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  // 2x kernels
  { NULL, NULL, NULL }, // 1
  // 3x kernels
  { NULL, canonical_pb_boundary_surfvz_3x3v_ser_p1, NULL }, // 2
};

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim].kernels[poly_order]

/**
 * Free vlasov eqn object.
 *
 * @param ref Reference counter for vlasov eqn
 */
void gkyl_canonical_pb_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);

  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<canonical_pb->pdim-canonical_pb->cdim; ++i)
    idx_vel[i] = idxC[canonical_pb->cdim+i];
  long vidx = gkyl_range_idx(&canonical_pb->vel_range, idx_vel);

  // TODO: Fix idx, remove qmem

  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);

  // Iterate over phase space
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx_vel);

  if (dir < canonical_pb->cdim) {
    return canonical_pb->stream_surf[dir]
      (xcC, dxC, 
        (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx), 
        qInL, qInC, qInR, qRhsOut);
  }
  else {
    long cidx = gkyl_range_idx(&canonical_pb->conf_range, idxC);
    return canonical_pb->accel_surf[dir-canonical_pb->cdim]
      (xcC, dxC,
        (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx),
        canonical_pb->auxfields.qmem ? (const double*) gkyl_array_cfetch(canonical_pb->auxfields.qmem, cidx) : 0,
        qInL, qInC, qInR, qRhsOut);
  }
}

GKYL_CU_D
static double
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_canonical_pb *canonical_pb = container_of(eqn, struct dg_canonical_pb, eqn);
  long pidx = gkyl_range_idx(&canonical_pb->phase_range, idx);  

  if (dir >= canonical_pb->cdim) {
    long cidx = gkyl_range_idx(&canonical_pb->conf_range, idxSkin);
    return canonical_pb->accel_boundary_surf[dir-canonical_pb->cdim]
      (xcSkin, dxSkin,
        (const double*) gkyl_array_cfetch(canonical_pb->auxfields.hamil, pidx),
        canonical_pb->auxfields.qmem ? (const double*) gkyl_array_cfetch(canonical_pb->auxfields.qmem, cidx) : 0,
        edge, qInEdge, qInSkin, qRhsOut);
  }
  return 0.;
}
