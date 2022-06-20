#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_kernels.h>
#include <gkyl_vlasov.h>

// Types for various kernels
typedef double (*vlasov_vol_t)(const double *w, const double *dxv,
  const double *qmem, const double *f, double* GKYL_RESTRICT out);

typedef void (*vlasov_stream_surf_t)(const double *w, const double *dxv,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef void (*vlasov_accel_surf_t)(const double *w, const double *dxv,
  const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef void (*vlasov_accel_boundary_surf_t)(const double *w, const double *dxv,
  const double *qmem, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out);

typedef void (*bc_funcf_t)(size_t nc, double *out, const double *inp, void *ctx);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { vlasov_vol_t kernels[3]; } gkyl_dg_vlasov_stream_vol_kern_list;
typedef struct { vlasov_vol_t kernels[3]; } gkyl_dg_vlasov_vol_kern_list;
typedef struct { vlasov_stream_surf_t kernels[3]; } gkyl_dg_vlasov_stream_surf_kern_list;
typedef struct { vlasov_accel_surf_t kernels[3]; } gkyl_dg_vlasov_accel_surf_kern_list;
typedef struct { vlasov_accel_boundary_surf_t kernels[3]; } gkyl_dg_vlasov_accel_boundary_surf_kern_list;

//
// Serendipity basis kernels
//

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_vlasov_stream_vol_kern_list ser_stream_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_stream_vol_1x1v_ser_p1, vlasov_stream_vol_1x1v_ser_p2 }, // 0
  { NULL, vlasov_stream_vol_1x2v_ser_p1, vlasov_stream_vol_1x2v_ser_p2 }, // 1
  { NULL, vlasov_stream_vol_1x3v_ser_p1, vlasov_stream_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_stream_vol_2x2v_ser_p1, vlasov_stream_vol_2x2v_ser_p2 }, // 3
  { NULL, vlasov_stream_vol_2x3v_ser_p1, vlasov_stream_vol_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_stream_vol_3x3v_ser_p1, NULL               }, // 5
};

GKYL_CU_D
static const gkyl_dg_vlasov_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_vol_1x1v_ser_p1, vlasov_vol_1x1v_ser_p2 }, // 0
  { NULL, vlasov_vol_1x2v_ser_p1, vlasov_vol_1x2v_ser_p2 }, // 1
  { NULL, vlasov_vol_1x3v_ser_p1, vlasov_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_vol_2x2v_ser_p1, vlasov_vol_2x2v_ser_p2 }, // 3
  { NULL, vlasov_vol_2x3v_ser_p1, vlasov_vol_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_vol_3x3v_ser_p1, NULL               }, // 5
};

// Streaming surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ser_stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfx_1x1v_ser_p1, vlasov_surfx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_surfx_1x2v_ser_p1, vlasov_surfx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_surfx_1x3v_ser_p1, vlasov_surfx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surfx_2x2v_ser_p1, vlasov_surfx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_surfx_2x3v_ser_p1, vlasov_surfx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfx_3x3v_ser_p1, NULL                  }, // 5
};

// Streaming surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ser_stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, vlasov_surfy_2x2v_ser_p1, vlasov_surfy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_surfy_2x3v_ser_p1, vlasov_surfy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfy_3x3v_ser_p1, NULL                  }, // 5
};

// Streaming surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ser_stream_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_surfz_3x3v_ser_p1, NULL }, // 5
};

// Acceleration surface kernel list: vx-direction
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ser_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfvx_1x1v_ser_p1, vlasov_surfvx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_surfvx_1x2v_ser_p1, vlasov_surfvx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_surfvx_1x3v_ser_p1, vlasov_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surfvx_2x2v_ser_p1, vlasov_surfvx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_surfvx_2x3v_ser_p1, vlasov_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ser_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_surfvy_1x2v_ser_p1, vlasov_surfvy_1x2v_ser_p2 }, // 1
  { NULL, vlasov_surfvy_1x3v_ser_p1, vlasov_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surfvy_2x2v_ser_p1, vlasov_surfvy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_surfvy_2x3v_ser_p1, vlasov_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ser_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_surfvz_1x3v_ser_p1, vlasov_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_surfvz_2x3v_ser_p1, vlasov_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_surfvz_3x3v_ser_p1, NULL }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vx-direction
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ser_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_surfvx_1x1v_ser_p1, vlasov_boundary_surfvx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_boundary_surfvx_1x2v_ser_p1, vlasov_boundary_surfvx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_boundary_surfvx_1x3v_ser_p1, vlasov_boundary_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfvx_2x2v_ser_p1, vlasov_boundary_surfvx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_boundary_surfvx_2x3v_ser_p1, vlasov_boundary_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ser_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_boundary_surfvy_1x2v_ser_p1, vlasov_boundary_surfvy_1x2v_ser_p2 }, // 1
  { NULL, vlasov_boundary_surfvy_1x3v_ser_p1, vlasov_boundary_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfvy_2x2v_ser_p1, vlasov_boundary_surfvy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_boundary_surfvy_2x3v_ser_p1, vlasov_boundary_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ser_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_boundary_surfvz_1x3v_ser_p1, vlasov_boundary_surfvz_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_boundary_surfvz_2x3v_ser_p1, vlasov_boundary_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfvz_3x3v_ser_p1, NULL }, // 5
};

//
// Tensor-product basis kernels
//

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_vlasov_stream_vol_kern_list ten_stream_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_stream_vol_1x1v_ser_p1, vlasov_stream_vol_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_stream_vol_1x2v_ser_p1, vlasov_stream_vol_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_stream_vol_1x3v_ser_p1, vlasov_stream_vol_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_stream_vol_2x2v_ser_p1, vlasov_stream_vol_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_stream_vol_2x3v_ser_p1, NULL               }, // 4
  // 3x kernels
  { NULL, vlasov_stream_vol_3x3v_ser_p1, NULL               }, // 5
};

GKYL_CU_D
static const gkyl_dg_vlasov_vol_kern_list ten_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_vol_1x1v_ser_p1, vlasov_vol_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_vol_1x2v_ser_p1, vlasov_vol_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_vol_1x3v_ser_p1, vlasov_vol_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_vol_2x2v_ser_p1, vlasov_vol_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_vol_2x3v_ser_p1, NULL               }, // 4
  // 3x kernels
  { NULL, vlasov_vol_3x3v_ser_p1, NULL               }, // 5
};

// Streaming surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ten_stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfx_1x1v_ser_p1, vlasov_surfx_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_surfx_1x2v_ser_p1, vlasov_surfx_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_surfx_1x3v_ser_p1, vlasov_surfx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surfx_2x2v_ser_p1, vlasov_surfx_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_surfx_2x3v_ser_p1, NULL                  }, // 4
  // 3x kernels
  { NULL, vlasov_surfx_3x3v_ser_p1, NULL                  }, // 5
};

// Streaming surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ten_stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, vlasov_surfy_2x2v_ser_p1, vlasov_surfy_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_surfy_2x3v_ser_p1, NULL                  }, // 4
  // 3x kernels
  { NULL, vlasov_surfy_3x3v_ser_p1, NULL                  }, // 5
};

// Streaming surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_vlasov_stream_surf_kern_list ten_stream_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_surfz_3x3v_ser_p1, NULL }, // 5
};

// Acceleration surface kernel list: vx-direction
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ten_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfvx_1x1v_ser_p1, vlasov_surfvx_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_surfvx_1x2v_ser_p1, vlasov_surfvx_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_surfvx_1x3v_ser_p1, vlasov_surfvx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surfvx_2x2v_ser_p1, vlasov_surfvx_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_surfvx_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ten_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_surfvy_1x2v_ser_p1, vlasov_surfvy_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_surfvy_1x3v_ser_p1, vlasov_surfvy_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surfvy_2x2v_ser_p1, vlasov_surfvy_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_surfvy_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_accel_surf_kern_list ten_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_surfvz_1x3v_ser_p1, vlasov_surfvz_1x3v_tensor_p2}, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_surfvz_2x3v_ser_p1, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_surfvz_3x3v_ser_p1, NULL }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vx-direction
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ten_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_surfvx_1x1v_ser_p1, vlasov_boundary_surfvx_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_boundary_surfvx_1x2v_ser_p1, vlasov_boundary_surfvx_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_boundary_surfvx_1x3v_ser_p1, vlasov_boundary_surfvx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfvx_2x2v_ser_p1, vlasov_boundary_surfvx_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_boundary_surfvx_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ten_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, vlasov_boundary_surfvy_1x2v_ser_p1, vlasov_boundary_surfvy_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_boundary_surfvy_1x3v_ser_p1, vlasov_boundary_surfvy_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_surfvy_2x2v_ser_p1, vlasov_boundary_surfvy_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_boundary_surfvy_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_accel_boundary_surf_kern_list ten_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, vlasov_boundary_surfvz_1x3v_ser_p1, vlasov_boundary_surfvz_1x3v_tensor_p2}, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, vlasov_boundary_surfvz_2x3v_ser_p1, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_surfvz_3x3v_ser_p1, NULL }, // 5
};

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,cdim,vd,poly_order) lst[cv_index[cdim].vdim[vd]].kernels[poly_order]

struct dg_vlasov {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  vlasov_vol_t vol; // Volume kernel
  vlasov_stream_surf_t stream_surf[3]; // Surface terms for streaming
  vlasov_accel_surf_t accel_surf[3]; // Surface terms for acceleration
  vlasov_accel_boundary_surf_t accel_boundary_surf[3]; // Surface terms for acceleration
  bc_funcf_t wall_bc; // wall BCs function
  bc_funcf_t absorb_bc; // Absorbing BCs function
  struct gkyl_range conf_range; // configuration space range
  struct gkyl_dg_vlasov_auxfields auxfields; // Auxiliary fields.
};

/**
 * Free vlasov eqn object.
 *
 * @param ref Reference counter for vlasov eqn
 */
void gkyl_vlasov_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  long cidx = gkyl_range_idx(&vlasov->conf_range, idx);
  return vlasov->vol(xc, dx,
    vlasov->auxfields.qmem ? (const double*) gkyl_array_cfetch(vlasov->auxfields.qmem, cidx) : 0,
    qIn, qRhsOut);
}

GKYL_CU_D
static void
surf(const struct gkyl_dg_eqn *eqn, 
  int dir,
  const double*  xcL, const double*  xcC, const double*  xcR, 
  const double*  dxL, const double* dxC, const double* dxR,
  const int*  idxL, const int*  idxC, const int*  idxR,
  const double* qInL, const double*  qInC, const double*  qInR, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  if (dir < vlasov->cdim) {
    vlasov->stream_surf[dir]
      (xcC, dxC, qInL, qInC, qInR, qRhsOut);
  }
  else {
    long cidx = gkyl_range_idx(&vlasov->conf_range, idxC);
    vlasov->accel_surf[dir-vlasov->cdim]
      (xcC, dxC,
        vlasov->auxfields.qmem ? (const double*) gkyl_array_cfetch(vlasov->auxfields.qmem, cidx) : 0,
        qInL, qInC, qInR, qRhsOut);
  }
}

GKYL_CU_D
static void
boundary_surf(const struct gkyl_dg_eqn *eqn,
  int dir,
  const double*  xcEdge, const double*  xcSkin,
  const double*  dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);

  if (dir >= vlasov->cdim) {
    long cidx = gkyl_range_idx(&vlasov->conf_range, idxSkin);
    vlasov->accel_boundary_surf[dir-vlasov->cdim]
      (xcSkin, dxSkin,
        vlasov->auxfields.qmem ? (const double*) gkyl_array_cfetch(vlasov->auxfields.qmem, cidx) : 0,
        edge, qInEdge, qInSkin, qRhsOut);
  }
}

GKYL_CU_D
static void
species_wall_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int dir = mc->dir, cdim = mc->cdim;

  mc->basis->flip_odd_sign(dir, inp, out);
  mc->basis->flip_odd_sign(dir+cdim, out, out);
}

GKYL_CU_D
static void
species_absorb_bc(size_t nc, double *out, const double *inp, void *ctx)
{
  struct dg_bc_ctx *mc = (struct dg_bc_ctx*) ctx;
  int nbasis = mc->basis->num_basis;
  for (int c=0; c<nbasis; ++c) out[c] = 0.0;
}
