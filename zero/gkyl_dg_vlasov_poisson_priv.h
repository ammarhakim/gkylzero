#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_dg_eqn.h>
#include <gkyl_eqn_type.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_kernels.h>

// Types for various kernels
typedef double (*vlasov_poisson_vol_t)(const double *w, const double *dxv,
  const double *fac_phi, const double *vecA, const double *f, double* GKYL_RESTRICT out);

typedef void (*vlasov_poisson_stream_surf_t)(const double *w, const double *dxv,
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef void (*vlasov_poisson_accel_surf_t)(const double *w, const double *dxv,
  const double *fac_phi, const double *vecA, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out);

typedef void (*vlasov_poisson_accel_boundary_surf_t)(const double *w, const double *dxv,
  const double *fac_phi, const double *vecA, const int edge, const double *fEdge, const double *fSkin, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { vlasov_poisson_vol_t kernels[3]; } gkyl_dg_vlasov_poisson_vol_kern_list;
typedef struct { vlasov_poisson_stream_surf_t kernels[3]; } gkyl_dg_vlasov_poisson_stream_surf_kern_list;
typedef struct { vlasov_poisson_accel_surf_t kernels[3]; } gkyl_dg_vlasov_poisson_accel_surf_kern_list;
typedef struct { vlasov_poisson_accel_boundary_surf_t kernels[3]; } gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list;

//
// Serendipity basis kernels
//

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_vol_kern_list ser_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_poisson_vol_1x1v_ser_p1, vlasov_poisson_vol_1x1v_ser_p2 }, // 0
  { NULL, vlasov_poisson_vol_1x2v_ser_p1, vlasov_poisson_vol_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_vol_1x3v_ser_p1, vlasov_poisson_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_vol_2x2v_ser_p1, vlasov_poisson_vol_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_vol_2x3v_ser_p1, NULL               }, // 4
  // 3x kernels
  { NULL, vlasov_poisson_vol_3x3v_ser_p1, NULL               }, // 5
};

// Streaming surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_stream_surf_kern_list ser_stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, vlasov_surfx_1x1v_ser_p1, vlasov_surfx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_surfx_1x2v_ser_p1, vlasov_surfx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_surfx_1x3v_ser_p1, vlasov_surfx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_surfx_2x2v_ser_p1, vlasov_surfx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_surfx_2x3v_ser_p1, NULL                  }, // 4
  // 3x kernels
  { NULL, vlasov_surfx_3x3v_ser_p1, NULL                  }, // 5
};

// Streaming surface kernel list: y-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_stream_surf_kern_list ser_stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, vlasov_surfy_2x2v_ser_p1, vlasov_surfy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_surfy_2x3v_ser_p1, NULL                  }, // 4
  // 3x kernels
  { NULL, vlasov_surfy_3x3v_ser_p1, NULL                  }, // 5
};

// Streaming surface kernel list: z-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_stream_surf_kern_list ser_stream_surf_z_kernels[] = {
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
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ser_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_poisson_surfvx_1x1v_ser_p1, vlasov_poisson_surfvx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_poisson_surfvx_1x2v_ser_p1, vlasov_poisson_surfvx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_surfvx_1x3v_ser_p1, vlasov_poisson_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_surfvx_2x2v_ser_p1, vlasov_poisson_surfvx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_surfvx_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_poisson_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ser_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_surfvy_2x2v_ser_p1, vlasov_poisson_surfvy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_surfvy_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_poisson_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ser_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_poisson_surfvz_3x3v_ser_p1, NULL }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vx-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ser_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_poisson_boundary_surfvx_1x1v_ser_p1, vlasov_poisson_boundary_surfvx_1x1v_ser_p2 }, // 0
  { NULL, vlasov_poisson_boundary_surfvx_1x2v_ser_p1, vlasov_poisson_boundary_surfvx_1x2v_ser_p2 }, // 1
  { NULL, vlasov_poisson_boundary_surfvx_1x3v_ser_p1, vlasov_poisson_boundary_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_boundary_surfvx_2x2v_ser_p1, vlasov_poisson_boundary_surfvx_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_boundary_surfvx_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_poisson_boundary_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ser_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_boundary_surfvy_2x2v_ser_p1, vlasov_poisson_boundary_surfvy_2x2v_ser_p2 }, // 3
  { NULL, vlasov_poisson_boundary_surfvy_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_poisson_boundary_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ser_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_poisson_boundary_surfvz_3x3v_ser_p1, NULL }, // 5
};

//
// Tensor-product basis kernels
//

// Volume kernel list
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_vol_kern_list ten_vol_kernels[] = {
  // 1x kernels
  { NULL, vlasov_poisson_vol_1x1v_ser_p1, vlasov_poisson_vol_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_poisson_vol_1x2v_ser_p1, vlasov_poisson_vol_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_poisson_vol_1x3v_ser_p1, vlasov_poisson_vol_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_vol_2x2v_ser_p1, vlasov_poisson_vol_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_poisson_vol_2x3v_ser_p1, NULL               }, // 4
  // 3x kernels
  { NULL, vlasov_poisson_vol_3x3v_ser_p1, NULL               }, // 5
};

// Streaming surface kernel list: x-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_stream_surf_kern_list ten_stream_surf_x_kernels[] = {
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
static const gkyl_dg_vlasov_poisson_stream_surf_kern_list ten_stream_surf_y_kernels[] = {
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
static const gkyl_dg_vlasov_poisson_stream_surf_kern_list ten_stream_surf_z_kernels[] = {
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
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ten_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_poisson_surfvx_1x1v_ser_p1, vlasov_poisson_surfvx_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_poisson_surfvx_1x2v_ser_p1, vlasov_poisson_surfvx_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_poisson_surfvx_1x3v_ser_p1, vlasov_poisson_surfvx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_surfvx_2x2v_ser_p1, vlasov_poisson_surfvx_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_poisson_surfvx_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_poisson_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration surface kernel list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ten_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_surfvy_2x2v_ser_p1, vlasov_poisson_surfvy_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_poisson_surfvy_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_poisson_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration surface kernel list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_surf_kern_list ten_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_poisson_surfvz_3x3v_ser_p1, NULL }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vx-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ten_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, vlasov_poisson_boundary_surfvx_1x1v_ser_p1, vlasov_poisson_boundary_surfvx_1x1v_tensor_p2 }, // 0
  { NULL, vlasov_poisson_boundary_surfvx_1x2v_ser_p1, vlasov_poisson_boundary_surfvx_1x2v_tensor_p2 }, // 1
  { NULL, vlasov_poisson_boundary_surfvx_1x3v_ser_p1, vlasov_poisson_boundary_surfvx_1x3v_tensor_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_boundary_surfvx_2x2v_ser_p1, vlasov_poisson_boundary_surfvx_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_poisson_boundary_surfvx_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_poisson_boundary_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vy-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ten_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, vlasov_poisson_boundary_surfvy_2x2v_ser_p1, vlasov_poisson_boundary_surfvy_2x2v_tensor_p2 }, // 3
  { NULL, vlasov_poisson_boundary_surfvy_2x3v_ser_p1, NULL                   }, // 4
  // 3x kernels
  { NULL, vlasov_poisson_boundary_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vz-direction
GKYL_CU_D
static const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list ten_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_poisson_boundary_surfvz_3x3v_ser_p1, NULL }, // 5
};

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,cdim,vd,poly_order) lst[cv_index[cdim].vdim[vd]].kernels[poly_order]

struct dg_vlasov_poisson {
  struct gkyl_dg_eqn eqn; // Base object
  int cdim; // Config-space dimensions
  int pdim; // Phase-space dimensions
  vlasov_poisson_vol_t vol; // Volume kernel
  vlasov_poisson_stream_surf_t stream_surf[3]; // Surface terms for streaming
  vlasov_poisson_accel_surf_t accel_surf[3]; // Surface terms for acceleration
  vlasov_poisson_accel_boundary_surf_t accel_boundary_surf[3]; // Surface terms for acceleration
  struct gkyl_range conf_range; // configuration space range
  const struct gkyl_array *fac_phi; // Pointer to fac*phi, where phi is the potential
  const struct gkyl_array *vecA; // Pointer to q/m*A, where A is the vector potential
};

GKYL_CU_D
static double
vol(const struct gkyl_dg_eqn *eqn, const double*  xc, const double*  dx, 
  const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_vlasov_poisson *vlasov_poisson = container_of(eqn, struct dg_vlasov_poisson, eqn);

  long cidx = gkyl_range_idx(&vlasov_poisson->conf_range, idx);
  return vlasov_poisson->vol(xc, dx, (const double*) gkyl_array_cfetch(vlasov_poisson->fac_phi, cidx), NULL,
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
  struct dg_vlasov_poisson *vlasov_poisson = container_of(eqn, struct dg_vlasov_poisson, eqn);

  if (dir < vlasov_poisson->cdim) {
    vlasov_poisson->stream_surf[dir]
      (xcC, dxC, qInL, qInC, qInR, qRhsOut);
  }
  else {
    long cidx = gkyl_range_idx(&vlasov_poisson->conf_range, idxC);
    vlasov_poisson->accel_surf[dir-vlasov_poisson->cdim]
      (xcC, dxC, (const double*) gkyl_array_cfetch(vlasov_poisson->fac_phi, cidx), NULL,
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
  struct dg_vlasov_poisson *vlasov_poisson = container_of(eqn, struct dg_vlasov_poisson, eqn);

  if (dir >= vlasov_poisson->cdim) {
    long cidx = gkyl_range_idx(&vlasov_poisson->conf_range, idxSkin);
    vlasov_poisson->accel_boundary_surf[dir-vlasov_poisson->cdim]
      (xcSkin, dxSkin, (const double*) gkyl_array_cfetch(vlasov_poisson->fac_phi, cidx), NULL,
        edge, qInEdge, qInSkin, qRhsOut);
  }
}

