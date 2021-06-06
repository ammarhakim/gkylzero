/* -*- c -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_dg_vlasov.h>    
#include <gkyl_dg_vlasov_priv.h>
}

// various pointers to functions on device
__device__ static vol_termf_t d_vol = &vol;
__device__ static surf_termf_t d_surf = &surf;
__device__ static boundary_surf_termf_t d_boundary_surf = &boundary_surf;

// Volume kernel list
__device__ static struct { vlasov_stream_vol_t kernels[3]; } d_stream_vol_kernels[] = {
  // 1x kernels
  { NULL, &vlasov_stream_vol_1x1v_ser_p1, &vlasov_stream_vol_1x1v_ser_p2 }, // 0
  { NULL, &vlasov_stream_vol_1x2v_ser_p1, &vlasov_stream_vol_1x2v_ser_p2 }, // 1
  { NULL, &vlasov_stream_vol_1x3v_ser_p1, &vlasov_stream_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, &vlasov_stream_vol_2x2v_ser_p1, &vlasov_stream_vol_2x2v_ser_p2 }, // 3
  { NULL, &vlasov_stream_vol_2x3v_ser_p1, &vlasov_stream_vol_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_stream_vol_3x3v_ser_p1, NULL               }, // 5
};

__device__ static struct { vlasov_vol_t kernels[3]; } d_vol_kernels[] = {
  // 1x kernels
  { NULL, &vlasov_vol_1x1v_ser_p1, &vlasov_vol_1x1v_ser_p2 }, // 0
  { NULL, &vlasov_vol_1x2v_ser_p1, &vlasov_vol_1x2v_ser_p2 }, // 1
  { NULL, &vlasov_vol_1x3v_ser_p1, &vlasov_vol_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, &vlasov_vol_2x2v_ser_p1, &vlasov_vol_2x2v_ser_p2 }, // 3
  { NULL, &vlasov_vol_2x3v_ser_p1, &vlasov_vol_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_vol_3x3v_ser_p1, NULL               }, // 5
};

// Streaming surface kernel list: x-direction
__device__ static struct { vlasov_stream_surf_t kernels[3]; } d_stream_surf_x_kernels[] = {
  // 1x kernels
  { NULL, &vlasov_surfx_1x1v_ser_p1, &vlasov_surfx_1x1v_ser_p2 }, // 0
  { NULL, &vlasov_surfx_1x2v_ser_p1, &vlasov_surfx_1x2v_ser_p2 }, // 1
  { NULL, &vlasov_surfx_1x3v_ser_p1, &vlasov_surfx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, &vlasov_surfx_2x2v_ser_p1, &vlasov_surfx_2x2v_ser_p2 }, // 3
  { NULL, &vlasov_surfx_2x3v_ser_p1, &vlasov_surfx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_surfx_3x3v_ser_p1, NULL                  }, // 5
};

// Streaming surface kernel list: y-direction
__device__ static struct { vlasov_stream_surf_t kernels[3]; } d_stream_surf_y_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, &vlasov_surfy_2x2v_ser_p1, &vlasov_surfy_2x2v_ser_p2 }, // 3
  { NULL, &vlasov_surfy_2x3v_ser_p1, &vlasov_surfy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_surfy_3x3v_ser_p1, NULL                  }, // 5
};

// Streaming surface kernel list: z-direction
__device__ static struct { vlasov_stream_surf_t kernels[3]; } d_stream_surf_z_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, NULL, NULL }, // 2  
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, &vlasov_surfz_3x3v_ser_p1, NULL }, // 5
};

// Acceleration surface kernel list: vx-direction
__device__ static struct { vlasov_accel_surf_t kernels[3]; } d_accel_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, &vlasov_surfvx_1x1v_ser_p1, &vlasov_surfvx_1x1v_ser_p2 }, // 0
  { NULL, &vlasov_surfvx_1x2v_ser_p1, &vlasov_surfvx_1x2v_ser_p2 }, // 1
  { NULL, &vlasov_surfvx_1x3v_ser_p1, &vlasov_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, &vlasov_surfvx_2x2v_ser_p1, &vlasov_surfvx_2x2v_ser_p2 }, // 3
  { NULL, &vlasov_surfvx_2x3v_ser_p1, &vlasov_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration surface kernel list: vy-direction
__device__ static struct { vlasov_accel_surf_t kernels[3]; } d_accel_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, &vlasov_surfvy_1x2v_ser_p1, &vlasov_surfvy_1x2v_ser_p2 }, // 1
  { NULL, &vlasov_surfvy_1x3v_ser_p1, &vlasov_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, &vlasov_surfvy_2x2v_ser_p1, &vlasov_surfvy_2x2v_ser_p2 }, // 3
  { NULL, &vlasov_surfvy_2x3v_ser_p1, &vlasov_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration surface kernel list: vz-direction
__device__ static struct { vlasov_accel_surf_t kernels[3]; } d_accel_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, &vlasov_surfvz_1x3v_ser_p1, &vlasov_surfvz_1x3v_ser_p2}, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, &vlasov_surfvz_2x3v_ser_p1, &vlasov_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_surfvz_3x3v_ser_p1, NULL }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vx-direction
__device__ static struct { vlasov_accel_boundary_surf_t kernels[3]; } d_accel_boundary_surf_vx_kernels[] = {
  // 1x kernels
  { NULL, &vlasov_boundary_surfvx_1x1v_ser_p1, &vlasov_boundary_surfvx_1x1v_ser_p2 }, // 0
  { NULL, &vlasov_boundary_surfvx_1x2v_ser_p1, &vlasov_boundary_surfvx_1x2v_ser_p2 }, // 1
  { NULL, &vlasov_boundary_surfvx_1x3v_ser_p1, &vlasov_boundary_surfvx_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, &vlasov_boundary_surfvx_2x2v_ser_p1, &vlasov_boundary_surfvx_2x2v_ser_p2 }, // 3
  { NULL, &vlasov_boundary_surfvx_2x3v_ser_p1, &vlasov_boundary_surfvx_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_boundary_surfvx_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vy-direction
__device__ static struct { vlasov_accel_boundary_surf_t kernels[3]; } d_accel_boundary_surf_vy_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, &vlasov_boundary_surfvy_1x2v_ser_p1, &vlasov_boundary_surfvy_1x2v_ser_p2 }, // 1
  { NULL, &vlasov_boundary_surfvy_1x3v_ser_p1, &vlasov_boundary_surfvy_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, &vlasov_boundary_surfvy_2x2v_ser_p1, &vlasov_boundary_surfvy_2x2v_ser_p2 }, // 3
  { NULL, &vlasov_boundary_surfvy_2x3v_ser_p1, &vlasov_boundary_surfvy_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_boundary_surfvy_3x3v_ser_p1, NULL                   }, // 5
};

// Acceleration boundary surface kernel (zero-flux BCs) list: vz-direction
__device__ static struct { vlasov_accel_boundary_surf_t kernels[3]; } d_accel_boundary_surf_vz_kernels[] = {
  // 1x kernels
  { NULL, NULL, NULL }, // 0
  { NULL, NULL, NULL }, // 1
  { NULL, &vlasov_boundary_surfvz_1x3v_ser_p1, &vlasov_boundary_surfvz_1x3v_ser_p2}, // 2
  // 2x kernels
  { NULL, NULL, NULL }, // 3
  { NULL, &vlasov_boundary_surfvz_2x3v_ser_p1, &vlasov_boundary_surfvz_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_boundary_surfvz_3x3v_ser_p1, NULL }, // 5
};

// CUDA kernel to set pointer to qmem = q/m*EM
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ void
gkyl_vlasov_set_qmem_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *qmem)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);
  vlasov->qmem = qmem;
}

// Host-side wrapper for set_qmem_cu_kernel
void
gkyl_vlasov_set_qmem_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *qmem)
{
  gkyl_vlasov_set_qmem_cu_kernel<<<1,1>>>(eqn, qmem);
}

// CUDA kernel to set device pointers to range object and vlasov kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ void
dg_vlasov_set_cu_dev_ptrs(struct dg_vlasov *vlasov, const struct gkyl_range *conf_range, int cv_index, int cdim, int vdim, int polyOrder)
{
  vlasov->qmem = 0; 
  vlasov->conf_range = *conf_range;

  vlasov->eqn.num_equations = 1;
  vlasov->eqn.vol_term = d_vol;
  vlasov->eqn.surf_term = d_surf;
  vlasov->eqn.boundary_surf_term = d_boundary_surf;
 
  vlasov->vol = d_vol_kernels[cv_index].kernels[polyOrder];

  vlasov->stream_surf[0] = d_stream_surf_x_kernels[cv_index].kernels[polyOrder];
  if (cdim>1)
    vlasov->stream_surf[1] = d_stream_surf_y_kernels[cv_index].kernels[polyOrder];
  if (cdim>2)
    vlasov->stream_surf[2] = d_stream_surf_z_kernels[cv_index].kernels[polyOrder];

  vlasov->accel_surf[0] = d_accel_surf_vx_kernels[cv_index].kernels[polyOrder];
  if (vdim>1)
    vlasov->accel_surf[1] = d_accel_surf_vy_kernels[cv_index].kernels[polyOrder];
  if (vdim>2)
    vlasov->accel_surf[2] = d_accel_surf_vz_kernels[cv_index].kernels[polyOrder];

  vlasov->accel_boundary_surf[0] = d_accel_boundary_surf_vx_kernels[cv_index].kernels[polyOrder];
  if (vdim>1)
    vlasov->accel_boundary_surf[1] = d_accel_boundary_surf_vy_kernels[cv_index].kernels[polyOrder];
  if (vdim>2)
    vlasov->accel_boundary_surf[2] = d_accel_boundary_surf_vz_kernels[cv_index].kernels[polyOrder];
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range_cu)
{
  struct dg_vlasov *vlasov = (struct dg_vlasov*) gkyl_malloc(sizeof(struct dg_vlasov));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int polyOrder = cbasis->polyOrder;

  vlasov->cdim = cdim;
  vlasov->pdim = pdim;

  vlasov->eqn.num_equations = 1;

  // copy the host struct to device struct
  struct dg_vlasov *vlasov_cu = (struct dg_vlasov*) gkyl_cu_malloc(sizeof(struct dg_vlasov));
  gkyl_cu_memcpy(vlasov_cu, vlasov, sizeof(struct dg_vlasov), GKYL_CU_MEMCPY_H2D);

  dg_vlasov_set_cu_dev_ptrs<<<1,1>>>(vlasov_cu, conf_range_cu, cv_index[cdim].vdim[vdim], cdim, vdim, polyOrder);

  gkyl_free(vlasov);  
  
  return &vlasov_cu->eqn;
}
