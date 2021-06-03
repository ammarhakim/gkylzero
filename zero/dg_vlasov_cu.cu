/* -*- c -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_dg_vlasov.h>    
#include <gkyl_dg_vlasov_priv.h>
}

// various pointers to functions on device
__device__ static vol_termf_t p_vol = &vol;
__device__ static surf_termf_t p_surf = &surf;
__device__ static boundary_surf_termf_t p_boundary_surf = &boundary_surf;

// Volume kernel list
__device__ static struct { vlasov_stream_vol_t kernels[3]; } p_stream_vol_kernels[] = {
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

__device__ static struct { vlasov_vol_t kernels[3]; } p_vol_kernels[] = {
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
__device__ static struct { vlasov_stream_surf_t kernels[3]; } p_stream_surf_x_kernels[] = {
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
__device__ static struct { vlasov_stream_surf_t kernels[3]; } p_stream_surf_y_kernels[] = {
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
__device__ static struct { vlasov_stream_surf_t kernels[3]; } p_stream_surf_z_kernels[] = {
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
__device__ static struct { vlasov_accel_surf_t kernels[3]; } p_accel_surf_vx_kernels[] = {
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
__device__ static struct { vlasov_accel_surf_t kernels[3]; } p_accel_surf_vy_kernels[] = {
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
__device__ static struct { vlasov_accel_surf_t kernels[3]; } p_accel_surf_vz_kernels[] = {
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
__device__ static struct { vlasov_accel_boundary_surf_t kernels[3]; } p_accel_boundary_surf_vx_kernels[] = {
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
__device__ static struct { vlasov_accel_boundary_surf_t kernels[3]; } p_accel_boundary_surf_vy_kernels[] = {
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
__device__ static struct { vlasov_accel_boundary_surf_t kernels[3]; } p_accel_boundary_surf_vz_kernels[] = {
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

__device__
void
gkyl_vlasov_set_qmem_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *qmem)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);
  vlasov->qmem = qmem;
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range)
{
  struct dg_vlasov *vlasov = (struct dg_vlasov*) gkyl_malloc(sizeof(struct dg_vlasov));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int polyOrder = cbasis->polyOrder;

  vlasov->cdim = cdim;
  vlasov->pdim = pdim;

  vlasov->eqn.num_equations = 1;
  
  cudaMemcpyFromSymbol(&vlasov->eqn.vol_term, p_vol, sizeof(vol_termf_t));
  cudaMemcpyFromSymbol(&vlasov->eqn.surf_term, p_surf, sizeof(surf_termf_t));
  cudaMemcpyFromSymbol(&vlasov->eqn.boundary_surf_term, p_boundary_surf, sizeof(boundary_surf_termf_t));

  cudaMemcpyFromSymbol(&vlasov->vol, CK(p_vol_kernels,cdim,vdim,polyOrder), sizeof(vlasov_vol_t));

  cudaMemcpyFromSymbol(&vlasov->stream_surf[0],
    CK(p_stream_surf_x_kernels,cdim,vdim,polyOrder), sizeof(vlasov_stream_surf_t));
  if (cdim>1)
    cudaMemcpyFromSymbol(&vlasov->stream_surf[1],
      CK(p_stream_surf_y_kernels,cdim,vdim,polyOrder), sizeof(vlasov_stream_surf_t));
  if (cdim>2)
    cudaMemcpyFromSymbol(&vlasov->stream_surf[2],
      CK(p_stream_surf_z_kernels,cdim,vdim,polyOrder), sizeof(vlasov_stream_surf_t));

  cudaMemcpyFromSymbol(&vlasov->accel_surf[0],
    CK(p_accel_surf_vx_kernels,cdim,vdim,polyOrder), sizeof(vlasov_accel_surf_t));
  if (vdim>1)
    cudaMemcpyFromSymbol(&vlasov->accel_surf[1],
      CK(p_accel_surf_vy_kernels,cdim,vdim,polyOrder), sizeof(vlasov_accel_surf_t));
  if (vdim>2)
    cudaMemcpyFromSymbol(&vlasov->accel_surf[2],
      CK(p_accel_surf_vz_kernels,cdim,vdim,polyOrder), sizeof(vlasov_accel_surf_t));

  cudaMemcpyFromSymbol(&vlasov->accel_boundary_surf[0],
    CK(p_accel_boundary_surf_vx_kernels,cdim,vdim,polyOrder), sizeof(vlasov_accel_boundary_surf_t));
  if (vdim>1)
    cudaMemcpyFromSymbol(&vlasov->accel_boundary_surf[1],
      CK(p_accel_boundary_surf_vy_kernels,cdim,vdim,polyOrder), sizeof(vlasov_accel_boundary_surf_t));
  if (vdim>2)
    cudaMemcpyFromSymbol(&vlasov->accel_boundary_surf[2],
      CK(p_accel_boundary_surf_vz_kernels,cdim,vdim,polyOrder), sizeof(vlasov_accel_boundary_surf_t));

  vlasov->qmem = 0;

  // copy the host struct to device struct
  struct dg_vlasov *vlasov_cu = (struct dg_vlasov*) gkyl_cu_malloc(sizeof(struct dg_vlasov));
  gkyl_cu_memcpy(vlasov_cu, vlasov, sizeof(struct dg_vlasov), GKYL_CU_MEMCPY_H2D);

  // we need to copy conf_range seperately as it is already on device
  gkyl_cu_memcpy(&vlasov_cu->conf_range, (void*) conf_range, sizeof(struct gkyl_range), GKYL_CU_MEMCPY_D2D);

  gkyl_free(vlasov);  
  
  return &vlasov_cu->eqn;
}
