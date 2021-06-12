/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_dg_vlasov.h>    
#include <gkyl_dg_vlasov_priv.h>
}

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
  gkyl_vlasov_set_qmem_cu_kernel<<<1,1>>>(eqn, qmem->on_device);
}

// CUDA kernel to set device pointers to range object and vlasov kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ void
dg_vlasov_set_cu_dev_ptrs(struct dg_vlasov *vlasov, int cv_index, int cdim, int vdim, int polyOrder)
{
  vlasov->qmem = 0; 

  vlasov->eqn.vol_term = vol;
  vlasov->eqn.surf_term = surf;
  vlasov->eqn.boundary_surf_term = boundary_surf;
 
  vlasov->vol = vol_kernels[cv_index].kernels[polyOrder];

  vlasov->stream_surf[0] = stream_surf_x_kernels[cv_index].kernels[polyOrder];
  if (cdim>1)
    vlasov->stream_surf[1] = stream_surf_y_kernels[cv_index].kernels[polyOrder];
  if (cdim>2)
    vlasov->stream_surf[2] = stream_surf_z_kernels[cv_index].kernels[polyOrder];

  vlasov->accel_surf[0] = accel_surf_vx_kernels[cv_index].kernels[polyOrder];
  if (vdim>1)
    vlasov->accel_surf[1] = accel_surf_vy_kernels[cv_index].kernels[polyOrder];
  if (vdim>2)
    vlasov->accel_surf[2] = accel_surf_vz_kernels[cv_index].kernels[polyOrder];

  vlasov->accel_boundary_surf[0] = accel_boundary_surf_vx_kernels[cv_index].kernels[polyOrder];
  if (vdim>1)
    vlasov->accel_boundary_surf[1] = accel_boundary_surf_vy_kernels[cv_index].kernels[polyOrder];
  if (vdim>2)
    vlasov->accel_boundary_surf[2] = accel_boundary_surf_vz_kernels[cv_index].kernels[polyOrder];
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
  vlasov->conf_range = *conf_range;

  // copy the host struct to device struct
  struct dg_vlasov *vlasov_cu = (struct dg_vlasov*) gkyl_cu_malloc(sizeof(struct dg_vlasov));
  gkyl_cu_memcpy(vlasov_cu, vlasov, sizeof(struct dg_vlasov), GKYL_CU_MEMCPY_H2D);

  dg_vlasov_set_cu_dev_ptrs<<<1,1>>>(vlasov_cu, cv_index[cdim].vdim[vdim], cdim, vdim, polyOrder);

  gkyl_free(vlasov);  
  
  return &vlasov_cu->eqn;
}
