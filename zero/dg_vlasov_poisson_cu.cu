/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_dg_vlasov_poisson_priv.h>    
#include <gkyl_dg_vlasov_poisson.h>
}

#include <cassert>

// CUDA kernel to set pointer to fac_phi = factor*phi
// This factor is q/m for plasmas and G*m for self-gravitating systems
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_vlasov_poisson_set_fac_phi_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *fac_phi)
{
  struct dg_vlasov_poisson *vlasov_poisson = container_of(eqn, struct dg_vlasov_poisson, eqn);
  vlasov_poisson->fac_phi = fac_phi;
}

// Host-side wrapper for set_fac_phi_cu_kernel
void
gkyl_vlasov_poisson_set_fac_phi_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *fac_phi)
{
  gkyl_vlasov_poisson_set_fac_phi_cu_kernel<<<1,1>>>(eqn, fac_phi->on_dev);
}

// CUDA kernel to set pointer to vecA = q/m*A, where A is the vector potential
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_vlasov_poisson_set_vecA_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *vecA)
{
  struct dg_vlasov_poisson *vlasov_poisson = container_of(eqn, struct dg_vlasov_poisson, eqn);
  vlasov_poisson->vecA = vecA;
}

// Host-side wrapper for set_vecA_cu_kernel
void
gkyl_vlasov_poisson_set_vecA_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *vecA)
{
  gkyl_vlasov_poisson_set_vecA_cu_kernel<<<1,1>>>(eqn, vecA->on_dev);
}

// CUDA kernel to set device pointers to range object and vlasov kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_vlasov_poisson_set_cu_dev_ptrs(struct dg_vlasov_poisson *vlasov_poisson, enum gkyl_basis_type b_type,
  int cv_index, int cdim, int vdim, int poly_order)
{
  vlasov_poisson->fac_phi = 0; 
  vlasov_poisson->vecA = 0; 

  vlasov_poisson->eqn.vol_term = vol;
  vlasov_poisson->eqn.surf_term = surf;
  vlasov_poisson->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_vlasov_poisson_vol_kern_list *vol_kernels;
  const gkyl_dg_vlasov_poisson_stream_surf_kern_list *stream_surf_x_kernels, *stream_surf_y_kernels, *stream_surf_z_kernels;
  const gkyl_dg_vlasov_poisson_accel_surf_kern_list *accel_surf_vx_kernels, *accel_surf_vy_kernels, *accel_surf_vz_kernels;
  const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list *accel_boundary_surf_vx_kernels, *accel_boundary_surf_vy_kernels,
    *accel_boundary_surf_vz_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      stream_surf_x_kernels = ser_stream_surf_x_kernels;
      stream_surf_y_kernels = ser_stream_surf_y_kernels;
      stream_surf_z_kernels = ser_stream_surf_z_kernels;
      accel_surf_vx_kernels = ser_accel_surf_vx_kernels;
      accel_surf_vy_kernels = ser_accel_surf_vy_kernels;
      accel_surf_vz_kernels = ser_accel_surf_vz_kernels;
      accel_boundary_surf_vx_kernels = ser_accel_boundary_surf_vx_kernels;
      accel_boundary_surf_vy_kernels = ser_accel_boundary_surf_vy_kernels;
      accel_boundary_surf_vz_kernels = ser_accel_boundary_surf_vz_kernels;
      
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      vol_kernels = ten_vol_kernels;
      stream_surf_x_kernels = ten_stream_surf_x_kernels;
      stream_surf_y_kernels = ten_stream_surf_y_kernels;
      stream_surf_z_kernels = ten_stream_surf_z_kernels;
      accel_surf_vx_kernels = ten_accel_surf_vx_kernels;
      accel_surf_vy_kernels = ten_accel_surf_vy_kernels;
      accel_surf_vz_kernels = ten_accel_surf_vz_kernels;
      accel_boundary_surf_vx_kernels = ten_accel_boundary_surf_vx_kernels;
      accel_boundary_surf_vy_kernels = ten_accel_boundary_surf_vy_kernels;
      accel_boundary_surf_vz_kernels = ten_accel_boundary_surf_vz_kernels;
      break;

    default:
      assert(false);
      break;    
  }  
 
  vlasov_poisson->vol = vol_kernels[cv_index].kernels[poly_order];

  vlasov_poisson->stream_surf[0] = stream_surf_x_kernels[cv_index].kernels[poly_order];
  if (cdim>1)
    vlasov_poisson->stream_surf[1] = stream_surf_y_kernels[cv_index].kernels[poly_order];
  if (cdim>2)
    vlasov_poisson->stream_surf[2] = stream_surf_z_kernels[cv_index].kernels[poly_order];

  vlasov_poisson->accel_surf[0] = accel_surf_vx_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    vlasov_poisson->accel_surf[1] = accel_surf_vy_kernels[cv_index].kernels[poly_order];
  if (vdim>2)
    vlasov_poisson->accel_surf[2] = accel_surf_vz_kernels[cv_index].kernels[poly_order];

  vlasov_poisson->accel_boundary_surf[0] = accel_boundary_surf_vx_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    vlasov_poisson->accel_boundary_surf[1] = accel_boundary_surf_vy_kernels[cv_index].kernels[poly_order];
  if (vdim>2)
    vlasov_poisson->accel_boundary_surf[2] = accel_boundary_surf_vz_kernels[cv_index].kernels[poly_order];
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_poisson_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range)
{
  struct dg_vlasov_poisson *vlasov_poisson = (struct dg_vlasov_poisson*) gkyl_malloc(sizeof(struct dg_vlasov_poisson));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  vlasov_poisson->cdim = cdim;
  vlasov_poisson->pdim = pdim;

  vlasov_poisson->eqn.num_equations = 1;
  vlasov_poisson->conf_range = *conf_range;

  // copy the host struct to device struct
  struct dg_vlasov_poisson *vlasov_poisson_cu = (struct dg_vlasov_poisson*) gkyl_cu_malloc(sizeof(struct dg_vlasov_poisson));
  gkyl_cu_memcpy(vlasov_poisson_cu, vlasov_poisson, sizeof(struct dg_vlasov_poisson), GKYL_CU_MEMCPY_H2D);

  dg_vlasov_poisson_set_cu_dev_ptrs<<<1,1>>>(vlasov_poisson_cu, cbasis->b_type, cv_index[cdim].vdim[vdim], cdim, vdim, poly_order);

  gkyl_free(vlasov_poisson);  
  
  return &vlasov_poisson_cu->eqn;
}
