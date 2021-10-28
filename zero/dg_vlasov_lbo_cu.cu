/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_dg_vlasov_lbo.h>    
#include <gkyl_dg_vlasov_lbo_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to nuSum, sum of collisionalities
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ void
gkyl_vlasov_lbo_set_nuSum_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuSum)
{
  struct dg_vlasov_lbo *vlasov_lbo = container_of(eqn, struct dg_vlasov_lbo, eqn);
  vlasov_lbo->nuSum = nuSum;
}

// Host-side wrapper for set_nuSum_cu_kernel
void
gkyl_vlasov_lbo_set_nuSum_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuSum)
{
  gkyl_vlasov_lbo_set_nuSum_cu_kernel<<<1,1>>>(eqn, nuSum->on_dev);
}

// CUDA kernel to set pointer to nuUSum, sum of nu*u for updating the drag flux term
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ void
gkyl_vlasov_lbo_set_nuUSum_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuUSum)
{
  struct dg_vlasov_lbo *vlasov_lbo = container_of(eqn, struct dg_vlasov_lbo, eqn);
  vlasov_lbo->nuUSum = nuUSum;
}

// Host-side wrapper for set_nuUSum_cu_kernel
void
gkyl_vlasov_lbo_set_nuUSum_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuUSum)
{
  gkyl_vlasov_lbo_set_nuUSum_cu_kernel<<<1,1>>>(eqn, nuUSum->on_dev);
}

// CUDA kernel to set pointer to nuVtSqSum, sum of nu*vth^2 for updating the diffusion flux term.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ void
gkyl_vlasov_lbo_set_nuVtSqSum_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuVtSqSum)
{
  struct dg_vlasov_lbo *vlasov_lbo = container_of(eqn, struct dg_vlasov_lbo, eqn);
  vlasov_lbo->nuVtSqSum = nuVtSqSum;
}

// Host-side wrapper for set_nuVtSqSum_cu_kernel
void
gkyl_vlasov_lbo_set_nuVtSqSum_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuVtSqSum)
{
  gkyl_vlasov_lbo_set_nuVtSqSum_cu_kernel<<<1,1>>>(eqn, nuVtSqSum->on_dev);
}

// CUDA kernel to set device pointers to range object and vlasov LBO kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ void static
dg_vlasov_lbo_set_cu_dev_ptrs(struct dg_vlasov_lbo *vlasov_lbo, enum gkyl_basis_type b_type,
  int cv_index, int cdim, int vdim, int poly_order)
{
  vlasov_lbo->nuSum = 0; 
  vlasov_lbo->nuUSum = 0; 
  vlasov_lbo->nuVtSqSum = 0; 

  vlasov_lbo->eqn.vol_term = vol;
  vlasov_lbo->eqn.surf_term = surf;
  vlasov_lbo->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_vlasov_lbo_vol_kern_list *vol_kernels;
  const gkyl_dg_vlasov_lbo_surf_kern_list *surf_vx_kernels, *surf_vy_kernels, *surf_vz_kernels;
  const gkyl_dg_vlasov_lbo_boundary_surf_kern_list *boundary_surf_vx_kernels, *boundary_surf_vy_kernels,
    *boundary_surf_vz_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_vx_kernels = ser_surf_vx_kernels;
      surf_vy_kernels = ser_surf_vy_kernels;
      surf_vz_kernels = ser_surf_vz_kernels;
      boundary_surf_vx_kernels = ser_boundary_surf_vx_kernels;
      boundary_surf_vy_kernels = ser_boundary_surf_vy_kernels;
      boundary_surf_vz_kernels = ser_boundary_surf_vz_kernels;
      
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      vol_kernels = ten_vol_kernels;
      surf_vx_kernels = ten_surf_vx_kernels;
      surf_vy_kernels = ten_surf_vy_kernels;
      surf_vz_kernels = ten_surf_vz_kernels;
      boundary_surf_vx_kernels = ten_boundary_surf_vx_kernels;
      boundary_surf_vy_kernels = ten_boundary_surf_vy_kernels;
      boundary_surf_vz_kernels = ten_boundary_surf_vz_kernels;
      break;

    default:
      assert(false);
      break;    
  }  
 
  vlasov_lbo->vol = vol_kernels[cv_index].kernels[poly_order];

  vlasov_lbo->surf[0] = surf_vx_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    vlasov_lbo->surf[1] = surf_vy_kernels[cv_index].kernels[poly_order];
  if (vdim>2)
    vlasov_lbo->surf[2] = surf_vz_kernels[cv_index].kernels[poly_order];

  vlasov_lbo->boundary_surf[0] = boundary_surf_vx_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    vlasov_lbo->boundary_surf[1] = boundary_surf_vy_kernels[cv_index].kernels[poly_order];
  if (vdim>2)
    vlasov_lbo->boundary_surf[2] = boundary_surf_vz_kernels[cv_index].kernels[poly_order];
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_lbo_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range)
{
  struct dg_vlasov_lbo *vlasov_lbo = (struct dg_vlasov_lbo*) gkyl_malloc(sizeof(struct dg_vlasov_lbo));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  vlasov_lbo->cdim = cdim;
  vlasov_lbo->pdim = pdim;

  vlasov_lbo->eqn.num_equations = 1;
  vlasov_lbo->conf_range = *conf_range;

  // copy the host struct to device struct
  struct dg_vlasov_lbo *vlasov_lbo_cu = (struct dg_vlasov_lbo*) gkyl_cu_malloc(sizeof(struct dg_vlasov_lbo));
  gkyl_cu_memcpy(vlasov_lbo_cu, vlasov_lbo, sizeof(struct dg_vlasov_lbo), GKYL_CU_MEMCPY_H2D);

  dg_vlasov_lbo_set_cu_dev_ptrs<<<1,1>>>(vlasov_lbo_cu, cbasis->b_type, cv_index[cdim].vdim[vdim], cdim, vdim, poly_order);

  gkyl_free(vlasov_lbo);  
  
  return &vlasov_lbo_cu->eqn;
}
