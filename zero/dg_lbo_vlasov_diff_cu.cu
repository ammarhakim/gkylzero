/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_lbo_vlasov_diff.h>    
#include <gkyl_dg_lbo_vlasov_diff_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to nuSum, sum of collisionalities
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_lbo_vlasov_diff_set_nuSum_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuSum)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  lbo_vlasov_diff->nuSum = nuSum;
}

// Host-side wrapper for set_nuSum_cu_kernel
void
gkyl_lbo_vlasov_diff_set_nuSum_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuSum)
{
  gkyl_lbo_vlasov_diff_set_nuSum_cu_kernel<<<1,1>>>(eqn, nuSum->on_dev);
}

// CUDA kernel to set pointer to nuUSum, sum of nu*u for updating the drag flux term
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_lbo_vlasov_diff_set_nuUSum_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuUSum)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  lbo_vlasov_diff->nuUSum = nuUSum;
}

// Host-side wrapper for set_nuUSum_cu_kernel
void
gkyl_lbo_vlasov_diff_set_nuUSum_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuUSum)
{
  gkyl_lbo_vlasov_diff_set_nuUSum_cu_kernel<<<1,1>>>(eqn, nuUSum->on_dev);
}

// CUDA kernel to set pointer to nuVtSqSum, sum of nu*vth^2 for updating the diffusion flux term.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_lbo_vlasov_diff_set_nuVtSqSum_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuVtSqSum)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff = container_of(eqn, struct dg_lbo_vlasov_diff, eqn);
  lbo_vlasov_diff->nuVtSqSum = nuVtSqSum;
}

// Host-side wrapper for set_nuVtSqSum_cu_kernel
void
gkyl_lbo_vlasov_diff_set_nuVtSqSum_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuVtSqSum)
{
  gkyl_lbo_vlasov_diff_set_nuVtSqSum_cu_kernel<<<1,1>>>(eqn, nuVtSqSum->on_dev);
}

// CUDA kernel to set device pointers to range object and vlasov LBO kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_lbo_vlasov_diff_set_cu_dev_ptrs(struct dg_lbo_vlasov_diff *lbo_vlasov_diff, enum gkyl_basis_type b_type,
  int cv_index, int cdim, int vdim, int poly_order)
{
  lbo_vlasov_diff->nuSum = 0; 
  lbo_vlasov_diff->nuUSum = 0; 
  lbo_vlasov_diff->nuVtSqSum = 0; 

  lbo_vlasov_diff->eqn.vol_term = vol;
  lbo_vlasov_diff->eqn.surf_term = surf;
  lbo_vlasov_diff->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_lbo_vlasov_diff_vol_kern_list *vol_kernels;
  const gkyl_dg_lbo_vlasov_diff_surf_kern_list *surf_vx_kernels, *surf_vy_kernels, *surf_vz_kernels;
  const gkyl_dg_lbo_vlasov_diff_boundary_surf_kern_list *boundary_surf_vx_kernels, *boundary_surf_vy_kernels,
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

    // case GKYL_BASIS_MODAL_TENSOR:
    //   vol_kernels = ten_vol_kernels;
    //   surf_vx_kernels = ten_surf_vx_kernels;
    //   surf_vy_kernels = ten_surf_vy_kernels;
    //   surf_vz_kernels = ten_surf_vz_kernels;
    //   boundary_surf_vx_kernels = ten_boundary_surf_vx_kernels;
    //   boundary_surf_vy_kernels = ten_boundary_surf_vy_kernels;
    //   boundary_surf_vz_kernels = ten_boundary_surf_vz_kernels;
    //   break;

    default:
      assert(false);
      break;    
  }  
 
  lbo_vlasov_diff->vol = vol_kernels[cv_index].kernels[poly_order];

  lbo_vlasov_diff->surf[0] = surf_vx_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    lbo_vlasov_diff->surf[1] = surf_vy_kernels[cv_index].kernels[poly_order];
  if (vdim>2)
    lbo_vlasov_diff->surf[2] = surf_vz_kernels[cv_index].kernels[poly_order];

  lbo_vlasov_diff->boundary_surf[0] = boundary_surf_vx_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    lbo_vlasov_diff->boundary_surf[1] = boundary_surf_vy_kernels[cv_index].kernels[poly_order];
  if (vdim>2)
    lbo_vlasov_diff->boundary_surf[2] = boundary_surf_vz_kernels[cv_index].kernels[poly_order];
}

struct gkyl_dg_eqn*
gkyl_dg_lbo_vlasov_diff_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range)
{
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff =
    (struct dg_lbo_vlasov_diff*) gkyl_malloc(sizeof(struct dg_lbo_vlasov_diff));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  lbo_vlasov_diff->cdim = cdim;
  lbo_vlasov_diff->pdim = pdim;

  lbo_vlasov_diff->eqn.num_equations = 1;
  lbo_vlasov_diff->conf_range = *conf_range;

  lbo_vlasov_diff->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(lbo_vlasov_diff->eqn.flags);
  lbo_vlasov_diff->eqn.ref_count = gkyl_ref_count_init(gkyl_lbo_vlasov_diff_free);

  // copy the host struct to device struct
  struct dg_lbo_vlasov_diff *lbo_vlasov_diff_cu =
    (struct dg_lbo_vlasov_diff*) gkyl_cu_malloc(sizeof(struct dg_lbo_vlasov_diff));

  gkyl_cu_memcpy(lbo_vlasov_diff_cu, lbo_vlasov_diff,
    sizeof(struct dg_lbo_vlasov_diff), GKYL_CU_MEMCPY_H2D);

  dg_lbo_vlasov_diff_set_cu_dev_ptrs<<<1,1>>>(lbo_vlasov_diff_cu,
    cbasis->b_type, cv_index[cdim].vdim[vdim], cdim, vdim, poly_order);

  lbo_vlasov_diff->eqn.on_dev = &lbo_vlasov_diff_cu->eqn;
  
  return &lbo_vlasov_diff->eqn;
}
