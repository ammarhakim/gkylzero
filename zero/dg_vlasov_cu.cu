/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_vlasov.h>    
#include <gkyl_dg_vlasov_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to qmem = q/m*EM
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_vlasov_set_qmem_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *qmem)
{
  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);
  vlasov->qmem = qmem;
}

// Host-side wrapper for set_qmem_cu_kernel
void
gkyl_vlasov_set_qmem_cu(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *qmem)
{
  gkyl_vlasov_set_qmem_cu_kernel<<<1,1>>>(eqn, qmem->on_dev);
}

// CUDA kernel to set device pointers to range object and vlasov kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_vlasov_set_cu_dev_ptrs(struct dg_vlasov *vlasov, enum gkyl_basis_type b_type,
  int cv_index, int cdim, int vdim, int poly_order, enum gkyl_field_id field_id)
{
  vlasov->qmem = 0; 

  vlasov->eqn.vol_term = vol;
  vlasov->eqn.surf_term = surf;
  vlasov->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_vlasov_stream_vol_kern_list *stream_vol_kernels;
  const gkyl_dg_vlasov_vol_kern_list *vol_kernels;
  const gkyl_dg_vlasov_stream_surf_kern_list *stream_surf_x_kernels, *stream_surf_y_kernels, *stream_surf_z_kernels;
  const gkyl_dg_vlasov_accel_surf_kern_list *accel_surf_vx_kernels, *accel_surf_vy_kernels, *accel_surf_vz_kernels;
  const gkyl_dg_vlasov_accel_boundary_surf_kern_list *accel_boundary_surf_vx_kernels, *accel_boundary_surf_vy_kernels,
    *accel_boundary_surf_vz_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      stream_vol_kernels = ser_stream_vol_kernels;
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
      stream_vol_kernels = ten_stream_vol_kernels;
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
  if (field_id == GKYL_FIELD_NULL)
    vlasov->vol = stream_vol_kernels[cv_index].kernels[poly_order];
  else
    vlasov->vol = vol_kernels[cv_index].kernels[poly_order];

  vlasov->stream_surf[0] = stream_surf_x_kernels[cv_index].kernels[poly_order];
  if (cdim>1)
    vlasov->stream_surf[1] = stream_surf_y_kernels[cv_index].kernels[poly_order];
  if (cdim>2)
    vlasov->stream_surf[2] = stream_surf_z_kernels[cv_index].kernels[poly_order];

  vlasov->accel_surf[0] = accel_surf_vx_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    vlasov->accel_surf[1] = accel_surf_vy_kernels[cv_index].kernels[poly_order];
  if (vdim>2)
    vlasov->accel_surf[2] = accel_surf_vz_kernels[cv_index].kernels[poly_order];

  vlasov->accel_boundary_surf[0] = accel_boundary_surf_vx_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    vlasov->accel_boundary_surf[1] = accel_boundary_surf_vy_kernels[cv_index].kernels[poly_order];
  if (vdim>2)
    vlasov->accel_boundary_surf[2] = accel_boundary_surf_vz_kernels[cv_index].kernels[poly_order];
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, enum gkyl_field_id field_id)
{
  struct dg_vlasov *vlasov = (struct dg_vlasov*) gkyl_malloc(sizeof(struct dg_vlasov));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  vlasov->cdim = cdim;
  vlasov->pdim = pdim;

  vlasov->eqn.num_equations = 1;
  vlasov->conf_range = *conf_range;

  vlasov->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(vlasov->eqn.flags);
  vlasov->eqn.ref_count = gkyl_ref_count_init(gkyl_vlasov_free);

  // copy the host struct to device struct
  struct dg_vlasov *vlasov_cu = (struct dg_vlasov*) gkyl_cu_malloc(sizeof(struct dg_vlasov));
  gkyl_cu_memcpy(vlasov_cu, vlasov, sizeof(struct dg_vlasov), GKYL_CU_MEMCPY_H2D);

  dg_vlasov_set_cu_dev_ptrs<<<1,1>>>(vlasov_cu, cbasis->b_type, cv_index[cdim].vdim[vdim],
    cdim, vdim, poly_order, field_id);

  // set parent on_dev pointer
  vlasov->eqn.on_dev = &vlasov_cu->eqn;
  
  return &vlasov->eqn;
}
