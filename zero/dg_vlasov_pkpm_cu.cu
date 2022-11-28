/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_vlasov_pkpm.h>    
#include <gkyl_dg_vlasov_pkpm_priv.h>
}

#include <cassert>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_vlasov_pkpm_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, 
  const struct gkyl_array *bvar, const struct gkyl_array *u_i, 
  const struct gkyl_array *bb_grad_u, const struct gkyl_array *p_force, const struct gkyl_array *vth_sq)
{
  struct dg_vlasov_pkpm *vlasov_pkpm = container_of(eqn, struct dg_vlasov_pkpm, eqn);
  vlasov_pkpm->auxfields.bvar = bvar;
  vlasov_pkpm->auxfields.u_i = u_i;
  vlasov_pkpm->auxfields.bb_grad_u = bb_grad_u;
  vlasov_pkpm->auxfields.p_force = p_force;
  vlasov_pkpm->auxfields.vth_sq = vth_sq;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_vlasov_pkpm_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_vlasov_pkpm_auxfields auxin)
{
  gkyl_vlasov_pkpm_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.bvar->on_dev, auxin.u_i->on_dev, 
    auxin.bb_grad_u->on_dev, auxin.p_force->on_dev, auxin.vth_sq->on_dev);
}

// CUDA kernel to set device pointers to range object and vlasov_pkpm kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_vlasov_pkpm_set_cu_dev_ptrs(struct dg_vlasov_pkpm *vlasov_pkpm, enum gkyl_basis_type b_type,
  int cdim, int poly_order)
{
  vlasov_pkpm->auxfields.bvar = 0;  
  vlasov_pkpm->auxfields.u_i = 0;
  vlasov_pkpm->auxfields.bb_grad_u = 0;
  vlasov_pkpm->auxfields.p_force = 0;  
  vlasov_pkpm->auxfields.vth_sq = 0;  
  
  vlasov_pkpm->eqn.surf_term = surf;
  vlasov_pkpm->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_vlasov_pkpm_vol_kern_list *vol_kernels;
  const gkyl_dg_vlasov_pkpm_stream_surf_kern_list *stream_surf_x_kernels, *stream_surf_y_kernels, *stream_surf_z_kernels;
  const gkyl_dg_vlasov_pkpm_accel_surf_kern_list *accel_surf_vpar_kernels;
  const gkyl_dg_vlasov_pkpm_accel_boundary_surf_kern_list *accel_boundary_surf_vpar_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      stream_surf_x_kernels = ser_stream_surf_x_kernels;
      stream_surf_y_kernels = ser_stream_surf_y_kernels;
      stream_surf_z_kernels = ser_stream_surf_z_kernels;
      accel_surf_vpar_kernels = ser_accel_surf_vpar_kernels;
      accel_boundary_surf_vpar_kernels = ser_accel_boundary_surf_vpar_kernels;
      
      break;
      
    default:
      assert(false);
      break;    
  }  

  vlasov_pkpm->eqn.vol_term = CK(vol_kernels,cdim,poly_order);

  vlasov_pkpm->stream_surf[0] = CK(stream_surf_x_kernels,cdim,poly_order);
  if (cdim>1)
    vlasov_pkpm->stream_surf[1] = CK(stream_surf_y_kernels,cdim,poly_order);
  if (cdim>2)
    vlasov_pkpm->stream_surf[2] = CK(stream_surf_z_kernels,cdim,poly_order);

  vlasov_pkpm->accel_surf = CK(accel_surf_vpar_kernels,cdim,poly_order);

  vlasov_pkpm->accel_boundary_surf = CK(accel_boundary_surf_vpar_kernels,cdim,poly_order);
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_pkpm_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range)
{
  struct dg_vlasov_pkpm *vlasov_pkpm = (struct dg_vlasov_pkpm*) gkyl_malloc(sizeof(struct dg_vlasov_pkpm));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  vlasov_pkpm->cdim = cdim;
  vlasov_pkpm->pdim = pdim;

  vlasov_pkpm->eqn.num_equations = 1;
  vlasov_pkpm->conf_range = *conf_range;

  vlasov_pkpm->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(vlasov_pkpm->eqn.flags);
  vlasov_pkpm->eqn.ref_count = gkyl_ref_count_init(gkyl_vlasov_pkpm_free);

  // copy the host struct to device struct
  struct dg_vlasov_pkpm *vlasov_pkpm_cu = (struct dg_vlasov_pkpm*) gkyl_cu_malloc(sizeof(struct dg_vlasov_pkpm));
  gkyl_cu_memcpy(vlasov_pkpm_cu, vlasov_pkpm, sizeof(struct dg_vlasov_pkpm), GKYL_CU_MEMCPY_H2D);

  dg_vlasov_pkpm_set_cu_dev_ptrs<<<1,1>>>(vlasov_pkpm_cu, cbasis->b_type, cdim, poly_order);

  // set parent on_dev pointer
  vlasov_pkpm->eqn.on_dev = &vlasov_pkpm_cu->eqn;
  
  return &vlasov_pkpm->eqn;
}
