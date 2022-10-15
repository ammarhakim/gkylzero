/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_diffusion_euler_iso.h>    
#include <gkyl_dg_diffusion_euler_iso_priv.h>
}

#include <cassert>

#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_diffusion_euler_iso_set_auxfields_cu_kernel(const struct gkyl_dg_eqn* eqn, const struct gkyl_array* D, const struct gkyl_array* u_i)
{
  struct dg_diffusion_euler_iso* diffusion_euler_iso = container_of(eqn, struct dg_diffusion_euler_iso, eqn);
  diffusion_euler_iso->auxfields.D = D;
  diffusion_euler_iso->auxfields.u_i = u_i;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_diffusion_euler_iso_set_auxfields_cu(const struct gkyl_dg_eqn* eqn, struct gkyl_dg_diffusion_euler_iso_auxfields auxin)
{
  gkyl_diffusion_euler_iso_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.D->on_dev, auxin.u_i->on_dev);
}

__global__ void static
dg_diffusion_euler_iso_set_cu_dev_ptrs(struct dg_diffusion_euler_iso* diffusion_euler_iso, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  diffusion_euler_iso->auxfields.D = 0; 
  diffusion_euler_iso->auxfields.u_i = 0; 

  const gkyl_dg_diffusion_euler_iso_vol_kern_list* vol_kernels;
  const gkyl_dg_diffusion_euler_iso_surf_kern_list* surf_x_kernels;
  const gkyl_dg_diffusion_euler_iso_surf_kern_list* surf_y_kernels;
  const gkyl_dg_diffusion_euler_iso_surf_kern_list* surf_z_kernels; 

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_x_kernels = ser_surf_x_kernels;
      surf_y_kernels = ser_surf_y_kernels;
      surf_z_kernels = ser_surf_z_kernels;
      break;

    default:
      assert(false);
      break;    
  } 
  
  diffusion_euler_iso->eqn.vol_term = vol;
  diffusion_euler_iso->eqn.surf_term = surf;
  //advection->eqn.boundary_surf_term = boundary_surf;

  diffusion_euler_iso->vol = CK(vol_kernels, cdim, poly_order);

  diffusion_euler_iso->surf[0] = CK(surf_x_kernels, cdim, poly_order);
  if (cdim>1)
    diffusion_euler_iso->surf[1] = CK(surf_y_kernels, cdim, poly_order);
  if (cdim>2)
    diffusion_euler_iso->surf[2] = CK(surf_z_kernels, cdim, poly_order);
}

struct gkyl_dg_eqn*
gkyl_dg_diffusion_euler_iso_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range)
{
  struct dg_diffusion_euler_iso* diffusion_euler_iso = (struct dg_diffusion_euler_iso*) gkyl_malloc(sizeof(struct dg_diffusion_euler_iso));

  // set basic parameters
  diffusion_euler_iso->eqn.num_equations = 4;
  diffusion_euler_iso->conf_range = *conf_range;

  diffusion_euler_iso->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(diffusion_euler_iso->eqn.flags);
  diffusion_euler_iso->eqn.ref_count = gkyl_ref_count_init(gkyl_diffusion_euler_iso_free);

  // copy the host struct to device struct
  struct dg_diffusion_euler_iso* diffusion_euler_iso_cu = (struct dg_diffusion_euler_iso*) gkyl_cu_malloc(sizeof(struct dg_diffusion_euler_iso));
  gkyl_cu_memcpy(diffusion_euler_iso_cu, diffusion_euler_iso, sizeof(struct dg_diffusion_euler_iso), GKYL_CU_MEMCPY_H2D);
  dg_diffusion_euler_iso_set_cu_dev_ptrs<<<1,1>>>(diffusion_euler_iso_cu, cbasis->b_type, cbasis->ndim, cbasis->poly_order);

  // set parent on_dev pointer
  diffusion_euler_iso->eqn.on_dev = &diffusion_euler_iso_cu->eqn;

  return &diffusion_euler_iso->eqn;
}
