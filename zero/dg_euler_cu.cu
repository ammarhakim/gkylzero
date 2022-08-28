/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_euler.h>    
#include <gkyl_dg_euler_priv.h>
}

#include <cassert>

#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_euler_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *u_i, const struct gkyl_array *p_ij)
{
  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  euler->auxfields.u_i = uvar;
  euler->auxfields.p_ij = pvar;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_euler_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_euler_auxfields auxin)
{
  gkyl_euler_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.u_i->on_dev, auxin.p_ij->on_dev);
}

__global__ void static
dg_euler_set_cu_dev_ptrs(struct dg_euler* euler, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  euler->auxfields.u_i = 0; 
  euler->auxfields.p_ij = 0; 

  const gkyl_dg_euler_vol_kern_list *vol_kernels;
  const gkyl_dg_euler_surf_kern_list *surf_x_kernels, *surf_y_kernels, *surf_z_kernels;  
  
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
  
  euler->eqn.vol_term = vol;
  euler->eqn.surf_term = surf;
  euler->eqn.boundary_surf_term = boundary_surf;

  euler->vol =  CK(vol_kernels, cdim, poly_order);

  euler->surf[0] = CK(surf_x_kernels, cdim, poly_order);
  if (cdim>1)
    euler->surf[1] = CK(surf_y_kernels, cdim, poly_order);
  if (cdim>2)
    euler->surf[2] = CK(surf_z_kernels, cdim, poly_order);
}

struct gkyl_dg_eqn*
gkyl_dg_euler_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range, 
  double gas_gamma)
{
  struct dg_euler *euler = (struct dg_euler*) gkyl_malloc(sizeof(struct dg_euler));

  // set basic parameters
  euler->eqn.num_equations = 5;
  euler->gas_gamma = gas_gamma;

  euler->conf_range = *conf_range;

  euler->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(euler->eqn.flags);
  euler->eqn.ref_count = gkyl_ref_count_init(gkyl_euler_free);

  // copy the host struct to device struct
  struct dg_euler *euler_cu = (struct dg_euler*) gkyl_cu_malloc(sizeof(struct dg_euler));
  gkyl_cu_memcpy(euler_cu, euler, sizeof(struct dg_euler), GKYL_CU_MEMCPY_H2D);
  dg_euler_set_cu_dev_ptrs<<<1,1>>>(euler_cu, cbasis->b_type, cbasis->ndim, cbasis->poly_order);

  // set parent on_dev pointer
  euler->eqn.on_dev = &euler_cu->eqn;

  return &euler->eqn;
}
