/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_advection.h>    
#include <gkyl_dg_advection_priv.h>
}

#include <cassert>

#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

__global__ static void
gkyl_advection_absorb_bc_create_set_cu_dev_ptrs(const struct gkyl_dg_eqn *eqn, int dir,
  const struct gkyl_basis* cbasis, struct dg_bc_ctx *ctx, struct gkyl_array_copy_func *bc)
{
  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  
  ctx->basis = cbasis;

  bc->func = advection->absorb_bc;
  bc->ctx = ctx;
}

struct gkyl_array_copy_func*
gkyl_advection_absorb_bc_create_cu(const struct gkyl_dg_eqn *eqn, int dir, const struct gkyl_basis* cbasis)
{
  // create host context and bc func structs
  struct dg_bc_ctx *ctx = (struct dg_bc_ctx*) gkyl_malloc(sizeof(struct dg_bc_ctx));
  struct gkyl_array_copy_func *bc = (struct gkyl_array_copy_func*) gkyl_malloc(sizeof(struct gkyl_array_copy_func));
  bc->ctx = ctx;

  bc->flags = 0;
  GKYL_SET_CU_ALLOC(bc->flags);

  // create device context and bc func structs
  struct dg_bc_ctx *ctx_cu = (struct dg_bc_ctx*) gkyl_cu_malloc(sizeof(struct dg_bc_ctx));
  struct gkyl_array_copy_func *bc_cu = (struct gkyl_array_copy_func*) gkyl_cu_malloc(sizeof(struct gkyl_array_copy_func));

  gkyl_cu_memcpy(ctx_cu, ctx, sizeof(struct dg_bc_ctx), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(bc_cu, bc, sizeof(struct gkyl_array_copy_func), GKYL_CU_MEMCPY_H2D);

  bc->ctx_on_dev = ctx_cu;

  gkyl_advection_absorb_bc_create_set_cu_dev_ptrs<<<1,1>>>(eqn, dir, cbasis, ctx_cu, bc_cu);

  // set parent on_dev pointer 
  bc->on_dev = bc_cu;  
  return bc;
}

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_advection_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *u)
{
  struct dg_advection *advection = container_of(eqn, struct dg_advection, eqn);
  advection->auxfields.u = u;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_advection_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_advection_auxfields auxin)
{
  gkyl_advection_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.u->on_dev);
}

__global__ void static
dg_advection_set_cu_dev_ptrs(struct dg_advection* advection, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  advection->auxfields.u = 0; 

  const gkyl_dg_advection_vol_kern_list *vol_kernels;
  const gkyl_dg_advection_surf_kern_list *surf_x_kernels, *surf_y_kernels, *surf_z_kernels;  
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_x_kernels = ser_surf_x_kernels;
      surf_y_kernels = ser_surf_y_kernels;
      surf_z_kernels = ser_surf_z_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      vol_kernels = ten_vol_kernels;
      surf_x_kernels = ten_surf_x_kernels;
      surf_y_kernels = ten_surf_y_kernels;
      surf_z_kernels = ten_surf_z_kernels;
      break;

    default:
      assert(false);
      break;    
  }  
  
  advection->eqn.vol_term = vol;
  advection->eqn.surf_term = surf;
  advection->eqn.boundary_surf_term = boundary_surf;

  advection->vol =  CK(vol_kernels, cdim, poly_order);

  advection->surf[0] = CK(surf_x_kernels, cdim, poly_order);
  if (cdim>1)
    advection->surf[1] = CK(surf_y_kernels, cdim, poly_order);
  if (cdim>2)
    advection->surf[2] = CK(surf_z_kernels, cdim, poly_order);

  // setup pointer for absorbing BC function
  advection->absorb_bc = advection_absorb_bc;
}

struct gkyl_dg_eqn*
gkyl_dg_advection_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range)
{
  struct dg_advection *advection = (struct dg_advection*) gkyl_malloc(sizeof(struct dg_advection));

  // set basic parameters
  advection->eqn.num_equations = 1;
  advection->conf_range = *conf_range;

  advection->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(advection->eqn.flags);
  advection->eqn.ref_count = gkyl_ref_count_init(gkyl_advection_free);

  // copy the host struct to device struct
  struct dg_advection *advection_cu = (struct dg_advection*) gkyl_cu_malloc(sizeof(struct dg_advection));
  gkyl_cu_memcpy(advection_cu, advection, sizeof(struct dg_advection), GKYL_CU_MEMCPY_H2D);
  dg_advection_set_cu_dev_ptrs<<<1,1>>>(advection_cu, cbasis->b_type, cbasis->ndim, cbasis->poly_order);

  // set parent on_dev pointer
  advection->eqn.on_dev = &advection_cu->eqn;

  return &advection->eqn;
}
