/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_maxwell.h>    
#include <gkyl_dg_maxwell_priv.h>
}

#include <cassert>

#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

__global__ static void
gkyl_maxwell_wall_bc_create_set_cu_dev_ptrs(const struct gkyl_dg_eqn *eqn, int dir,
  const struct gkyl_basis* cbasis, struct maxwell_wall_bc_ctx *ctx, struct gkyl_array_copy_func *bc)
{
  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);

  ctx->dir = dir;
  ctx->basis = cbasis;

  bc->func = maxwell->wall_bc;
  bc->ctx = ctx;
}

struct gkyl_array_copy_func*
gkyl_maxwell_wall_bc_create_cu(const struct gkyl_dg_eqn *eqn, int dir, const struct gkyl_basis* cbasis)
{
  // create host context and bc func structs
  struct maxwell_wall_bc_ctx *ctx = (struct maxwell_wall_bc_ctx*) gkyl_malloc(sizeof(struct maxwell_wall_bc_ctx));
  struct gkyl_array_copy_func *bc = (struct gkyl_array_copy_func*) gkyl_malloc(sizeof(struct gkyl_array_copy_func));
  bc->ctx = ctx;

  bc->flags = 0;
  GKYL_SET_CU_ALLOC(bc->flags);

  // create device context and bc func structs
  struct maxwell_wall_bc_ctx *ctx_cu = (struct maxwell_wall_bc_ctx*) gkyl_cu_malloc(sizeof(struct maxwell_wall_bc_ctx));
  struct gkyl_array_copy_func *bc_cu = (struct gkyl_array_copy_func*) gkyl_cu_malloc(sizeof(struct gkyl_array_copy_func));

  gkyl_cu_memcpy(ctx_cu, ctx, sizeof(struct maxwell_wall_bc_ctx), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(bc_cu, bc, sizeof(struct gkyl_array_copy_func), GKYL_CU_MEMCPY_H2D);

  gkyl_maxwell_wall_bc_create_set_cu_dev_ptrs<<<1,1>>>(eqn, dir, cbasis, ctx_cu, bc_cu);

  // set parent on_dev pointer 
  bc->on_dev = bc_cu;  
  return bc;
}

__global__ void static
dg_maxwell_set_cu_dev_ptrs(struct dg_maxwell* maxwell, enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  const gkyl_dg_maxwell_vol_kern_list *vol_kernels;
  const gkyl_dg_maxwell_surf_kern_list *surf_x_kernels, *surf_y_kernels, *surf_z_kernels;  
  
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
  
  maxwell->eqn.vol_term = vol;
  maxwell->eqn.surf_term = surf;
  maxwell->eqn.boundary_surf_term = boundary_surf;

  maxwell->vol =  CK(vol_kernels, cdim, poly_order);

  maxwell->surf[0] = CK(surf_x_kernels, cdim, poly_order);
  if (cdim>1)
    maxwell->surf[1] = CK(surf_y_kernels, cdim, poly_order);
  if (cdim>2)
    maxwell->surf[2] = CK(surf_z_kernels, cdim, poly_order);

  // setup pointer for wall BC function
  maxwell->wall_bc = maxwell_wall_bc;
}

struct gkyl_dg_eqn*
gkyl_dg_maxwell_cu_dev_new(const struct gkyl_basis* cbasis,
  double lightSpeed, double elcErrorSpeedFactor, double mgnErrorSpeedFactor)
{
  struct dg_maxwell *maxwell = (struct dg_maxwell*) gkyl_malloc(sizeof(struct dg_maxwell));

  // set basic parameters
  maxwell->eqn.num_equations = 8;
  maxwell->maxwell_data.c = lightSpeed;
  maxwell->maxwell_data.chi = lightSpeed*elcErrorSpeedFactor;
  maxwell->maxwell_data.gamma = lightSpeed*mgnErrorSpeedFactor;

  maxwell->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(maxwell->eqn.flags);
  maxwell->eqn.ref_count = gkyl_ref_count_init(gkyl_maxwell_free);

  // copy the host struct to device struct
  struct dg_maxwell *maxwell_cu = (struct dg_maxwell*) gkyl_cu_malloc(sizeof(struct dg_maxwell));
  gkyl_cu_memcpy(maxwell_cu, maxwell, sizeof(struct dg_maxwell), GKYL_CU_MEMCPY_H2D);
  dg_maxwell_set_cu_dev_ptrs<<<1,1>>>(maxwell_cu, cbasis->b_type, cbasis->ndim, cbasis->poly_order);

  // set parent on_dev pointer
  maxwell->eqn.on_dev = &maxwell_cu->eqn;

  return &maxwell->eqn;
}
