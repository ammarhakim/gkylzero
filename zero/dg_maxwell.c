#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_maxwell_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

void 
gkyl_maxwell_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_maxwell *maxwell = container_of(base->on_dev, struct dg_maxwell, eqn);
    gkyl_cu_free(maxwell);
  }  
  
  struct dg_maxwell *maxwell = container_of(base, struct dg_maxwell, eqn);
  gkyl_free(maxwell);
}

struct gkyl_array_copy_func*
gkyl_maxwell_wall_bc_create(const struct gkyl_dg_eqn *eqn, int dir, const struct gkyl_basis* cbasis)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_dg_eqn_is_cu_dev(eqn)) {
    return gkyl_maxwell_wall_bc_create_cu(eqn->on_dev, dir, cbasis);
  }
#endif

  struct dg_maxwell *maxwell = container_of(eqn, struct dg_maxwell, eqn);

  struct maxwell_wall_bc_ctx *ctx = (struct maxwell_wall_bc_ctx*) gkyl_malloc(sizeof(struct maxwell_wall_bc_ctx));
  ctx->dir = dir;
  ctx->basis = cbasis;

  struct gkyl_array_copy_func *bc = (struct gkyl_array_copy_func*) gkyl_malloc(sizeof(struct gkyl_array_copy_func));
  bc->func = maxwell->wall_bc;
  bc->ctx = ctx;

  bc->flags = 0;
  GKYL_CLEAR_CU_ALLOC(bc->flags);
  bc->on_dev = bc; // CPU eqn obj points to itself
  return bc;
}

void
gkyl_maxwell_wall_bc_release(struct gkyl_array_copy_func* bc)
{
  if (gkyl_array_copy_func_is_cu_dev(bc)) {
    //I think the context also needs to be freed but this produces
    //a segmentation faults (JJ 05/14/2022)
    //gkyl_cu_free(bc->on_dev->ctx);
    gkyl_cu_free(bc->on_dev);
  }
  gkyl_free(bc->ctx);
  gkyl_free(bc);
}

struct gkyl_dg_eqn*
gkyl_dg_maxwell_new(const struct gkyl_basis* cbasis,
  double lightSpeed, double elcErrorSpeedFactor, double mgnErrorSpeedFactor, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_maxwell_cu_dev_new(cbasis, lightSpeed, elcErrorSpeedFactor, mgnErrorSpeedFactor);
  } 
#endif
  struct dg_maxwell *maxwell = gkyl_malloc(sizeof(struct dg_maxwell));

  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;

  const gkyl_dg_maxwell_vol_kern_list *vol_kernels;
  const gkyl_dg_maxwell_surf_kern_list *surf_x_kernels, *surf_y_kernels, *surf_z_kernels;

  switch (cbasis->b_type) {
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
    
  maxwell->eqn.num_equations = 8;
  maxwell->eqn.vol_term = vol;
  maxwell->eqn.surf_term = surf;
  maxwell->eqn.boundary_surf_term = boundary_surf;

  maxwell->maxwell_data.c = lightSpeed;
  maxwell->maxwell_data.chi = lightSpeed*elcErrorSpeedFactor;
  maxwell->maxwell_data.gamma = lightSpeed*mgnErrorSpeedFactor;

  maxwell->vol =  CK(vol_kernels, cdim, poly_order);
  assert(maxwell->vol);

  maxwell->surf[0] = CK(surf_x_kernels, cdim, poly_order);
  if (cdim>1)
    maxwell->surf[1] = CK(surf_y_kernels, cdim, poly_order);
  if (cdim>2)
    maxwell->surf[2] = CK(surf_z_kernels, cdim, poly_order);

  // setup pointer for wall BC function
  maxwell->wall_bc = maxwell_wall_bc;

  // ensure non-NULL pointers 
  for (int i=0; i<cdim; ++i) assert(maxwell->surf[i]);

  maxwell->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(maxwell->eqn.flags);
  maxwell->eqn.ref_count = gkyl_ref_count_init(gkyl_maxwell_free);
  maxwell->eqn.on_dev = &maxwell->eqn; // CPU eqn obj points to itself
  
  return &maxwell->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_maxwell_cu_dev_new(const struct gkyl_basis* cbasis,
  double lightSpeed, double elcErrorSpeedFactor, double mgnErrorSpeedFactor)
{
  assert(false);
  return 0;
}

#endif
