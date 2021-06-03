/* -*- c -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_dg_maxwell.h>    
#include <gkyl_dg_maxwell_priv.h>
}

#include <stdio.h>

// various pointers to functions on device
__device__ static vol_termf_t p_vol = &vol;
__device__ static surf_termf_t p_surf = &surf;
__device__ static boundary_surf_termf_t p_boundary_surf = &boundary_surf;

struct gkyl_dg_eqn*
gkyl_dg_maxwell_cu_dev_new(const struct gkyl_basis* cbasis,
  double lightSpeed, double elcErrorSpeedFactor, double mgnErrorSpeedFactor)
{
  struct dg_maxwell *maxwell = (struct dg_maxwell*) gkyl_malloc(sizeof(struct dg_maxwell));

  int cdim = cbasis->ndim;
  int polyOrder = cbasis->polyOrder;

  maxwell->eqn.num_equations = 8;

  cudaMemcpyFromSymbol(&maxwell->eqn.vol_term, p_vol, sizeof(vol_termf_t));
  cudaMemcpyFromSymbol(&maxwell->eqn.surf_term, p_surf, sizeof(surf_termf_t));
  cudaMemcpyFromSymbol(&maxwell->eqn.boundary_surf_term, p_boundary_surf, sizeof(boundary_surf_termf_t));

  maxwell->maxwell_data.c = lightSpeed;
  maxwell->maxwell_data.chi = lightSpeed*elcErrorSpeedFactor;
  maxwell->maxwell_data.gamma = lightSpeed*mgnErrorSpeedFactor;

  maxwell->vol = CK(vol_kernels, cdim, polyOrder);

  maxwell->surf[0] = CK(surf_x_kernels, cdim, polyOrder);
  if (cdim>1)
    maxwell->surf[1] = CK(surf_y_kernels, cdim, polyOrder);
  if (cdim>2)
    maxwell->surf[2] = CK(surf_z_kernels, cdim, polyOrder);

  // copy to device
  struct dg_maxwell *maxwell_cu = (struct dg_maxwell*) gkyl_cu_malloc(sizeof(struct dg_maxwell));
  gkyl_cu_memcpy(maxwell_cu, maxwell, sizeof(struct dg_maxwell), GKYL_CU_MEMCPY_H2D);

  // cleanup
  gkyl_free(maxwell);

  return &maxwell_cu->eqn;
}
