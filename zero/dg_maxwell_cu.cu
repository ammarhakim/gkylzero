/* -*- c -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_dg_maxwell.h>    
#include <gkyl_dg_maxwell_priv.h>
}

#define CK(lst,cdim,polyOrder) lst[cdim-1].kernels[polyOrder]

__global__
void dg_maxwell_set_cu_dev_ptrs(struct dg_maxwell* maxwell, int cdim, int polyOrder)
{
  maxwell->eqn.vol_term = vol;
  maxwell->eqn.surf_term = surf;
  maxwell->eqn.boundary_surf_term = boundary_surf;

  maxwell->vol =  CK(vol_kernels, cdim, polyOrder);

  maxwell->surf[0] = CK(surf_x_kernels, cdim, polyOrder);
  if (cdim>1)
    maxwell->surf[1] = CK(surf_y_kernels, cdim, polyOrder);
  if (cdim>2)
    maxwell->surf[2] = CK(surf_z_kernels, cdim, polyOrder);
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

  // copy the host struct to device struct
  struct dg_maxwell *maxwell_cu = (struct dg_maxwell*) gkyl_cu_malloc(sizeof(struct dg_maxwell));
  gkyl_cu_memcpy(maxwell_cu, maxwell, sizeof(struct dg_maxwell), GKYL_CU_MEMCPY_H2D);

  dg_maxwell_set_cu_dev_ptrs<<<1,1>>>(maxwell_cu, cbasis->ndim, cbasis->polyOrder);

  gkyl_free(maxwell);
  
  return &maxwell_cu->eqn;
}
