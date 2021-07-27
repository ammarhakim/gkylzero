#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_dg_maxwell.h>
#include <gkyl_dg_maxwell_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

static void
maxwell_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_maxwell *maxwell = container_of(base, struct dg_maxwell, eqn);
  free(maxwell);
}

struct gkyl_dg_eqn*
gkyl_dg_maxwell_new(const struct gkyl_basis* cbasis,
  double lightSpeed, double elcErrorSpeedFactor, double mgnErrorSpeedFactor)
{
  struct dg_maxwell *maxwell = gkyl_malloc(sizeof(struct dg_maxwell));

  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;

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

  // ensure non-NULL pointers 
  for (int i=0; i<cdim; ++i) assert(maxwell->surf[i]);

  // set reference counter
  maxwell->eqn.ref_count = (struct gkyl_ref_count) { maxwell_free, 1 };
  
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
