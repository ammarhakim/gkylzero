#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_diffusion.h>
#include <gkyl_dg_diffusion_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polynomial order
#define CK(lst, cdim, poly_order) lst[cdim-1].kernels[poly_order]

void
gkyl_diffusion_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_diffusion* diffusion = container_of(base, struct dg_diffusion, eqn);
  gkyl_free(diffusion);
}

void
gkyl_diffusion_set_auxfields(const struct gkyl_dg_eqn* eqn, struct gkyl_dg_diffusion_auxfields auxin)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  diffusion->auxfields.D = auxin.D;
}

struct gkyl_dg_eqn*
gkyl_dg_diffusion_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range, bool use_gpu)
{
  struct dg_diffusion* diffusion = gkyl_malloc(sizeof(struct dg_diffusion));

  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;

  diffusion->eqn.num_equations = 1;
  diffusion->eqn.vol_term = vol;
  diffusion->eqn.surf_term = surf;
  //diffusion->eqn.boundary_surf_term = boundary_surf;

  diffusion->vol = CK(vol_kernels, cdim, poly_order);

  diffusion->surf[0] = CK(surf_x_kernels, cdim, poly_order);
  if (cdim>1)
    diffusion->surf[1] = CK(surf_y_kernels, cdim, poly_order);
  /* if (cdim>2) */
  /*   diffusion->surf[2] = CK(const_surf_z_kernels, cdim, poly_order); */

  // ensure non-NULL pointers
  assert(diffusion->vol);
  for (int i=0; i<cdim; ++i) assert(diffusion->surf[i]);

  diffusion->auxfields.D = 0;
  diffusion->conf_range = *conf_range;

  diffusion->eqn.flags = 0;
  diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_diffusion_free);
  diffusion->eqn.on_dev = &diffusion->eqn;
  
  return &diffusion->eqn;
}
