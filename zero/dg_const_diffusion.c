#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_dg_const_diffusion.h>
#include <gkyl_dg_const_diffusion_priv.h>
#include <gkyl_util.h>

static void
dg_const_diffusion_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_const_diffusion* const_diffusion = container_of(base, struct dg_const_diffusion, eqn);
  gkyl_free(const_diffusion);
}

struct gkyl_dg_eqn*
gkyl_dg_const_diffusion_new(const struct gkyl_basis* basis, const double* D)
{
  struct dg_const_diffusion* const_diffusion = gkyl_malloc(sizeof(struct dg_const_diffusion));

  int dim = basis->ndim;
  int poly_order = basis->poly_order;

  const_diffusion->dim = dim;

  const_diffusion->eqn.num_equations = 1;
  const_diffusion->eqn.vol_term = vol;
  const_diffusion->eqn.surf_term = surf;
  //const_diffusion->eqn.boundary_surf_term = boundary_surf;

  const_diffusion->vol = CK(vol_kernels, dim, poly_order);

  const_diffusion->surf[0] = CK(surf_x_kernels, dim, poly_order);
  /* if (cdim>1) */
  /*   const_diffusion->surf[1] = CK(const_surf_y_kernels, dim, poly_order); */
  /* if (cdim>2) */
  /*   const_diffusion->surf[2] = CK(const_surf_z_kernels, dim, poly_order); */

  // ensure non-NULL pointers
  assert(const_diffusion->vol);
  for (int i=0; i<dim; ++i) assert(const_diffusion->surf[i]);

  const_diffusion->D = D; 

  // set reference counter
  const_diffusion->eqn.ref_count = (struct gkyl_ref_count) { dg_const_diffusion_free, 1 };
  
  return &const_diffusion->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_const_diffusion_cu_dev_new(const struct gkyl_basis* basis, const double* D)
{
  assert(false);
  return 0;
}

#endif
