#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_diffusion_euler_iso.h>
#include <gkyl_dg_diffusion_euler_iso_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim, vdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

void
gkyl_diffusion_euler_iso_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_diffusion_euler_iso* diffusion_euler_iso = container_of(base->on_dev, struct dg_diffusion_euler_iso, eqn);
    gkyl_cu_free(diffusion_euler_iso);
  }

  struct dg_diffusion_euler_iso* diffusion_euler_iso = container_of(base, struct dg_diffusion_euler_iso, eqn);
  gkyl_free(diffusion_euler_iso);
}

void
gkyl_diffusion_euler_iso_set_auxfields(const struct gkyl_dg_eqn* eqn, struct gkyl_dg_diffusion_euler_iso_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_dg_eqn_is_cu_dev(eqn)) {
    gkyl_diffusion_euler_iso_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif

  struct dg_diffusion_euler_iso* diffusion_euler_iso = container_of(eqn, struct dg_diffusion_euler_iso, eqn);
  diffusion_euler_iso->auxfields.D = auxin.D;
  diffusion_euler_iso->auxfields.u_i = auxin.u_i;
}

struct gkyl_dg_eqn*
gkyl_dg_diffusion_euler_iso_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu)
    return gkyl_dg_diffusion_euler_iso_cu_dev_new(cbasis, conf_range);
#endif

  struct dg_diffusion_euler_iso* diffusion_euler_iso = gkyl_malloc(sizeof(struct dg_diffusion_euler_iso));

  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;

  const gkyl_dg_diffusion_euler_iso_vol_kern_list* vol_kernels;
  const gkyl_dg_diffusion_euler_iso_surf_kern_list* surf_x_kernels;
  const gkyl_dg_diffusion_euler_iso_surf_kern_list* surf_y_kernels;
  const gkyl_dg_diffusion_euler_iso_surf_kern_list* surf_z_kernels;

  switch (cbasis->b_type) {
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

  diffusion_euler_iso->eqn.num_equations = 4;
  diffusion_euler_iso->eqn.surf_term = surf;
  //diffusion_euler_iso->eqn.boundary_surf_term = boundary_surf;

  diffusion_euler_iso->eqn.vol_term = CK(vol_kernels, cdim, poly_order);

  diffusion_euler_iso->surf[0] = CK(surf_x_kernels, cdim, poly_order);
  if (cdim>1)
    diffusion_euler_iso->surf[1] = CK(surf_y_kernels, cdim, poly_order);
  if (cdim>2)
    diffusion_euler_iso->surf[2] = CK(surf_z_kernels, cdim, poly_order);

  // ensure non-NULL pointers
  for (int i=0; i<cdim; ++i) assert(diffusion_euler_iso->surf[i]);

  diffusion_euler_iso->auxfields.D = 0;
  diffusion_euler_iso->auxfields.u_i = 0;
  diffusion_euler_iso->conf_range = *conf_range;

  diffusion_euler_iso->eqn.flags = 0;
  diffusion_euler_iso->eqn.ref_count = gkyl_ref_count_init(gkyl_diffusion_euler_iso_free);
  diffusion_euler_iso->eqn.on_dev = &diffusion_euler_iso->eqn;

  return &diffusion_euler_iso->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_diffusion_euler_iso_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range)
{
  assert(false);
  return 0;
}

#endif
