#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_gen_diffusion.h>
#include <gkyl_dg_gen_diffusion_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polynomial order
#define CK(lst, cdim, poly_order) lst[cdim-1].kernels[poly_order]

void
gkyl_gen_diffusion_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_gen_diffusion* gen_diffusion = container_of(base->on_dev, struct dg_gen_diffusion, eqn);
    gkyl_cu_free(gen_diffusion);
  }
  
  struct dg_gen_diffusion* gen_diffusion = container_of(base, struct dg_gen_diffusion, eqn);
  gkyl_free(gen_diffusion);
}

void
gkyl_gen_diffusion_set_auxfields(const struct gkyl_dg_eqn* eqn, struct gkyl_dg_gen_diffusion_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.Dij)) {
    gkyl_gen_diffusion_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif
  
  struct dg_gen_diffusion* gen_diffusion = container_of(eqn, struct dg_gen_diffusion, eqn);
  gen_diffusion->auxfields.Dij = auxin.Dij;
}

struct gkyl_dg_eqn*
gkyl_dg_gen_diffusion_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu)
    return gkyl_dg_gen_diffusion_cu_dev_new(cbasis, conf_range);
#endif
  
  struct dg_gen_diffusion* gen_diffusion = gkyl_malloc(sizeof(struct dg_gen_diffusion));

  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;

  const gkyl_dg_gen_diffusion_vol_kern_list* vol_kernels;
  const gkyl_dg_gen_diffusion_surf_kern_list* surf_xx_kernels;
  const gkyl_dg_gen_diffusion_surf_kern_list* surf_xy_kernels;
  const gkyl_dg_gen_diffusion_surf_kern_list* surf_xz_kernels;
  const gkyl_dg_gen_diffusion_surf_kern_list* surf_yx_kernels;
  const gkyl_dg_gen_diffusion_surf_kern_list* surf_yy_kernels;
  const gkyl_dg_gen_diffusion_surf_kern_list* surf_yz_kernels;
  const gkyl_dg_gen_diffusion_surf_kern_list* surf_zx_kernels;
  const gkyl_dg_gen_diffusion_surf_kern_list* surf_zy_kernels;
  const gkyl_dg_gen_diffusion_surf_kern_list* surf_zz_kernels; 

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_xx_kernels = ser_surf_xx_kernels;
      surf_xy_kernels = ser_surf_xy_kernels;
      surf_xz_kernels = ser_surf_xz_kernels;
      surf_yx_kernels = ser_surf_yx_kernels;
      surf_yy_kernels = ser_surf_yy_kernels;
      surf_yz_kernels = ser_surf_yz_kernels;
      surf_zx_kernels = ser_surf_zx_kernels;
      surf_zy_kernels = ser_surf_zy_kernels;
      surf_zz_kernels = ser_surf_zz_kernels;
      break;

    default:
      assert(false);
      break;    
  } 

  gen_diffusion->eqn.num_equations = 1;
  gen_diffusion->eqn.vol_term = vol;
  gen_diffusion->eqn.gen_surf_term = surf;

  gen_diffusion->vol = CK(vol_kernels, cdim, poly_order);

  gen_diffusion->surf[0][0] = CK(surf_xx_kernels, cdim, poly_order);
  if (cdim>1) {
    gen_diffusion->surf[0][1] = CK(surf_xy_kernels, cdim, poly_order);
    gen_diffusion->surf[1][0] = CK(surf_yx_kernels, cdim, poly_order);
    gen_diffusion->surf[1][1] = CK(surf_yy_kernels, cdim, poly_order);
  }
  if (cdim>2) {
    gen_diffusion->surf[0][2] = CK(surf_xz_kernels, cdim, poly_order);
    gen_diffusion->surf[1][2] = CK(surf_yz_kernels, cdim, poly_order);
    gen_diffusion->surf[2][0] = CK(surf_zx_kernels, cdim, poly_order);
    gen_diffusion->surf[2][1] = CK(surf_zy_kernels, cdim, poly_order);
    gen_diffusion->surf[2][2] = CK(surf_zz_kernels, cdim, poly_order);
  }

  // ensure non-NULL pointers
  assert(gen_diffusion->vol);
  for (int i=0; i<cdim; ++i) assert(gen_diffusion->surf[i]);

  gen_diffusion->auxfields.Dij = 0;
  gen_diffusion->conf_range = *conf_range;

  gen_diffusion->eqn.flags = 0;
  gen_diffusion->eqn.ref_count = gkyl_ref_count_init(gkyl_gen_diffusion_free);
  gen_diffusion->eqn.on_dev = &gen_diffusion->eqn;
  
  return &gen_diffusion->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_gen_diffusion_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range)
{
  assert(false);
  return 0;
}

#endif
