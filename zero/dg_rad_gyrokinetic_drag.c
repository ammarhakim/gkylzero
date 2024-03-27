#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_dg_rad_gyrokinetic_drag.h>
#include <gkyl_dg_rad_gyrokinetic_drag_priv.h>
#include <gkyl_util.h>

void
gkyl_rad_gyrokinetic_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag  = container_of(base, struct dg_rad_gyrokinetic_drag, eqn);

  if (GKYL_IS_CU_ALLOC(rad_gyrokinetic_drag->eqn.flags))
    gkyl_cu_free(rad_gyrokinetic_drag->eqn.on_dev);
  
  gkyl_free(rad_gyrokinetic_drag);
}

void
gkyl_rad_gyrokinetic_drag_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_rad_gyrokinetic_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
 if (gkyl_dg_eqn_is_cu_dev(eqn)) {
   gkyl_rad_gyrokinetic_drag_set_auxfields_cu(eqn->on_dev, auxin);
   return;
}
#endif

  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  rad_gyrokinetic_drag->auxfields.nvnu_surf = auxin.nvnu_surf;
  rad_gyrokinetic_drag->auxfields.nvnu = auxin.nvnu;
  rad_gyrokinetic_drag->auxfields.nvsqnu_surf = auxin.nvsqnu_surf;
  rad_gyrokinetic_drag->auxfields.nvsqnu = auxin.nvsqnu;
  rad_gyrokinetic_drag->auxfields.vtsq = auxin.vtsq;
  rad_gyrokinetic_drag->auxfields.vtsq_min = auxin.vtsq_min;
}

struct gkyl_dg_eqn*
gkyl_dg_rad_gyrokinetic_drag_new(const struct gkyl_basis *conf_basis, 
  const struct gkyl_basis *phase_basis, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
 if (use_gpu)    
   return gkyl_dg_rad_gyrokinetic_drag_cu_dev_new(conf_basis, phase_basis, phase_range, conf_range);
#endif
    
  struct dg_rad_gyrokinetic_drag* rad_gyrokinetic_drag = gkyl_malloc(sizeof(struct dg_rad_gyrokinetic_drag));

  int cdim = conf_basis->ndim, pdim = phase_basis->ndim, vdim = pdim-cdim;
  int poly_order = conf_basis->poly_order;
  
  rad_gyrokinetic_drag->cdim = cdim;
  rad_gyrokinetic_drag->pdim = pdim;

  rad_gyrokinetic_drag->eqn.num_equations = 1;
  rad_gyrokinetic_drag->eqn.surf_term = surf;
  rad_gyrokinetic_drag->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_rad_gyrokinetic_vol_kern_list *vol_kernels;
  const gkyl_dg_rad_gyrokinetic_surf_kern_list *surf_vpar_kernels, *surf_mu_kernels;
  const gkyl_dg_rad_gyrokinetic_boundary_surf_kern_list *boundary_surf_vpar_kernels, *boundary_surf_mu_kernels;
  
  switch (conf_basis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_vpar_kernels = ser_surf_vpar_kernels;
      surf_mu_kernels = ser_surf_mu_kernels;
      boundary_surf_vpar_kernels = ser_boundary_surf_vpar_kernels;
      boundary_surf_mu_kernels = ser_boundary_surf_mu_kernels;
      break;

    default:
      assert(false);
      break;    
  }  

  rad_gyrokinetic_drag->eqn.vol_term = CK(vol_kernels, cdim, vdim, poly_order);

  rad_gyrokinetic_drag->surf[0] = CK(surf_vpar_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    rad_gyrokinetic_drag->surf[1] = CK(surf_mu_kernels, cdim, vdim, poly_order);

  rad_gyrokinetic_drag->boundary_surf[0] = CK(boundary_surf_vpar_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    rad_gyrokinetic_drag->boundary_surf[1] = CK(boundary_surf_mu_kernels, cdim, vdim, poly_order);
  
  // ensure non-NULL pointers
  for (int i=0; i<vdim; ++i) assert(rad_gyrokinetic_drag->surf[i]);
  for (int i=0; i<vdim; ++i) assert(rad_gyrokinetic_drag->boundary_surf[i]);

  rad_gyrokinetic_drag->auxfields.nvnu_surf = 0;
  rad_gyrokinetic_drag->auxfields.nvnu = 0;
  rad_gyrokinetic_drag->auxfields.nvsqnu_surf = 0;
  rad_gyrokinetic_drag->auxfields.nvsqnu = 0;
  rad_gyrokinetic_drag->auxfields.vtsq = 0;
  rad_gyrokinetic_drag->phase_range = *phase_range;
  rad_gyrokinetic_drag->conf_range = *conf_range;

  rad_gyrokinetic_drag->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(rad_gyrokinetic_drag->eqn.flags);
  rad_gyrokinetic_drag->eqn.ref_count = gkyl_ref_count_init(gkyl_rad_gyrokinetic_free);
  rad_gyrokinetic_drag->eqn.on_dev = &rad_gyrokinetic_drag->eqn;

  return &rad_gyrokinetic_drag->eqn;
}


#ifndef GKYL_HAVE_CUDA
struct gkyl_dg_eqn*
gkyl_dg_rad_gyrokinetic_drag_cu_dev_new(const struct gkyl_basis* conf_basis, 
  const struct gkyl_basis* phase_basis, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range)
{
  assert(false);
  return 0;
}

#endif
