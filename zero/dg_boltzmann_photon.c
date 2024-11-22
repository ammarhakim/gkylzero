#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_boltzmann_photon.h>
#include <gkyl_dg_boltzmann_photon_priv.h>
#include <gkyl_util.h>

void
gkyl_boltzmann_photon_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  
  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_boltzmann_photon *boltzmann_photon = container_of(base->on_dev, struct dg_boltzmann_photon, eqn);
    gkyl_cu_free(boltzmann_photon);
  }
  
  struct dg_boltzmann_photon *boltzmann_photon = container_of(base, struct dg_boltzmann_photon, eqn);
  gkyl_free(boltzmann_photon);
}

void
gkyl_boltzmann_photon_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_boltzmann_photon_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_dg_eqn_is_cu_dev(eqn)) {
    gkyl_boltzmann_photon_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif

  struct dg_boltzmann_photon *boltzmann_photon = container_of(eqn, struct dg_boltzmann_photon, eqn);
  boltzmann_photon->auxfields.kpar_abs = auxin.kpar_abs;
  boltzmann_photon->auxfields.jacob_vel_inv = auxin.jacob_vel_inv;
}

struct gkyl_dg_eqn*
gkyl_dg_boltzmann_photon_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* vel_range, double light_speed, double rho_curv, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_boltzmann_photon_cu_dev_new(cbasis, pbasis, vel_range, 
      light_speed, rho_curv);
  } 
#endif
  struct dg_boltzmann_photon *boltzmann_photon = gkyl_malloc(sizeof(struct dg_boltzmann_photon));


  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  // Equation system only works for 1x2v for now
  assert(cdim == 1);
  assert(vdim == 2);

  boltzmann_photon->cdim = cdim;
  boltzmann_photon->pdim = pdim;

  boltzmann_photon->eqn.num_equations = 1;
  boltzmann_photon->eqn.surf_term = surf;
  boltzmann_photon->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_boltzmann_photon_vol_kern_list *vol_kernels;
  const gkyl_dg_boltzmann_photon_surf_kern_list *surf_x_kernels, *surf_vx_kernels, *surf_vy_kernels;
  const gkyl_dg_boltzmann_photon_boundary_surf_kern_list *boundary_surf_x_kernels, *boundary_surf_vx_kernels, *boundary_surf_vy_kernels;
  
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_x_kernels = ser_surf_x_kernels;
      surf_vx_kernels = ser_surf_vx_kernels;
      surf_vy_kernels = ser_surf_vy_kernels;
      boundary_surf_x_kernels = ser_boundary_surf_x_kernels;
      boundary_surf_vx_kernels = ser_boundary_surf_vx_kernels;
      boundary_surf_vy_kernels = ser_boundary_surf_vy_kernels;
      
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      vol_kernels = tensor_vol_kernels;
      surf_x_kernels = tensor_surf_x_kernels;
      surf_vx_kernels = tensor_surf_vx_kernels;
      surf_vy_kernels = tensor_surf_vy_kernels;
      boundary_surf_x_kernels = tensor_boundary_surf_x_kernels;
      boundary_surf_vx_kernels = tensor_boundary_surf_vx_kernels;
      boundary_surf_vy_kernels = tensor_boundary_surf_vy_kernels;
      
      break;

    default:
      assert(false);
      break;    
  }  
  boltzmann_photon->eqn.vol_term = CK(vol_kernels,cdim,vdim,poly_order);

  boltzmann_photon->surf[0] = CK(surf_x_kernels,cdim,vdim,poly_order);
  boltzmann_photon->surf[1] = CK(surf_vx_kernels,cdim,vdim,poly_order);
  boltzmann_photon->surf[2] = CK(surf_vy_kernels,cdim,vdim,poly_order);

  boltzmann_photon->boundary_surf[0] = CK(boundary_surf_x_kernels,cdim,vdim,poly_order);
  boltzmann_photon->boundary_surf[1] = CK(boundary_surf_vx_kernels,cdim,vdim,poly_order);
  boltzmann_photon->boundary_surf[2] = CK(boundary_surf_vy_kernels,cdim,vdim,poly_order);

  // ensure non-NULL pointers
  for (int i=0; i<cdim+vdim; ++i) assert(boltzmann_photon->surf[i]);
  for (int i=0; i<vdim; ++i) assert(boltzmann_photon->boundary_surf[i]);

  boltzmann_photon->vel_range = *vel_range;
  boltzmann_photon->auxfields.kpar_abs = 0;
  boltzmann_photon->auxfields.jacob_vel_inv = 0;

  boltzmann_photon->light_speed = light_speed;
  boltzmann_photon->rho_curv = rho_curv;  

  boltzmann_photon->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(boltzmann_photon->eqn.flags);

  boltzmann_photon->eqn.ref_count = gkyl_ref_count_init(gkyl_boltzmann_photon_free);
  boltzmann_photon->eqn.on_dev = &boltzmann_photon->eqn; // CPU eqn obj points to itself
  
  return &boltzmann_photon->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_boltzmann_photon_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* vel_range, double light_speed, double rho_curv)
{
  assert(false);
  return 0;
}

#endif
