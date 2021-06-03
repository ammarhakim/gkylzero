/* -*- c -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_dg_vlasov.h>    
#include <gkyl_dg_vlasov_priv.h>
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range)
{
  struct dg_vlasov *vlasov = (struct dg_vlasov*) gkyl_malloc(sizeof(struct dg_vlasov));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int polyOrder = cbasis->polyOrder;

  vlasov->cdim = cdim;
  vlasov->pdim = pdim;

  vlasov->eqn.num_equations = 1;
  vlasov->eqn.vol_term = vol;
  vlasov->eqn.surf_term = surf;
  vlasov->eqn.boundary_surf_term = boundary_surf;

  vlasov->vol = CK(vol_kernels,cdim,vdim,polyOrder);

  vlasov->stream_surf[0] = CK(stream_surf_x_kernels,cdim,vdim,polyOrder);
  if (cdim>1)
    vlasov->stream_surf[1] = CK(stream_surf_y_kernels,cdim,vdim,polyOrder);
  if (cdim>2)
    vlasov->stream_surf[2] = CK(stream_surf_z_kernels,cdim,vdim,polyOrder);

  vlasov->accel_surf[0] = CK(accel_surf_vx_kernels,cdim,vdim,polyOrder);
  if (vdim>1)
    vlasov->accel_surf[1] = CK(accel_surf_vy_kernels,cdim,vdim,polyOrder);
  if (vdim>2)
    vlasov->accel_surf[2] = CK(accel_surf_vz_kernels,cdim,vdim,polyOrder);

  vlasov->accel_boundary_surf[0] = CK(accel_boundary_surf_vx_kernels,cdim,vdim,polyOrder);
  if (vdim>1)
    vlasov->accel_boundary_surf[1] = CK(accel_boundary_surf_vy_kernels,cdim,vdim,polyOrder);
  if (vdim>2)
    vlasov->accel_boundary_surf[2] = CK(accel_boundary_surf_vz_kernels,cdim,vdim,polyOrder);

  vlasov->qmem = 0; 
  vlasov->conf_range = *conf_range;
  
  return &vlasov->eqn;
}
