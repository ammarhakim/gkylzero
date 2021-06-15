#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_dg_vlasov.h>
#include <gkyl_dg_vlasov_priv.h>
#include <gkyl_util.h>

static void
vlasov_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_vlasov *vlasov = container_of(base, struct dg_vlasov, eqn);
  gkyl_free(vlasov);
}

void
gkyl_vlasov_set_qmem(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *qmem)
{
#ifdef GKYL_HAVE_CUDA
  if(gkyl_array_is_cu_dev(qmem)) {gkyl_vlasov_set_qmem_cu(eqn, qmem); return;}
#endif

  struct dg_vlasov *vlasov = container_of(eqn, struct dg_vlasov, eqn);
  vlasov->qmem = qmem;
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range)
{
  struct dg_vlasov *vlasov = gkyl_malloc(sizeof(struct dg_vlasov));

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

  // ensure non-NULL pointers
  assert(vlasov->vol);
  for (int i=0; i<cdim; ++i) assert(vlasov->stream_surf[i]);
  for (int i=0; i<vdim; ++i) assert(vlasov->accel_surf[i]);
  for (int i=0; i<vdim; ++i) assert(vlasov->accel_boundary_surf[i]);

  vlasov->qmem = 0; 
  vlasov->conf_range = *conf_range;

  // set reference counter
  vlasov->eqn.ref_count = (struct gkyl_ref_count) { vlasov_free, 1 };
  
  return &vlasov->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range)
{
  assert(false);
  return 0;
}

#endif
