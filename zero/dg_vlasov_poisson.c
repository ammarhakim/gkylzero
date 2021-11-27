#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_dg_vlasov_poisson.h>
#include <gkyl_dg_vlasov_poisson_priv.h>
#include <gkyl_util.h>

static void
vlasov_poisson_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_vlasov_poisson *vlasov_poisson = container_of(base, struct dg_vlasov_poisson, eqn);
  gkyl_free(vlasov_poisson);
}

void
gkyl_vlasov_poisson_set_fac_phi(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *fac_phi)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(fac_phi)) {gkyl_vlasov_poisson_set_fac_phi_cu(eqn, fac_phi); return;}
#endif

  struct dg_vlasov_poisson *vlasov_poisson = container_of(eqn, struct dg_vlasov_poisson, eqn);
  vlasov_poisson->fac_phi = fac_phi;
}

void
gkyl_vlasov_poisson_set_vecA(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *vecA)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(vecA)) {gkyl_vlasov_poisson_set_vecA_cu(eqn, vecA); return;}
#endif

  struct dg_vlasov_poisson *vlasov_poisson = container_of(eqn, struct dg_vlasov_poisson, eqn);
  vlasov_poisson->vecA = vecA;
}

struct gkyl_dg_eqn*
gkyl_dg_vlasov_poisson_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, enum gkyl_field_id field_id)
{
  struct dg_vlasov_poisson *vlasov_poisson = gkyl_malloc(sizeof(struct dg_vlasov_poisson));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  vlasov_poisson->cdim = cdim;
  vlasov_poisson->pdim = pdim;

  vlasov_poisson->eqn.num_equations = 1;
  vlasov_poisson->eqn.vol_term = vol;
  vlasov_poisson->eqn.surf_term = surf;
  vlasov_poisson->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_vlasov_poisson_vol_kern_list *vol_kernels;
  const gkyl_dg_vlasov_poisson_extem_vol_kern_list *extem_vol_kernels;

  const gkyl_dg_vlasov_poisson_steam_surf_kern_list *stream_surf_x_kernels, *stream_surf_y_kernels, *stream_surf_z_kernels;

  const gkyl_dg_vlasov_poisson_accel_surf_kern_list *accel_surf_vx_kernels, *accel_surf_vy_kernels, *accel_surf_vz_kernels;
  const gkyl_dg_vlasov_poisson_extem_accel_surf_kern_list *extem_accel_surf_vx_kernels, *extem_accel_surf_vy_kernels, *extem_accel_surf_vz_kernels;

  const gkyl_dg_vlasov_poisson_accel_boundary_surf_kern_list *accel_boundary_surf_vx_kernels, *accel_boundary_surf_vy_kernels,
    *accel_boundary_surf_vz_kernels;
  const gkyl_dg_vlasov_poisson_extem_accel_boundary_surf_kern_list *extem_accel_boundary_surf_vx_kernels, *extem_accel_boundary_surf_vy_kernels,
    *extem_accel_boundary_surf_vz_kernels;
  
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      extem_vol_kernels = ser_extem_vol_kernels;
      stream_surf_x_kernels = ser_stream_surf_x_kernels;
      stream_surf_y_kernels = ser_stream_surf_y_kernels;
      stream_surf_z_kernels = ser_stream_surf_z_kernels;
      accel_surf_vx_kernels = ser_accel_surf_vx_kernels;
      extem_accel_surf_vx_kernels = ser_extem_accel_surf_vx_kernels;
      accel_surf_vy_kernels = ser_accel_surf_vy_kernels;
      extem_accel_surf_vy_kernels = ser_extem_accel_surf_vy_kernels;
      accel_surf_vz_kernels = ser_accel_surf_vz_kernels;
      extem_accel_surf_vz_kernels = ser_extem_accel_surf_vz_kernels;
      accel_boundary_surf_vx_kernels = ser_accel_boundary_surf_vx_kernels;
      extem_accel_boundary_surf_vx_kernels = ser_extem_accel_boundary_surf_vx_kernels;
      accel_boundary_surf_vy_kernels = ser_accel_boundary_surf_vy_kernels;
      extem_accel_boundary_surf_vy_kernels = ser_extem_accel_boundary_surf_vy_kernels;
      accel_boundary_surf_vz_kernels = ser_accel_boundary_surf_vz_kernels;
      extem_accel_boundary_surf_vz_kernels = ser_extem_accel_boundary_surf_vz_kernels;
      
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      vol_kernels = ten_vol_kernels;
      extem_vol_kernels = ten_extem_vol_kernels;
      stream_surf_x_kernels = ten_stream_surf_x_kernels;
      stream_surf_y_kernels = ten_stream_surf_y_kernels;
      stream_surf_z_kernels = ten_stream_surf_z_kernels;
      accel_surf_vx_kernels = ten_accel_surf_vx_kernels;
      extem_accel_surf_vx_kernels = ten_extem_accel_surf_vx_kernels;
      accel_surf_vy_kernels = ten_accel_surf_vy_kernels;
      extem_accel_surf_vy_kernels = ten_extem_accel_surf_vy_kernels;
      accel_surf_vz_kernels = ten_accel_surf_vz_kernels;
      extem_accel_surf_vz_kernels = ten_extem_accel_surf_vz_kernels;
      accel_boundary_surf_vx_kernels = ten_accel_boundary_surf_vx_kernels;
      extem_accel_boundary_surf_vx_kernels = ten_extem_accel_boundary_surf_vx_kernels;
      accel_boundary_surf_vy_kernels = ten_accel_boundary_surf_vy_kernels;
      extem_accel_boundary_surf_vy_kernels = ten_extem_accel_boundary_surf_vy_kernels;
      accel_boundary_surf_vz_kernels = ten_accel_boundary_surf_vz_kernels;
      extem_accel_boundary_surf_vz_kernels = ten_extem_accel_boundary_surf_vz_kernels;

      break;

    default:
      assert(false);
      break;    
  }  

  if (field_id == GKYL_FIELD_PHI_A) {
    vlasov_poisson->vol = extem_vol_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    vlasov_poisson->accel_surf[0] = extem_accel_surf_vx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    if (vdim>1)
      vlasov_poisson->accel_surf[1] = extem_accel_surf_vy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    if (vdim>2)
      vlasov_poisson->accel_surf[2] = extem_accel_surf_vz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];

    vlasov_poisson->accel_boundary_surf[0] = extem_accel_boundary_surf_vx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    if (vdim>1)
      vlasov_poisson->accel_boundary_surf[1] = extem_accel_boundary_surf_vy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    if (vdim>2)
      vlasov_poisson->accel_boundary_surf[2] = extem_accel_boundary_surf_vz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  }
  else {
    vlasov_poisson->vol = extem_vol_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    vlasov_poisson->accel_surf[0] = accel_surf_vx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    if (vdim>1)
      vlasov_poisson->accel_surf[1] = accel_surf_vy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    if (vdim>2)
      vlasov_poisson->accel_surf[2] = accel_surf_vz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];

    vlasov_poisson->accel_boundary_surf[0] = accel_boundary_surf_vx_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    if (vdim>1)
      vlasov_poisson->accel_boundary_surf[1] = accel_boundary_surf_vy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    if (vdim>2)
      vlasov_poisson->accel_boundary_surf[2] = accel_boundary_surf_vz_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  }

  vlasov_poisson->stream_surf[0] = stream_surf_x_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  if (cdim>1)
    vlasov_poisson->stream_surf[1] = stream_surf_y_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  if (cdim>2)
    vlasov_poisson->stream_surf[2] = stream_surf_z_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];

  // ensure non-NULL pointers
  assert(vlasov_poisson->vol);
  for (int i=0; i<cdim; ++i) assert(vlasov_poisson->stream_surf[i]);
  for (int i=0; i<vdim; ++i) assert(vlasov_poisson->accel_surf[i]);
  for (int i=0; i<vdim; ++i) assert(vlasov_poisson->accel_boundary_surf[i]);

  vlasov_poisson->fac_phi = 0;
  vlasov_poisson->vecA = 0; 
  vlasov_poisson->conf_range = *conf_range;

  // set reference counter
  vlasov_poisson->eqn.ref_count = (struct gkyl_ref_count) { vlasov_poisson_free, 1 };
  
  return &vlasov_poisson->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_vlasov_poisson_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, enum gkyl_field_id field_id)
{
  assert(false);
  return 0;
}

#endif
