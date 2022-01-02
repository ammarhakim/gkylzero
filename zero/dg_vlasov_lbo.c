#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_dg_vlasov_lbo.h>
#include <gkyl_dg_vlasov_lbo_priv.h>
#include <gkyl_util.h>

static void
dg_vlasov_lbo_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_vlasov_lbo *vlasov_lbo  = container_of(base, struct dg_vlasov_lbo, eqn);
  gkyl_free(vlasov_lbo);
}

void
gkyl_vlasov_lbo_set_nuSum(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuSum)
{
//#ifdef GKYL_HAVE_CUDA
//  if (gkyl_array_is_cu_dev(nuSum)) {gkyl_vlasov_lbo_set_nuSum_cu(eqn, nuSum); return;}
//#endif

  struct dg_vlasov_lbo *vlasov_lbo = container_of(eqn, struct dg_vlasov_lbo, eqn);
  vlasov_lbo->nuSum = nuSum;
}

void
gkyl_vlasov_lbo_set_nuUSum(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuUSum)
{
//#ifdef GKYL_HAVE_CUDA
//  if (gkyl_array_is_cu_dev(nuUSum)) {gkyl_vlasov_lbo_set_nuUSum_cu(eqn, nuUSum); return;}
//#endif

  struct dg_vlasov_lbo *vlasov_lbo = container_of(eqn, struct dg_vlasov_lbo, eqn);
  vlasov_lbo->nuUSum = nuUSum;
}


void
gkyl_vlasov_lbo_set_nuVtSqSum(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nuVtSqSum)
{
//#ifdef GKYL_HAVE_CUDA
//  if (gkyl_array_is_cu_dev(nuVtSqSum)) {gkyl_vlasov_lbo_set_nuVtSqSum_cu(eqn, nuVtSqSum); return;}
//#endif

  struct dg_vlasov_lbo *vlasov_lbo = container_of(eqn, struct dg_vlasov_lbo, eqn);
  vlasov_lbo->nuVtSqSum = nuVtSqSum;
}


struct gkyl_dg_eqn*
gkyl_dg_vlasov_lbo_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range)
{
  struct dg_vlasov_lbo* vlasov_lbo = gkyl_malloc(sizeof(struct dg_vlasov_lbo));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  vlasov_lbo->cdim = cdim;
  vlasov_lbo->pdim = pdim;

  vlasov_lbo->eqn.num_equations = 1;
  vlasov_lbo->eqn.vol_term = vol;
  vlasov_lbo->eqn.surf_term = surf;
  vlasov_lbo->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_vlasov_lbo_vol_kern_list *vol_kernels;
  const gkyl_dg_vlasov_lbo_surf_kern_list *surf_vx_kernels, *surf_vy_kernels, *surf_vz_kernels;
  const gkyl_dg_vlasov_lbo_boundary_surf_kern_list *boundary_surf_vx_kernels, *boundary_surf_vy_kernels,
    *boundary_surf_vz_kernels;
  
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_vx_kernels = ser_surf_vx_kernels;
      surf_vy_kernels = ser_surf_vy_kernels;
      surf_vz_kernels = ser_surf_vz_kernels;
      boundary_surf_vx_kernels = ser_boundary_surf_vx_kernels;
      boundary_surf_vy_kernels = ser_boundary_surf_vy_kernels;
      boundary_surf_vz_kernels = ser_boundary_surf_vz_kernels;
      break;

    /* case GKYL_BASIS_MODAL_TENSOR: */
    /*   vol_kernels = ten_vol_kernels; */
    /*   surf_vx_kernels = ten_surf_vx_kernels; */
    /*   surf_vy_kernels = ten_surf_vy_kernels; */
    /*   surf_vz_kernels = ten_surf_vz_kernels; */
    /*   boundary_surf_vx_kernels = ten_boundary_surf_vx_kernels; */
    /*   boundary_surf_vy_kernels = ten_boundary_surf_vy_kernels; */
    /*   boundary_surf_vz_kernels = ten_boundary_surf_vz_kernels; */
    /*   break; */

    default:
      assert(false);
      break;    
  }  

  vlasov_lbo->vol = CK(vol_kernels, cdim, vdim, poly_order);

  vlasov_lbo->surf[0] = CK(surf_vx_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    vlasov_lbo->surf[1] = CK(surf_vy_kernels, cdim, vdim, poly_order);
  if (vdim>2)
    vlasov_lbo->surf[2] = CK(surf_vz_kernels, cdim, vdim, poly_order);

  vlasov_lbo->boundary_surf[0] = CK(boundary_surf_vx_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    vlasov_lbo->boundary_surf[1] = CK(boundary_surf_vy_kernels, cdim, vdim, poly_order);
  if (vdim>2)
    vlasov_lbo->boundary_surf[2] = CK(boundary_surf_vz_kernels, cdim, vdim, poly_order);

  // ensure non-NULL pointers
  assert(vlasov_lbo->vol);
  for (int i=0; i<vdim; ++i) assert(vlasov_lbo->surf[i]);
  for (int i=0; i<vdim; ++i) assert(vlasov_lbo->boundary_surf[i]);

  vlasov_lbo->nuSum = 0;
  vlasov_lbo->nuUSum = 0;
  vlasov_lbo->nuVtSqSum = 0;
  vlasov_lbo->conf_range = *conf_range;

  // set reference counter
  vlasov_lbo->eqn.ref_count = gkyl_ref_count_init(dg_vlasov_lbo_free);
  
  return &vlasov_lbo->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_vlasov_lbo_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range)
{
  assert(false);
  return 0;
}

#endif
