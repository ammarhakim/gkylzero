#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_dg_lbo.h>
#include <gkyl_dg_lbo_priv.h>
#include <gkyl_util.h>

static void
dg_lbo_free(const struct gkyl_ref_count* ref)
{
  struct gkyl_dg_eqn* base = container_of(ref, struct gkyl_dg_eqn, ref_count);
  struct dg_lbo *lbo  = container_of(base, struct dg_lbo, eqn);
  gkyl_free(lbo);
}

struct gkyl_dg_eqn*
gkyl_dg_lbo_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const double nuSum,
  const double* nuUSum_l, const double* nuUSum_r, const double* nuVtSqSum_l, const double* nuVtSqSum_r)
{
  struct dg_lbo* lbo = gkyl_malloc(sizeof(struct dg_lbo));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  lbo->cdim = cdim;
  lbo->pdim = pdim;

  lbo->eqn.num_equations = 1;
  //lbo->eqn.vol_term = vol;
  lbo->eqn.surf_term = surf;
  lbo->eqn.boundary_surf_term = boundary_surf;

  //lbo->vol = CK(vol_kernels, cdim, vdim, poly_order);

  lbo->surf[0] = CK(constNu_surf_vx_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    lbo->surf[1] = CK(constNu_surf_vy_kernels, cdim, vdim, poly_order);
  if (vdim>2)
    lbo->surf[2] = CK(constNu_surf_vz_kernels, cdim, vdim, poly_order);

  lbo->boundary_surf[0] = CK(constNu_boundary_surf_vx_kernels, cdim, vdim, poly_order);
  if (vdim>1)
    lbo->boundary_surf[1] = CK(constNu_boundary_surf_vy_kernels, cdim, vdim, poly_order);
  if (vdim>2)
    lbo->boundary_surf[2] = CK(constNu_boundary_surf_vz_kernels, cdim, vdim, poly_order);

  // ensure non-NULL pointers
  //assert(lbo->vol);
  for (int i=0; i<vdim; ++i) assert(lbo->surf[i]);

  lbo->nuSum = nuSum;
  lbo->nuUSum_l = nuUSum_l;
  lbo->nuUSum_r = nuUSum_r;
  lbo->nuVtSqSum_l = nuVtSqSum_l;
  lbo->nuVtSqSum_r = nuVtSqSum_r;

  // set reference counter
  lbo->eqn.ref_count = (struct gkyl_ref_count) { dg_lbo_free, 1 };
  
  return &lbo->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_lbo_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const double nuUSum,
  const double* nuUSum_l, const double* nuUSum_r, const double* nuVtSqSum_l, const double* nuVtSqSum_r)
{
  assert(false);
  return 0;
}

#endif
