#include "gkyl_dg_eqn.h"
#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_gyrokinetic.h>
#include <gkyl_dg_gyrokinetic_priv.h>
#include <gkyl_util.h>

void
gkyl_gyrokinetic_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_gyrokinetic *gyrokinetic = container_of(base->on_dev, struct dg_gyrokinetic, eqn);
    gkyl_cu_free(gyrokinetic);
  }

  struct dg_gyrokinetic *gyrokinetic = container_of(base, struct dg_gyrokinetic, eqn);
  gkyl_free(gyrokinetic);
}

void
gkyl_gyrokinetic_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_gyrokinetic_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.bmag) && gkyl_array_is_cu_dev(auxin.jacobtot_inv) &&
      gkyl_array_is_cu_dev(auxin.cmag) && gkyl_array_is_cu_dev(auxin.b_i) &&
      gkyl_array_is_cu_dev(auxin.phi) && gkyl_array_is_cu_dev(auxin.apar) &&
      gkyl_array_is_cu_dev(auxin.apardot) ) {
    gkyl_gyrokinetic_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif

  struct dg_gyrokinetic *gyrokinetic = container_of(eqn, struct dg_gyrokinetic, eqn);
  gyrokinetic->auxfields.bmag         = auxin.bmag;
  gyrokinetic->auxfields.jacobtot_inv = auxin.jacobtot_inv;
  gyrokinetic->auxfields.cmag         = auxin.cmag;
  gyrokinetic->auxfields.b_i          = auxin.b_i;
  gyrokinetic->auxfields.phi     = auxin.phi;
  gyrokinetic->auxfields.apar    = auxin.apar;
  gyrokinetic->auxfields.apardot = auxin.apardot;
}

struct gkyl_dg_eqn*
gkyl_dg_gyrokinetic_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, const double charge, const double mass, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    return gkyl_dg_gyrokinetic_cu_dev_new(cbasis, pbasis, conf_range, charge, mass);
#endif

  struct dg_gyrokinetic *gyrokinetic = gkyl_malloc(sizeof(struct dg_gyrokinetic));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  gyrokinetic->cdim = cdim;
  gyrokinetic->pdim = pdim;

  gyrokinetic->charge = charge;
  gyrokinetic->mass = mass;

  gyrokinetic->eqn.num_equations = 1;
  gyrokinetic->eqn.surf_term = surf;
  gyrokinetic->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_gyrokinetic_vol_kern_list *vol_kernels;
  const gkyl_dg_gyrokinetic_surf_kern_list *surf_x_kernels, *surf_y_kernels, *surf_z_kernels;
  const gkyl_dg_gyrokinetic_surf_kern_list *surf_vpar_kernels;
  const gkyl_dg_gyrokinetic_boundary_surf_kern_list *boundary_surf_x_kernels,*boundary_surf_y_kernels,*boundary_surf_z_kernels;
  const gkyl_dg_gyrokinetic_boundary_surf_kern_list *boundary_surf_vpar_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_x_kernels = ser_surf_x_kernels;
      surf_y_kernels = ser_surf_y_kernels;
      surf_z_kernels = ser_surf_z_kernels;
      surf_vpar_kernels = ser_surf_vpar_kernels;
      boundary_surf_x_kernels = ser_boundary_surf_x_kernels;
      boundary_surf_y_kernels = ser_boundary_surf_y_kernels;
      boundary_surf_z_kernels = ser_boundary_surf_z_kernels;
      boundary_surf_vpar_kernels = ser_boundary_surf_vpar_kernels;

      break;

    default:
      assert(false);
      break;
  }

  gyrokinetic->eqn.vol_term = CK(vol_kernels,cdim,vdim,poly_order);

  gyrokinetic->surf[0] = CK(surf_x_kernels,cdim,vdim,poly_order);
  if (cdim>1)
    gyrokinetic->surf[1] = CK(surf_y_kernels,cdim,vdim,poly_order);
  if (cdim>2)
    gyrokinetic->surf[2] = CK(surf_z_kernels,cdim,vdim,poly_order);
  gyrokinetic->surf[cdim] = CK(surf_vpar_kernels,cdim,vdim,poly_order);

  gyrokinetic->boundary_surf[0] = CK(boundary_surf_x_kernels,cdim,vdim,poly_order);
  if (cdim>1)
    gyrokinetic->boundary_surf[1] = CK(boundary_surf_y_kernels,cdim,vdim,poly_order);
  if (cdim>2)
    gyrokinetic->boundary_surf[2] = CK(boundary_surf_z_kernels,cdim,vdim,poly_order);
  gyrokinetic->boundary_surf[cdim] = CK(boundary_surf_vpar_kernels,cdim,vdim,poly_order);

  // Ensure non-NULL pointers.
  for (int i=0; i<cdim; ++i) assert(gyrokinetic->surf[i]);
  assert(gyrokinetic->surf[cdim]);
  for (int i=0; i<cdim; ++i) assert(gyrokinetic->boundary_surf[i]);
  assert(gyrokinetic->boundary_surf[cdim]);

  gyrokinetic->auxfields.bmag         = 0;
  gyrokinetic->auxfields.jacobtot_inv = 0;
  gyrokinetic->auxfields.cmag         = 0;
  gyrokinetic->auxfields.b_i          = 0;
  gyrokinetic->auxfields.phi     = 0;
  gyrokinetic->auxfields.apar    = 0;
  gyrokinetic->auxfields.apardot = 0;
  gyrokinetic->conf_range = *conf_range;

  gyrokinetic->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(gyrokinetic->eqn.flags);

  gyrokinetic->eqn.ref_count = gkyl_ref_count_init(gkyl_gyrokinetic_free);
  gyrokinetic->eqn.on_dev = &gyrokinetic->eqn; // CPU eqn obj points to itself

  return &gyrokinetic->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const double charge, const double mass) 
{
  assert(false);
  return 0;
}

#endif
