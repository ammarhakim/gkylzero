#include "gkyl_dg_eqn.h"
#include <assert.h>
#include <stdio.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array.h>
#include <gkyl_dg_isoeuler.h>
#include <gkyl_dg_isoeuler_priv.h>
#include <gkyl_util.h>

void
gkyl_isoeuler_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_isoeuler *isoeuler = container_of(base->on_dev, struct dg_isoeuler, eqn);
    gkyl_cu_free(isoeuler);
  }

  struct dg_isoeuler *isoeuler = container_of(base, struct dg_isoeuler, eqn);
  gkyl_free(isoeuler);
}

void
gkyl_isoeuler_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_isoeuler_auxfields auxin)
{
  struct dg_isoeuler *isoeuler = container_of(eqn, struct dg_isoeuler, eqn);
  isoeuler->auxfields.uvar = auxin.uvar;
}

struct gkyl_dg_eqn*
gkyl_dg_isoeuler_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, //TODO: remove pbasis
  const struct gkyl_range* conf_range, const double vth, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    print("TODO: GKYL_HAVE_CUDA in dg_isoeuler.c....\n");
    return 0;
  }
#endif
  struct dg_isoeuler *isoeuler = gkyl_malloc(sizeof(struct dg_isoeuler));

  int cdim = cbasis->ndim-1, pdim = pbasis->ndim-1, vdim = pdim-cdim; //have cdim+1 number of dimensions in statevec
  int poly_order = cbasis->poly_order;

  isoeuler->cdim = cdim;
  isoeuler->pdim = pdim;

  isoeuler->vth = vth;

  isoeuler->eqn.num_equations = cdim+1;
  isoeuler->eqn.vol_term = vol;
  isoeuler->eqn.surf_term = surf;

  const gkyl_dg_isoeuler_vol_kern_list *vol_kernels;
  const gkyl_dg_isoeuler_surf_kern_list *surf_x_kernels, *surf_y_kernels, *surf_z_kernels;

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

  printf("Choosing kernel, TODO: double check that correct one is chosen...\n");
  printf("Choosing kernel with cdim %i and polyorder %i \n",cdim, poly_order);

  isoeuler->vol = CK(vol_kernels,cdim,poly_order); //TODO: clean up passing cdim+1 and then -2 in

  isoeuler->surf[0] = CK(surf_x_kernels,cdim,poly_order);
  if (cdim>1)
    isoeuler->surf[1] = CK(surf_y_kernels,cdim,poly_order);
  if (cdim>2)
    isoeuler->surf[2] = CK(surf_z_kernels,cdim,poly_order);

  // Ensure non-NULL pointers.
  assert(isoeuler->vol);
  for (int i=0; i<cdim; ++i) assert(isoeuler->surf[i]);
  assert(isoeuler->surf[cdim]); //TODO: not sure why this is needed, check if it is

  isoeuler->conf_range = *conf_range;

  isoeuler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(isoeuler->eqn.flags);

  isoeuler->eqn.ref_count = gkyl_ref_count_init(gkyl_isoeuler_free);
  isoeuler->eqn.on_dev = &isoeuler->eqn; // CPU eqn obj points to itself

  return &isoeuler->eqn;
}
