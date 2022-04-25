#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_util.h>
#include <gkyl_prim_lbo_vlasov_with_fluid.h>
#include <gkyl_prim_lbo_vlasov_with_fluid_priv.h>

void
prim_lbo_vlasov_with_fluid_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_prim_lbo_type *prim = container_of(ref, struct gkyl_prim_lbo_type, ref_count);
  if (GKYL_IS_CU_ALLOC(prim->flag))
    gkyl_cu_free(prim->on_dev);
  gkyl_free(prim);
}

void
gkyl_prim_lbo_vlasov_with_fluid_set_auxfields(const struct gkyl_prim_lbo_type *prim,
  struct gkyl_prim_lbo_vlasov_with_fluid_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.fluid)) {
    gkyl_prim_lbo_vlasov_with_fluid_set_auxfields_cu(prim->on_dev, auxin);
    return;
  }
#endif

  struct prim_lbo_type_vlasov_with_fluid *prim_vlasov_with_fluid = container_of(prim, struct prim_lbo_type_vlasov_with_fluid, prim);
  prim_vlasov_with_fluid->auxfields.fluid = auxin.fluid;
}


struct gkyl_prim_lbo_type*
gkyl_prim_lbo_vlasov_with_fluid_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct prim_lbo_type_vlasov_with_fluid *prim_vlasov_with_fluid = gkyl_malloc(sizeof(struct prim_lbo_type_vlasov_with_fluid));
  int cdim = prim_vlasov_with_fluid->prim.cdim = cbasis->ndim;
  int pdim = prim_vlasov_with_fluid->prim.pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int poly_order = prim_vlasov_with_fluid->prim.poly_order = cbasis->poly_order;
  prim_vlasov_with_fluid->prim.num_config = cbasis->num_basis;
  prim_vlasov_with_fluid->prim.num_phase = pbasis->num_basis;
  prim_vlasov_with_fluid->prim.self_prim = self_prim;
  prim_vlasov_with_fluid->prim.cross_prim = cross_prim;

  // choose kernel tables based on basis-function type
  const gkyl_prim_lbo_vlasov_with_fluid_self_kern_list *self_prim_kernels;
  const gkyl_prim_lbo_vlasov_with_fluid_cross_kern_list *cross_prim_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      self_prim_kernels = ser_self_prim_kernels;
      cross_prim_kernels = ser_cross_prim_kernels;
      break;

    default:
      assert(false);
      break;    
  }
  assert(cv_index[cdim].vdim[vdim] != -1);
  assert(NULL != self_prim_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
  assert(NULL != cross_prim_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
  prim_vlasov_with_fluid->self_prim = self_prim_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  prim_vlasov_with_fluid->cross_prim = cross_prim_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];

  prim_vlasov_with_fluid->conf_range = *conf_range;
  
  prim_vlasov_with_fluid->prim.flag = 0;
  GKYL_CLEAR_CU_ALLOC(prim_vlasov_with_fluid->prim.flag);
  prim_vlasov_with_fluid->prim.ref_count = gkyl_ref_count_init(prim_lbo_vlasov_with_fluid_free);

  prim_vlasov_with_fluid->prim.on_dev = &prim_vlasov_with_fluid->prim;
    
  return &prim_vlasov_with_fluid->prim;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_prim_lbo_type*
gkyl_prim_lbo_vlasov_with_fluid_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range)
{
  assert(false);
  return 0;
}

#endif
