#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_util.h>
#include <gkyl_prim_lbo_gyrokinetic.h>
#include <gkyl_prim_lbo_gyrokinetic_priv.h>

void
prim_lbo_gyrokinetic_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_prim_lbo_type *prim = container_of(ref, struct gkyl_prim_lbo_type, ref_count);
  if (GKYL_IS_CU_ALLOC(prim->flag))
    gkyl_cu_free(prim->on_dev);
  gkyl_free(prim);
}

struct gkyl_prim_lbo_type*
gkyl_prim_lbo_gyrokinetic_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct prim_lbo_type_gyrokinetic *prim_gyrokinetic = gkyl_malloc(sizeof(struct prim_lbo_type_gyrokinetic));
  int cdim = prim_gyrokinetic->prim.cdim = cbasis->ndim;
  int pdim = prim_gyrokinetic->prim.pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int poly_order = prim_gyrokinetic->prim.poly_order = cbasis->poly_order;
  prim_gyrokinetic->prim.num_config = cbasis->num_basis;
  prim_gyrokinetic->prim.num_phase = pbasis->num_basis;
  prim_gyrokinetic->prim.self_prim = self_prim;

  // choose kernel tables based on basis-function type
  const gkyl_prim_lbo_gyrokinetic_kern_list *self_prim_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      self_prim_kernels = ser_self_prim_kernels;
      break;

    default:
      assert(false);
      break;    
  }
  assert(cv_index[cdim].vdim[vdim] != -1);
  assert(NULL != self_prim_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
  prim_gyrokinetic->self_prim = self_prim_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];

  prim_gyrokinetic->prim.flag = 0;
  GKYL_CLEAR_CU_ALLOC(prim_gyrokinetic->prim.flag);
  prim_gyrokinetic->prim.ref_count = gkyl_ref_count_init(prim_lbo_gyrokinetic_free);

  prim_gyrokinetic->prim.on_dev = &prim_gyrokinetic->prim;
    
  return &prim_gyrokinetic->prim;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_prim_lbo_type*
gkyl_prim_lbo_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis)
{
  assert(false);
  return 0;
}

#endif
