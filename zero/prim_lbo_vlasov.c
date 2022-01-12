#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_prim_lbo_vlasov.h>
#include <gkyl_prim_lbo_vlasov_priv.h>

static void
prim_lbo_vlasov_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_prim_lbo *prim = container_of(ref, struct gkyl_prim_lbo, ref_count);
  gkyl_free(prim);
}

struct gkyl_prim_lbo*
gkyl_prim_lbo_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct prim_lbo_vlasov *prim_vlasov = gkyl_malloc(sizeof(struct prim_lbo_vlasov));
  int cdim = prim_vlasov->prim.cdim = cbasis->ndim;
  int pdim = prim_vlasov->prim.pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int poly_order = prim_vlasov->prim.poly_order = cbasis->poly_order;
  prim_vlasov->prim.num_config = cbasis->num_basis;
  prim_vlasov->prim.num_phase = pbasis->num_basis;
  prim_vlasov->prim.self_prim = self_prim;

  // choose kernel tables based on basis-function type
  const gkyl_prim_lbo_vlasov_kern_list *self_prim_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      self_prim_kernels = ser_self_prim_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      break;

    default:
      assert(false);
      break;    
  }
  assert(cv_index[cdim].vdim[vdim] != -1);
  assert(NULL != self_prim_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
  prim_vlasov->self_prim = self_prim_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];

  // set reference counter
  prim_vlasov->prim.ref_count = gkyl_ref_count_init(prim_lbo_vlasov_free);
    
  return &prim_vlasov->prim;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_prim_lbo*
gkyl_prim_lbo_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis)
{
  assert(false);
  return 0;
}

#endif
