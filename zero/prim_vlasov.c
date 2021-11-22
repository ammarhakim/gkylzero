#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_prim_vlasov.h>
#include <gkyl_prim_vlasov_priv.h>

struct gkyl_prim_vlasov*
gkyl_prim_vlasov_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct gkyl_prim_vlasov *prim_vlasov = gkyl_malloc(sizeof(struct gkyl_prim_vlasov));
  int cdim = prim_vlasov->cdim = cbasis->ndim;
  int pdim = prim_vlasov->pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int poly_order = prim_vlasov->poly_order = cbasis->poly_order;
  prim_vlasov->num_config = cbasis->num_basis;
  prim_vlasov->num_phase = pbasis->num_basis;

  // choose kernel tables based on basis-function type
  const gkyl_prim_vlasov_kern_list *self_prim_kernels;

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
    
  return prim_vlasov;
}

void
gkyl_prim_vlasov_release(struct gkyl_prim_vlasov* prim)
{
  gkyl_free(prim);
}
