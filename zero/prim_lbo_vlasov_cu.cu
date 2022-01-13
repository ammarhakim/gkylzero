/* -*- c++ -*- */

#include <assert.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_prim_lbo_vlasov.h>
#include <gkyl_prim_lbo_vlasov_priv.h>
}

__global__ static void
prim_lbo_vlasov_set_cu_dev_ptrs(struct prim_lbo_vlasov *prim_vlasov, int cdim, int vdim, int poly_order, enum gkyl_basis_type b_type, int cv_index)
{
  // choose kernel tables based on basis-function type
  const gkyl_prim_lbo_vlasov_kern_list *self_prim_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      self_prim_kernels = ser_self_prim_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      break;

    default:
      assert(false);
      break;    
  }

  prim_vlasov->self_prim = self_prim_kernels[cv_index].kernels[poly_order];
}

struct gkyl_prim_lbo*
gkyl_prim_lbo_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct prim_lbo_vlasov *prim_vlasov = (struct prim_lbo_vlasov*) gkyl_malloc(sizeof(struct prim_lbo_vlasov));
  int cdim = prim_vlasov->prim.cdim = cbasis->ndim;
  int pdim = prim_vlasov->prim.pdim = pbasis->ndim;
  int vdim = pdim - cdim;
  int poly_order = prim_vlasov->prim.poly_order = cbasis->poly_order;
  
  // copy the host struct to device struct
  struct prim_lbo_vlasov *prim_vlasov_cu = (struct prim_lbo_vlasov*) gkyl_cu_malloc(sizeof(struct prim_lbo_vlasov));
  gkyl_cu_memcpy(prim_vlasov_cu, prim_vlasov, sizeof(struct prim_lbo_vlasov), GKYL_CU_MEMCPY_H2D);

  prim_lbo_vlasov_set_cu_dev_ptrs<<<1,1>>>(prim_vlasov_cu, cdim, vdim, poly_order, cbasis->b_type, cv_index[cdim].vdim[vdim]);

  gkyl_free(prim_vlasov);  
    
  return &prim_vlasov_cu->prim;
}
