/* -*- c++ -*- */

#include <assert.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_util.h>
#include <gkyl_prim_lbo_vlasov_pkpm.h>
#include <gkyl_prim_lbo_vlasov_pkpm_priv.h>
}

#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_prim_lbo_vlasov_pkpm_set_auxfields_cu_kernel(const struct gkyl_prim_lbo_type *prim, const struct gkyl_array *pvar)
{
  struct prim_lbo_type_vlasov_pkpm *prim_vlasov_pkpm = container_of(prim, struct prim_lbo_type_vlasov_pkpm, prim);
  prim_vlasov_pkpm->auxfields.pvar = pvar;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_prim_lbo_vlasov_pkpm_set_auxfields_cu(const struct gkyl_prim_lbo_type *prim,
  struct gkyl_prim_lbo_vlasov_pkpm_auxfields auxin)
{
  gkyl_prim_lbo_vlasov_pkpm_set_auxfields_cu_kernel<<<1,1>>>(prim, auxin.pvar->on_dev);
}

__global__ static void
gkyl_prim_lbo_vlasov_pkpm_set_cu_dev_ptrs(struct prim_lbo_type_vlasov_pkpm *prim_vlasov_pkpm, 
  enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  prim_vlasov_pkpm->prim.self_prim = self_prim;
  
  // choose kernel tables based on basis-function type
  const gkyl_prim_lbo_vlasov_pkpm_self_kern_list *self_prim_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      self_prim_kernels = ser_self_prim_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  prim_vlasov_pkpm->self_prim = CK(self_prim_kernels, cdim, poly_order);
}

struct gkyl_prim_lbo_type*
gkyl_prim_lbo_vlasov_pkpm_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct prim_lbo_type_vlasov_pkpm *prim_vlasov_pkpm =
    (struct prim_lbo_type_vlasov_pkpm*) gkyl_malloc(sizeof(struct prim_lbo_type_vlasov_pkpm));
  
  int cdim = prim_vlasov_pkpm->prim.cdim = cbasis->ndim;
  int pdim = prim_vlasov_pkpm->prim.pdim = pbasis->ndim;

  int poly_order = prim_vlasov_pkpm->prim.poly_order = cbasis->poly_order;
  prim_vlasov_pkpm->prim.udim = 0;

  prim_vlasov_pkpm->conf_range = *conf_range;

  prim_vlasov_pkpm->prim.flag = 0;
  GKYL_SET_CU_ALLOC(prim_vlasov_pkpm->prim.flag);
  prim_vlasov_pkpm->prim.ref_count = gkyl_ref_count_init(prim_lbo_vlasov_pkpm_free);
  
  // copy the host struct to device struct
  struct prim_lbo_type_vlasov_pkpm *prim_vlasov_pkpm_cu = (struct prim_lbo_type_vlasov_pkpm*)
    gkyl_cu_malloc(sizeof(struct prim_lbo_type_vlasov_pkpm));
  gkyl_cu_memcpy(prim_vlasov_pkpm_cu, prim_vlasov_pkpm, sizeof(struct prim_lbo_type_vlasov_pkpm), GKYL_CU_MEMCPY_H2D);
  
  gkyl_prim_lbo_vlasov_pkpm_set_cu_dev_ptrs<<<1,1>>>(prim_vlasov_pkpm_cu, cbasis->b_type, cdim, poly_order);

  prim_vlasov_pkpm->prim.on_dev = &prim_vlasov_pkpm_cu->prim;
    
  return &prim_vlasov_pkpm->prim;
}
