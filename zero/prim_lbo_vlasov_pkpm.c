#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_util.h>
#include <gkyl_prim_lbo_vlasov_pkpm.h>
#include <gkyl_prim_lbo_vlasov_pkpm_priv.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

void
prim_lbo_vlasov_pkpm_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_prim_lbo_type *prim = container_of(ref, struct gkyl_prim_lbo_type, ref_count);
  if (GKYL_IS_CU_ALLOC(prim->flag))
    gkyl_cu_free(prim->on_dev);
  gkyl_free(prim);
}

void
gkyl_prim_lbo_vlasov_pkpm_set_auxfields(const struct gkyl_prim_lbo_type *prim,
  struct gkyl_prim_lbo_vlasov_pkpm_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.pvar)) {
    gkyl_prim_lbo_vlasov_pkpm_set_auxfields_cu(prim->on_dev, auxin);
    return;
  }
#endif

  struct prim_lbo_type_vlasov_pkpm *prim_vlasov_pkpm = container_of(prim, struct prim_lbo_type_vlasov_pkpm, prim);
  prim_vlasov_pkpm->auxfields.pvar = auxin.pvar;
}


struct gkyl_prim_lbo_type*
gkyl_prim_lbo_vlasov_pkpm_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_prim_lbo_vlasov_pkpm_cu_dev_new(cbasis, pbasis, conf_range);
  } 
#endif     
  struct prim_lbo_type_vlasov_pkpm *prim_vlasov_pkpm = gkyl_malloc(sizeof(struct prim_lbo_type_vlasov_pkpm));
  int cdim = prim_vlasov_pkpm->prim.cdim = cbasis->ndim;
  int pdim = prim_vlasov_pkpm->prim.pdim = pbasis->ndim;

  int poly_order = prim_vlasov_pkpm->prim.poly_order = cbasis->poly_order;
  prim_vlasov_pkpm->prim.num_config = cbasis->num_basis;
  prim_vlasov_pkpm->prim.num_phase = pbasis->num_basis;
  prim_vlasov_pkpm->prim.self_prim = self_prim;

  // kinetic equation is in local rest frame, so no flow
  prim_vlasov_pkpm->prim.udim = 0;

  // choose kernel tables based on basis-function type
  const gkyl_prim_lbo_vlasov_pkpm_self_kern_list *self_prim_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      self_prim_kernels = ser_self_prim_kernels;

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      self_prim_kernels = ten_self_prim_kernels;
      
      break;

    default:
      assert(false);
      break;    
  }
    
  prim_vlasov_pkpm->self_prim = CK(self_prim_kernels, cdim, poly_order);

  prim_vlasov_pkpm->conf_range = *conf_range;
  
  prim_vlasov_pkpm->prim.flag = 0;
  GKYL_CLEAR_CU_ALLOC(prim_vlasov_pkpm->prim.flag);
  prim_vlasov_pkpm->prim.ref_count = gkyl_ref_count_init(prim_lbo_vlasov_pkpm_free);

  prim_vlasov_pkpm->prim.on_dev = &prim_vlasov_pkpm->prim;
    
  return &prim_vlasov_pkpm->prim;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_prim_lbo_type*
gkyl_prim_lbo_vlasov_pkpm_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range)
{
  assert(false);
  return 0;
}

#endif
