/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_vlasov_pkpm.h>
#include <gkyl_mom_vlasov_pkpm_priv.h>
#include <gkyl_util.h>
}

#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

__global__
static void
set_cu_ptrs(struct mom_type_vlasov_pkpm *mom_vlasov_pkpm,
  enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  // choose kernel tables based on basis-function type
  const gkyl_mom_vlasov_pkpm_kern_list *mom_vlasov_pkpm_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      mom_vlasov_pkpm_kernels = ser_mom_vlasov_pkpm_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  mom_vlasov_pkpm->momt.kernel = CK(mom_vlasov_pkpm_kernels, cdim, poly_order);
}

struct gkyl_mom_type*
gkyl_mom_vlasov_pkpm_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  double mass)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_vlasov_pkpm *mom_vlasov_pkpm = (struct mom_type_vlasov_pkpm*)
    gkyl_malloc(sizeof(struct mom_type_vlasov_pkpm));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim;
  int poly_order = cbasis->poly_order;

  mom_vlasov_pkpm->momt.cdim = cdim;
  mom_vlasov_pkpm->momt.pdim = pdim;
  mom_vlasov_pkpm->momt.poly_order = poly_order;
  mom_vlasov_pkpm->momt.num_config = cbasis->num_basis;
  mom_vlasov_pkpm->momt.num_phase = pbasis->num_basis;

  mom_vlasov_pkpm->momt.num_mom = 3; // rho, p_parallel, q_parallel 

  mom_vlasov_pkpm->mass = mass;

  mom_vlasov_pkpm->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_vlasov_pkpm->momt.flags);
  mom_vlasov_pkpm->momt.ref_count = gkyl_ref_count_init(gkyl_gk_mom_free);
  
  // copy struct to device
  struct mom_type_vlasov_pkpm *mom_vlasov_pkpm_cu = (struct mom_type_vlasov_pkpm*)
    gkyl_cu_malloc(sizeof(struct mom_type_vlasov_pkpm));
  gkyl_cu_memcpy(mom_vlasov_pkpm_cu, mom_vlasov_pkpm, sizeof(struct mom_type_vlasov_pkpm), GKYL_CU_MEMCPY_H2D);

  set_cu_ptrs<<<1,1>>>(mom_vlasov_pkpm_cu, cbasis->b_type, cdim, poly_order);

  mom_vlasov_pkpm->momt.on_dev = &mom_vlasov_pkpm_cu->momt;
  
  return &mom_vlasov_pkpm->momt;
}
