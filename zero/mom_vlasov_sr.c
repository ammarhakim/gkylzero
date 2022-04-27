#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_vlasov_sr.h>
#include <gkyl_mom_vlasov_sr_priv.h>
#include <gkyl_util.h>

void
gkyl_mom_vm_sr_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *momt = container_of(ref, struct gkyl_mom_type, ref_count);
  if (GKYL_IS_CU_ALLOC(momt->flags))
    gkyl_cu_free(momt->on_dev);
  gkyl_free(momt);
}

void
gkyl_mom_vlasov_sr_set_auxfields(const struct gkyl_mom_type *momt, struct gkyl_mom_vlasov_sr_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_mom_type_is_cu_dev(momt)) {
    gkyl_mom_vlasov_sr_set_auxfields_cu(momt->on_dev, auxin);
    return;
  }
#endif

  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);
  mom_vm_sr->auxfields.p_over_gamma = auxin.p_over_gamma;
}

struct gkyl_mom_type*
gkyl_mom_vlasov_sr_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* vel_range,
  const char *mom)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct mom_type_vlasov_sr *mom_vm_sr = gkyl_malloc(sizeof(struct mom_type_vlasov_sr));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_vm_sr->momt.cdim = cdim;
  mom_vm_sr->momt.pdim = pdim;
  mom_vm_sr->momt.poly_order = poly_order;
  mom_vm_sr->momt.num_config = cbasis->num_basis;
  mom_vm_sr->momt.num_phase = pbasis->num_basis;
  mom_vm_sr->momt.kernel = kernel;

  // choose kernel tables based on basis-function type
  const gkyl_mom_kern_list *m0_kernels, *m1i_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m1i_kernels = ser_m1i_kernels;
      break;

    // case GKYL_BASIS_MODAL_TENSOR:
    //   m0_kernels = ten_m0_kernels;
    //   m1i_kernels = ten_m1i_kernels;
    //   break;

    default:
      assert(false);
      break;    
  }

  if (strcmp(mom, "M0") == 0) { // density
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_vm_sr->kernel = m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_vm_sr->momt.num_mom = 1;
  }
  else if (strcmp(mom, "M1i") == 0) { // momentum
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m1i_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_vm_sr->kernel = m1i_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_vm_sr->momt.num_mom = vdim;
  }
  else {
    // string not recognized
    gkyl_exit("gkyl_mom_type_vlasov_sr: Unrecognized moment requested!");
  }

  mom_vm_sr->vel_range = *vel_range;

  mom_vm_sr->auxfields.p_over_gamma = 0;

  mom_vm_sr->momt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mom_vm_sr->momt.flags);
  mom_vm_sr->momt.ref_count = gkyl_ref_count_init(gkyl_mom_vm_sr_free);
  
  mom_vm_sr->momt.on_dev = &mom_vm_sr->momt; // on host, self-reference
    
  return &mom_vm_sr->momt;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_mom_type*
gkyl_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* vel_range,
  const char *mom)
{
  assert(false);
  return 0;
}

#endif
