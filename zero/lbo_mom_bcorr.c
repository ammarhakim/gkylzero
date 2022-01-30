#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_lbo_mom_bcorr.h>
#include <gkyl_lbo_mom_bcorr_priv.h>
#include <gkyl_util.h>

void
mom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *momt = container_of(ref, struct gkyl_mom_type, ref_count);
  if (GKYL_IS_CU_ALLOC(momt->flag))
    gkyl_cu_free(momt->on_dev);
  gkyl_free(momt);
}


struct gkyl_mom_type*
gkyl_vlasov_lbo_mom_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const char *mom, const double* vBoundary)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct lbo_mom_type *mom_bcorr = gkyl_malloc(sizeof(struct lbo_mom_type));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_bcorr->momt.cdim = cdim;
  mom_bcorr->momt.pdim = pdim;
  mom_bcorr->momt.poly_order = poly_order;
  mom_bcorr->momt.num_config = cbasis->num_basis;
  mom_bcorr->momt.num_phase = pbasis->num_basis;
  mom_bcorr->momt.kernel = kernel;
  for (int d=0; d<vdim; ++d) {
    mom_bcorr->vBoundary[d] = vBoundary[d];
    mom_bcorr->vBoundary[d + vdim] = vBoundary[d + vdim];
  }

  // choose kernel tables based on basis-function type
  const gkyl_lbo_mom_kern_list *boundary_integral_f_kernels, *boundary_integral_vf_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      boundary_integral_f_kernels = ser_boundary_integral_f_kernels;
      boundary_integral_vf_kernels = ser_boundary_integral_vf_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      break;

    default:
      assert(false);
      break;    
  }
  if (strcmp(mom, "f") == 0) { // momentum correction
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != boundary_integral_f_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_bcorr->kernel = boundary_integral_f_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_bcorr->momt.num_mom = vdim;
  }
  else if (strcmp(mom, "vf") == 0) { // energy correction
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != boundary_integral_vf_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_bcorr->kernel = boundary_integral_vf_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_bcorr->momt.num_mom = vdim;
  }

  mom_bcorr->momt.flag = 0;
  GKYL_CLEAR_CU_ALLOC(mom_bcorr->momt.flag);
  mom_bcorr->momt.ref_count = gkyl_ref_count_init(mom_free);

  mom_bcorr->momt.on_dev = &mom_bcorr->momt;
    
  return &mom_bcorr->momt;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_mom_type*
gkyl_vlasov_lbo_mom_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const char *mom, const double* vBoundary)
{
  assert(false);
}

#endif
