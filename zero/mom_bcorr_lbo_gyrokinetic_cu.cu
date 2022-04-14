/* -*- c++ -*- */

#include <assert.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_util.h>
#include <gkyl_mom_bcorr_lbo_gyrokinetic.h>
#include <gkyl_mom_bcorr_lbo_gyrokinetic_priv.h>
}

__global__
static void
gkyl_mom_bcorr_lbo_gyrokinetic_set_cu_dev_ptrs(struct mom_type_bcorr_lbo_gyrokinetic* mom_bcorr, enum gkyl_basis_type b_type, int vdim, int poly_order, int tblidx)
{
  mom_bcorr->momt.kernel = kernel;

  // choose kernel tables based on basis-function type
  const gkyl_mom_bcorr_lbo_gyrokinetic_kern_list *mom_bcorr_lbo_gyrokinetic_kernels;

  switch (b_type) {
  case GKYL_BASIS_MODAL_SERENDIPITY:
    mom_bcorr_lbo_gyrokinetic_kernels = ser_mom_bcorr_lbo_gyrokinetic_kernels;
    break;

  // case GKYL_BASIS_MODAL_TENSOR:
  //   mom_bcorr_lbo_gyrokinetic_kernels = ten_mom_bcorr_lbo_gyrokinetic_kernels;
  //   break;

  default:
    assert(false);
    break;
  }
  mom_bcorr->kernel = mom_bcorr_lbo_gyrokinetic_kernels[tblidx].kernels[poly_order];
  mom_bcorr->momt.num_mom = 2;
}

struct gkyl_mom_type*
gkyl_mom_bcorr_lbo_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const double *vBoundary, double mass)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_bcorr_lbo_gyrokinetic *mom_bcorr = (struct mom_type_bcorr_lbo_gyrokinetic*) gkyl_malloc(sizeof(struct mom_type_bcorr_lbo_gyrokinetic));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_bcorr->momt.cdim = cdim;
  mom_bcorr->momt.pdim = pdim;
  mom_bcorr->momt.poly_order = poly_order;
  mom_bcorr->momt.num_config = cbasis->num_basis;
  mom_bcorr->momt.num_phase = pbasis->num_basis;
  for (int d=0; d<vdim; ++d) {
    mom_bcorr->vBoundary[d] = vBoundary[d];
    mom_bcorr->vBoundary[d + vdim] = vBoundary[d + vdim];
  }
  mom_bcorr->momt.num_mom = 2; // number of moments

  mom_bcorr->_m = mass;

  mom_bcorr->momt.flag = 0;
  GKYL_SET_CU_ALLOC(mom_bcorr->momt.flag);
  mom_bcorr->momt.ref_count = gkyl_ref_count_init(gk_mom_free);

  // copy struct to device
  struct mom_type_bcorr_lbo_gyrokinetic *mom_bcorr_cu = (struct mom_type_bcorr_lbo_gyrokinetic*)
    gkyl_cu_malloc(sizeof(struct mom_type_bcorr_lbo_gyrokinetic));
  gkyl_cu_memcpy(mom_bcorr_cu, mom_bcorr, sizeof(struct mom_type_bcorr_lbo_gyrokinetic), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);

  gkyl_mom_bcorr_lbo_gyrokinetic_set_cu_dev_ptrs<<<1,1>>>(mom_bcorr_cu, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  mom_bcorr->momt.on_dev = &mom_bcorr_cu->momt;

  return &mom_bcorr->momt;
}
