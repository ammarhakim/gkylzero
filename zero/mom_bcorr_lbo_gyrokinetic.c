#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_bcorr_lbo_gyrokinetic.h>
#include <gkyl_mom_bcorr_lbo_gyrokinetic_priv.h>
#include <gkyl_util.h>

void
gk_mom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *momt = container_of(ref, struct gkyl_mom_type, ref_count);

  struct mom_type_bcorr_lbo_gyrokinetic *mom_bcorr = container_of(momt, struct mom_type_bcorr_lbo_gyrokinetic, momt);
  gkyl_velocity_map_release(mom_bcorr->vel_map);

  if (GKYL_IS_CU_ALLOC(momt->flags))
    gkyl_cu_free(momt->on_dev);
  gkyl_free(momt);
}


struct gkyl_mom_type*
gkyl_mom_bcorr_lbo_gyrokinetic_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  double mass, const struct gkyl_velocity_map *vel_map, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    return gkyl_mom_bcorr_lbo_gyrokinetic_cu_dev_new(cbasis, pbasis, mass, vel_map);
  } 
#endif    
  struct mom_type_bcorr_lbo_gyrokinetic *mom_bcorr = gkyl_malloc(sizeof(struct mom_type_bcorr_lbo_gyrokinetic));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_bcorr->momt.cdim = cdim;
  mom_bcorr->momt.pdim = pdim;
  mom_bcorr->momt.poly_order = poly_order;
  mom_bcorr->momt.num_config = cbasis->num_basis;
  mom_bcorr->momt.num_phase = pbasis->num_basis;
  mom_bcorr->momt.kernel = kernel;

  // choose kernel tables based on basis-function type
  const gkyl_mom_bcorr_lbo_gyrokinetic_kern_list *mom_bcorr_lbo_gyrokinetic_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      mom_bcorr_lbo_gyrokinetic_kernels = ser_mom_bcorr_lbo_gyrokinetic_kernels;
      break;

    /* case GKYL_BASIS_MODAL_TENSOR: */
    /*   mom_bcorr_lbo_gyrokinetic_kernels = ten_mom_bcorr_lbo_gyrokinetic_kernels; */
    /*   break; */

    default:
      assert(false);
      break;    
  }
  assert(cv_index[cdim].vdim[vdim] != -1);
  assert(NULL != mom_bcorr_lbo_gyrokinetic_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);

  mom_bcorr->kernel = mom_bcorr_lbo_gyrokinetic_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  mom_bcorr->momt.num_mom = 2;

  mom_bcorr->_m = mass;
  mom_bcorr->vel_map = gkyl_velocity_map_acquire(vel_map);
  
  mom_bcorr->momt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mom_bcorr->momt.flags);
  mom_bcorr->momt.ref_count = gkyl_ref_count_init(gk_mom_free);

  mom_bcorr->momt.on_dev = &mom_bcorr->momt;
    
  return &mom_bcorr->momt;
}
