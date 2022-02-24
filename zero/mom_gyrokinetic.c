#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_mom_gyrokinetic_priv.h>
#include <gkyl_util.h>

void
gkyl_gk_mom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *momt = container_of(ref, struct gkyl_mom_type, ref_count);
  if (GKYL_IS_CU_ALLOC(momt->flag))
    gkyl_cu_free(momt->on_dev);
  gkyl_free(momt);
}

void
gkyl_vlasov_set_bmag(const struct gkyl_mom_type *momt, const struct gkyl_array *bmag)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(bmag)) {
    //gkyl_vlasov_set_bmag_cu(momt->on_dev, bmag);
    return;
  }
#endif

  struct mom_type_gyrokinetic *mom_gyrokinetic = container_of(momt, struct mom_type_gyrokinetic, momt);
  mom_gyrokinetic->bmag = bmag;
}

struct gkyl_mom_type*
gkyl_mom_gyrokinetic_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, double mass, const char *mom)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct mom_type_gyrokinetic *mom_gyrokinetic = gkyl_malloc(sizeof(struct mom_type_gyrokinetic));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_gyrokinetic->momt.cdim = cdim;
  mom_gyrokinetic->momt.pdim = pdim;
  mom_gyrokinetic->momt.poly_order = poly_order;
  mom_gyrokinetic->momt.num_config = cbasis->num_basis;
  mom_gyrokinetic->momt.num_phase = pbasis->num_basis;
  mom_gyrokinetic->momt.kernel = kernel;

  // choose kernel tables based on basis-function type
  const gkyl_gyrokinetic_mom_kern_list *m0_kernels, *m1_kernels, *m2_kernels, 
    *m2_par_kernels, *m2_perp_kernels, *m3_par_kernels, *m3_perp_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m1_kernels = ser_m1_kernels;
      m2_kernels = ser_m2_kernels;
      m2_par_kernels = ser_m2_par_kernels;
      m2_perp_kernels = ser_m2_perp_kernels;
      m3_par_kernels = ser_m3_par_kernels;
      m3_perp_kernels = ser_m3_perp_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      m0_kernels = ten_m0_kernels;
      m1_kernels = ten_m1_kernels;
      m2_kernels = ten_m2_kernels;
      m2_par_kernels = ten_m2_par_kernels;
      m2_perp_kernels = ten_m2_perp_kernels;
      m3_par_kernels = ten_m3_par_kernels;
      m3_perp_kernels = ten_m3_perp_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  if (strcmp(mom, "GkM0") == 0) { // density
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gyrokinetic->kernel = m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gyrokinetic->momt.num_mom = 1;
  }
  else if (strcmp(mom, "GkM1") == 0) { // parallel momentum
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m1_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gyrokinetic->kernel = m1_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gyrokinetic->momt.num_mom = 1;
  }
  else if (strcmp(mom, "GkM2") == 0) { // energy
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gyrokinetic->kernel = m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gyrokinetic->momt.num_mom = 1;
  }
  else if (strcmp(mom, "GkM2par") == 0) { // parallel energy
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2_par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gyrokinetic->kernel = m2_par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gyrokinetic->momt.num_mom = 1;
  }
  else if (strcmp(mom, "GkM2perp") == 0) { // perpendicular energy
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2_perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gyrokinetic->kernel = m2_perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gyrokinetic->momt.num_mom = 1;
  }
  else if (strcmp(mom, "GkM3par") == 0) { // parallel heat flux
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3_par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gyrokinetic->kernel = m3_par_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gyrokinetic->momt.num_mom = 1;
  }
  else if (strcmp(mom, "GkM3perp") == 0) { // perpendicular heat flux
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3_perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_gyrokinetic->kernel = m3_perp_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_gyrokinetic->momt.num_mom = 1;
  }
  else {
    // string not recognized
    gkyl_exit("gkyl_mom_type_gyrokinetic: Unrecognized moment requested!");
  }

  mom_gyrokinetic->_m = mass;
  mom_gyrokinetic->bmag = 0;  
  mom_gyrokinetic->conf_range = *conf_range;
  
  mom_gyrokinetic->momt.flag = 0;
  GKYL_CLEAR_CU_ALLOC(mom_gyrokinetic->momt.flag);
  mom_gyrokinetic->momt.ref_count = gkyl_ref_count_init(gkyl_gk_mom_free);
  
  mom_gyrokinetic->momt.on_dev = &mom_gyrokinetic->momt; // on host, self-reference
    
  return &mom_gyrokinetic->momt;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_mom_type*
gkyl_mom_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, double mass, const char *mom)
{
  assert(false);
  return 0;
}

#endif
