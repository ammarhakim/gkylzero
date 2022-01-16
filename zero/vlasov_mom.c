#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_mom.h>
#include <gkyl_vlasov_mom_priv.h>

void
gkyl_mom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *momt = container_of(ref, struct gkyl_mom_type, ref_count);
  if (GKYL_IS_CU_ALLOC(momt->flag))
    gkyl_cu_free(momt->on_dev);
  gkyl_free(momt);
}

struct gkyl_mom_type*
gkyl_vlasov_mom_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct vlasov_mom_type *vlasov_mom = gkyl_malloc(sizeof(struct vlasov_mom_type));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  vlasov_mom->momt.cdim = cdim;
  vlasov_mom->momt.pdim = pdim;
  vlasov_mom->momt.poly_order = poly_order;
  vlasov_mom->momt.num_config = cbasis->num_basis;
  vlasov_mom->momt.num_phase = pbasis->num_basis;
  vlasov_mom->momt.kernel = kernel;

  // choose kernel tables based on basis-function type
  const gkyl_mom_kern_list *m0_kernels, *m1i_kernels,
    *m2_kernels, *m2ij_kernels, *m3i_kernels, *m3ijk_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m1i_kernels = ser_m1i_kernels;
      m2_kernels = ser_m2_kernels;
      m2ij_kernels = ser_m2ij_kernels;
      m3i_kernels = ser_m3i_kernels;
      m3ijk_kernels = ser_m3ijk_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      m0_kernels = ten_m0_kernels;
      m1i_kernels = ten_m1i_kernels;
      m2_kernels = ten_m2_kernels;
      m2ij_kernels = ten_m2ij_kernels;
      m3i_kernels = ten_m3i_kernels;
      m3ijk_kernels = ten_m3ijk_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  if (strcmp(mom, "M0") == 0) { // density
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    vlasov_mom->kernel = m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    vlasov_mom->momt.num_mom = 1;
  }
  else if (strcmp(mom, "M1i") == 0) { // momentum
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m1i_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    vlasov_mom->kernel = m1i_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    vlasov_mom->momt.num_mom = vdim;
  }
  else if (strcmp(mom, "M2") == 0) { // energy
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    vlasov_mom->kernel = m2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    vlasov_mom->momt.num_mom = 1;
  }
  else if (strcmp(mom, "M2ij") == 0) { // pressure tensor in lab-frame
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2ij_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    vlasov_mom->kernel = m2ij_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    vlasov_mom->momt.num_mom = vdim*(vdim+1)/2;
  }
  else if (strcmp(mom, "M3i") == 0) { // heat-flux vector in lab-frame
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3i_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    vlasov_mom->kernel = m3i_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    vlasov_mom->momt.num_mom = vdim;
  }
  else if (strcmp(mom, "M3ijk") == 0) { // heat-flux tensor in lab-frame
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3ijk_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    vlasov_mom->kernel = m3ijk_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];

    int m3ijk_count[] = { 1, 4, 10 };
    vlasov_mom->momt.num_mom = m3ijk_count[vdim-1];
  }
  else {
    // string not recognized
    gkyl_exit("gkyl_vlasov_mom_type: Unrecognized moment requested!");
  }

  vlasov_mom->momt.flag = 0;
  GKYL_CLEAR_CU_ALLOC(vlasov_mom->momt.flag);
  vlasov_mom->momt.ref_count = gkyl_ref_count_init(gkyl_mom_free);
  
  vlasov_mom->momt.on_dev = &vlasov_mom->momt; // on host, self-reference
    
  return &vlasov_mom->momt;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_mom_type*
gkyl_vlasov_mom_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom)
{
  assert(false);
  return 0;
}

#endif
