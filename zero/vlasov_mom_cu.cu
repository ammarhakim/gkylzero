#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_mom.h>
#include <gkyl_vlasov_mom_priv.h>
}

// M0 kernel list
static struct { momf_t kernels[3]; } p_m0_kernels[] = {
  // 1x kernels
  { NULL, &vlasov_M0_1x1v_ser_p1, &vlasov_M0_1x1v_ser_p2 }, // 0
  { NULL, &vlasov_M0_1x2v_ser_p1, &vlasov_M0_1x2v_ser_p2 }, // 1
  { NULL, &vlasov_M0_1x3v_ser_p1, &vlasov_M0_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, &vlasov_M0_2x2v_ser_p1, &vlasov_M0_2x2v_ser_p2 }, // 3
  { NULL, &vlasov_M0_2x3v_ser_p1, &vlasov_M0_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_M0_3x3v_ser_p1, NULL                  }, // 5
};

// M1i kernel list
static struct { momf_t kernels[3]; } p_m1i_kernels[] = {
  // 1x kernels
  { NULL, &vlasov_M1i_1x1v_ser_p1, &vlasov_M1i_1x1v_ser_p2 }, // 0
  { NULL, &vlasov_M1i_1x2v_ser_p1, &vlasov_M1i_1x2v_ser_p2 }, // 1
  { NULL, &vlasov_M1i_1x3v_ser_p1, &vlasov_M1i_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, &vlasov_M1i_2x2v_ser_p1, &vlasov_M1i_2x2v_ser_p2 }, // 3
  { NULL, &vlasov_M1i_2x3v_ser_p1, &vlasov_M1i_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_M1i_3x3v_ser_p1, NULL                   }, // 5
};

// M2 kernel list
static struct { momf_t kernels[3]; } p_m2_kernels[] = {
  // 1x kernels
  { NULL, &vlasov_M2_1x1v_ser_p1, &vlasov_M2_1x1v_ser_p2 }, // 0
  { NULL, &vlasov_M2_1x2v_ser_p1, &vlasov_M2_1x2v_ser_p2 }, // 1
  { NULL, &vlasov_M2_1x3v_ser_p1, &vlasov_M2_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, &vlasov_M2_2x2v_ser_p1, &vlasov_M2_2x2v_ser_p2 }, // 3
  { NULL, &vlasov_M2_2x3v_ser_p1, &vlasov_M2_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_M2_3x3v_ser_p1, NULL                   }, // 5
};

// M2ij kernel list
static struct { momf_t kernels[3]; } p_m2ij_kernels[] = {
  // 1x kernels
  { NULL, &vlasov_M2ij_1x1v_ser_p1, &vlasov_M2ij_1x1v_ser_p2 }, // 0
  { NULL, &vlasov_M2ij_1x2v_ser_p1, &vlasov_M2ij_1x2v_ser_p2 }, // 1
  { NULL, &vlasov_M2ij_1x3v_ser_p1, &vlasov_M2ij_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, &vlasov_M2ij_2x2v_ser_p1, &vlasov_M2ij_2x2v_ser_p2 }, // 3
  { NULL, &vlasov_M2ij_2x3v_ser_p1, &vlasov_M2ij_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_M2ij_3x3v_ser_p1, NULL                   }, // 5
};

// M3i kernel list
static struct { momf_t kernels[3]; } p_m3i_kernels[] = {
  // 1x kernels
  { NULL, &vlasov_M3i_1x1v_ser_p1, &vlasov_M3i_1x1v_ser_p2 }, // 0
  { NULL, &vlasov_M3i_1x2v_ser_p1, &vlasov_M3i_1x2v_ser_p2 }, // 1
  { NULL, &vlasov_M3i_1x3v_ser_p1, &vlasov_M3i_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, &vlasov_M3i_2x2v_ser_p1, &vlasov_M3i_2x2v_ser_p2 }, // 3
  { NULL, &vlasov_M3i_2x3v_ser_p1, &vlasov_M3i_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_M3i_3x3v_ser_p1, NULL                   }, // 5
};

// M3ijk kernel list
static struct { momf_t kernels[3]; } p_m3ijk_kernels[] = {
  // 1x kernels
  { NULL, &vlasov_M3ijk_1x1v_ser_p1, &vlasov_M3ijk_1x1v_ser_p2 }, // 0
  { NULL, &vlasov_M3ijk_1x2v_ser_p1, &vlasov_M3ijk_1x2v_ser_p2 }, // 1
  { NULL, &vlasov_M3ijk_1x3v_ser_p1, &vlasov_M3ijk_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, &vlasov_M3ijk_2x2v_ser_p1, &vlasov_M3ijk_2x2v_ser_p2 }, // 3
  { NULL, &vlasov_M3ijk_2x3v_ser_p1, &vlasov_M3ijk_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, &vlasov_M3ijk_3x3v_ser_p1, NULL                   }, // 5
};

struct gkyl_mom_type*
gkyl_vlasov_mom_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom)
{
  assert(cbasis->polyOrder == pbasis->polyOrder);
  
  struct gkyl_mom_type *momt = (struct gkyl_mom_type*) gkyl_malloc(sizeof(struct gkyl_mom_type));
  int cdim = momt->cdim = cbasis->ndim;
  int pdim = momt->pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int polyOrder = momt->polyOrder = cbasis->polyOrder;
  momt->num_config = cbasis->numBasis;
  momt->num_phase = pbasis->numBasis;

  if (strcmp(mom, "M0") == 0) { // density
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    cudaMemcpyFromSymbol(&momt->kernel,
      p_m0_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder], sizeof(momf_t));
    momt->num_mom = 1;
  }
  else if (strcmp(mom, "M1i") == 0) { // momentum
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m1i_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    cudaMemcpyFromSymbol(&momt->kernel,
      p_m1i_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder], sizeof(momf_t));
    momt->num_mom = vdim;
  }
  else if (strcmp(mom, "M2") == 0) { // energy
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    cudaMemcpyFromSymbol(&momt->kernel,
      p_m2_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder], sizeof(momf_t));
    momt->num_mom = 1;
  }
  else if (strcmp(mom, "M2ij") == 0) { // pressure tensor in lab-frame
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2ij_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    cudaMemcpyFromSymbol(&momt->kernel,
      p_m2ij_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder], sizeof(momf_t));
    momt->num_mom = vdim*(vdim+1)/2;
  }
  else if (strcmp(mom, "M3i") == 0) { // heat-flux vector in lab-frame
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3i_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    cudaMemcpyFromSymbol(&momt->kernel,
      p_m3i_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder], sizeof(momf_t));
    momt->num_mom = vdim;
  }
  else if (strcmp(mom, "M3ijk") == 0) { // heat-flux tensor in lab-frame
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3ijk_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    cudaMemcpyFromSymbol(&momt->kernel,
      p_m3ijk_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder], sizeof(momf_t));

    int m3ijk_count[] = { 1, 4, 10 };
    momt->num_mom = m3ijk_count[vdim-1];
  }
  else {
    // string not recognized
    gkyl_exit("gkyl_vlasov_mom_type: Unrecognized moment requested!");
  }

  // copy struct to device
  struct gkyl_mom_type *momt_cu = (struct gkyl_mom_type*) gkyl_cu_malloc(sizeof(struct gkyl_mom_type));
  gkyl_cu_memcpy(momt_cu, momt, sizeof(struct gkyl_mom_type), GKYL_CU_MEMCPY_H2D);
  gkyl_free(momt);
    
  return momt_cu;
}
