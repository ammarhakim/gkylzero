#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_mom.h>
#include <gkyl_vlasov_mom_priv.h>

static void
mom_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *momt = container_of(ref, struct gkyl_mom_type, ref_count);
  gkyl_free(momt);
}

struct gkyl_mom_type*
gkyl_vlasov_mom_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom)
{
  assert(cbasis->polyOrder == pbasis->polyOrder);
  
  struct gkyl_mom_type *momt = gkyl_malloc(sizeof(struct gkyl_mom_type));
  int cdim = momt->cdim = cbasis->ndim;
  int pdim = momt->pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int polyOrder = momt->polyOrder = cbasis->polyOrder;
  momt->num_config = cbasis->numBasis;
  momt->num_phase = pbasis->numBasis;

  if (strcmp(mom, "M0") == 0) { // density
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    momt->kernel = m0_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder];
    momt->num_mom = 1;
  }
  else if (strcmp(mom, "M1i") == 0) { // momentum
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m1i_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    momt->kernel = m1i_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder];
    momt->num_mom = vdim;
  }
  else if (strcmp(mom, "M2") == 0) { // energy
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    momt->kernel = m2_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder];
    momt->num_mom = 1;
  }
  else if (strcmp(mom, "M2ij") == 0) { // pressure tensor in lab-frame
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2ij_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    momt->kernel = m2ij_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder];
    momt->num_mom = 1/2*vdim*(vdim+1);
  }
  else if (strcmp(mom, "M3i") == 0) { // heat-flux vector in lab-frame
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3i_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    momt->kernel = m3i_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder];
    momt->num_mom = vdim;
  }
  else if (strcmp(mom, "M3ijk") == 0) { // heat-flux tensor in lab-frame
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3ijk_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    momt->kernel = m3ijk_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder];
    momt->num_mom = vdim;
  }
  else {
    // string not recognized
    gkyl_exit("gkyl_vlasov_mom_type: Unrecognized moment requested!");
  }

  // set reference counter
  momt->ref_count = (struct gkyl_ref_count) { mom_free, 1 };
    
  return momt;
}
