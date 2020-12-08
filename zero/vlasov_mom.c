#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_mom.h>
#include <kernels/vlasov/gkyl_vlasov_mom_kernels.h>

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// M0 kernel list
static struct { momf_t kernels[3]; } m0_kernels[] = {
  // 1x kernels
  { NULL, vlasov_mom_1x1v_m0_p1, vlasov_mom_1x1v_m0_p2 }, // 0
  { NULL, vlasov_mom_1x2v_m0_p1, vlasov_mom_1x2v_m0_p2 }, // 1
  { NULL, vlasov_mom_1x3v_m0_p1, vlasov_mom_1x3v_m0_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_mom_2x2v_m0_p1, vlasov_mom_2x2v_m0_p2 }, // 3
  { NULL, vlasov_mom_2x3v_m0_p1, vlasov_mom_2x3v_m0_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_mom_3x3v_m0_p1, NULL                  }, // 5
};

// M1i kernel list
static struct { momf_t kernels[3]; } m1i_kernels[] = {
  // 1x kernels
  { NULL, vlasov_mom_1x1v_m1i_p1, vlasov_mom_1x1v_m1i_p2 }, // 0
  { NULL, vlasov_mom_1x2v_m1i_p1, vlasov_mom_1x2v_m1i_p2 }, // 1
  { NULL, vlasov_mom_1x3v_m1i_p1, vlasov_mom_1x3v_m1i_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_mom_2x2v_m1i_p1, vlasov_mom_2x2v_m1i_p2 }, // 3
  { NULL, vlasov_mom_2x3v_m1i_p1, vlasov_mom_2x3v_m1i_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_mom_3x3v_m1i_p1, NULL                   }, // 5
};

// M2ij kernel list
static struct { momf_t kernels[3]; } m2ij_kernels[] = {
  // 1x kernels
  { NULL, vlasov_mom_1x1v_m2ij_p1, vlasov_mom_1x1v_m2ij_p2 }, // 0
  { NULL, vlasov_mom_1x2v_m2ij_p1, vlasov_mom_1x2v_m2ij_p2 }, // 1
  { NULL, vlasov_mom_1x3v_m2ij_p1, vlasov_mom_1x3v_m2ij_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_mom_2x2v_m2ij_p1, vlasov_mom_2x2v_m2ij_p2 }, // 3
  { NULL, vlasov_mom_2x3v_m2ij_p1, vlasov_mom_2x3v_m2ij_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_mom_3x3v_m2ij_p1, NULL                   }, // 5
};

// M2 kernel list
static struct { momf_t kernels[3]; } m2_kernels[] = {
  // 1x kernels
  { NULL, vlasov_mom_1x1v_m2_p1, vlasov_mom_1x1v_m2_p2 }, // 0
  { NULL, vlasov_mom_1x2v_m2_p1, vlasov_mom_1x2v_m2_p2 }, // 1
  { NULL, vlasov_mom_1x3v_m2_p1, vlasov_mom_1x3v_m2_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_mom_2x2v_m2_p1, vlasov_mom_2x2v_m2_p2 }, // 3
  { NULL, vlasov_mom_2x3v_m2_p1, vlasov_mom_2x3v_m2_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_mom_3x3v_m2_p1, NULL                   }, // 5
};

// M3i kernel list
static struct { momf_t kernels[3]; } m3i_kernels[] = {
  // 1x kernels
  { NULL, vlasov_mom_1x1v_m3i_p1, vlasov_mom_1x1v_m3i_p2 }, // 0
  { NULL, vlasov_mom_1x2v_m3i_p1, vlasov_mom_1x2v_m3i_p2 }, // 1
  { NULL, vlasov_mom_1x3v_m3i_p1, vlasov_mom_1x3v_m3i_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_mom_2x2v_m3i_p1, vlasov_mom_2x2v_m3i_p2 }, // 3
  { NULL, vlasov_mom_2x3v_m3i_p1, vlasov_mom_2x3v_m3i_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_mom_3x3v_m3i_p1, NULL                   }, // 5
};

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
  else if (strcmp(mom, "M2ij") == 0) { // pressure tensor in lab-frame
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2ij_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    momt->kernel = m2ij_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder];
    momt->num_mom = 1/2*vdim*(vdim+1);
  }
  else if (strcmp(mom, "M2") == 0) { // trace of M2ij
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m2_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    momt->kernel = m2_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder];
    momt->num_mom = 1;
  }
  else if (strcmp(mom, "M3i") == 0) { // heat-flux vector in lab-frame
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m3i_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder]);
    
    momt->kernel = m3i_kernels[cv_index[cdim].vdim[vdim]].kernels[polyOrder];
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
