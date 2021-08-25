/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_mom.h>
#include <gkyl_vlasov_mom_priv.h>
}

enum { M0, M1i, M2, M2ij, M3i, M3ijk, BAD };

static int
get_mom_id(const char *mom)
{
  int mom_idx = BAD;

  if (strcmp(mom, "M0") == 0) { // density
    mom_idx = M0;
  }
  else if (strcmp(mom, "M1i") == 0) { // momentum
    mom_idx = M1i;
  }
  else if (strcmp(mom, "M2") == 0) { // energy
    mom_idx = M2;
  }
  else if (strcmp(mom, "M2ij") == 0) { // pressure tensor in lab-frame
    mom_idx = M2ij;    
  }
  else if (strcmp(mom, "M3i") == 0) { // heat-flux vector in lab-frame
    mom_idx = M3i;
  }
  else if (strcmp(mom, "M3ijk") == 0) { // heat-flux tensor in lab-frame
    mom_idx = M3ijk;
  }
  else {
    mom_idx = BAD;
  }    

  return mom_idx;
}

__global__
static void
vlasov_mom_set_cu_dev_ptrs(struct gkyl_mom_type* momt, int mom_id, enum gkyl_basis_type b_type, int vdim,
  int poly_order, int tblidx)
{
  int m3ijk_count[] = { 1, 4, 10 };

  // choose kernel tables based on basis-function type
  const mom_kern_list *m0_kernels, *m1i_kernels,
    *m2_kernels, *m2ij_kernels, *m3i_kernels, *m3ijk_kernels;

  switch (b_type) {
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
  
  switch (mom_id) {
    case M0:
      momt->kernel = m0_kernels[tblidx].kernels[poly_order];
      momt->num_mom = 1;
      break;

    case M1i:
      momt->kernel = m1i_kernels[tblidx].kernels[poly_order];
      momt->num_mom = vdim;
      break;

    case M2:
      momt->kernel = m2_kernels[tblidx].kernels[poly_order];
      momt->num_mom = 1;
      break;

    case M2ij:
      momt->kernel = m2ij_kernels[tblidx].kernels[poly_order];
      momt->num_mom = vdim*(vdim+1)/2;
      break;

    case M3i:
      momt->kernel = m3i_kernels[tblidx].kernels[poly_order];
      momt->num_mom = vdim;
      break;

    case M3ijk:
      momt->kernel = m3ijk_kernels[tblidx].kernels[poly_order];
      momt->num_mom = m3ijk_count[vdim-1];
      break;
      
    default: // can't happen
      break;
  }
}

struct gkyl_mom_type*
gkyl_vlasov_mom_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct gkyl_mom_type *momt = (struct gkyl_mom_type*) gkyl_malloc(sizeof(struct gkyl_mom_type));
  int cdim = momt->cdim = cbasis->ndim;
  int pdim = momt->pdim = pbasis->ndim;
  int vdim = pdim-cdim;
  int poly_order = momt->poly_order = cbasis->poly_order;
  momt->num_config = cbasis->num_basis;
  momt->num_phase = pbasis->num_basis;

  // copy struct to device
  struct gkyl_mom_type *momt_cu = (struct gkyl_mom_type*) gkyl_cu_malloc(sizeof(struct gkyl_mom_type));
  gkyl_cu_memcpy(momt_cu, momt, sizeof(struct gkyl_mom_type), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);
  int mom_id = get_mom_id(mom);
  assert(mom_id != BAD);

  vlasov_mom_set_cu_dev_ptrs<<<1,1>>>(momt_cu, mom_id, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);
  
  gkyl_free(momt);
    
  return momt_cu;
}
