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
vlasov_mom_set_cu_dev_ptrs(struct vlasov_mom_type* vlasov_mom, int mom_id, enum gkyl_basis_type b_type, int vdim,
  int poly_order, int tblidx)
{
  int m3ijk_count[] = { 1, 4, 10 };

  vlasov_mom->momt.kernel = kernel;  
  
  // choose kernel tables based on basis-function type
  const gkyl_mom_kern_list *m0_kernels, *m1i_kernels,
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
      vlasov_mom->kernel = m0_kernels[tblidx].kernels[poly_order];
      vlasov_mom->momt.num_mom = 1;
      break;

    case M1i:
      vlasov_mom->kernel = m1i_kernels[tblidx].kernels[poly_order];
      vlasov_mom->momt.num_mom = vdim;
      break;

    case M2:
      vlasov_mom->kernel = m2_kernels[tblidx].kernels[poly_order];
      vlasov_mom->momt.num_mom = 1;
      break;

    case M2ij:
      vlasov_mom->kernel = m2ij_kernels[tblidx].kernels[poly_order];
      vlasov_mom->momt.num_mom = vdim*(vdim+1)/2;
      break;

    case M3i:
      vlasov_mom->kernel = m3i_kernels[tblidx].kernels[poly_order];
      vlasov_mom->momt.num_mom = vdim;
      break;

    case M3ijk:
      vlasov_mom->kernel = m3ijk_kernels[tblidx].kernels[poly_order];
      vlasov_mom->momt.num_mom = m3ijk_count[vdim-1];
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

  struct vlasov_mom_type *vlasov_mom = (struct vlasov_mom_type*) gkyl_malloc(sizeof(struct vlasov_mom_type));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  vlasov_mom->momt.cdim = cdim;
  vlasov_mom->momt.pdim = pdim;
  vlasov_mom->momt.poly_order = poly_order;
  vlasov_mom->momt.num_config = cbasis->num_basis;
  vlasov_mom->momt.num_phase = pbasis->num_basis;
  
  // copy struct to device
  struct vlasov_mom_type *vlasov_mom_cu = (struct vlasov_mom_type*) gkyl_cu_malloc(sizeof(struct vlasov_mom_type));
  gkyl_cu_memcpy(vlasov_mom_cu, vlasov_mom, sizeof(struct vlasov_mom_type), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);
  int mom_id = get_mom_id(mom);
  assert(mom_id != BAD);

  vlasov_mom_set_cu_dev_ptrs<<<1,1>>>(vlasov_mom_cu, mom_id, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);
  
  gkyl_free(vlasov_mom);
    
  return &vlasov_mom_cu->momt;
}
