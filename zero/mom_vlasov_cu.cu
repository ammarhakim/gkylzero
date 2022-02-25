/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_vlasov.h>
#include <gkyl_mom_vlasov_priv.h>
#include <gkyl_util.h>
}

enum { M0, M1i, M2, M2ij, M3i, M3ijk, FiveMoments, BAD };

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
  else if (strcmp(mom, "FiveMoments") == 0) { // heat-flux tensor in lab-frame
    mom_idx = FiveMoments;
  }
  else {
    mom_idx = BAD;
  }    

  return mom_idx;
}

static int
v_num_mom(int vdim, int mom_id)
{
  int m3ijk_count[] = { 1, 4, 10 };
  int num_mom = 0;
  
  switch (mom_id) {
    case M0:
      num_mom = 1;
      break;

    case M1i:
      num_mom = vdim;
      break;

    case M2:
      num_mom = 1;
      break;

    case M2ij:
      num_mom = vdim*(vdim+1)/2;
      break;

    case M3i:
      num_mom = vdim;
      break;

    case M3ijk:
      num_mom = m3ijk_count[vdim-1];
      break;

    case FiveMoments:
      num_mom = vdim+2;
      break;      
      
    default: // can't happen
      break;
  }

  return num_mom;
}

__global__
static void
gkyl_mom_vlasov_set_cu_dev_ptrs(struct mom_type_vlasov* mom_vlasov, int mom_id, enum gkyl_basis_type b_type, int vdim,
  int poly_order, int tblidx)
{
  int m3ijk_count[] = { 1, 4, 10 };

  mom_vlasov->momt.kernel = kernel;  
  
  // choose kernel tables based on basis-function type
  const gkyl_mom_kern_list *m0_kernels, *m1i_kernels,
    *m2_kernels, *m2ij_kernels, *m3i_kernels, *m3ijk_kernels, *five_moments_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m1i_kernels = ser_m1i_kernels;
      m2_kernels = ser_m2_kernels;
      m2ij_kernels = ser_m2ij_kernels;
      m3i_kernels = ser_m3i_kernels;
      m3ijk_kernels = ser_m3ijk_kernels;
      five_moments_kernels = ser_five_moments_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      m0_kernels = ten_m0_kernels;
      m1i_kernels = ten_m1i_kernels;
      m2_kernels = ten_m2_kernels;
      m2ij_kernels = ten_m2ij_kernels;
      m3i_kernels = ten_m3i_kernels;
      m3ijk_kernels = ten_m3ijk_kernels;
      five_moments_kernels = ten_five_moments_kernels;
      break;

    default:
      assert(false);
      break;    
  }  
  
  switch (mom_id) {
    case M0:
      mom_vlasov->kernel = m0_kernels[tblidx].kernels[poly_order];
      mom_vlasov->momt.num_mom = 1;
      break;

    case M1i:
      mom_vlasov->kernel = m1i_kernels[tblidx].kernels[poly_order];
      mom_vlasov->momt.num_mom = vdim;
      break;

    case M2:
      mom_vlasov->kernel = m2_kernels[tblidx].kernels[poly_order];
      mom_vlasov->momt.num_mom = 1;
      break;

    case M2ij:
      mom_vlasov->kernel = m2ij_kernels[tblidx].kernels[poly_order];
      mom_vlasov->momt.num_mom = vdim*(vdim+1)/2;
      break;

    case M3i:
      mom_vlasov->kernel = m3i_kernels[tblidx].kernels[poly_order];
      mom_vlasov->momt.num_mom = vdim;
      break;

    case M3ijk:
      mom_vlasov->kernel = m3ijk_kernels[tblidx].kernels[poly_order];
      mom_vlasov->momt.num_mom = m3ijk_count[vdim-1];
      break;

    case FiveMoments:
      mom_vlasov->kernel = five_moments_kernels[tblidx].kernels[poly_order];
      mom_vlasov->momt.num_mom = vdim+2;
      break;
      
    default: // can't happen
      break;
  }
}

struct gkyl_mom_type*
gkyl_mom_vlasov_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *mom)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_vlasov *mom_vlasov = (struct mom_type_vlasov*)
    gkyl_malloc(sizeof(struct mom_type_vlasov));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_vlasov->momt.cdim = cdim;
  mom_vlasov->momt.pdim = pdim;
  mom_vlasov->momt.poly_order = poly_order;
  mom_vlasov->momt.num_config = cbasis->num_basis;
  mom_vlasov->momt.num_phase = pbasis->num_basis;

  int mom_id = get_mom_id(mom);
  assert(mom_id != BAD);
  mom_vlasov->momt.num_mom = v_num_mom(vdim, mom_id); // number of moments

  mom_vlasov->momt.flag = 0;
  GKYL_SET_CU_ALLOC(mom_vlasov->momt.flag);
  mom_vlasov->momt.ref_count = gkyl_ref_count_init(gkyl_mom_free);
  
  // copy struct to device
  struct mom_type_vlasov *mom_vlasov_cu = (struct mom_type_vlasov*)
    gkyl_cu_malloc(sizeof(struct mom_type_vlasov));
  gkyl_cu_memcpy(mom_vlasov_cu, mom_vlasov, sizeof(struct mom_type_vlasov), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);

  gkyl_mom_vlasov_set_cu_dev_ptrs<<<1,1>>>(mom_vlasov_cu, mom_id, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  mom_vlasov->momt.on_dev = &mom_vlasov_cu->momt;
  
  return &mom_vlasov->momt;
}
