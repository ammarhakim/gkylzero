/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_vlasov_sr.h>
#include <gkyl_mom_vlasov_sr_priv.h>
#include <gkyl_util.h>
}

enum { M0, M1i, BAD };

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
  else {
    mom_idx = BAD;
  }    

  return mom_idx;
}

static int
v_num_mom(int vdim, int mom_id)
{
  int num_mom = 0;
  
  switch (mom_id) {
    case M0:
      num_mom = 1;
      break;

    case M1i:
      num_mom = vdim;
      break;   
      
    default: // can't happen
      break;
  }

  return num_mom;
}

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_mom_vlasov_sr_set_auxfields_cu_kernel(const struct gkyl_mom_type *momt, const struct gkyl_array *p_over_gamma)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);
  mom_vm_sr->auxfields.p_over_gamma = p_over_gamma;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_mom_vlasov_sr_set_auxfields_cu(const struct gkyl_mom_type *momt, struct gkyl_mom_vlasov_sr_auxfields auxin)
{
  gkyl_mom_vlasov_sr_set_auxfields_cu_kernel<<<1,1>>>(momt, auxin.p_over_gamma->on_dev);
}


__global__
static void
set_cu_ptrs(struct mom_type_vlasov_sr* momt, int mom_id, enum gkyl_basis_type b_type, int vdim,
  int poly_order, int tblidx)
{
  momt->auxfields.p_over_gamma= 0; 
  momt->momt.kernel = kernel;
  
  // choose kernel tables based on basis-function type
  const gkyl_mom_kern_list *m0_kernels, *m1i_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m1i_kernels = ser_m1i_kernels;
      break;

    // case GKYL_BASIS_MODAL_TENSOR:
    //   m0_kernels = ten_m0_kernels;
    //   m1i_kernels = ten_m1i_kernels;
    //   break;

    default:
      assert(false);
      break;    
  }  
  
  switch (mom_id) {
    case M0:
      momt->kernel = m0_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = 1;
      break;

    case M1i:
      momt->kernel = m1i_kernels[tblidx].kernels[poly_order];
      momt->momt.num_mom = vdim;
      break;
      
    default: // can't happen
      break;
  }
}

struct gkyl_mom_type*
gkyl_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* vel_range,
  const char *mom)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_vlasov_sr *momt = (struct mom_type_vlasov_sr*)
    gkyl_malloc(sizeof(struct mom_type_vlasov_sr));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  momt->momt.cdim = cdim;
  momt->momt.pdim = pdim;
  momt->momt.poly_order = poly_order;
  momt->momt.num_config = cbasis->num_basis;
  momt->momt.num_phase = pbasis->num_basis;

  momt->vel_range = *vel_range;

  int mom_id = get_mom_id(mom);
  assert(mom_id != BAD);
  momt->momt.num_mom = v_num_mom(vdim, mom_id); // number of moments

  momt->momt.flags = 0;
  GKYL_SET_CU_ALLOC(momt->momt.flags);
  momt->momt.ref_count = gkyl_ref_count_init(gkyl_mom_vm_sr_free);
  
  // copy struct to device
  struct mom_type_vlasov_sr *momt_cu = (struct mom_type_vlasov_sr*)
    gkyl_cu_malloc(sizeof(struct mom_type_vlasov_sr));
  gkyl_cu_memcpy(momt_cu, momt, sizeof(struct mom_type_vlasov_sr), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);

  set_cu_ptrs<<<1,1>>>(momt_cu, mom_id, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  momt->momt.on_dev = &momt_cu->momt;
  
  return &momt->momt;
}

__global__
static void
set_int_cu_ptrs(struct mom_type_vlasov* momt, enum gkyl_basis_type b_type, int vdim,
  int poly_order, int tblidx)
{
  momt->auxfields.p_over_gamma= 0; 
  momt->momt.kernel = kernel;

  // set kernel pointer
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      momt->kernel = ser_int_mom_kernels[tblidx].kernels[poly_order];
      break;

    // case GKYL_BASIS_MODAL_TENSOR:
    //   momt->kernel = ten_int_mom_kernels[tblidx].kernels[poly_order];
    //   break;

    default:
      assert(false);
      break;    
  }
}

struct gkyl_mom_type *
gkyl_int_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_vlasov_sr *momt = (struct mom_type_vlasov_sr*)
    gkyl_malloc(sizeof(struct mom_type_vlasov_sr));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  momt->momt.cdim = cdim;
  momt->momt.pdim = pdim;
  momt->momt.poly_order = poly_order;
  momt->momt.num_config = cbasis->num_basis;
  momt->momt.num_phase = pbasis->num_basis;

  momt->momt.num_mom = vdim+2;

  momt->vel_range = *vel_range;

  momt->momt.flags = 0;
  GKYL_SET_CU_ALLOC(momt->momt.flags);
  momt->momt.ref_count = gkyl_ref_count_init(gkyl_mom_vm_sr_free);
  
  // copy struct to device
  struct mom_type_vlasov *momt_cu = (struct mom_type_vlasov_sr*)
    gkyl_cu_malloc(sizeof(struct mom_type_vlasov_sr));
  gkyl_cu_memcpy(momt_cu, momt, sizeof(struct mom_type_vlasov_sr), GKYL_CU_MEMCPY_H2D);

  set_int_cu_ptrs<<<1,1>>>(momt_cu, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  momt->momt.on_dev = &momt_cu->momt;
  
  return &momt->momt;
}
