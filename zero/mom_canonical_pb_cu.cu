/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_canonical_pb.h>
#include <gkyl_mom_canonical_pb_priv.h>
#include <gkyl_util.h>
}

enum { MEnergy, BAD };

static int
get_mom_id(const char *mom)
{
  int mom_idx = BAD;

  if (strcmp(mom, "MEnergy") == 0) { // total energy = integral(hamil*f) velocity moment
    mom_idx = MEnergy;
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
    case MEnergy:
      num_mom = 1;
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
gkyl_mom_canonical_pb_set_auxfields_cu_kernel(const struct gkyl_mom_type *momt, 
  const struct gkyl_array *hamil)
{
  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  mom_can_pb->auxfields.hamil = hamil;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_mom_canonical_pb_set_auxfields_cu(const struct gkyl_mom_type *momt, struct gkyl_mom_canonical_pb_auxfields auxin)
{
  gkyl_mom_canonical_pb_set_auxfields_cu_kernel<<<1,1>>>(momt, auxin.hamil->on_dev);
}


__global__
static void
set_cu_ptrs(struct mom_type_canonical_pb* mom_can_pb, int mom_id, enum gkyl_basis_type b_type, int vdim,
  int poly_order, int tblidx)
{
  mom_can_pb->auxfields.hamil = 0;
  
  // choose kernel tables based on basis-function type
  const gkyl_canonical_pb_mom_kern_list *menergy_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      menergy_kernels = ser_menergy_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      menergy_kernels = tensor_menergy_kernels;
      break;

    default:
      assert(false);
      break;    
  }
  
  switch (mom_id) {
    case MEnergy:
      mom_can_pb->momt.kernel = menergy_kernels[tblidx].kernels[poly_order];
      mom_can_pb->momt.num_mom = 1;
      break;

    default: // can't happen
      break;
  }
}

struct gkyl_mom_type*
gkyl_mom_canonical_pb_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range, const char *mom)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_canonical_pb *mom_can_pb = (struct mom_type_canonical_pb*)
    gkyl_malloc(sizeof(struct mom_type_canonical_pb));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_can_pb->momt.cdim = cdim;
  mom_can_pb->momt.pdim = pdim;
  mom_can_pb->momt.poly_order = poly_order;
  mom_can_pb->momt.num_config = cbasis->num_basis;
  mom_can_pb->momt.num_phase = pbasis->num_basis;

  mom_can_pb->phase_range = *phase_range;

  int mom_id = get_mom_id(mom);
  assert(mom_id != BAD);
  mom_can_pb->momt.num_mom = v_num_mom(vdim, mom_id); // number of moments

  mom_can_pb->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_can_pb->momt.flags);
  mom_can_pb->momt.ref_count = gkyl_ref_count_init(gkyl_mom_can_pb_free);
  
  // copy struct to device
  struct mom_type_canonical_pb *momt_cu = (struct mom_type_canonical_pb*)
    gkyl_cu_malloc(sizeof(struct mom_type_canonical_pb));
  gkyl_cu_memcpy(momt_cu, mom_can_pb, sizeof(struct mom_type_canonical_pb), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);

  set_cu_ptrs<<<1,1>>>(momt_cu, mom_id, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  mom_can_pb->momt.on_dev = &momt_cu->momt;
  
  return &mom_can_pb->momt;
}

__global__
static void
set_int_cu_ptrs(struct mom_type_canonical_pb* mom_can_pb, enum gkyl_basis_type b_type, int vdim,
  int poly_order, int tblidx)
{
  mom_can_pb->auxfields.hamil = 0;

  // choose kernel tables based on basis-function type
  const gkyl_canonical_pb_mom_kern_list *int_mom_kernels;  
  
  // set kernel pointer
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      int_mom_kernels = ser_int_mom_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      int_mom_kernels = tensor_int_mom_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  mom_can_pb->momt.kernel = int_mom_kernels[tblidx].kernels[poly_order];
  mom_can_pb->momt.num_mom = 2+vdim;
}

struct gkyl_mom_type *
gkyl_int_mom_canonical_pb_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_canonical_pb *mom_can_pb = (struct mom_type_canonical_pb*)
    gkyl_malloc(sizeof(struct mom_type_canonical_pb));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_can_pb->momt.cdim = cdim;
  mom_can_pb->momt.pdim = pdim;
  mom_can_pb->momt.poly_order = poly_order;
  mom_can_pb->momt.num_config = cbasis->num_basis;
  mom_can_pb->momt.num_phase = pbasis->num_basis;

  mom_can_pb->momt.num_mom = vdim+2;

  mom_can_pb->phase_range = *phase_range;

  mom_can_pb->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_can_pb->momt.flags);
  mom_can_pb->momt.ref_count = gkyl_ref_count_init(gkyl_mom_can_pb_free);
  
  // copy struct to device
  struct mom_type_canonical_pb *momt_cu = (struct mom_type_canonical_pb*)
    gkyl_cu_malloc(sizeof(struct mom_type_canonical_pb));
  gkyl_cu_memcpy(momt_cu, mom_can_pb, sizeof(struct mom_type_canonical_pb), GKYL_CU_MEMCPY_H2D);

  set_int_cu_ptrs<<<1,1>>>(momt_cu, pbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  mom_can_pb->momt.on_dev = &momt_cu->momt;
  
  return &mom_can_pb->momt;
}
