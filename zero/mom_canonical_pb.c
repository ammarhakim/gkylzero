#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_canonical_pb.h>
#include <gkyl_mom_canonical_pb_priv.h>
#include <gkyl_util.h>

void
gkyl_mom_can_pb_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *momt = container_of(ref, struct gkyl_mom_type, ref_count);
  if (GKYL_IS_CU_ALLOC(momt->flags))
    gkyl_cu_free(momt->on_dev);
  gkyl_free(momt);
}

void
gkyl_mom_canonical_pb_set_auxfields(const struct gkyl_mom_type *momt, struct gkyl_mom_canonical_pb_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_mom_type_is_cu_dev(momt)) {
    gkyl_mom_canonical_pb_set_auxfields_cu(momt->on_dev, auxin);
    return;
  }
#endif

  struct mom_type_canonical_pb *mom_can_pb = container_of(momt, struct mom_type_canonical_pb, momt);
  mom_can_pb->auxfields.hamil = auxin.hamil;
}

struct gkyl_mom_type*
gkyl_mom_canonical_pb_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range, const char *mom, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_mom_canonical_pb_cu_dev_new(cbasis, pbasis, phase_range, mom);
  } 
#endif  
  struct mom_type_canonical_pb *mom_can_pb = gkyl_malloc(sizeof(struct mom_type_canonical_pb));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_can_pb->momt.cdim = cdim;
  mom_can_pb->momt.pdim = pdim;
  mom_can_pb->momt.poly_order = poly_order;
  mom_can_pb->momt.num_config = cbasis->num_basis;
  mom_can_pb->momt.num_phase = pbasis->num_basis;

  // choose kernel tables based on basis-function type
  const gkyl_canonical_pb_mom_kern_list *menergy_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      menergy_kernels = ser_menergy_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  if (strcmp(mom, "MEnergy") == 0) { // Energy int( f*H ) 
      assert(cdim == vdim); // Can-pb needs equal vdim and cdim
      assert(cv_index[cdim].vdim[vdim] != -1);
      assert(NULL != menergy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
      
      mom_can_pb->momt.kernel = menergy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
      mom_can_pb->momt.num_mom = 1;
  }
  else {
    // string not recognized
    gkyl_exit("gkyl_mom_type_canonical_pb: Unrecognized moment requested!");
  }

  mom_can_pb->phase_range = *phase_range;

  mom_can_pb->auxfields.hamil = 0;

  mom_can_pb->momt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mom_can_pb->momt.flags);
  mom_can_pb->momt.ref_count = gkyl_ref_count_init(gkyl_mom_can_pb_free);
  
  mom_can_pb->momt.on_dev = &mom_can_pb->momt; // on host, self-reference
    
  return &mom_can_pb->momt;
}

struct gkyl_mom_type*
gkyl_int_mom_canonical_pb_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_int_mom_canonical_pb_cu_dev_new(cbasis, pbasis, phase_range);
  } 
#endif  
  struct mom_type_canonical_pb *mom_can_pb = gkyl_malloc(sizeof(struct mom_type_canonical_pb));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_can_pb->momt.cdim = cdim;
  mom_can_pb->momt.pdim = pdim;
  mom_can_pb->momt.poly_order = poly_order;
  mom_can_pb->momt.num_config = cbasis->num_basis;
  mom_can_pb->momt.num_phase = pbasis->num_basis;

  // choose kernel tables based on basis-function type
  const gkyl_canonical_pb_mom_kern_list *int_mom_kernels;  
  
  // set kernel pointer
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      int_mom_kernels = ser_int_mom_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  assert(cv_index[cdim].vdim[vdim] != -1);
  assert(NULL != int_mom_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
  mom_can_pb->momt.kernel = int_mom_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  mom_can_pb->momt.num_mom = vdim+2;

  mom_can_pb->phase_range = *phase_range;

  mom_can_pb->auxfields.hamil = 0;

  mom_can_pb->momt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mom_can_pb->momt.flags);
  mom_can_pb->momt.ref_count = gkyl_ref_count_init(gkyl_mom_can_pb_free);
  
  mom_can_pb->momt.on_dev = &mom_can_pb->momt; // on host, self-reference
    
  return &mom_can_pb->momt;  
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_mom_type*
gkyl_mom_canonical_pb_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range, const char *mom)
{
  assert(false);
  return 0;
}

struct gkyl_mom_type *
gkyl_int_mom_canonical_pb_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* phase_range)
{
  assert(false);
  return 0;
}

#endif