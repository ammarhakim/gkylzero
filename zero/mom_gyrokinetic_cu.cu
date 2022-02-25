/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_gyrokinetic.h>
#include <gkyl_mom_gyrokinetic_priv.h>
#include <gkyl_util.h>
}

enum { GkM0, GkM1, GkM2, GkM2par, GkM2perp, GkM3par, GkM3perp, ThreeMoments, BAD };

static int
get_gk_mom_id(const char *mom)
{
  int mom_idx = BAD;

  if (strcmp(mom, "GkM0") == 0) { // density
    mom_idx = GkM0;
  }
  else if (strcmp(mom, "GkM1") == 0) { // parallel momentum
    mom_idx = GkM1;
  }
  else if (strcmp(mom, "GkM2") == 0) { // total energy
    mom_idx = GkM2;
  }
  else if (strcmp(mom, "GkM2par") == 0) { // parallel energy
    mom_idx = GkM2par;
  }
  else if (strcmp(mom, "GkM2perp") == 0) { // perpendicular energy
    mom_idx = GkM2perp;
  }
  else if (strcmp(mom, "GkM3par") == 0) { // parallel heat flux
    mom_idx = GkM3par;
  }
  else if (strcmp(mom, "GkM3perp") == 0) { // perpendicular heat flux
    mom_idx = GkM3perp;
  }
  else if (strcmp(mom, "ThreeMoments") == 0) { // Zeroth (density), First (parallel momentum), 
    mom_idx = ThreeMoments;                    // and Second (total energy) computed together
  }
  else {
    mom_idx = BAD;
  }    

  return mom_idx;
}

static int
gk_num_mom(int vdim, int mom_id)
{
  int num_mom = 0;
  
  switch (mom_id) {
    case GkM0:
      num_mom = 1;
      break;

    case GkM1:
      num_mom = 1;
      break;

    case GkM2:
      num_mom = 1;
      break;

    case GkM2par:
      num_mom = 1;
      break;

    case GkM2perp:
      num_mom = 1;
      break;

    case GkM3par:
      num_mom = 1;
      break;

    case GkM3perp:
      num_mom = 1;
      break;

    case ThreeMoments:
      num_mom = 3;
      break;      
      
    default: // can't happen
      break;
  }

  return num_mom;
}

// CUDA kernel to set pointer to bmag
// This is required because moment type object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_gyrokinetic_set_bmag_cu_kernel(const struct gkyl_mom_type *momt, const struct gkyl_array *bmag)
{
  struct mom_type_gyrokinetic *mom_gyrokinetic = container_of(momt, struct mom_type_gyrokinetic, momt);
  mom_gyrokinetic->bmag = bmag;
}

// Host-side wrapper for set_bmag_cu_kernel
void
gkyl_gyrokinetic_set_bmag_cu(const struct gkyl_mom_type *momt, const struct gkyl_array *bmag)
{
  gkyl_gyrokinetic_set_bmag_cu_kernel<<<1,1>>>(momt, bmag->on_dev);
}


__global__
static void
gkyl_mom_gyrokinetic_set_cu_dev_ptrs(struct mom_type_gyrokinetic* mom_gyrokinetic, int mom_id, enum gkyl_basis_type b_type, int vdim,
  int poly_order, int tblidx)
{

  printf("******** FIX BUG IN mom_gyrokinetic to enable it to run on GPUs!");    
  assert(false);
  // NOTE: FIX ME. the following line is a problem. However, the issue
  // appears in the priv header and not here, apparently. The problem
  // is the return statement exactly like the issue with vlasov_poisson volume
  
  //  mom_gyrokinetic->momt.kernel = kernel;  
  
  // choose kernel tables based on basis-function type
  const gkyl_gyrokinetic_mom_kern_list *m0_kernels, *m1_kernels, *m2_kernels, 
    *m2_par_kernels, *m2_perp_kernels, *m3_par_kernels, *m3_perp_kernels, *three_moments_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m1_kernels = ser_m1_kernels;
      m2_kernels = ser_m2_kernels;
      m2_par_kernels = ser_m2_par_kernels;
      m2_perp_kernels = ser_m2_perp_kernels;
      m3_par_kernels = ser_m3_par_kernels;
      m3_perp_kernels = ser_m3_perp_kernels;
      three_moments_kernels = ser_three_moments_kernels;
      break;

    case GKYL_BASIS_MODAL_TENSOR:
      m0_kernels = ten_m0_kernels;
      m1_kernels = ten_m1_kernels;
      m2_kernels = ten_m2_kernels;
      m2_par_kernels = ten_m2_par_kernels;
      m2_perp_kernels = ten_m2_perp_kernels;
      m3_par_kernels = ten_m3_par_kernels;
      m3_perp_kernels = ten_m3_perp_kernels;
      three_moments_kernels = ten_three_moments_kernels;
      break;

    default:
      assert(false);
      break;    
  }  
  
  switch (mom_id) {
    case GkM0:
      mom_gyrokinetic->kernel = m0_kernels[tblidx].kernels[poly_order];
      mom_gyrokinetic->momt.num_mom = 1;
      break;

    case GkM1:
      mom_gyrokinetic->kernel = m1_kernels[tblidx].kernels[poly_order];
      mom_gyrokinetic->momt.num_mom = 1;
      break;

    case GkM2:
      mom_gyrokinetic->kernel = m2_kernels[tblidx].kernels[poly_order];
      mom_gyrokinetic->momt.num_mom = 1;
      break;

    case GkM2par:
      mom_gyrokinetic->kernel = m2_par_kernels[tblidx].kernels[poly_order];
      mom_gyrokinetic->momt.num_mom = 1;
      break;

    case GkM2perp:
      mom_gyrokinetic->kernel = m2_perp_kernels[tblidx].kernels[poly_order];
      mom_gyrokinetic->momt.num_mom = 1;
      break;

    case GkM3par:
      mom_gyrokinetic->kernel = m3_par_kernels[tblidx].kernels[poly_order];
      mom_gyrokinetic->momt.num_mom = 1;
      break;

    case GkM3perp:
      mom_gyrokinetic->kernel = m3_perp_kernels[tblidx].kernels[poly_order];
      mom_gyrokinetic->momt.num_mom = 1;
      break;

    case ThreeMoments:
      mom_gyrokinetic->kernel = three_moments_kernels[tblidx].kernels[poly_order];
      mom_gyrokinetic->momt.num_mom = 3;
      break;
      
    default: // can't happen
      break;
  }
}

struct gkyl_mom_type*
gkyl_mom_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, double mass, const char *mom)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_gyrokinetic *mom_gyrokinetic = (struct mom_type_gyrokinetic*)
    gkyl_malloc(sizeof(struct mom_type_gyrokinetic));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_gyrokinetic->momt.cdim = cdim;
  mom_gyrokinetic->momt.pdim = pdim;
  mom_gyrokinetic->momt.poly_order = poly_order;
  mom_gyrokinetic->momt.num_config = cbasis->num_basis;
  mom_gyrokinetic->momt.num_phase = pbasis->num_basis;

  int mom_id = get_gk_mom_id(mom);
  assert(mom_id != BAD);
  mom_gyrokinetic->momt.num_mom = gk_num_mom(vdim, mom_id); // number of moments

  mom_gyrokinetic->mass = mass;
  mom_gyrokinetic->conf_range = *conf_range;

  mom_gyrokinetic->momt.flag = 0;
  GKYL_SET_CU_ALLOC(mom_gyrokinetic->momt.flag);
  mom_gyrokinetic->momt.ref_count = gkyl_ref_count_init(gkyl_gk_mom_free);
  
  // copy struct to device
  struct mom_type_gyrokinetic *mom_gyrokinetic_cu = (struct mom_type_gyrokinetic*)
    gkyl_cu_malloc(sizeof(struct mom_type_gyrokinetic));
  gkyl_cu_memcpy(mom_gyrokinetic_cu, mom_gyrokinetic, sizeof(struct mom_type_gyrokinetic), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);

  gkyl_mom_gyrokinetic_set_cu_dev_ptrs<<<1,1>>>(mom_gyrokinetic_cu, mom_id, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  mom_gyrokinetic->momt.on_dev = &mom_gyrokinetic_cu->momt;
  
  return &mom_gyrokinetic->momt;
}
