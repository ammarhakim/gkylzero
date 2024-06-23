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

enum { M0, M1i, Ni, Energy, Pressure, Tij, BAD };

static int
get_mom_id(const char *mom)
{
  int mom_idx = BAD;

  if (strcmp(mom, "M0") == 0) { // density (GammaV*n)
    mom_idx = M0;
  }
  else if (strcmp(mom, "M1i") == 0) { // momentum (GammaV*n*V)
    mom_idx = M1i;
  }
  else if (strcmp(mom, "Ni") == 0) { // 4-momentum (GammaV*n, GammaV*n*V)
    mom_idx = Ni;
  }
  else if (strcmp(mom, "Energy") == 0) { // total energy = gamma*mc^2 moment
    mom_idx = Energy;
  }
  else if (strcmp(mom, "Pressure") == 0) { // total fluid-frame pressure
                                           // P = n*T where n is the fluid-frame density
    mom_idx = Pressure;
  }
  else if (strcmp(mom, "Tij") == 0) { // Stress-energy tensor 
                                      // (Energy, Energy flux (vdim components), Stress tensor (vdim*(vdim+1))/2 components))
    mom_idx = Tij;
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

    case Ni:
      num_mom = vdim+1;
      break;   

    case Energy:
      num_mom = 1;
      break;   
      
    case Pressure:
      num_mom = 1;
      break;   

    case Tij:
      num_mom = 1+vdim+(vdim*(vdim+1))/2;
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
gkyl_mom_vlasov_sr_set_auxfields_cu_kernel(const struct gkyl_mom_type *momt, 
  const struct gkyl_array *gamma, const struct gkyl_array *gamma_inv, 
  const struct gkyl_array *V_drift, const struct gkyl_array *GammaV2)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);
  mom_vm_sr->auxfields.gamma = gamma;
  mom_vm_sr->auxfields.gamma_inv = gamma_inv;
  mom_vm_sr->auxfields.V_drift = V_drift;
  mom_vm_sr->auxfields.GammaV2 = GammaV2;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_mom_vlasov_sr_set_auxfields_cu(const struct gkyl_mom_type *momt, struct gkyl_mom_vlasov_sr_auxfields auxin)
{
  gkyl_mom_vlasov_sr_set_auxfields_cu_kernel<<<1,1>>>(momt, auxin.gamma->on_dev, auxin.gamma_inv->on_dev, 
    auxin.V_drift->on_dev, auxin.GammaV2->on_dev);
}


__global__
static void
set_cu_ptrs(struct mom_type_vlasov_sr* mom_vm_sr, int mom_id, enum gkyl_basis_type b_type, int vdim,
  int poly_order, int tblidx)
{
  mom_vm_sr->auxfields.gamma = 0;
  mom_vm_sr->auxfields.gamma_inv = 0;
  mom_vm_sr->auxfields.V_drift = 0;
  mom_vm_sr->auxfields.GammaV2 = 0;
  
  // choose kernel tables based on basis-function type
  const gkyl_vlasov_sr_mom_kern_list *m0_kernels, *m1i_kernels, 
    *Ni_kernels, *Energy_kernels, *Pressure_kernels, *Tij_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      m0_kernels = ser_m0_kernels;
      m1i_kernels = ser_m1i_kernels;
      Ni_kernels = ser_Ni_kernels;
      Energy_kernels = ser_Energy_kernels;
      Pressure_kernels = ser_Pressure_kernels;
      Tij_kernels = ser_Tij_kernels;
      break;

    default:
      assert(false);
      break;    
  }
  
  switch (mom_id) {
    case M0:
      mom_vm_sr->momt.kernel = m0_kernels[tblidx].kernels[poly_order];
      mom_vm_sr->momt.num_mom = 1;
      break;

    case M1i:
      mom_vm_sr->momt.kernel = m1i_kernels[tblidx].kernels[poly_order];
      mom_vm_sr->momt.num_mom = vdim;
      break;

    case Ni:
      mom_vm_sr->momt.kernel = Ni_kernels[tblidx].kernels[poly_order];
      mom_vm_sr->momt.num_mom = 1+vdim;
      break;

    case Energy:
      mom_vm_sr->momt.kernel = Energy_kernels[tblidx].kernels[poly_order];
      mom_vm_sr->momt.num_mom = 1;
      break;

    case Pressure:
      mom_vm_sr->momt.kernel = Pressure_kernels[tblidx].kernels[poly_order];
      mom_vm_sr->momt.num_mom = 1;
      break;

    case Tij:
      mom_vm_sr->momt.kernel = Tij_kernels[tblidx].kernels[poly_order];
      mom_vm_sr->momt.num_mom = 1+vdim+(vdim*(vdim+1))/2;
      break;

    default: // can't happen
      break;
  }
}

struct gkyl_mom_type*
gkyl_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, 
  const char *mom)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_vlasov_sr *mom_vm_sr = (struct mom_type_vlasov_sr*)
    gkyl_malloc(sizeof(struct mom_type_vlasov_sr));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_vm_sr->momt.cdim = cdim;
  mom_vm_sr->momt.pdim = pdim;
  mom_vm_sr->momt.poly_order = poly_order;
  mom_vm_sr->momt.num_config = cbasis->num_basis;
  mom_vm_sr->momt.num_phase = pbasis->num_basis;

  mom_vm_sr->conf_range = *conf_range;
  mom_vm_sr->vel_range = *vel_range;

  int mom_id = get_mom_id(mom);
  assert(mom_id != BAD);
  mom_vm_sr->momt.num_mom = v_num_mom(vdim, mom_id); // number of moments

  mom_vm_sr->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_vm_sr->momt.flags);
  mom_vm_sr->momt.ref_count = gkyl_ref_count_init(gkyl_mom_vm_sr_free);
  
  // copy struct to device
  struct mom_type_vlasov_sr *momt_cu = (struct mom_type_vlasov_sr*)
    gkyl_cu_malloc(sizeof(struct mom_type_vlasov_sr));
  gkyl_cu_memcpy(momt_cu, mom_vm_sr, sizeof(struct mom_type_vlasov_sr), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);

  set_cu_ptrs<<<1,1>>>(momt_cu, mom_id, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  mom_vm_sr->momt.on_dev = &momt_cu->momt;
  
  return &mom_vm_sr->momt;
}

__global__
static void
set_int_cu_ptrs(struct mom_type_vlasov_sr* mom_vm_sr, enum gkyl_basis_type b_type, int vdim,
  int poly_order, int tblidx)
{
  mom_vm_sr->auxfields.gamma = 0;
  mom_vm_sr->auxfields.gamma_inv = 0;
  mom_vm_sr->auxfields.V_drift = 0;
  mom_vm_sr->auxfields.GammaV2 = 0;

  // choose kernel tables based on basis-function type
  const gkyl_vlasov_sr_mom_kern_list *int_mom_kernels;  
  
  // set kernel pointer
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      int_mom_kernels = ser_int_mom_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  mom_vm_sr->momt.kernel = int_mom_kernels[tblidx].kernels[poly_order];
  mom_vm_sr->momt.num_mom = 2+vdim;
}

struct gkyl_mom_type *
gkyl_int_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct mom_type_vlasov_sr *mom_vm_sr = (struct mom_type_vlasov_sr*)
    gkyl_malloc(sizeof(struct mom_type_vlasov_sr));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_vm_sr->momt.cdim = cdim;
  mom_vm_sr->momt.pdim = pdim;
  mom_vm_sr->momt.poly_order = poly_order;
  mom_vm_sr->momt.num_config = cbasis->num_basis;
  mom_vm_sr->momt.num_phase = pbasis->num_basis;

  mom_vm_sr->momt.num_mom = vdim+2;

  mom_vm_sr->conf_range = *conf_range;
  mom_vm_sr->vel_range = *vel_range;

  mom_vm_sr->momt.flags = 0;
  GKYL_SET_CU_ALLOC(mom_vm_sr->momt.flags);
  mom_vm_sr->momt.ref_count = gkyl_ref_count_init(gkyl_mom_vm_sr_free);
  
  // copy struct to device
  struct mom_type_vlasov_sr *momt_cu = (struct mom_type_vlasov_sr*)
    gkyl_cu_malloc(sizeof(struct mom_type_vlasov_sr));
  gkyl_cu_memcpy(momt_cu, mom_vm_sr, sizeof(struct mom_type_vlasov_sr), GKYL_CU_MEMCPY_H2D);

  set_int_cu_ptrs<<<1,1>>>(momt_cu, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  mom_vm_sr->momt.on_dev = &momt_cu->momt;
  
  return &mom_vm_sr->momt;
}
