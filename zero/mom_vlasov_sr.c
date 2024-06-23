#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_mom_vlasov_sr.h>
#include <gkyl_mom_vlasov_sr_priv.h>
#include <gkyl_util.h>

void
gkyl_mom_vm_sr_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_mom_type *momt = container_of(ref, struct gkyl_mom_type, ref_count);
  if (GKYL_IS_CU_ALLOC(momt->flags))
    gkyl_cu_free(momt->on_dev);
  gkyl_free(momt);
}

void
gkyl_mom_vlasov_sr_set_auxfields(const struct gkyl_mom_type *momt, struct gkyl_mom_vlasov_sr_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_mom_type_is_cu_dev(momt)) {
    gkyl_mom_vlasov_sr_set_auxfields_cu(momt->on_dev, auxin);
    return;
  }
#endif

  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);
  mom_vm_sr->auxfields.gamma = auxin.gamma;
  mom_vm_sr->auxfields.gamma_inv = auxin.gamma_inv;
  mom_vm_sr->auxfields.V_drift = auxin.V_drift;
  mom_vm_sr->auxfields.GammaV2 = auxin.GammaV2;
}

struct gkyl_mom_type*
gkyl_mom_vlasov_sr_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, 
  const char *mom, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_mom_vlasov_sr_cu_dev_new(cbasis, pbasis, conf_range, vel_range, mom);
  } 
#endif  
  struct mom_type_vlasov_sr *mom_vm_sr = gkyl_malloc(sizeof(struct mom_type_vlasov_sr));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_vm_sr->momt.cdim = cdim;
  mom_vm_sr->momt.pdim = pdim;
  mom_vm_sr->momt.poly_order = poly_order;
  mom_vm_sr->momt.num_config = cbasis->num_basis;
  mom_vm_sr->momt.num_phase = pbasis->num_basis;

  // choose kernel tables based on basis-function type
  const gkyl_vlasov_sr_mom_kern_list *m0_kernels, *m1i_kernels, 
    *Ni_kernels, *Energy_kernels, *Pressure_kernels, *Tij_kernels;

  switch (cbasis->b_type) {
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

  if (strcmp(mom, "M0") == 0) { // density (GammaV*n)
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_vm_sr->momt.kernel = m0_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_vm_sr->momt.num_mom = 1;
  }
  else if (strcmp(mom, "M1i") == 0) { // momentum (GammaV*n*V)
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != m1i_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_vm_sr->momt.kernel = m1i_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_vm_sr->momt.num_mom = vdim;
  }
  else if (strcmp(mom, "Ni") == 0) { // 4-momentum (GammaV*n, GammaV*n*V)
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != Ni_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_vm_sr->momt.kernel = Ni_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_vm_sr->momt.num_mom = 1+vdim;
  }
  else if (strcmp(mom, "Energy") == 0) { // total energy = gamma*mc^2 moment
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != Energy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_vm_sr->momt.kernel = Energy_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_vm_sr->momt.num_mom = 1;
  }
  else if (strcmp(mom, "Pressure") == 0) { // total fluid-frame pressure = n*T where n is the fluid-frame density
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != Pressure_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_vm_sr->momt.kernel = Pressure_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_vm_sr->momt.num_mom = 1;
  }
  else if (strcmp(mom, "Tij") == 0) { // Stress-energy tensor (Energy, Energy flux (vdim components), Stress tensor (vdim*(vdim+1))/2 components))
    assert(cv_index[cdim].vdim[vdim] != -1);
    assert(NULL != Tij_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order]);
    
    mom_vm_sr->momt.kernel = Tij_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    mom_vm_sr->momt.num_mom = 1+vdim+(vdim*(vdim+1))/2;
  }
  else {
    // string not recognized
    gkyl_exit("gkyl_mom_type_vlasov_sr: Unrecognized moment requested!");
  }

  mom_vm_sr->conf_range = *conf_range;
  mom_vm_sr->vel_range = *vel_range;

  mom_vm_sr->auxfields.gamma = 0;
  mom_vm_sr->auxfields.gamma_inv = 0;
  mom_vm_sr->auxfields.V_drift = 0;
  mom_vm_sr->auxfields.GammaV2 = 0;

  mom_vm_sr->momt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mom_vm_sr->momt.flags);
  mom_vm_sr->momt.ref_count = gkyl_ref_count_init(gkyl_mom_vm_sr_free);
  
  mom_vm_sr->momt.on_dev = &mom_vm_sr->momt; // on host, self-reference
    
  return &mom_vm_sr->momt;
}

struct gkyl_mom_type*
gkyl_int_mom_vlasov_sr_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_int_mom_vlasov_sr_cu_dev_new(cbasis, pbasis, conf_range, vel_range);
  } 
#endif  
  struct mom_type_vlasov_sr *mom_vm_sr = gkyl_malloc(sizeof(struct mom_type_vlasov_sr));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  mom_vm_sr->momt.cdim = cdim;
  mom_vm_sr->momt.pdim = pdim;
  mom_vm_sr->momt.poly_order = poly_order;
  mom_vm_sr->momt.num_config = cbasis->num_basis;
  mom_vm_sr->momt.num_phase = pbasis->num_basis;

  // choose kernel tables based on basis-function type
  const gkyl_vlasov_sr_mom_kern_list *int_mom_kernels;  
  
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
  mom_vm_sr->momt.kernel = int_mom_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
  mom_vm_sr->momt.num_mom = 2+vdim;

  mom_vm_sr->conf_range = *conf_range;
  mom_vm_sr->vel_range = *vel_range;

  mom_vm_sr->auxfields.gamma = 0;
  mom_vm_sr->auxfields.gamma_inv = 0;
  mom_vm_sr->auxfields.V_drift = 0;
  mom_vm_sr->auxfields.GammaV2 = 0;

  mom_vm_sr->momt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(mom_vm_sr->momt.flags);
  mom_vm_sr->momt.ref_count = gkyl_ref_count_init(gkyl_mom_vm_sr_free);
  
  mom_vm_sr->momt.on_dev = &mom_vm_sr->momt; // on host, self-reference
    
  return &mom_vm_sr->momt;  
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_mom_type*
gkyl_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range, const char *mom)
{
  assert(false);
  return 0;
}

struct gkyl_mom_type *
gkyl_int_mom_vlasov_sr_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const struct gkyl_range* vel_range)
{
  assert(false);
  return 0;
}

#endif
