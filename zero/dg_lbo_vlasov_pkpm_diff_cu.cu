/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_lbo_vlasov_pkpm_diff.h>    
#include <gkyl_dg_lbo_vlasov_pkpm_diff_priv.h>
}

#include <cassert>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

// CUDA kernel to set pointer to nu*vt^2 (collision frequency * thermal velocity squared)
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_lbo_vlasov_pkpm_diff_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, const struct gkyl_array *nu, const struct gkyl_array *nuVtSq)
{
  struct dg_lbo_vlasov_pkpm_diff *lbo_vlasov_pkpm_diff = container_of(eqn, struct dg_lbo_vlasov_pkpm_diff, eqn);
  lbo_vlasov_pkpm_diff->auxfields.nu = nu;
  lbo_vlasov_pkpm_diff->auxfields.nuVtSq = nuVtSq;
}

//// Host-side wrapper for device kernels setting nuVtSq.
void
gkyl_lbo_vlasov_pkpm_diff_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_vlasov_pkpm_diff_auxfields auxin)
{
  gkyl_lbo_vlasov_pkpm_diff_set_auxfields_cu_kernel<<<1,1>>>(eqn, auxin.nu->on_dev, auxin.nuVtSq->on_dev);
}

// CUDA kernel to set device pointers to range object and Vlasov PKPM LBO diffusion kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_lbo_vlasov_pkpm_diff_set_cu_dev_ptrs(struct dg_lbo_vlasov_pkpm_diff *lbo_vlasov_pkpm_diff, enum gkyl_basis_type b_type,
  int cdim, int poly_order)
{
  lbo_vlasov_pkpm_diff->auxfields.nu = 0; 
  lbo_vlasov_pkpm_diff->auxfields.nuVtSq = 0; 

  lbo_vlasov_pkpm_diff->eqn.vol_term = vol;
  lbo_vlasov_pkpm_diff->eqn.surf_term = surf;
  lbo_vlasov_pkpm_diff->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_lbo_vlasov_pkpm_diff_vol_kern_list *vol_kernels;
  const gkyl_dg_lbo_vlasov_pkpm_diff_surf_kern_list *surf_vpar_kernels;
  const gkyl_dg_lbo_vlasov_pkpm_diff_boundary_surf_kern_list *boundary_surf_vpar_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_vpar_kernels = ser_surf_vpar_kernels;
      boundary_surf_vpar_kernels = ser_boundary_surf_vpar_kernels;
      
      break;

    default:
      assert(false);
      break;    
  }  
 
  lbo_vlasov_pkpm_diff->vol = CK(vol_kernels, cdim, poly_order);

  lbo_vlasov_pkpm_diff->surf = CK(surf_vpar_kernels, cdim, poly_order);

  lbo_vlasov_pkpm_diff->boundary_surf = CK(boundary_surf_vpar_kernels, cdim, poly_order);
}

struct gkyl_dg_eqn*
gkyl_dg_lbo_vlasov_pkpm_diff_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, const struct gkyl_rect_grid *pgrid)
{
  struct dg_lbo_vlasov_pkpm_diff *lbo_vlasov_pkpm_diff =
    (struct dg_lbo_vlasov_pkpm_diff*) gkyl_malloc(sizeof(struct dg_lbo_vlasov_pkpm_diff));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  lbo_vlasov_pkpm_diff->cdim = cdim;
  lbo_vlasov_pkpm_diff->pdim = pdim;

  lbo_vlasov_pkpm_diff->eqn.num_equations = 1;
  lbo_vlasov_pkpm_diff->conf_range = *conf_range;
  lbo_vlasov_pkpm_diff->vMaxSq = pow(pgrid->upper[cdim],2);

  lbo_vlasov_pkpm_diff->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(lbo_vlasov_pkpm_diff->eqn.flags);
  lbo_vlasov_pkpm_diff->eqn.ref_count = gkyl_ref_count_init(gkyl_lbo_vlasov_pkpm_diff_free);

  // copy the host struct to device struct
  struct dg_lbo_vlasov_pkpm_diff *lbo_vlasov_pkpm_diff_cu =
    (struct dg_lbo_vlasov_pkpm_diff*) gkyl_cu_malloc(sizeof(struct dg_lbo_vlasov_pkpm_diff));

  gkyl_cu_memcpy(lbo_vlasov_pkpm_diff_cu, lbo_vlasov_pkpm_diff,
    sizeof(struct dg_lbo_vlasov_pkpm_diff), GKYL_CU_MEMCPY_H2D);

  dg_lbo_vlasov_pkpm_diff_set_cu_dev_ptrs<<<1,1>>>(lbo_vlasov_pkpm_diff_cu,
    cbasis->b_type, cdim, poly_order);

  lbo_vlasov_pkpm_diff->eqn.on_dev = &lbo_vlasov_pkpm_diff_cu->eqn;  
  
  return &lbo_vlasov_pkpm_diff->eqn;
}
