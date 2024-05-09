/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_rad_gyrokinetic_drag.h>    
#include <gkyl_dg_rad_gyrokinetic_drag_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_rad_gyrokinetic_drag_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, 
  const struct gkyl_array* nvnu_surf, const struct gkyl_array* nvnu, 
  const struct gkyl_array* nvsqnu_surf, const struct gkyl_array* nvsqnu, const struct gkyl_array* vtsq,
  const double vtsq_min)
{
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = container_of(eqn, struct dg_rad_gyrokinetic_drag, eqn);
  rad_gyrokinetic_drag->auxfields.nvnu_surf = nvnu_surf;
  rad_gyrokinetic_drag->auxfields.nvnu = nvnu;
  rad_gyrokinetic_drag->auxfields.nvsqnu_surf = nvsqnu_surf;
  rad_gyrokinetic_drag->auxfields.nvsqnu = nvsqnu;
  rad_gyrokinetic_drag->auxfields.vtsq = vtsq;
  rad_gyrokinetic_drag->auxfields.vtsq_min = vtsq_min;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_rad_gyrokinetic_drag_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_rad_gyrokinetic_auxfields auxin)
{
  gkyl_rad_gyrokinetic_drag_set_auxfields_cu_kernel<<<1,1>>>(eqn, 
    auxin.nvnu_surf->on_dev, auxin.nvnu->on_dev, 
    auxin.nvsqnu_surf->on_dev, auxin.nvsqnu->on_dev, auxin.vtsq->on_dev, auxin.vtsq_min);
}

// CUDA kernel to set device pointers to range object and rad_gyrokinetic_drag kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_rad_gyrokinetic_drag_set_cu_dev_ptrs(struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag, enum gkyl_basis_type b_type,
  int cv_index, int cdim, int vdim, int poly_order)
{
  rad_gyrokinetic_drag->auxfields.nvnu_surf = 0; 
  rad_gyrokinetic_drag->auxfields.nvnu = 0; 
  rad_gyrokinetic_drag->auxfields.nvsqnu_surf = 0; 
  rad_gyrokinetic_drag->auxfields.nvsqnu = 0;
  rad_gyrokinetic_drag->auxfields.vtsq = 0; 

  rad_gyrokinetic_drag->eqn.surf_term = surf;
  rad_gyrokinetic_drag->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_rad_gyrokinetic_vol_kern_list *vol_kernels;
  const gkyl_dg_rad_gyrokinetic_surf_kern_list *surf_vpar_kernels, *surf_mu_kernels;
  const gkyl_dg_rad_gyrokinetic_boundary_surf_kern_list *boundary_surf_vpar_kernels, *boundary_surf_mu_kernels;
  
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_vpar_kernels = ser_surf_vpar_kernels;
      surf_mu_kernels = ser_surf_mu_kernels;
      boundary_surf_vpar_kernels = ser_boundary_surf_vpar_kernels;
      boundary_surf_mu_kernels = ser_boundary_surf_mu_kernels;
      break;

    default:
      assert(false);
      break;    
  }  
  rad_gyrokinetic_drag->eqn.vol_term = vol_kernels[cv_index].kernels[poly_order];

  rad_gyrokinetic_drag->surf[0] = surf_vpar_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    rad_gyrokinetic_drag->surf[1] = surf_mu_kernels[cv_index].kernels[poly_order];

  rad_gyrokinetic_drag->boundary_surf[0] = boundary_surf_vpar_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    rad_gyrokinetic_drag->boundary_surf[1] = boundary_surf_mu_kernels[cv_index].kernels[poly_order];
}

struct gkyl_dg_eqn*
gkyl_dg_rad_gyrokinetic_drag_cu_dev_new(const struct gkyl_basis* conf_basis, 
  const struct gkyl_basis* phase_basis, const struct gkyl_range *phase_range,
  const struct gkyl_range *conf_range, const struct gkyl_velocity_map *vel_map)
{
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag = (struct dg_rad_gyrokinetic_drag*) gkyl_malloc(sizeof(struct dg_rad_gyrokinetic_drag));

  int cdim = conf_basis->ndim, pdim = phase_basis->ndim, vdim = pdim-cdim;
  int poly_order = conf_basis->poly_order;
  
  rad_gyrokinetic_drag->cdim = cdim;
  rad_gyrokinetic_drag->pdim = pdim;
  rad_gyrokinetic_drag->phase_range = *phase_range;
  rad_gyrokinetic_drag->conf_range = *conf_range;

  // Acquire pointers to on_dev objects so memcpy below copies those too.
  struct gkyl_velocity_map *vel_map_ho = gkyl_velocity_map_acquire(vel_map);
  gyrokinetic->vel_map = vel_map_ho->on_dev;

  rad_gyrokinetic_drag->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(rad_gyrokinetic_drag->eqn.flags);
  rad_gyrokinetic_drag->eqn.ref_count = gkyl_ref_count_init(gkyl_rad_gyrokinetic_free);

  // copy the host struct to device struct
  struct dg_rad_gyrokinetic_drag *rad_gyrokinetic_drag_cu = (struct dg_rad_gyrokinetic_drag*) gkyl_cu_malloc(sizeof(struct dg_rad_gyrokinetic_drag));
  gkyl_cu_memcpy(rad_gyrokinetic_drag_cu, rad_gyrokinetic_drag, sizeof(struct dg_rad_gyrokinetic_drag), GKYL_CU_MEMCPY_H2D);

  dg_rad_gyrokinetic_drag_set_cu_dev_ptrs<<<1,1>>>(rad_gyrokinetic_drag_cu, conf_basis->b_type, cv_index[cdim].vdim[vdim],
    cdim, vdim, poly_order);

  // set parent on_dev pointer
  rad_gyrokinetic_drag->eqn.on_dev = &rad_gyrokinetic_drag_cu->eqn;

  // Updater should store host pointers.
  gyrokinetic->vel_map = vel_map_ho; 

  return &rad_gyrokinetic_drag->eqn;
}
