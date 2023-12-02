/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_lbo_gyrokinetic_drag.h>    
#include <gkyl_dg_lbo_gyrokinetic_drag_priv.h>
}

#include <cassert>

// CUDA kernel to set pointer to nuSum, sum of collisionalities
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_lbo_gyrokinetic_drag_set_auxfields_cu_kernel(const struct gkyl_dg_eqn *eqn, 
  const struct gkyl_array *nuSum, const struct gkyl_array *nuPrimMomsSum, const struct gkyl_array *m2self)
{
  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag = container_of(eqn, struct dg_lbo_gyrokinetic_drag, eqn);
  lbo_gyrokinetic_drag->auxfields.nuSum = nuSum;
  lbo_gyrokinetic_drag->auxfields.nuPrimMomsSum = nuPrimMomsSum;
  lbo_gyrokinetic_drag->auxfields.m2self = m2self;
}

//// Host-side wrapper for device kernels setting nuSum, nuUSum and nuVtSqSum.
void
gkyl_lbo_gyrokinetic_drag_set_auxfields_cu(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_lbo_gyrokinetic_drag_auxfields auxin)
{
  gkyl_lbo_gyrokinetic_drag_set_auxfields_cu_kernel<<<1,1>>>(eqn, 
    auxin.nuSum->on_dev, auxin.nuPrimMomsSum->on_dev, auxin.m2self->on_dev);
}

// CUDA kernel to set device pointers to range object and gyrokinetic LBO kernel function
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_lbo_gyrokinetic_drag_set_cu_dev_ptrs(struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag, enum gkyl_basis_type b_type,
  int cv_index, int cdim, int vdim, int poly_order)
{
  lbo_gyrokinetic_drag->auxfields.nuSum = 0; 
  lbo_gyrokinetic_drag->auxfields.nuPrimMomsSum = 0; 
  lbo_gyrokinetic_drag->auxfields.m2self = 0; 

  lbo_gyrokinetic_drag->eqn.surf_term = surf;
  lbo_gyrokinetic_drag->eqn.boundary_surf_term = boundary_surf;

  const gkyl_dg_lbo_gyrokinetic_drag_vol_kern_list *vol_kernels;
  const gkyl_dg_lbo_gyrokinetic_drag_surf_kern_list *surf_vpar_kernels, *surf_mu_kernels;
  const gkyl_dg_lbo_gyrokinetic_drag_boundary_surf_kern_list *boundary_surf_vpar_kernels, *boundary_surf_mu_kernels;
  
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

  lbo_gyrokinetic_drag->eqn.vol_term = vol_kernels[cv_index].kernels[poly_order];

  lbo_gyrokinetic_drag->surf[0] = surf_vpar_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    lbo_gyrokinetic_drag->surf[1] = surf_mu_kernels[cv_index].kernels[poly_order];

  lbo_gyrokinetic_drag->boundary_surf[0] = boundary_surf_vpar_kernels[cv_index].kernels[poly_order];
  if (vdim>1)
    lbo_gyrokinetic_drag->boundary_surf[1] = boundary_surf_mu_kernels[cv_index].kernels[poly_order];

}

struct gkyl_dg_eqn*
gkyl_dg_lbo_gyrokinetic_drag_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, const struct gkyl_rect_grid *pgrid, double mass, const struct gk_geometry *gk_geom)
{
  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag =
    (struct dg_lbo_gyrokinetic_drag*) gkyl_malloc(sizeof(struct dg_lbo_gyrokinetic_drag));

  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  lbo_gyrokinetic_drag->cdim = cdim;
  lbo_gyrokinetic_drag->pdim = pdim;

  lbo_gyrokinetic_drag->eqn.num_equations = 1;
  lbo_gyrokinetic_drag->mass = mass;
  // acquire pointer to geometry object
  struct gk_geometry *geom = gkyl_gk_geometry_acquire(gk_geom);
  lbo_gyrokinetic_drag->gk_geom = geom->on_dev; // this is so the memcpy below has geometry on_dev
  lbo_gyrokinetic_drag->conf_range = *conf_range;

  lbo_gyrokinetic_drag->vparMax = pgrid->upper[cdim];
  lbo_gyrokinetic_drag->vparMaxSq = pow(pgrid->upper[cdim],2);
  lbo_gyrokinetic_drag->num_cbasis = cbasis->num_basis;

  lbo_gyrokinetic_drag->eqn.flags = 0;
  GKYL_SET_CU_ALLOC(lbo_gyrokinetic_drag->eqn.flags);
  lbo_gyrokinetic_drag->eqn.ref_count = gkyl_ref_count_init(gkyl_lbo_gyrokinetic_drag_free);

  // copy the host struct to device struct
  struct dg_lbo_gyrokinetic_drag *lbo_gyrokinetic_drag_cu =
    (struct dg_lbo_gyrokinetic_drag*) gkyl_cu_malloc(sizeof(struct dg_lbo_gyrokinetic_drag));

  gkyl_cu_memcpy(lbo_gyrokinetic_drag_cu, lbo_gyrokinetic_drag,
    sizeof(struct dg_lbo_gyrokinetic_drag), GKYL_CU_MEMCPY_H2D);

  dg_lbo_gyrokinetic_drag_set_cu_dev_ptrs<<<1,1>>>(lbo_gyrokinetic_drag_cu,
    cbasis->b_type, cv_index[cdim].vdim[vdim], cdim, vdim, poly_order);

  lbo_gyrokinetic_drag->eqn.on_dev = &lbo_gyrokinetic_drag_cu->eqn;  
  // updater should store host pointers
  lbo_gyrokinetic_diff->gk_geom = geom;
  
  return &lbo_gyrokinetic_drag->eqn;
}
