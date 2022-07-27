/* -*- c++ -*- */

#include <assert.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_util.h>
#include <gkyl_prim_lbo_vlasov_with_fluid.h>
#include <gkyl_prim_lbo_vlasov_with_fluid_priv.h>
}

// CUDA kernel to set pointer to auxiliary fields.
// This is required because eqn object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_prim_lbo_vlasov_with_fluid_set_auxfields_cu_kernel(const struct gkyl_prim_lbo_type *prim, const struct gkyl_array *fluid)
{
  struct prim_lbo_type_vlasov_with_fluid *prim_vlasov_with_fluid = container_of(prim, struct prim_lbo_type_vlasov_with_fluid, prim);
  prim_vlasov_with_fluid->auxfields.fluid = fluid;
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_prim_lbo_vlasov_with_fluid_set_auxfields_cu(const struct gkyl_prim_lbo_type *prim,
  struct gkyl_prim_lbo_vlasov_with_fluid_auxfields auxin)
{
  gkyl_prim_lbo_vlasov_with_fluid_set_auxfields_cu_kernel<<<1,1>>>(prim, auxin.fluid->on_dev);
}

__global__ static void
gkyl_prim_lbo_vlasov_with_fluid_set_cu_dev_ptrs(struct prim_lbo_type_vlasov_with_fluid *prim_vlasov_with_fluid, int cdim, int vdim, int poly_order, enum gkyl_basis_type b_type, int tblidx)
{
  prim_vlasov_with_fluid->prim.self_prim = self_prim;
  prim_vlasov_with_fluid->prim.cross_prim = cross_prim;
  
  // choose kernel tables based on basis-function type
  const gkyl_prim_lbo_vlasov_with_fluid_self_kern_list *self_prim_kernels;
  const gkyl_prim_lbo_vlasov_with_fluid_cross_kern_list *cross_prim_kernels;

  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      self_prim_kernels = ser_self_prim_kernels;
      cross_prim_kernels = ser_cross_prim_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  prim_vlasov_with_fluid->self_prim = self_prim_kernels[tblidx].kernels[poly_order];
  prim_vlasov_with_fluid->cross_prim = cross_prim_kernels[tblidx].kernels[poly_order];
}

struct gkyl_prim_lbo_type*
gkyl_prim_lbo_vlasov_with_fluid_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range)
{
  assert(cbasis->poly_order == pbasis->poly_order);
  
  struct prim_lbo_type_vlasov_with_fluid *prim_vlasov_with_fluid =
    (struct prim_lbo_type_vlasov_with_fluid*) gkyl_malloc(sizeof(struct prim_lbo_type_vlasov_with_fluid));
  
  int cdim = prim_vlasov_with_fluid->prim.cdim = cbasis->ndim;
  int pdim = prim_vlasov_with_fluid->prim.pdim = pbasis->ndim;
  int vdim = pdim - cdim;
  int poly_order = prim_vlasov_with_fluid->prim.poly_order = cbasis->poly_order;
  prim_vlasov_with_fluid->prim.udim = vdim;

  prim_vlasov_with_fluid->conf_range = *conf_range;

  prim_vlasov_with_fluid->prim.flag = 0;
  GKYL_SET_CU_ALLOC(prim_vlasov_with_fluid->prim.flag);
  prim_vlasov_with_fluid->prim.ref_count = gkyl_ref_count_init(prim_lbo_vlasov_with_fluid_free);
  
  // copy the host struct to device struct
  struct prim_lbo_type_vlasov_with_fluid *prim_vlasov_with_fluid_cu = (struct prim_lbo_type_vlasov_with_fluid*)
    gkyl_cu_malloc(sizeof(struct prim_lbo_type_vlasov_with_fluid));
  gkyl_cu_memcpy(prim_vlasov_with_fluid_cu, prim_vlasov_with_fluid, sizeof(struct prim_lbo_type_vlasov_with_fluid), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);
  
  gkyl_prim_lbo_vlasov_with_fluid_set_cu_dev_ptrs<<<1,1>>>(prim_vlasov_with_fluid_cu, cdim, vdim, poly_order,
    cbasis->b_type, cv_index[cdim].vdim[vdim]);

  prim_vlasov_with_fluid->prim.on_dev = &prim_vlasov_with_fluid_cu->prim;
    
  return &prim_vlasov_with_fluid->prim;
}
