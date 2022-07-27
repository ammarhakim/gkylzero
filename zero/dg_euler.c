#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_euler.h>
#include <gkyl_dg_euler_priv.h>
#include <gkyl_util.h>

// "Choose Kernel" based on cdim and polyorder
#define CK(lst,cdim,poly_order) lst[cdim-1].kernels[poly_order]

void 
gkyl_euler_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_eqn *base = container_of(ref, struct gkyl_dg_eqn, ref_count);

  if (gkyl_dg_eqn_is_cu_dev(base)) {
    // free inner on_dev object
    struct dg_euler *euler = container_of(base->on_dev, struct dg_euler, eqn);
    gkyl_cu_free(euler);
  }  
  
  struct dg_euler *euler = container_of(base, struct dg_euler, eqn);
  gkyl_free(euler);
}

void
gkyl_euler_set_auxfields(const struct gkyl_dg_eqn *eqn, struct gkyl_dg_euler_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.u)) {
    gkyl_euler_set_auxfields_cu(eqn->on_dev, auxin);
    return;
  }
#endif

  struct dg_euler *euler = container_of(eqn, struct dg_euler, eqn);
  euler->auxfields.u = auxin.u;
}

struct gkyl_dg_eqn*
gkyl_dg_euler_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range,
  double gas_gamma, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_euler_cu_dev_new(cbasis, conf_range, gas_gamma);
  } 
#endif
  struct dg_euler *euler = gkyl_malloc(sizeof(struct dg_euler));

  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;

  const gkyl_dg_euler_vol_kern_list *vol_kernels;
  const gkyl_dg_euler_surf_kern_list *surf_x_kernels, *surf_y_kernels, *surf_z_kernels;

  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      vol_kernels = ser_vol_kernels;
      surf_x_kernels = ser_surf_x_kernels;
      surf_y_kernels = ser_surf_y_kernels;
      surf_z_kernels = ser_surf_z_kernels;
      break;

    default:
      assert(false);
      break;    
  }  
    
  euler->eqn.num_equations = 5;
  euler->eqn.vol_term = vol;
  euler->eqn.surf_term = surf;
  euler->eqn.boundary_surf_term = boundary_surf;

  euler->gas_gamma = gas_gamma;

  euler->vol =  CK(vol_kernels, cdim, poly_order);
  assert(euler->vol);

  euler->surf[0] = CK(surf_x_kernels, cdim, poly_order);
  if (cdim>1)
    euler->surf[1] = CK(surf_y_kernels, cdim, poly_order);
  if (cdim>2)
    euler->surf[2] = CK(surf_z_kernels, cdim, poly_order);

  // ensure non-NULL pointers 
  for (int i=0; i<cdim; ++i) assert(euler->surf[i]);

  euler->auxfields.u = 0;  
  euler->conf_range = *conf_range;
  
  euler->eqn.flags = 0;
  GKYL_CLEAR_CU_ALLOC(euler->eqn.flags);
  euler->eqn.ref_count = gkyl_ref_count_init(gkyl_euler_free);
  euler->eqn.on_dev = &euler->eqn; // CPU eqn obj points to itself
  
  return &euler->eqn;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_eqn*
gkyl_dg_euler_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_range* conf_range,
  double gas_gamma)
{
  assert(false);
  return 0;
}

#endif
