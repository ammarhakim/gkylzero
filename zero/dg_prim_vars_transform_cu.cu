/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_prim_vars_transform.h>
#include <gkyl_dg_prim_vars_transform_priv.h>  
#include <gkyl_util.h>
}

enum { u_par_i, u_par, prim_gk, prim_vlasov, BAD };

static int
get_prim_id(const char *prim_nm)
{
  int prim_idx = BAD;

  if (strcmp(prim_nm, "u_par_i") == 0) { // u_par b_i
    prim_idx = u_par_i;
  }
  else if (strcmp(prim_nm, "u_par") == 0) { // u_i . b_i
    prim_idx = u_par;
  }
  else if (strcmp(prim_nm, "prim_gk") == 0) { // combined (u_par = u_i . b_i, vth_GK^2 = (M2/M0 - upar^2))
    prim_idx = prim;
  }
  else if (strcmp(prim_nm, "prim_vlasov") == 0) { // combined (u_par_ and vth^2)
    prim_idx = prim;
  }
  else {
    prim_idx = BAD;
  }    

  return prim_idx;
}

static int
v_num_prim(int vdim, int prim_id)
{
  int num_prim = 0;
  
  switch (prim_id) {
    case u_par_i:
      num_prim = vdim;
      break;

    case u_par:
      num_prim = 1;
      break;

    case prim_gk:
      num_prim = 2;
      break;

    case prim_vlasov:
      num_prim = 3;
      break; 
      
    default: // can't happen
      break;
  }

  return num_prim;
}

// CUDA kernel to set pointer to auxiliary fields.
// This is required because prim vars object lives on device,
// and so its members cannot be modified without a full __global__ kernel on device.
__global__ static void
gkyl_dg_prim_vars_transform_set_auxfields_cu_kernel(const struct gkyl_dg_prim_vars_type *pvt, 
  const struct gkyl_array *b_i)
{
  struct dg_prim_vars_type_transform *prim = container_of(pvt, struct dg_prim_vars_type_transform, pvt);
  prim->auxfields.b_i = b_i; // covariant components of the field aligned unit vector.
}

// Host-side wrapper for set_auxfields_cu_kernel
void
gkyl_dg_prim_vars_transform_set_auxfields_cu(const struct gkyl_dg_prim_vars_type *pvt, 
  struct gkyl_dg_prim_vars_auxfields auxin)
{
  gkyl_dg_prim_vars_transform_set_auxfields_cu_kernel<<<1,1>>>(pvt, auxin.b_i->on_dev);
}

__global__
static void
set_cu_ptrs(struct dg_prim_vars_type_transform* pvt, int prim_id, enum gkyl_basis_type b_type, int vdim,
  int poly_order, int tblidx)
{
  pvt->auxfields.b_i = 0;

  const gkyl_dg_prim_vars_transform_kern_list *dg_prim_vars_transform_u_par_i_kernels, 
    *dg_prim_vars_transform_u_par_kernels, *dg_prim_vars_transform_kernels;

  // choose kernel tables based on basis-function type
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      dg_prim_vars_transform_u_par_i_kernels = ser_dg_prim_vars_transform_u_par_i_kernels;
      dg_prim_vars_transform_u_par_kernels = ser_dg_prim_vars_transform_u_par_kernels;
      dg_prim_vars_transform_gk_kernels = ser_dg_prim_vars_transform_gk_kernels;
      dg_prim_vars_transform_vlasov_kernels = ser_dg_prim_vars_transform_vlasov_kernels;
      break;

    default:
      assert(false);
      break;    
  }
  
  switch (prim_id) {
    case u_par_i:
      pvt->pvt.kernel = dg_prim_vars_transform_u_par_i_kernels[tblidx].kernels[poly_order];
      break;

    case u_par:
      pvt->pvt.kernel = dg_prim_vars_transform_u_par_kernels[tblidx].kernels[poly_order];
      break;

    case prim_gk:
      pvt->pvt.kernel = dg_prim_vars_transform_gk_kernels[tblidx].kernels[poly_order];
      break;

    case prim_vlasov:
      pvt->pvt.kernel = dg_prim_vars_transform_vlasov_kernels[tblidx].kernels[poly_order];
      break;
      
    default: // can't happen
      break;
  }
}

struct gkyl_dg_prim_vars_type*
gkyl_dg_prim_vars_transform_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const struct gkyl_range* conf_range, const char *prim_nm)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct dg_prim_vars_type_transform_vlasov_gk *pvt = (struct dg_prim_vars_type_transform_vlasov_gk*)
    gkyl_malloc(sizeof(struct dg_prim_vars_type_transform_vlasov_gk));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  pvt->pvt.cdim = cdim;
  pvt->pvt.vdim = vdim;
  pvt->pvt.poly_order = poly_order;
  pvt->pvt.num_config = cbasis->num_basis;

  pvt->conf_range = *conf_range;

  int prim_id = get_prim_id(prim_nm);
  assert(prim_id != BAD);
  pvt->pvt.num_mom = v_num_prim(vdim, prim_id); // number of primitive variables

  pvt->pvt.flags = 0;
  GKYL_SET_CU_ALLOC(pvt->pvt.flags);
  pvt->pvt.ref_count = gkyl_ref_count_init(gkyl_dg_prim_vars_transform_free);  
  
  // copy struct to device
  struct dg_prim_vars_type_transform *pvt_cu = (struct dg_prim_vars_type_transform*)
    gkyl_cu_malloc(sizeof(struct dg_prim_vars_type_transform));
  gkyl_cu_memcpy(pvt_cu, pvt, sizeof(struct dg_prim_vars_type_transform), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);

  if ((strcmp(prim_nm, "u_par_i") == 0) || (strcmp(prim_nm, "prim_vlasov") == 0)) {
    set_cu_ptrs<<<1,1>>>(pvt_cu, prim_id, cbasis->b_type,
      vdim, poly_order, cv_vlasov_index[cdim].vdim[vdim]);
  }
  else if ((strcmp(prim_nm, "u_par") == 0) || (strcmp(prim_nm, "prim_gk") == 0)) {
    set_cu_ptrs<<<1,1>>>(pvt_cu, prim_id, cbasis->b_type,
      vdim, poly_order, cv_gk_index[cdim].vdim[vdim]);
  }
  pvt->pvt.on_dev = &pvt_cu->pvt;
  
  return &pvt->pvt;
}
