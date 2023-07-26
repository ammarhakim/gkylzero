/* -*- c++ -*- */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_prim_vars_gyrokinetic.h>
#include <gkyl_dg_prim_vars_gyrokinetic_priv.h>
#include <gkyl_util.h>
}

enum { upar, vtSq, prim, BAD };

static int
get_prim_id(const char *prim_nm)
{
  int prim_idx = BAD;

  if (strcmp(prim_nm, "upar") == 0) { // flow velocity
    prim_idx = upar;
  }
  else if (strcmp(prim_nm, "vtSq") == 0) { // vth^2
    prim_idx = vtSq;
  }
  else if (strcmp(prim_nm, "prim") == 0) { // combined (upar, vth^2)
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
    case upar:
      num_prim = 1;
      break;

    case vtSq:
      num_prim = 1;
      break;

    case prim:
      num_prim = 2;
      break;    
      
    default: // can't happen
      break;
  }

  return num_prim;
}

__global__
static void
set_cu_ptrs(struct dg_prim_vars_type_gyrokinetic* pvt, int prim_id, enum gkyl_basis_type b_type, int vdim,
  int poly_order, int tblidx)
{
  const gkyl_dg_prim_vars_gyrokinetic_kern_list *dg_prim_vars_gyrokinetic_upar_kernels, 
    *dg_prim_vars_gyrokinetic_vtSq_kernels, *dg_prim_vars_gyrokinetic_kernels;

  // choose kernel tables based on basis-function type
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      dg_prim_vars_gyrokinetic_upar_kernels = ser_dg_prim_vars_gyrokinetic_upar_kernels;
      dg_prim_vars_gyrokinetic_vtSq_kernels = ser_dg_prim_vars_gyrokinetic_vtSq_kernels;
      dg_prim_vars_gyrokinetic_kernels = ser_dg_prim_vars_gyrokinetic_kernels;
      break;

    default:
      assert(false);
      break;    
  }
  
  switch (prim_id) {
    case upar:
      pvt->pvt.kernel = dg_prim_vars_gyrokinetic_upar_kernels[tblidx].kernels[poly_order];
      break;

    case vtSq:
      pvt->pvt.kernel = dg_prim_vars_gyrokinetic_vtSq_kernels[tblidx].kernels[poly_order];
      break;

    case prim:
      pvt->pvt.kernel = dg_prim_vars_gyrokinetic_kernels[tblidx].kernels[poly_order];
      break;
      
    default: // can't happen
      break;
  }
}

struct gkyl_dg_prim_vars_type*
gkyl_dg_prim_vars_gyrokinetic_cu_dev_new(const struct gkyl_basis* cbasis,
  const struct gkyl_basis* pbasis, const char *prim_nm)
{
  assert(cbasis->poly_order == pbasis->poly_order);

  struct dg_prim_vars_type_gyrokinetic *pvt = (struct dg_prim_vars_type_gyrokinetic*)
    gkyl_malloc(sizeof(struct dg_prim_vars_type_gyrokinetic));
  
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  pvt->pvt.cdim = cdim;
  pvt->pvt.vdim = vdim;
  pvt->pvt.poly_order = poly_order;
  pvt->pvt.num_config = cbasis->num_basis;

  int prim_id = get_prim_id(prim_nm);
  assert(prim_id != BAD);
  pvt->pvt.num_mom = v_num_prim(vdim, prim_id); // number of primitive variables

  pvt->pvt.flags = 0;
  GKYL_SET_CU_ALLOC(pvt->pvt.flags);
  pvt->pvt.ref_count = gkyl_ref_count_init(gkyl_dg_prim_vars_gyrokinetic_free);
  
  // copy struct to device
  struct dg_prim_vars_type_gyrokinetic *pvt_cu = (struct dg_prim_vars_type_gyrokinetic*)
    gkyl_cu_malloc(sizeof(struct dg_prim_vars_type_gyrokinetic));
  gkyl_cu_memcpy(pvt_cu, pvt, sizeof(struct dg_prim_vars_type_gyrokinetic), GKYL_CU_MEMCPY_H2D);

  assert(cv_index[cdim].vdim[vdim] != -1);

  set_cu_ptrs<<<1,1>>>(pvt_cu, prim_id, cbasis->b_type,
    vdim, poly_order, cv_index[cdim].vdim[vdim]);

  pvt->pvt.on_dev = &pvt_cu->pvt;
  
  return &pvt->pvt;
}
