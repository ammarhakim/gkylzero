#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_prim_vars_vlasov.h>
#include <gkyl_dg_prim_vars_vlasov_priv.h>
#include <gkyl_util.h>

void
gkyl_dg_prim_vars_vlasov_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_prim_vars_type *pvt = container_of(ref, struct gkyl_dg_prim_vars_type, ref_count);
  if (GKYL_IS_CU_ALLOC(pvt->flags))
    gkyl_cu_free(pvt->on_dev);
  gkyl_free(pvt);
}

struct gkyl_dg_prim_vars_type*
gkyl_dg_prim_vars_vlasov_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const char *prim_nm, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_prim_vars_vlasov_cu_dev_new(cbasis, pbasis);
  } 
#endif    
  struct dg_prim_vars_type_vlasov *dg_prim_vars_vlasov = gkyl_malloc(sizeof(struct dg_prim_vars_type_vlasov));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  dg_prim_vars_vlasov->pvt.cdim = cdim;
  dg_prim_vars_vlasov->pvt.vdim = vdim;
  dg_prim_vars_vlasov->pvt.poly_order = poly_order;
  dg_prim_vars_vlasov->pvt.num_config = cbasis->num_basis;

  const gkyl_dg_prim_vars_vlasov_kern_list *dg_prim_vars_vlasov_u_i_kernels, 
    *dg_prim_vars_vlasov_vth2_kernels, *dg_prim_vars_vlasov_kernels;

  // choose kernel tables based on basis-function type
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      dg_prim_vars_vlasov_u_i_kernels = ser_dg_prim_vars_vlasov_u_i_kernels;
      dg_prim_vars_vlasov_vth2_kernels = ser_dg_prim_vars_vlasov_vth2_kernels;
      dg_prim_vars_vlasov_kernels = ser_dg_prim_vars_vlasov_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  if (strcmp(prim_nm, "u_i") == 0) { // flow velocity
    dg_prim_vars_vlasov->pvt.kernel = dg_prim_vars_vlasov_u_i_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    dg_prim_vars_vlasov->pvt.num_mom = vdim; 
  }
  else if (strcmp(prim_nm, "vth2") == 0) { // thermal velocity squared
    dg_prim_vars_vlasov->pvt.kernel = dg_prim_vars_vlasov_vth2_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    dg_prim_vars_vlasov->pvt.num_mom = 1; 
  }
  else if (strcmp(prim_nm, "prim") == 0) { // combined (u_i, vth^2)
    dg_prim_vars_vlasov->pvt.kernel = dg_prim_vars_vlasov_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    dg_prim_vars_vlasov->pvt.num_mom = vdim+1; 
  }
  else {
    // string not recognized
    gkyl_exit("gkyl_dg_prim_vars_type_vlasov: Unrecognized primitive variable requested!");
  }

  dg_prim_vars_vlasov->pvt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(dg_prim_vars_vlasov->pvt.flags);
  dg_prim_vars_vlasov->pvt.ref_count = gkyl_ref_count_init(gkyl_dg_prim_vars_vlasov_free);
  
  dg_prim_vars_vlasov->pvt.on_dev = &dg_prim_vars_vlasov->pvt; // on host, self-reference
    
  return &dg_prim_vars_vlasov->pvt;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_prim_vars_type*
gkyl_dg_prim_vars_vlasov_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const char *prim_nm)
{
  assert(false);
  return 0;
}

#endif
