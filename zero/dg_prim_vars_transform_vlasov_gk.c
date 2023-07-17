#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_prim_vars_transform_vlasov_gk.h>
#include <gkyl_dg_prim_vars_transform_vlasov_gk_priv.h>
#include <gkyl_util.h>

void
gkyl_dg_prim_vars_transform_vlasov_gk_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_dg_prim_vars_type *pvt = container_of(ref, struct gkyl_dg_prim_vars_type, ref_count);
  if (GKYL_IS_CU_ALLOC(pvt->flags))
    gkyl_cu_free(pvt->on_dev);
  gkyl_free(pvt);
}

void
gkyl_dg_prim_vars_transform_vlasov_gk_set_auxfields(const struct gkyl_dg_prim_vars_type *pvt, struct gkyl_dg_prim_vars_auxfields auxin)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(auxin.b_i)) {
    gkyl_dg_prim_vars_transform_vlasov_gk_set_auxfields_cu(pvt->on_dev, auxin);
    return;
  }
#endif

  struct dg_prim_vars_type_transform_vlasov_gk *prim = container_of(pvt, struct dg_prim_vars_type_transform_vlasov_gk, pvt);
  prim->auxfields.b_i = auxin.b_i; // covariant components of the field aligned unit vector.
}

struct gkyl_dg_prim_vars_type*
gkyl_dg_prim_vars_transform_vlasov_gk_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis,
  const struct gkyl_range* conf_range, const char *prim_nm, bool use_gpu)
{
  assert(cbasis->poly_order == pbasis->poly_order);

#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_prim_vars_transform_vlasov_gk_cu_dev_new(cbasis, pbasis, prim_nm);
  } 
#endif    
  struct dg_prim_vars_type_transform_vlasov_gk *pvt = gkyl_malloc(sizeof(struct dg_prim_vars_type_transform_vlasov_gk));
  int cdim = cbasis->ndim, pdim = pbasis->ndim, vdim = pdim-cdim;
  int poly_order = cbasis->poly_order;

  pvt->pvt.cdim = cdim;
  pvt->pvt.vdim = vdim;
  pvt->pvt.poly_order = poly_order;
  pvt->pvt.num_config = cbasis->num_basis;

  const gkyl_dg_prim_vars_transform_vlasov_gk_kern_list *dg_prim_vars_transform_vlasov_gk_u_par_i_kernels, 
    *dg_prim_vars_transform_vlasov_gk_kernels;

  // choose kernel tables based on basis-function type
  switch (cbasis->b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      dg_prim_vars_transform_vlasov_gk_u_par_i_kernels = ser_dg_prim_vars_transform_vlasov_gk_u_par_i_kernels;
      dg_prim_vars_transform_vlasov_gk_kernels = ser_dg_prim_vars_transform_vlasov_gk_kernels;
      break;

    default:
      assert(false);
      break;    
  }

  if (strcmp(prim_nm, "u_par_i") == 0) { // projection of parallel velocity from GK to Vlasov u_par b_i
    pvt->pvt.kernel = dg_prim_vars_transform_vlasov_gk_u_par_i_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    pvt->pvt.num_mom = vdim; 
  }
  else if (strcmp(prim_nm, "prim") == 0) { // projection of Vlasov u_i, vth^2 to upar (u_i . b_i) and vth_GK^2 (M2/M0 - upar^2)
    pvt->pvt.kernel = dg_prim_vars_transform_vlasov_gk_kernels[cv_index[cdim].vdim[vdim]].kernels[poly_order];
    pvt->pvt.num_mom = 2; 
  }
  else {
    // string not recognized
    gkyl_exit("gkyl_dg_prim_vars_type_transform_vlasov_gk: Unrecognized primitive variable requested!");
  }

  pvt->auxfields.b_i = 0;
  pvt->conf_range = *conf_range;

  pvt->pvt.flags = 0;
  GKYL_CLEAR_CU_ALLOC(pvt->pvt.flags);
  pvt->pvt.ref_count = gkyl_ref_count_init(gkyl_dg_prim_vars_transform_vlasov_gk_free);
  
  pvt->pvt.on_dev = &pvt->pvt; // on host, self-reference
    
  return &pvt->pvt;
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_prim_vars_type*
gkyl_dg_prim_vars_transform_vlasov_gk_cu_dev_new(const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, 
  const struct gkyl_range* conf_range, const char *prim_nm)
{
  assert(false);
  return 0;
}

#endif
