#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_prim_vars.h>
#include <gkyl_dg_calc_prim_vars_priv.h>
#include <gkyl_util.h>

void gkyl_calc_prim_vars_u_from_statevec(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* statevec, struct gkyl_array* u_i)
{
  // Find number of components of flow vector
  int num_comp = u_i->ncomp/basis.num_basis;
  for (int i = 0; i<num_comp; ++i)
    gkyl_dg_div_op_range(mem, basis, 
      i, u_i, i+1, statevec, 0, statevec, range);  
}

void gkyl_calc_prim_vars_u_from_rhou(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* rho, const struct gkyl_array* rhou, struct gkyl_array* u_i)
{
  // Find number of components of flow vector
  int num_comp = u_i->ncomp/basis.num_basis;
  for (int i = 0; i<num_comp; ++i)
    gkyl_dg_div_op_range(mem, basis, 
      i, u_i, i, rhou, 0, rho, range);  
}

void gkyl_calc_prim_vars_p_from_statevec(struct gkyl_basis basis, const struct gkyl_range *range,
  const double p_fac,  const struct gkyl_array* u_i, const struct gkyl_array* statevec,
  struct gkyl_array* p_ij)
{
  int cdim = basis.ndim;
  int poly_order = basis.poly_order;
  euler_pressure_t pressure = choose_ser_euler_pressure_kern(cdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *statevec_d = gkyl_array_cfetch(statevec, loc);
    const double *u_i_d = gkyl_array_cfetch(u_i, loc);
    double *p_ij_d = gkyl_array_fetch(p_ij, loc);
    pressure(p_fac, u_i_d, statevec_d, p_ij_d);
  }
}

void gkyl_calc_prim_vars_p_pkpm(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* p_ij)
{
  int cdim = basis.ndim;
  int poly_order = basis.poly_order;
  euler_pkpm_pressure_t pkpm_pressure = choose_ser_euler_pkpm_pressure_kern(cdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *bvar_d = gkyl_array_cfetch(bvar, loc);
    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);
    
    double *p_ij_d = gkyl_array_fetch(p_ij, loc);
    pkpm_pressure(bvar_d, vlasov_pkpm_moms_d, euler_pkpm_d, p_ij_d);
  }
}

void gkyl_calc_prim_vars_p_pkpm_source(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* nu, const struct gkyl_array* nu_vthsq,
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* u_i, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* rhs)
{
  int cdim = basis.ndim;
  int poly_order = basis.poly_order;
  euler_pkpm_pressure_source_t pkpm_pressure_source = choose_ser_euler_pkpm_pressure_source_kern(cdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *nu_d = gkyl_array_cfetch(nu, loc);
    const double *nu_vthsq_d = gkyl_array_cfetch(nu_vthsq, loc);
    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *u_i_d = gkyl_array_cfetch(u_i, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);

    double *rhs_d = gkyl_array_fetch(rhs, loc);
    pkpm_pressure_source(nu_d, nu_vthsq_d, vlasov_pkpm_moms_d, u_i_d, euler_pkpm_d, rhs_d);
  }
}