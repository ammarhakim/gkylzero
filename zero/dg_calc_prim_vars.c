#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
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

void gkyl_calc_prim_vars_pkpm(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* u_i, struct gkyl_array* p_ij)
{
// Check if more than one of the output arrays is on device? 
// Probably a better way to do this (JJ: 11/16/22)
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(u_i)) {
    return gkyl_calc_prim_vars_pkpm_cu(basis, range, bvar, vlasov_pkpm_moms, euler_pkpm, u_i, p_ij);
  }
#endif

  int cdim = basis.ndim;
  int poly_order = basis.poly_order;
  euler_pkpm_prim_vars_t pkpm_prim_vars = choose_ser_euler_pkpm_prim_vars_kern(cdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *bvar_d = gkyl_array_cfetch(bvar, loc);
    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);
    
    double *u_i_d = gkyl_array_fetch(u_i, loc);
    double *p_ij_d = gkyl_array_fetch(p_ij, loc);
    pkpm_prim_vars(bvar_d, vlasov_pkpm_moms_d, euler_pkpm_d, u_i_d, p_ij_d);
  }
}

void gkyl_calc_prim_vars_pkpm_source(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* qmem, const struct gkyl_array* nu, const struct gkyl_array* nu_vthsq,
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  const struct gkyl_array* p_perp_source, struct gkyl_array* rhs)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(rhs)) {
    return gkyl_calc_prim_vars_pkpm_source_cu(basis, range, 
      qmem, nu, nu_vthsq, vlasov_pkpm_moms, euler_pkpm, p_perp_source, rhs);
  }
#endif

  int cdim = basis.ndim;
  int poly_order = basis.poly_order;
  euler_pkpm_source_t pkpm_source = choose_ser_euler_pkpm_source_kern(cdim, poly_order);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(range, iter.idx);

    const double *qmem_d = gkyl_array_cfetch(qmem, loc);
    const double *nu_d = gkyl_array_cfetch(nu, loc);
    const double *nu_vthsq_d = gkyl_array_cfetch(nu_vthsq, loc);
    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);
    const double *p_perp_source_d = gkyl_array_cfetch(p_perp_source, loc);

    double *rhs_d = gkyl_array_fetch(rhs, loc);
    pkpm_source(qmem_d, nu_d, nu_vthsq_d, vlasov_pkpm_moms_d, euler_pkpm_d, p_perp_source_d, rhs_d);
  }
}

void gkyl_calc_prim_vars_pkpm_recovery(const struct gkyl_rect_grid *grid, 
  struct gkyl_basis basis, const struct gkyl_range *range, double nuHyp, 
  const struct gkyl_array* bvar, const struct gkyl_array* u_i, 
  const struct gkyl_array* p_ij, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* div_b, struct gkyl_array* bb_grad_u, 
  struct gkyl_array* div_p, struct gkyl_array* p_force, struct gkyl_array* p_perp_source)
{
// Check if more than one of the output arrays is on device? 
// Probably a better way to do this (JJ: 11/16/22)
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(div_p)) {
    return gkyl_calc_prim_vars_pkpm_recovery_cu(grid, basis, range, nuHyp, 
      bvar, u_i, p_ij, vlasov_pkpm_moms, euler_pkpm, div_b, bb_grad_u, div_p, p_force, p_perp_source);
  }
#endif

  int cdim = basis.ndim;
  int poly_order = basis.poly_order;
  euler_pkpm_recovery_t pkpm_recovery[3];
  // Fetch the kernels in each direction
  for (int d=0; d<cdim; ++d) 
    pkpm_recovery[d] = choose_ser_euler_pkpm_recovery_kern(d, cdim, poly_order);

  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idxc);
    long linc = gkyl_range_idx(range, idxc);

    const double *bvar_c = gkyl_array_cfetch(bvar, linc);
    const double *u_i_c = gkyl_array_cfetch(u_i, linc);
    const double *p_ij_c = gkyl_array_cfetch(p_ij, linc);
    const double *vlasov_pkpm_moms_c = gkyl_array_cfetch(vlasov_pkpm_moms, linc);
    const double *euler_pkpm_c = gkyl_array_cfetch(euler_pkpm, linc);

    double *div_b_d = gkyl_array_fetch(div_b, linc);
    double *bb_grad_u_d = gkyl_array_fetch(bb_grad_u, linc);
    double *div_p_d = gkyl_array_fetch(div_p, linc);
    double *p_force_d = gkyl_array_fetch(p_force, linc);
    double *p_perp_source_d = gkyl_array_fetch(p_perp_source, linc);

    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, iter.idx, idxl);
      gkyl_copy_int_arr(cdim, iter.idx, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(range, idxl); 
      long linr = gkyl_range_idx(range, idxr);

      const double *bvar_l = gkyl_array_cfetch(bvar, linl);
      const double *bvar_r = gkyl_array_cfetch(bvar, linr);

      const double *u_i_l = gkyl_array_cfetch(u_i, linl);
      const double *u_i_r = gkyl_array_cfetch(u_i, linr);

      const double *p_ij_l = gkyl_array_cfetch(p_ij, linl);
      const double *p_ij_r = gkyl_array_cfetch(p_ij, linr);

      const double *vlasov_pkpm_moms_l = gkyl_array_cfetch(vlasov_pkpm_moms, linl);
      const double *vlasov_pkpm_moms_r = gkyl_array_cfetch(vlasov_pkpm_moms, linr);

      const double *euler_pkpm_l = gkyl_array_cfetch(euler_pkpm, linl);
      const double *euler_pkpm_r = gkyl_array_cfetch(euler_pkpm, linr);

      pkpm_recovery[dir](grid->dx, nuHyp, 
        bvar_l, bvar_c, bvar_r, u_i_l, u_i_c, u_i_r, 
        p_ij_l, p_ij_c, p_ij_r, vlasov_pkpm_moms_l, vlasov_pkpm_moms_c, vlasov_pkpm_moms_r, 
        euler_pkpm_l, euler_pkpm_c, euler_pkpm_r, 
        div_b_d, bb_grad_u_d, div_p_d, p_force_d, p_perp_source_d);
    }
  }
}
