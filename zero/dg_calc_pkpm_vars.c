#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_calc_pkpm_vars_priv.h>
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_util.h>

gkyl_dg_calc_pkpm_vars*
gkyl_dg_calc_pkpm_vars_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *cbasis_2p, 
  const struct gkyl_range *mem_range, 
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *geom, 
  double limiter_fac, bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_calc_pkpm_vars_cu_dev_new(conf_grid, cbasis, cbasis_2p, 
      mem_range, wv_eqn, geom, limiter_fac);
  } 
#endif     
  gkyl_dg_calc_pkpm_vars *up = gkyl_malloc(sizeof(gkyl_dg_calc_pkpm_vars));

  up->conf_grid = *conf_grid;
  int cdim = cbasis->ndim;
  up->cdim = cdim;
  up->mem_range = *mem_range;

  up->wv_eqn = gkyl_wv_eqn_acquire(wv_eqn);
  up->geom = gkyl_wave_geom_acquire(geom);
  // Limiter factor for relationship between slopes and cell average differences
  // By default, this factor is 1/sqrt(3) because cell_avg(f) = f0/sqrt(2^cdim)
  // and a cell slope estimate from two adjacent cells is (for the x variation): 
  // integral(psi_1 [cell_avg(f_{i+1}) - cell_avg(f_{i})]*x) = sqrt(2^cdim)/sqrt(3)*[cell_avg(f_{i+1}) - cell_avg(f_{i})]
  // where psi_1 is the x cell slope basis in our orthonormal expansion psi_1 = sqrt(3)/sqrt(2^cdim)*x
  // This factor can be made smaller (larger) to increase (decrease) the diffusion from the slope limiter
  if (limiter_fac == 0.0) {
    up->limiter_fac = 0.5773502691896258;
  }
  else {
    up->limiter_fac = limiter_fac;
  }

  // Set the function pointers for the order p quantities we solve for 
  int poly_order = cbasis->poly_order;
  up->pkpm_u_set = choose_pkpm_u_set_kern(cdim, poly_order);
  up->pkpm_u_copy = choose_pkpm_u_copy_kern(cdim, poly_order);  
  up->pkpm_u_surf = choose_pkpm_u_surf_kern(cdim, poly_order);
  up->pkpm_explicit_source = choose_pkpm_explicit_source_kern(cdim, poly_order);
  for (int d=0; d<cdim; ++d) {
    up->pkpm_limiter[d] = choose_pkpm_limiter_kern(d, cdim, poly_order);
  }
  // Linear system for solving for ux, uy, uz
  up->Ncomp_u = 3;
  int nc = cbasis->num_basis;
  up->As_u = gkyl_nmat_new(up->Ncomp_u*mem_range->volume, nc, nc);
  up->xs_u = gkyl_nmat_new(up->Ncomp_u*mem_range->volume, nc, 1);
  up->mem_u = gkyl_nmat_linsolve_lu_new(up->As_u->num, up->As_u->nr);

  // Set the function pointers for the order 2*p quantities we solve for
  int poly_order_2p = cbasis_2p->poly_order;
  up->pkpm_set = choose_pkpm_set_kern(cdim, poly_order_2p);
  up->pkpm_copy = choose_pkpm_copy_kern(cdim, poly_order_2p);
  up->pkpm_surf_set = choose_pkpm_surf_set_kern(cdim, poly_order_2p);
  up->pkpm_surf_copy = choose_pkpm_surf_copy_kern(cdim, poly_order_2p);
  up->pkpm_pressure = choose_pkpm_pressure_kern(cdim, poly_order_2p);
  up->pkpm_int = choose_pkpm_int_kern(cdim, poly_order_2p);
  up->pkpm_io = choose_pkpm_io_kern(cdim, poly_order_2p);
  // Fetch the kernels in each direction
  for (int d=0; d<cdim; ++d) {
    up->pkpm_accel[d] = choose_pkpm_accel_kern(d, cdim, poly_order_2p);
  }
  // 6 components: 1/rho*div(p_par b), p_perp/rho, rho/p_perp, 3*Txx/m, 3*Tyy/m, 3*Tzz/m, p_perp/rho*div(b)
  up->Ncomp_prim = 7;
  int nc_2p = cbasis_2p->num_basis;
  up->As = gkyl_nmat_new(up->Ncomp_prim*mem_range->volume, nc_2p, nc_2p);
  up->xs = gkyl_nmat_new(up->Ncomp_prim*mem_range->volume, nc_2p, 1);
  up->mem = gkyl_nmat_linsolve_lu_new(up->As->num, up->As->nr);

  // 2*cdim components: 
  // 3*Txx/m at the left and right x surfaces 
  // 3*Tyy/m at the left and right y surfaces 
  // 3*Tzz/m at the left and right z surfaces
  up->Ncomp_surf = 2*cdim;
  int nc_surf_2p = cbasis_2p->num_basis/(poly_order_2p+1); 
  up->As_surf = gkyl_nmat_new(up->Ncomp_surf*mem_range->volume, nc_surf_2p, nc_surf_2p);
  up->xs_surf = gkyl_nmat_new(up->Ncomp_surf*mem_range->volume, nc_surf_2p, 1);
  up->mem_surf = gkyl_nmat_linsolve_lu_new(up->As_surf->num, up->As_surf->nr);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_dg_calc_pkpm_vars_advance(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* p_ij, 
  const struct gkyl_array* pkpm_div_ppar, const struct gkyl_array* div_b, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* prim, struct gkyl_array* pkpm_accel)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(prim)) {
    return gkyl_dg_calc_pkpm_vars_advance_cu(up, 
      vlasov_pkpm_moms, p_ij, pkpm_div_ppar, div_b, 
      cell_avg_prim, prim, pkpm_accel);
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *p_ij_d = gkyl_array_cfetch(p_ij, loc);
    const double *pkpm_div_ppar_d = gkyl_array_cfetch(pkpm_div_ppar, loc);
    const double *div_b_d = gkyl_array_cfetch(div_b, loc);

    int* cell_avg_prim_d = gkyl_array_fetch(cell_avg_prim, loc);
    // First index of cell_avg_prim is whether p = p_par + 2 p_perp < 0.0 at control points
    cell_avg_prim_d[1] = up->pkpm_set(count, up->As, up->xs, 
      vlasov_pkpm_moms_d, p_ij_d, pkpm_div_ppar_d, div_b_d);

    count += up->Ncomp_prim;
  }

  bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
  assert(status);

  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    double* prim_d = gkyl_array_fetch(prim, loc);
    double* pkpm_accel_d = gkyl_array_fetch(pkpm_accel, loc);

    up->pkpm_copy(count, up->xs, prim_d, pkpm_accel_d);

    count += up->Ncomp_prim;
  }
}

void gkyl_dg_calc_pkpm_vars_surf_advance(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* p_ij, 
  struct gkyl_array* prim_surf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(prim_surf)) {
    return gkyl_dg_calc_pkpm_vars_surf_advance_cu(up, 
      vlasov_pkpm_moms, p_ij, prim_surf);
  }
#endif

  // First loop over mem_range for solving linear systems to compute surface primitive moments
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *p_ij_d = gkyl_array_cfetch(p_ij, loc);

    up->pkpm_surf_set(count, up->As_surf, up->xs_surf, 
      vlasov_pkpm_moms_d, p_ij_d);

    count += up->Ncomp_surf;
  }

  // Only need to solve the linear system if cdim > 1
  if (up->cdim > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem_surf, up->As_surf, up->xs_surf);
    assert(status);
  }

  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    double* prim_surf_d = gkyl_array_fetch(prim_surf, loc);

    up->pkpm_surf_copy(count, up->xs_surf, prim_surf_d);

    count += up->Ncomp_surf;
  }
}

void gkyl_dg_calc_pkpm_vars_u(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* pkpm_u)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pkpm_u)) {
    return gkyl_dg_calc_pkpm_vars_u_cu(up, 
      vlasov_pkpm_moms, euler_pkpm, 
      cell_avg_prim, pkpm_u);
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);

    int* cell_avg_prim_d = gkyl_array_fetch(cell_avg_prim, loc);
    // Zeroth index of cell_avg_prim is whether rho < 0.0 at control points
    cell_avg_prim_d[0] = up->pkpm_u_set(count, up->As_u, up->xs_u, 
      vlasov_pkpm_moms_d, euler_pkpm_d);

    count += up->Ncomp_u;
  }

  bool status = gkyl_nmat_linsolve_lu_pa(up->mem_u, up->As_u, up->xs_u);
  assert(status);

  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    double* pkpm_u_d = gkyl_array_fetch(pkpm_u, loc);

    up->pkpm_u_copy(count, up->xs_u, pkpm_u_d);

    count += up->Ncomp_u;
  }
}

void gkyl_dg_calc_pkpm_vars_u_surf(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* pkpm_u, struct gkyl_array* pkpm_u_surf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pkpm_u)) {
    return gkyl_dg_calc_pkpm_vars_u_surf_cu(up, 
      pkpm_u, pkpm_u_surf);
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    const double* pkpm_u_d = gkyl_array_cfetch(pkpm_u, loc);
    double* pkpm_u_surf_d = gkyl_array_fetch(pkpm_u_surf, loc);

    up->pkpm_u_surf(pkpm_u_d, pkpm_u_surf_d);
  }
}

void gkyl_dg_calc_pkpm_vars_pressure(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bb, const struct gkyl_array* vlasov_pkpm_moms, struct gkyl_array* p_ij)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(p_ij)) {
    return gkyl_dg_calc_pkpm_vars_pressure_cu(up, conf_range, 
      bb, vlasov_pkpm_moms, p_ij);
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *bb_d = gkyl_array_cfetch(bb, loc);
    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);

    double* p_ij_d = gkyl_array_fetch(p_ij, loc);

    up->pkpm_pressure(bb_d, vlasov_pkpm_moms_d, p_ij_d);
  }
}

void gkyl_dg_calc_pkpm_vars_accel(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* pkpm_u_surf, const struct gkyl_array* pkpm_u, 
  const struct gkyl_array* pkpm_prim_surf, 
  const struct gkyl_array* bb, const struct gkyl_array* div_b, const struct gkyl_array* nu, 
  struct gkyl_array* pkpm_lax, struct gkyl_array* pkpm_accel)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pkpm_accel)) {
    return gkyl_dg_calc_pkpm_vars_accel_cu(up, conf_range, 
      pkpm_u_surf, pkpm_u, pkpm_prim_surf, 
      bb, div_b, nu, 
      pkpm_lax, pkpm_accel);
  }
#endif

  // Loop over configuration space range to compute gradients and pkpm acceleration variables
  int cdim = up->cdim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idxc);
    long linc = gkyl_range_idx(conf_range, idxc);

    const double *pkpm_u_surf_c = gkyl_array_cfetch(pkpm_u_surf, linc);
    const double *pkpm_prim_surf_c = gkyl_array_cfetch(pkpm_prim_surf, linc);

    const double *pkpm_u_d = gkyl_array_cfetch(pkpm_u, linc);
    const double *bb_d = gkyl_array_cfetch(bb, linc);
    const double *div_b_d = gkyl_array_cfetch(div_b, linc);
    const double *nu_d = gkyl_array_cfetch(nu, linc);

    double *pkpm_lax_d = gkyl_array_fetch(pkpm_lax, linc);
    double *pkpm_accel_d = gkyl_array_fetch(pkpm_accel, linc);

    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, iter.idx, idxl);
      gkyl_copy_int_arr(cdim, iter.idx, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(conf_range, idxl); 
      long linr = gkyl_range_idx(conf_range, idxr);

      const double *pkpm_u_surf_l = gkyl_array_cfetch(pkpm_u_surf, linl);
      const double *pkpm_u_surf_r = gkyl_array_cfetch(pkpm_u_surf, linr);
      const double *pkpm_prim_surf_l = gkyl_array_cfetch(pkpm_prim_surf, linl);
      const double *pkpm_prim_surf_r = gkyl_array_cfetch(pkpm_prim_surf, linr);

      up->pkpm_accel[dir](up->conf_grid.dx, 
        pkpm_u_surf_l, pkpm_u_surf_c, pkpm_u_surf_r, 
        pkpm_prim_surf_l, pkpm_prim_surf_c, pkpm_prim_surf_r, 
        pkpm_u_d, bb_d, nu_d,
        pkpm_lax_d, pkpm_accel_d);
    }
  }
}

void gkyl_dg_calc_pkpm_integrated_vars(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, 
  struct gkyl_array* pkpm_int_vars)
{
// Check if more than one of the output arrays is on device? 
// Probably a better way to do this (JJ: 11/16/22)
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pkpm_int_vars)) {
    return gkyl_dg_calc_pkpm_integrated_vars_cu(up, conf_range, 
      vlasov_pkpm_moms, pkpm_u, pkpm_int_vars);
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *pkpm_u_d = gkyl_array_cfetch(pkpm_u, loc);
    
    double *pkpm_int_vars_d = gkyl_array_fetch(pkpm_int_vars, loc);
    up->pkpm_int(vlasov_pkpm_moms_d, pkpm_u_d, pkpm_int_vars_d);
  }
}

void gkyl_dg_calc_pkpm_vars_explicit_source(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* qmem, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* rhs)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(rhs)) {
    return gkyl_dg_calc_pkpm_vars_explicit_source_cu(up, conf_range, 
      qmem, vlasov_pkpm_moms, euler_pkpm, rhs);
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *qmem_d = gkyl_array_cfetch(qmem, loc);
    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);

    double *rhs_d = gkyl_array_fetch(rhs, loc);
    up->pkpm_explicit_source(qmem_d, vlasov_pkpm_moms_d, euler_pkpm_d, rhs_d);
  }
}

void gkyl_dg_calc_pkpm_vars_io(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, const struct gkyl_array* p_ij, 
  const struct gkyl_array* prim, const struct gkyl_array* pkpm_accel, 
  struct gkyl_array* fluid_io, struct gkyl_array* pkpm_vars_io)
{
// Check if more than one of the output arrays is on device? 
// Probably a better way to do this (JJ: 11/16/22)
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pkpm_vars_io)) {
    return gkyl_dg_calc_pkpm_vars_io_cu(up, conf_range, 
      vlasov_pkpm_moms, pkpm_u, p_ij, prim, pkpm_accel, 
      fluid_io, pkpm_vars_io);
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *pkpm_u_d = gkyl_array_cfetch(pkpm_u, loc);
    const double *p_ij_d = gkyl_array_cfetch(p_ij, loc);
    const double *prim_d = gkyl_array_cfetch(prim, loc);
    const double *pkpm_accel_d = gkyl_array_cfetch(pkpm_accel, loc);
    
    double *fluid_io_d = gkyl_array_fetch(fluid_io, loc);
    double *pkpm_vars_io_d = gkyl_array_fetch(pkpm_vars_io, loc);
    up->pkpm_io(vlasov_pkpm_moms_d, pkpm_u_d, p_ij_d, prim_d, pkpm_accel_d, 
      fluid_io_d, pkpm_vars_io_d);
  }
}

void gkyl_dg_calc_pkpm_vars_limiter(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, const struct gkyl_array* p_ij, 
  struct gkyl_array* euler_pkpm)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(euler_pkpm)) {
    return gkyl_dg_calc_pkpm_vars_limiter_cu(up, conf_range, 
      vlasov_pkpm_moms, pkpm_u, p_ij, euler_pkpm);
  }
#endif
  int cdim = up->cdim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idxc);
    long linc = gkyl_range_idx(conf_range, idxc);
    const struct gkyl_wave_cell_geom *geom = gkyl_wave_geom_get(up->geom, idxc);

    const double *vlasov_pkpm_moms_c = gkyl_array_cfetch(vlasov_pkpm_moms, linc);
    const double *pkpm_u_c = gkyl_array_cfetch(pkpm_u, linc);
    const double *p_ij_c = gkyl_array_cfetch(p_ij, linc);

    double *euler_pkpm_c = gkyl_array_fetch(euler_pkpm, linc);
    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, iter.idx, idxl);
      gkyl_copy_int_arr(cdim, iter.idx, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(conf_range, idxl); 
      long linr = gkyl_range_idx(conf_range, idxr);

      const double *vlasov_pkpm_moms_l = gkyl_array_cfetch(vlasov_pkpm_moms, linl);
      const double *vlasov_pkpm_moms_r = gkyl_array_cfetch(vlasov_pkpm_moms, linr);
      const double *pkpm_u_l = gkyl_array_cfetch(pkpm_u, linl);
      const double *pkpm_u_r = gkyl_array_cfetch(pkpm_u, linr);
      const double *p_ij_l = gkyl_array_cfetch(p_ij, linl);
      const double *p_ij_r = gkyl_array_cfetch(p_ij, linr);

      double *euler_pkpm_l = gkyl_array_fetch(euler_pkpm, linl);
      double *euler_pkpm_r = gkyl_array_fetch(euler_pkpm, linr);

      up->pkpm_limiter[dir](up->limiter_fac, up->wv_eqn, geom, 
        vlasov_pkpm_moms_l, vlasov_pkpm_moms_c, vlasov_pkpm_moms_r, 
        pkpm_u_l, pkpm_u_c, pkpm_u_r, 
        p_ij_l, p_ij_c, p_ij_r, 
        euler_pkpm_l, euler_pkpm_c, euler_pkpm_r); 
    }
  }
}

void gkyl_dg_calc_pkpm_vars_release(gkyl_dg_calc_pkpm_vars *up)
{
  gkyl_wv_eqn_release(up->wv_eqn);
  gkyl_wave_geom_release(up->geom);

  gkyl_nmat_release(up->As);
  gkyl_nmat_release(up->xs);
  gkyl_nmat_linsolve_lu_release(up->mem);

  gkyl_nmat_release(up->As_surf);
  gkyl_nmat_release(up->xs_surf);
  gkyl_nmat_linsolve_lu_release(up->mem_surf);

  gkyl_nmat_release(up->As_u);
  gkyl_nmat_release(up->xs_u);
  gkyl_nmat_linsolve_lu_release(up->mem_u);
  
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
