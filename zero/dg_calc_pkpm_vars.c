#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_bin_ops_priv.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_calc_pkpm_vars_priv.h>
#include <gkyl_util.h>

gkyl_dg_calc_pkpm_vars*
gkyl_dg_calc_pkpm_vars_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_range *mem_range, 
  bool use_gpu)
{
#ifdef GKYL_HAVE_CUDA
  if(use_gpu) {
    return gkyl_dg_calc_pkpm_vars_cu_dev_new(conf_grid, cbasis, mem_range);
  } 
#endif     
  gkyl_dg_calc_pkpm_vars *up = gkyl_malloc(sizeof(gkyl_dg_calc_pkpm_vars));

  up->conf_grid = *conf_grid;
  int nc = cbasis->num_basis;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  int nc_surf = cbasis->num_basis/(poly_order+1); // *only valid for tensor bases for cdim > 1*
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->cdim = cdim;
  up->poly_order = poly_order;
  up->Ncomp = 6;
  up->Ncomp_surf = 2*cdim*3+2*cdim;

  up->pkpm_set = choose_pkpm_set_kern(b_type, cdim, poly_order);
  up->pkpm_surf_set = choose_pkpm_surf_set_kern(b_type, cdim, poly_order);
  up->pkpm_copy = choose_pkpm_copy_kern(b_type, cdim, poly_order);
  up->pkpm_surf_copy = choose_pkpm_surf_copy_kern(b_type, cdim, poly_order);
  up->pkpm_pressure = choose_pkpm_pressure_kern(b_type, cdim, poly_order);
  up->pkpm_source = choose_pkpm_source_kern(b_type, cdim, poly_order);
  up->pkpm_int = choose_pkpm_int_kern(b_type, cdim, poly_order);
  up->pkpm_io = choose_pkpm_io_kern(b_type, cdim, poly_order);
  // Fetch the kernels in each direction
  for (int d=0; d<cdim; ++d) 
    up->pkpm_accel[d] = choose_pkpm_accel_kern(d, b_type, cdim, poly_order);

  // There are Ncomp*range->volume linear systems to be solved 
  // 6 components: ux, uy, uz, div(p_par b)/rho, p_perp/rho, rho/p_perp
  up->As = gkyl_nmat_new(up->Ncomp*mem_range->volume, nc, nc);
  up->xs = gkyl_nmat_new(up->Ncomp*mem_range->volume, nc, 1);
  up->mem = gkyl_nmat_linsolve_lu_new(up->As->num, up->As->nr);

  // There are Ncomp_surf*range->volume linear systems to be solved 
  // Each linear system is nc_surf x nc_surf (only solved over the surface basis and only when poly_order and cdim > 1)
  // 2*cdim*3+2*cdim components: ux, uy, uz (3 components) at the left and right of the cell (2 components) in each dimension (cdim components)
  // Also solves for 3*Txx/m at the left and right x surfaces, 3*Tyy/m at the left and right y surfaces, and 3*Tzz/m at the left and right z surfaces
  up->As_surf = gkyl_nmat_new(up->Ncomp_surf*mem_range->volume, nc_surf, nc_surf);
  up->xs_surf = gkyl_nmat_new(up->Ncomp_surf*mem_range->volume, nc_surf, 1);
  up->mem_surf = gkyl_nmat_linsolve_lu_new(up->As_surf->num, up->As_surf->nr);

  up->mem_range = *mem_range;

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  
  return up;
}

void gkyl_dg_calc_pkpm_vars_advance(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* pkpm_div_ppar, struct gkyl_array* cell_avg_prim, 
  struct gkyl_array* prim)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(prim)) {
    return gkyl_dg_calc_pkpm_vars_advance_cu(up, 
      vlasov_pkpm_moms, euler_pkpm, 
      pkpm_div_ppar, 
      cell_avg_prim, prim);
  }
#endif

  // First loop over mem_range for solving linear systems to compute primitive moments
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);
    const double *pkpm_div_ppar_d = gkyl_array_cfetch(pkpm_div_ppar, loc);

    int* cell_avg_prim_d = gkyl_array_fetch(cell_avg_prim, loc);

    cell_avg_prim_d[0] = up->pkpm_set(count, up->As, up->xs, 
      vlasov_pkpm_moms_d, euler_pkpm_d, pkpm_div_ppar_d);

    count += up->Ncomp;
  }

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_range_iter_init(&iter, &up->mem_range);
  count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    double* prim_d = gkyl_array_fetch(prim, loc);

    up->pkpm_copy(count, up->xs, prim_d);

    count += up->Ncomp;
  }
}

void gkyl_dg_calc_pkpm_vars_surf_advance(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* p_ij, const struct gkyl_array* cell_avg_prim, 
  struct gkyl_array* prim_surf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(prim_surf)) {
    return gkyl_dg_calc_pkpm_vars_surf_advance_cu(up, 
      vlasov_pkpm_moms, euler_pkpm, 
      p_ij, cell_avg_prim, prim_surf);
  }
#endif

  // First loop over mem_range for solving linear systems to compute surface primitive moments
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->mem_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->mem_range, iter.idx);

    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);
    const double *p_ij_d = gkyl_array_cfetch(p_ij, loc);
    const int* cell_avg_prim_d = gkyl_array_cfetch(cell_avg_prim, loc);

    up->pkpm_surf_set(count, up->As_surf, up->xs_surf, 
      vlasov_pkpm_moms_d, euler_pkpm_d, p_ij_d, cell_avg_prim_d);

    count += up->Ncomp_surf;
  }

  // Only need to solve the linear system if both poly_order and cdim > 1
  if (up->poly_order > 1 && up->cdim > 1) {
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

void gkyl_dg_calc_pkpm_vars_pressure(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, struct gkyl_array* p_ij)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(p_ij)) {
    return gkyl_dg_calc_pkpm_vars_pressure_cu(up, conf_range, 
      bvar, vlasov_pkpm_moms, p_ij);
  }
#endif
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *bvar_d = gkyl_array_cfetch(bvar, loc);
    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);

    double* p_ij_d = gkyl_array_fetch(p_ij, loc);

    up->pkpm_pressure(bvar_d, vlasov_pkpm_moms_d, p_ij_d);
  }
}

void gkyl_dg_calc_pkpm_vars_accel(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bvar, const struct gkyl_array* prim_surf, 
  const struct gkyl_array* prim, const struct gkyl_array* nu, 
  struct gkyl_array* pkpm_lax, struct gkyl_array* pkpm_accel)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pkpm_accel)) {
    return gkyl_dg_calc_pkpm_vars_accel_cu(up, conf_range, 
      bvar, prim_surf, prim, nu, pkpm_lax, pkpm_accel);
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

    const double *bvar_c = gkyl_array_cfetch(bvar, linc);
    const double *prim_surf_c = gkyl_array_cfetch(prim_surf, linc);
    // Only need nu and the volume expansion of the primitive moments in center cell
    const double *prim_d = gkyl_array_cfetch(prim, linc);
    const double *nu_d = gkyl_array_cfetch(nu, linc);

    double *pkpm_lax_d = gkyl_array_fetch(pkpm_lax, linc);
    double *pkpm_accel_d = gkyl_array_fetch(pkpm_accel, linc);

    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, iter.idx, idxl);
      gkyl_copy_int_arr(cdim, iter.idx, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(conf_range, idxl); 
      long linr = gkyl_range_idx(conf_range, idxr);

      const double *bvar_l = gkyl_array_cfetch(bvar, linl);
      const double *bvar_r = gkyl_array_cfetch(bvar, linr);

      const double *prim_surf_l = gkyl_array_cfetch(prim_surf, linl);
      const double *prim_surf_r = gkyl_array_cfetch(prim_surf, linr);

      up->pkpm_accel[dir](up->conf_grid.dx, 
        bvar_l, bvar_c, bvar_r, 
        prim_surf_l, prim_surf_c, prim_surf_r, 
        prim_d, nu_d,
        pkpm_lax_d, pkpm_accel_d);
    }
  }
}

void gkyl_dg_calc_pkpm_integrated_vars(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* vlasov_pkpm_moms, 
  const struct gkyl_array* euler_pkpm, const struct gkyl_array* prim, 
  struct gkyl_array* pkpm_int_vars)
{
// Check if more than one of the output arrays is on device? 
// Probably a better way to do this (JJ: 11/16/22)
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pkpm_int_vars)) {
    return gkyl_dg_calc_pkpm_integrated_vars_cu(up, conf_range, 
      vlasov_pkpm_moms, euler_pkpm, prim, pkpm_int_vars);
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);
    const double *prim_d = gkyl_array_cfetch(prim, loc);
    
    double *pkpm_int_vars_d = gkyl_array_fetch(pkpm_int_vars, loc);
    up->pkpm_int(vlasov_pkpm_moms_d, euler_pkpm_d, prim_d, pkpm_int_vars_d);
  }
}

void gkyl_dg_calc_pkpm_vars_source(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* qmem, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* rhs)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(rhs)) {
    return gkyl_dg_calc_pkpm_vars_source_cu(up, conf_range, 
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
    up->pkpm_source(qmem_d, vlasov_pkpm_moms_d, euler_pkpm_d, rhs_d);
  }
}

void gkyl_dg_calc_pkpm_vars_io(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* vlasov_pkpm_moms, 
  const struct gkyl_array* euler_pkpm, const struct gkyl_array* p_ij, 
  const struct gkyl_array* prim, const struct gkyl_array* pkpm_accel, 
  struct gkyl_array* fluid_io, struct gkyl_array* pkpm_vars_io)
{
// Check if more than one of the output arrays is on device? 
// Probably a better way to do this (JJ: 11/16/22)
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(pkpm_vars_io)) {
    return gkyl_dg_calc_pkpm_vars_io_cu(up, conf_range, 
      vlasov_pkpm_moms, euler_pkpm, p_ij, prim, pkpm_accel, 
      fluid_io, pkpm_vars_io);
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(conf_range, iter.idx);

    const double *vlasov_pkpm_moms_d = gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = gkyl_array_cfetch(euler_pkpm, loc);
    const double *p_ij_d = gkyl_array_cfetch(p_ij, loc);
    const double *prim_d = gkyl_array_cfetch(prim, loc);
    const double *pkpm_accel_d = gkyl_array_cfetch(pkpm_accel, loc);
    
    double *fluid_io_d = gkyl_array_fetch(fluid_io, loc);
    double *pkpm_vars_io_d = gkyl_array_fetch(pkpm_vars_io, loc);
    up->pkpm_io(vlasov_pkpm_moms_d, euler_pkpm_d, p_ij_d, prim_d, pkpm_accel_d, 
      fluid_io_d, pkpm_vars_io_d);
  }
}

void gkyl_dg_calc_pkpm_vars_release(gkyl_dg_calc_pkpm_vars *up)
{
  gkyl_nmat_release(up->As);
  gkyl_nmat_release(up->xs);
  gkyl_nmat_linsolve_lu_release(up->mem);
  
  if (GKYL_IS_CU_ALLOC(up->flags))
    gkyl_cu_free(up->on_dev);
  gkyl_free(up);
}
