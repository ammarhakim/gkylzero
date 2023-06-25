/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_calc_pkpm_vars_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_calc_pkpm_vars_prim_cu_kernel(struct gkyl_basis basis, struct gkyl_range range, 
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* u_i, struct gkyl_array* p_ij, struct gkyl_array* T_ij, 
  struct gkyl_array* rho_inv, struct gkyl_array* T_perp_over_m, struct gkyl_array* T_perp_over_m_inv)
{
  int cdim = basis.ndim;
  int poly_order = basis.poly_order;

  pkpm_prim_t pkpm_prim = choose_ser_pkpm_prim_kern(cdim, poly_order);

  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&range, idx);

    const double *bvar_d = (const double*) gkyl_array_cfetch(bvar, start);
    const double *vlasov_pkpm_moms_d = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, start);
    const double *euler_pkpm_d = (const double*) gkyl_array_cfetch(euler_pkpm, start);
    
    double *u_i_d = (double*) gkyl_array_fetch(u_i, start);
    double *p_ij_d = (double*) gkyl_array_fetch(p_ij, start);
    double *T_ij_d = (double*) gkyl_array_fetch(T_ij, start);
    double *rho_inv_d = (double*) gkyl_array_fetch(rho_inv, start);
    double *T_perp_over_m_d = (double*) gkyl_array_fetch(T_perp_over_m, start);
    double *T_perp_over_m_inv_d = (double*) gkyl_array_fetch(T_perp_over_m_inv, start);

    pkpm_prim(bvar_d, vlasov_pkpm_moms_d, euler_pkpm_d, 
      u_i_d, p_ij_d, T_ij_d, 
      rho_inv_d, T_perp_over_m_d, T_perp_over_m_inv_d);
  }
}

// Host-side wrapper for pkpm primitive variable calculations
void
gkyl_calc_pkpm_vars_prim_cu(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* u_i, struct gkyl_array* p_ij, struct gkyl_array* T_ij, 
  struct gkyl_array* rho_inv, struct gkyl_array* T_perp_over_m, struct gkyl_array* T_perp_over_m_inv)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_calc_pkpm_vars_prim_cu_kernel<<<nblocks, nthreads>>>(basis, *range, 
    bvar->on_dev, vlasov_pkpm_moms->on_dev, euler_pkpm->on_dev, 
    u_i->on_dev, p_ij->on_dev, T_ij->on_dev, 
    rho_inv->on_dev, T_perp_over_m->on_dev, T_perp_over_m_inv->on_dev);
}

__global__ void
gkyl_calc_pkpm_vars_source_cu_kernel(struct gkyl_basis basis, struct gkyl_range range, 
  const struct gkyl_array* qmem, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* rhs)
{
  int cdim = basis.ndim;
  int poly_order = basis.poly_order;

  pkpm_source_t pkpm_source = choose_ser_pkpm_source_kern(cdim, poly_order);

  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start = gkyl_range_idx(&range, idx);

    const double *qmem_d = (const double*) gkyl_array_cfetch(qmem, start);
    const double *vlasov_pkpm_moms_d = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, start);
    const double *euler_pkpm_d = (const double*) gkyl_array_cfetch(euler_pkpm, start);

    double *rhs_d = (double*) gkyl_array_fetch(rhs, start);
    pkpm_source(qmem_d, vlasov_pkpm_moms_d, euler_pkpm_d, rhs_d);
  }
}

// Host-side wrapper for pkpm source term calculations
void
gkyl_calc_pkpm_vars_source_cu(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* qmem, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* rhs)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_calc_pkpm_vars_source_cu_kernel<<<nblocks, nthreads>>>(basis, *range, 
    qmem->on_dev, vlasov_pkpm_moms->on_dev, euler_pkpm->on_dev, 
    rhs->on_dev);
}

__global__ void
gkyl_calc_pkpm_vars_dist_mirror_force_cu_kernel(struct gkyl_rect_grid grid, struct gkyl_basis basis, 
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array* T_perp_over_m, const struct gkyl_array* T_perp_over_m_inv, 
  const struct gkyl_array* nu_vthsq, const struct gkyl_array* pkpm_accel_vars, 
  const struct gkyl_array* fIn, const struct gkyl_array* F_k_p_1,
  struct gkyl_array* g_dist_source, struct gkyl_array* F_k_m_1)
{
  double dx[GKYL_MAX_DIM] = {0.0};
  double xc[GKYL_MAX_DIM] = {0.0};
  int cdim = basis.ndim;
  int poly_order = basis.poly_order;
  for (int d=0; d<cdim+1; ++d) 
    dx[d] = grid.dx[d];

  pkpm_dist_mirror_force_t pkpm_dist_mirror_force;
  switch (basis.b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      pkpm_dist_mirror_force = choose_ser_pkpm_dist_mirror_force_kern(cdim, poly_order);

      break;

    case GKYL_BASIS_MODAL_TENSOR:
      pkpm_dist_mirror_force = choose_ten_pkpm_dist_mirror_force_kern(cdim, poly_order);
      
      break;

    default:
      assert(false);
      break;    
  }

  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    gkyl_rect_grid_cell_center(&grid, idx, xc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long start_conf = gkyl_range_idx(&conf_range, idx);
    long start_phase = gkyl_range_idx(&phase_range, idx);

    const double *T_perp_over_m_d = (const double*) gkyl_array_cfetch(T_perp_over_m, start_conf);
    const double *T_perp_over_m_inv_d = (const double*) gkyl_array_cfetch(T_perp_over_m_inv, start_conf);
    const double *nu_vthsq_d = (const double*) gkyl_array_cfetch(nu_vthsq, start_conf);
    const double *pkpm_accel_vars_d = (const double*) gkyl_array_cfetch(pkpm_accel_vars, start_conf);
    const double *fIn_d = (const double*) gkyl_array_cfetch(fIn, start_phase);
    const double *F_k_p_1_d = (const double*) gkyl_array_cfetch(F_k_p_1, start_phase);

    double *g_dist_source_d = (double*) gkyl_array_fetch(g_dist_source, start_phase);
    double *F_k_m_1_d = (double*) gkyl_array_fetch(F_k_m_1, start_phase);

    pkpm_dist_mirror_force(xc, dx, 
      T_perp_over_m_d, T_perp_over_m_inv_d, 
      nu_vthsq_d, pkpm_accel_vars_d, 
      fIn_d, F_k_p_1_d, g_dist_source_d, F_k_m_1_d);
  }  
}
// Host-side wrapper for pkpm mirror force source distribution function calculation
void 
gkyl_calc_pkpm_vars_dist_mirror_force_cu(const struct gkyl_rect_grid *grid, struct gkyl_basis basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_array* T_perp_over_m, const struct gkyl_array* T_perp_over_m_inv, 
  const struct gkyl_array* nu_vthsq, const struct gkyl_array* pkpm_accel_vars, 
  const struct gkyl_array* fIn, const struct gkyl_array* F_k_p_1,
  struct gkyl_array* g_dist_source, struct gkyl_array* F_k_m_1)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_calc_pkpm_vars_dist_mirror_force_cu_kernel<<<nblocks, nthreads>>>(*grid, basis, 
    *conf_range, *phase_range, 
    T_perp_over_m->on_dev, T_perp_over_m_inv->on_dev, 
    nu_vthsq->on_dev, pkpm_accel_vars->on_dev, 
    fIn->on_dev, F_k_p_1->on_dev, 
    g_dist_source->on_dev, F_k_m_1->on_dev);
}

__global__ void
gkyl_calc_pkpm_vars_recovery_cu_kernel(struct gkyl_rect_grid grid, struct gkyl_basis basis, struct gkyl_range range, 
  const struct gkyl_array* bvar, const struct gkyl_array* u_i, const struct gkyl_array* p_ij, 
  const struct gkyl_array* pkpm_div_ppar, const struct gkyl_array* rho_inv, const struct gkyl_array* T_perp_over_m, const struct gkyl_array* nu, 
  struct gkyl_array* div_p, struct gkyl_array* pkpm_accel_vars)
{
  int cdim = basis.ndim;
  int poly_order = basis.poly_order;
  double dx[GKYL_MAX_DIM] = {0.0};

  pkpm_recovery_t pkpm_recovery[3];
  // Fetch the kernels in each direction
  for (int d=0; d<cdim; ++d) {
    pkpm_recovery[d] = choose_ser_pkpm_recovery_kern(d, cdim, poly_order);
    dx[d] = grid.dx[d];
  }
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idxc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&range, idxc);

    const double *bvar_c = (const double*) gkyl_array_cfetch(bvar, linc);
    const double *u_i_c = (const double*) gkyl_array_cfetch(u_i, linc);
    const double *p_ij_c = (const double*) gkyl_array_cfetch(p_ij, linc);

    // Only need rho_inv, T_perp_over_m, T_perp_over_m_inv, nu, and nu_vthsq in center cell
    const double *pkpm_div_ppar_d = (const double *) gkyl_array_cfetch(pkpm_div_ppar, linc);
    const double *rho_inv_d = (const double*) gkyl_array_cfetch(rho_inv, linc);
    const double *T_perp_over_m_d = (const double*) gkyl_array_cfetch(T_perp_over_m, linc);
    const double *nu_d = (const double*) gkyl_array_cfetch(nu, linc);

    double *div_p_d = (double*) gkyl_array_fetch(div_p, linc);
    double *pkpm_accel_vars_d = (double*) gkyl_array_fetch(pkpm_accel_vars, linc);

    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, idxc, idxl);
      gkyl_copy_int_arr(cdim, idxc, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(&range, idxl); 
      long linr = gkyl_range_idx(&range, idxr);

      const double *bvar_l = (const double*) gkyl_array_cfetch(bvar, linl);
      const double *bvar_r = (const double*) gkyl_array_cfetch(bvar, linr);

      const double *u_i_l = (const double*) gkyl_array_cfetch(u_i, linl);
      const double *u_i_r = (const double*) gkyl_array_cfetch(u_i, linr);

      const double *p_ij_l = (const double*) gkyl_array_cfetch(p_ij, linl);
      const double *p_ij_r = (const double*) gkyl_array_cfetch(p_ij, linr);

      const double *euler_pkpm_l = (const double*) gkyl_array_cfetch(euler_pkpm, linl);
      const double *euler_pkpm_r = (const double*) gkyl_array_cfetch(euler_pkpm, linr);

      pkpm_recovery[dir](dx, 
        bvar_l, bvar_c, bvar_r, u_i_l, u_i_c, u_i_r, 
        p_ij_l, p_ij_c, p_ij_r, 
        pkpm_div_ppar_d, rho_inv_d, T_perp_over_m_d, nu_d, 
        div_p_d, pkpm_accel_vars_d);
    }
  }
}

// Host-side wrapper for pkpm derivative calculations with recovery
void
gkyl_calc_pkpm_vars_recovery_cu(const struct gkyl_rect_grid *grid, 
  struct gkyl_basis basis, const struct gkyl_range *range, 
  const struct gkyl_array* bvar, const struct gkyl_array* u_i, const struct gkyl_array* p_ij, 
  const struct gkyl_array* pkpm_div_ppar, const struct gkyl_array* rho_inv, const struct gkyl_array* T_perp_over_m, const struct gkyl_array* nu, 
  struct gkyl_array* div_p, struct gkyl_array* pkpm_accel_vars)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_calc_pkpm_vars_recovery_cu_kernel<<<nblocks, nthreads>>>(*grid, basis, *range, 
    bvar->on_dev, u_i->on_dev, p_ij->on_dev, 
    pkpm_div_ppar->on_dev, rho_inv->on_dev, T_perp_over_m->on_dev, nu->on_dev, 
    div_p->on_dev, pkpm_accel_vars->on_dev);
}

__global__ void
gkyl_calc_pkpm_vars_pressure_cu_kernel(struct gkyl_rect_grid grid, struct gkyl_basis basis, 
  struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array* bvar, const struct gkyl_array* fin, 
  struct gkyl_array* pkpm_div_ppar)
{
  double dx[GKYL_MAX_DIM] = {0.0};
  double xc[GKYL_MAX_DIM] = {0.0};
  int cdim = basis.ndim;
  int poly_order = basis.poly_order;
  for (int d=0; d<cdim+1; ++d) 
    dx[d] = grid.dx[d];

  pkpm_pressure_t pkpm_pressure[3];
  // Fetch the kernels in each direction
  for (int d=0; d<cdim; ++d) {
    switch (basis.b_type) {
      case GKYL_BASIS_MODAL_SERENDIPITY:
        pkpm_pressure[d] = choose_ser_pkpm_pressure_kern(d, cdim, poly_order);

        break;

      case GKYL_BASIS_MODAL_TENSOR:
        pkpm_pressure[d] = choose_ten_pkpm_pressure_kern(d, cdim, poly_order);
        
        break;

      default:
        assert(false);
        break;    
    }
  }

  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < phase_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&phase_range, linc1, idxc);
    gkyl_rect_grid_cell_center(&grid, idxc, xc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc_conf = gkyl_range_idx(&conf_range, idxc);
    long linc_phase = gkyl_range_idx(&phase_range, idxc);

    const double *bvar_c = (const double*) gkyl_array_cfetch(bvar, linc_conf);
    const double *f_c = (const double*) gkyl_array_cfetch(fin, linc_phase);

    double momLocal[96]; // hard-coded to 3 * max confBasis.num_basis (3x p=3 Ser) for now.
    for (unsigned int k=0; k<96; ++k)
      momLocal[k] = 0.0;

    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim+1, idxc, idxl);
      gkyl_copy_int_arr(cdim+1, idxc, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl_conf = gkyl_range_idx(&conf_range, idxl); 
      long linl_phase = gkyl_range_idx(&phase_range, idxl); 
      long linr_conf = gkyl_range_idx(&conf_range, idxr); 
      long linr_phase = gkyl_range_idx(&phase_range, idxr); 

      const double *bvar_l = (const double*) gkyl_array_cfetch(bvar, linl_conf);
      const double *f_l = (const double*) gkyl_array_cfetch(fin, linl_phase);
      const double *bvar_r = (const double*) gkyl_array_cfetch(bvar, linr_conf);
      const double *f_r = (const double*) gkyl_array_cfetch(fin, linr_phase);

      pkpm_pressure[dir](xc, dx, 
        bvar_l, bvar_c, bvar_r, f_l, f_c, f_r, 
        &momLocal[0]);
    }
    // Accumulate output to output array atomically to avoid race conditions
    double *pkpm_div_ppar_d = (double*) gkyl_array_fetch(pkpm_div_ppar, linc_conf);
    for (unsigned int k = 0; k < pkpm_div_ppar->ncomp; ++k) {
       if (linc1 < phase_range.volume)
         atomicAdd(&pkpm_div_ppar_d[k], momLocal[k]);
    }    
  }  
}
// Host-side wrapper for pkpm div(p_parallel b_hat) calculation
void 
gkyl_calc_pkpm_vars_pressure_cu(const struct gkyl_rect_grid *grid, struct gkyl_basis basis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_array* bvar, const struct gkyl_array* fin, 
  struct gkyl_array* pkpm_div_ppar)
{
  int nblocks = phase_range->nblocks;
  int nthreads = phase_range->nthreads;
  gkyl_calc_pkpm_vars_pressure_cu_kernel<<<nblocks, nthreads>>>(*grid, basis, 
    *conf_range, *phase_range, 
    bvar->on_dev, fin->on_dev, 
    pkpm_div_ppar->on_dev);
}
