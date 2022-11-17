/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_prim_vars.h>
#include <gkyl_dg_calc_prim_vars_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_calc_prim_vars_pkpm_cu_kernel(struct gkyl_basis basis, struct gkyl_range range, 
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* u_i, struct gkyl_array* u_perp_i, struct gkyl_array* rhou_perp_i,
  struct gkyl_array* p_perp, struct gkyl_array* p_ij)
{
  int cdim = basis.ndim;
  int poly_order = basis.poly_order;

  euler_pkpm_prim_vars_t pkpm_prim_vars = choose_ser_euler_pkpm_prim_vars_kern(cdim, poly_order);

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
    double *u_perp_i_d = (double*) gkyl_array_fetch(u_perp_i, start);
    double *rhou_perp_i_d = (double*) gkyl_array_fetch(rhou_perp_i, start);
    double *p_perp_d = (double*) gkyl_array_fetch(p_perp, start);
    double *p_ij_d = (double*) gkyl_array_fetch(p_ij, start);

    pkpm_prim_vars(bvar_d, vlasov_pkpm_moms_d, euler_pkpm_d, 
      u_i_d, u_perp_i_d, rhou_perp_i_d, p_perp_d, p_ij_d);
  }
}

// Host-side wrapper for pkpm primitive variable calculations
void
gkyl_calc_prim_vars_pkpm_cu(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* u_i, struct gkyl_array* u_perp_i, struct gkyl_array* rhou_perp_i,
  struct gkyl_array* p_perp, struct gkyl_array* p_ij)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_calc_prim_vars_pkpm_cu_kernel<<<nblocks, nthreads>>>(basis, *range, 
    bvar->on_dev, vlasov_pkpm_moms->on_dev, euler_pkpm->on_dev, 
    u_i->on_dev, u_perp_i->on_dev, rhou_perp_i->on_dev, p_perp->on_dev, p_ij->on_dev);
}

__global__ void
gkyl_calc_prim_vars_pkpm_source_cu_kernel(struct gkyl_basis basis, struct gkyl_range range, 
  const struct gkyl_array* qmem, const struct gkyl_array* nu, const struct gkyl_array* nu_vthsq,
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  const struct gkyl_array* rhou_perp_i,  const struct gkyl_array* p_perp, 
  struct gkyl_array* rhs)
{
  int cdim = basis.ndim;
  int poly_order = basis.poly_order;

  euler_pkpm_source_t pkpm_source = choose_ser_euler_pkpm_source_kern(cdim, poly_order);

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
    const double *nu_d = (const double*) gkyl_array_cfetch(nu, start);
    const double *nu_vthsq_d = (const double*) gkyl_array_cfetch(nu_vthsq, start);
    const double *vlasov_pkpm_moms_d = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, start);
    const double *euler_pkpm_d = (const double*) gkyl_array_cfetch(euler_pkpm, start);
    const double *rhou_perp_i_d = (const double*) gkyl_array_cfetch(rhou_perp_i, start);
    const double *p_perp_d = (const double*) gkyl_array_cfetch(p_perp, start);

    double *rhs_d = (double*) gkyl_array_fetch(rhs, start);
    pkpm_source(qmem_d, nu_d, nu_vthsq_d, vlasov_pkpm_moms_d, euler_pkpm_d, rhou_perp_i_d, p_perp_d, rhs_d);
  }
}

// Host-side wrapper for pkpm source term calculations
void
gkyl_calc_prim_vars_pkpm_source_cu(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* qmem, const struct gkyl_array* nu, const struct gkyl_array* nu_vthsq,
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  const struct gkyl_array* rhou_perp_i,  const struct gkyl_array* p_perp, 
  struct gkyl_array* rhs)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_calc_prim_vars_pkpm_source_cu_kernel<<<nblocks, nthreads>>>(basis, *range, 
    qmem->on_dev, nu->on_dev, nu_vthsq->on_dev, vlasov_pkpm_moms->on_dev, euler_pkpm->on_dev, 
    rhou_perp_i->on_dev, p_perp->on_dev, rhs->on_dev);
}

__global__ void
gkyl_calc_prim_vars_pkpm_recovery_cu_kernel(double dx, struct gkyl_basis basis, struct gkyl_range range, 
  const struct gkyl_array* bvar, const struct gkyl_array* u_i, const struct gkyl_array* p_ij, 
  struct gkyl_array* div_b, struct gkyl_array* bb_grad_u, struct gkyl_array* div_p)
{
  int cdim = basis.ndim;
  int poly_order = basis.poly_order;

  euler_pkpm_recovery_t pkpm_recovery[3];
  // Fetch the kernels in each direction
  for (int d=0; d<cdim; ++d) 
    pkpm_recovery[d] = choose_ser_euler_pkpm_recovery_kern(d, cdim, poly_order);

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

    double *div_b_d = (double*) gkyl_array_fetch(div_b, linc);
    double *bb_grad_u_d = (double*) gkyl_array_fetch(bb_grad_u, linc);
    double *div_p_d = (double*) gkyl_array_fetch(div_p, linc);

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

      pkpm_recovery[dir](&dx, bvar_l, bvar_c, bvar_r, u_i_l, u_i_c, u_i_r, p_ij_l, p_ij_c, p_ij_r, 
        div_b_d, bb_grad_u_d, div_p_d);
    }
  }
}

// Host-side wrapper for pkpm derivative calculations with recovery
void
gkyl_calc_prim_vars_pkpm_recovery_cu(const double *dx, 
  struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* u_i, const struct gkyl_array* p_ij, 
  struct gkyl_array* div_b, struct gkyl_array* bb_grad_u, struct gkyl_array* div_p)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_calc_prim_vars_pkpm_recovery_cu_kernel<<<nblocks, nthreads>>>(*dx, basis, *range, 
    bvar->on_dev, u_i->on_dev, p_ij->on_dev, 
    div_b->on_dev, bb_grad_u->on_dev, div_p->on_dev);
}

__global__ void
gkyl_calc_prim_vars_pkpm_p_force_cu_kernel(struct gkyl_basis basis, struct gkyl_range range, 
  const struct gkyl_array* bvar, const struct gkyl_array* div_p, const struct gkyl_array* vlasov_pkpm_moms, 
  const struct gkyl_array* euler_pkpm, const struct gkyl_array* div_b, struct gkyl_array* p_force)
{
  int cdim = basis.ndim;
  int poly_order = basis.poly_order;

  euler_pkpm_p_force_t pkpm_p_force = choose_ser_euler_pkpm_p_force_kern(cdim, poly_order);

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
    const double *div_p_d = (const double*) gkyl_array_cfetch(div_p, start);
    const double *vlasov_pkpm_moms_d = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, start);
    const double *euler_pkpm_d = (const double*) gkyl_array_cfetch(euler_pkpm, start);
    const double *div_b_d = (const double*) gkyl_array_cfetch(div_b, start);
    
    double *p_force_d = (double*) gkyl_array_fetch(p_force, start);
    pkpm_p_force(bvar_d, div_p_d, vlasov_pkpm_moms_d, euler_pkpm_d, div_b_d, p_force_d);
  }
}

// Host-side wrapper for pkpm primitive variable calculations
void
gkyl_calc_prim_vars_pkpm_p_force_cu(struct gkyl_basis basis, const struct gkyl_range *range,
  const struct gkyl_array* bvar, const struct gkyl_array* div_p, const struct gkyl_array* vlasov_pkpm_moms, 
  const struct gkyl_array* euler_pkpm, const struct gkyl_array* div_b, struct gkyl_array* p_force)
{
  int nblocks = range->nblocks;
  int nthreads = range->nthreads;
  gkyl_calc_prim_vars_pkpm_p_force_cu_kernel<<<nblocks, nthreads>>>(basis, *range, 
    bvar->on_dev, div_p->on_dev, vlasov_pkpm_moms->on_dev, 
    euler_pkpm->on_dev, div_b->on_dev, p_force->on_dev);
}
