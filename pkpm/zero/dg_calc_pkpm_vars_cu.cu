/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_pkpm_vars.h>
#include <gkyl_dg_calc_pkpm_vars_priv.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_dg_calc_pkpm_vars_u_set_cu_kernel(
  gkyl_dg_calc_pkpm_vars *up, struct gkyl_nmat *As, struct gkyl_nmat *xs,
  struct gkyl_range conf_range, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *euler_pkpm, struct gkyl_array *cell_avg_prim
)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1 * 3;

    const double *vlasov_pkpm_moms_d = (const double *)gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = (const double *)gkyl_array_cfetch(euler_pkpm, loc);

    int *cell_avg_prim_d = (int *)gkyl_array_fetch(cell_avg_prim, loc);

    cell_avg_prim_d[0] = up->pkpm_u_set(count, As, xs, vlasov_pkpm_moms_d, euler_pkpm_d);
  }
}

__global__ static void
gkyl_dg_calc_pkpm_vars_u_copy_cu_kernel(
  gkyl_dg_calc_pkpm_vars *up, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  struct gkyl_array *pkpm_u
)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1 * 3;

    double *pkpm_u_d = (double *)gkyl_array_fetch(pkpm_u, loc);

    up->pkpm_u_copy(count, xs, pkpm_u_d);
  }
}

// Host-side wrapper for pkpm flow velocity calculation
void
gkyl_dg_calc_pkpm_vars_u_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *euler_pkpm, struct gkyl_array *cell_avg_prim, struct gkyl_array *pkpm_u
)
{
  struct gkyl_range conf_range = up->mem_range;

  gkyl_dg_calc_pkpm_vars_u_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(
    up->on_dev, up->As_u->on_dev, up->xs_u->on_dev, conf_range, vlasov_pkpm_moms->on_dev,
    euler_pkpm->on_dev, cell_avg_prim->on_dev
  );

  gkyl_dg_calc_pkpm_vars_u_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(
    up->on_dev, up->xs_u->on_dev, conf_range, pkpm_u->on_dev
  );
}

__global__ void
gkyl_dg_calc_pkpm_vars_closure_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  double w_lo, double w_up, double dv, const struct gkyl_array *norm_moms, const struct gkyl_array *rhs_moms,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *grad_u,
  const struct gkyl_array *xi, const struct gkyl_array *grad_xi, const struct gkyl_array *theta,
  const struct gkyl_array *fin, struct gkyl_array *c_0, struct gkyl_array *alpha_0, struct gkyl_array *c_h,
  struct gkyl_array *vel_flux_surf, struct gkyl_array *cflrate, int *num_null
)
{
  int cdim = conf_range.ndim;
  int idx[GKYL_MAX_DIM];
  double cfl[2];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    idx[cdim] = phase_range.lower[cdim];
    long plo = gkyl_range_idx(&phase_range, idx);
    idx[cdim] = phase_range.upper[cdim];
    long pup = gkyl_range_idx(&phase_range, idx);
    idx[cdim] = phase_range.upper[cdim] + 1; // velocity ghost cell: its slot is the upper edge face
    long pgh = gkyl_range_idx(&phase_range, idx);
    int null = up->pkpm_closure(
      w_lo, w_up, dv, (const double *)gkyl_array_cfetch(norm_moms, loc), (const double *)gkyl_array_cfetch(rhs_moms, loc),
      (const double *)gkyl_array_cfetch(bvar, loc), (const double *)gkyl_array_cfetch(grad_b, loc),
      (const double *)gkyl_array_cfetch(grad_u, loc), (const double *)gkyl_array_cfetch(xi, loc),
      (const double *)gkyl_array_cfetch(grad_xi, loc), (const double *)gkyl_array_cfetch(theta, loc),
      (const double *)gkyl_array_cfetch(alpha_0, loc), (const double *)gkyl_array_cfetch(c_h, loc),
      (const double *)gkyl_array_cfetch(fin, plo), (const double *)gkyl_array_cfetch(fin, pup),
      (double *)gkyl_array_fetch(c_0, loc), (double *)gkyl_array_fetch(alpha_0, loc), (double *)gkyl_array_fetch(c_h, loc),
      (double *)gkyl_array_fetch(vel_flux_surf, plo), (double *)gkyl_array_fetch(vel_flux_surf, pgh), cfl
    );
    double *cfl_lo = (double *)gkyl_array_fetch(cflrate, plo), *cfl_up = (double *)gkyl_array_fetch(cflrate, pup);
    cfl_lo[0] += cfl[0];
    cfl_up[0] += cfl[1];
    if (null > 0) {
      atomicAdd(num_null, null);
    }
  }
}

__global__ void
gkyl_dg_calc_pkpm_vars_closure_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  double w_lo, double w_up, double dv, const struct gkyl_array *norm_moms, const struct gkyl_array *rhs_moms,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *grad_u,
  const struct gkyl_array *xi, const struct gkyl_array *grad_xi, const struct gkyl_array *theta,
  const struct gkyl_array *fin, struct gkyl_array *c_0, struct gkyl_array *alpha_0, struct gkyl_array *c_h,
  struct gkyl_array *vel_flux_surf, struct gkyl_array *cflrate, int *num_null
)
{
  int cdim = conf_range.ndim;
  int idx[GKYL_MAX_DIM];
  double cfl[2];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    idx[cdim] = phase_range.lower[cdim];
    long plo = gkyl_range_idx(&phase_range, idx);
    idx[cdim] = phase_range.upper[cdim];
    long pup = gkyl_range_idx(&phase_range, idx);
    idx[cdim] = phase_range.upper[cdim] + 1; // velocity ghost cell: its slot is the upper edge face
    long pgh = gkyl_range_idx(&phase_range, idx);
    int null = up->pkpm_closure(
      w_lo, w_up, dv, (const double *)gkyl_array_cfetch(norm_moms, loc), (const double *)gkyl_array_cfetch(rhs_moms, loc),
      (const double *)gkyl_array_cfetch(bvar, loc), (const double *)gkyl_array_cfetch(grad_b, loc),
      (const double *)gkyl_array_cfetch(grad_u, loc), (const double *)gkyl_array_cfetch(xi, loc),
      (const double *)gkyl_array_cfetch(grad_xi, loc), (const double *)gkyl_array_cfetch(theta, loc),
      (const double *)gkyl_array_cfetch(alpha_0, loc), (const double *)gkyl_array_cfetch(c_h, loc),
      (const double *)gkyl_array_cfetch(fin, plo), (const double *)gkyl_array_cfetch(fin, pup),
      (double *)gkyl_array_fetch(c_0, loc), (double *)gkyl_array_fetch(alpha_0, loc), (double *)gkyl_array_fetch(c_h, loc),
      (double *)gkyl_array_fetch(vel_flux_surf, plo), (double *)gkyl_array_fetch(vel_flux_surf, pgh), cfl
    );
    double *cfl_lo = (double *)gkyl_array_fetch(cflrate, plo), *cfl_up = (double *)gkyl_array_fetch(cflrate, pup);
    cfl_lo[0] += cfl[0];
    cfl_up[0] += cfl[1];
    if (null > 0) {
      atomicAdd(num_null, null);
    }
  }
}
// Host-side wrapper for the closures
int
gkyl_dg_calc_pkpm_vars_closure_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, const struct gkyl_range *phase_range,
  const struct gkyl_rect_grid *phase_grid, const struct gkyl_array *norm_moms, const struct gkyl_array *rhs_moms,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *grad_u,
  const struct gkyl_array *xi, const struct gkyl_array *grad_xi, const struct gkyl_array *theta,
  const struct gkyl_array *fin, struct gkyl_array *c_0, struct gkyl_array *alpha_0, struct gkyl_array *c_h,
  struct gkyl_array *vel_flux_surf, struct gkyl_array *cflrate
)
{
  int cdim = conf_range->ndim;
  gkyl_cu_memset(up->num_null_rho_dev, 0, sizeof(int));
  gkyl_dg_calc_pkpm_vars_closure_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, phase_grid->lower[cdim], phase_grid->upper[cdim], phase_grid->dx[cdim],
    norm_moms->on_dev, rhs_moms->on_dev, bvar->on_dev, grad_b->on_dev, grad_u->on_dev, xi->on_dev, grad_xi->on_dev,
    theta->on_dev, fin->on_dev, c_0->on_dev, alpha_0->on_dev, c_h->on_dev, vel_flux_surf->on_dev,
    cflrate->on_dev, up->num_null_rho_dev
  );
  int num_null = 0;
  gkyl_cu_memcpy(&num_null, up->num_null_rho_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  return num_null;
}

__global__ void
gkyl_dg_calc_pkpm_vars_prim_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, const struct gkyl_array *fluid,
  struct gkyl_array *u, struct gkyl_array *xi, int *num_null
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    int null = up->pkpm_prim(
      (const double *)gkyl_array_cfetch(fluid, loc), (double *)gkyl_array_fetch(u, loc), (double *)gkyl_array_fetch(xi, loc)
    );
    if (null > 0) {
      atomicAdd(num_null, null);
    }
  }
}

__global__ void
gkyl_dg_calc_pkpm_vars_prim_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, const struct gkyl_array *fluid,
  struct gkyl_array *u, struct gkyl_array *xi, int *num_null
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    int null = up->pkpm_prim(
      (const double *)gkyl_array_cfetch(fluid, loc), (double *)gkyl_array_fetch(u, loc), (double *)gkyl_array_fetch(xi, loc)
    );
    if (null > 0) {
      atomicAdd(num_null, null);
    }
  }
}
// Host-side wrapper for the nodal primitive variables
int
gkyl_dg_calc_pkpm_vars_prim_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, struct gkyl_array *u, struct gkyl_array *xi
)
{
  gkyl_cu_memset(up->num_null_rho_dev, 0, sizeof(int));
  gkyl_dg_calc_pkpm_vars_prim_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, fluid->on_dev, u->on_dev, xi->on_dev, up->num_null_rho_dev
  );
  int num_null = 0;
  gkyl_cu_memcpy(&num_null, up->num_null_rho_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  return num_null;
}

__global__ void
gkyl_dg_calc_pkpm_vars_int_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, const struct gkyl_array *fluid,
  const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *u, struct gkyl_array *int_pkpm_vars
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    up->pkpm_int(
      (const double *)gkyl_array_cfetch(fluid, loc), (const double *)gkyl_array_cfetch(vlasov_pkpm_moms, loc),
      (const double *)gkyl_array_cfetch(u, loc), (double *)gkyl_array_fetch(int_pkpm_vars, loc)
    );
  }
}

__global__ void
gkyl_dg_calc_pkpm_vars_int_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, const struct gkyl_array *fluid,
  const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *u, struct gkyl_array *int_pkpm_vars
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    up->pkpm_int(
      (const double *)gkyl_array_cfetch(fluid, loc), (const double *)gkyl_array_cfetch(vlasov_pkpm_moms, loc),
      (const double *)gkyl_array_cfetch(u, loc), (double *)gkyl_array_fetch(int_pkpm_vars, loc)
    );
  }
}
// Host-side wrapper for the integrated diagnostics
void
gkyl_dg_calc_pkpm_vars_int_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *fluid, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *u, struct gkyl_array *int_pkpm_vars
)
{
  gkyl_dg_calc_pkpm_vars_int_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, fluid->on_dev, vlasov_pkpm_moms->on_dev, u->on_dev, int_pkpm_vars->on_dev
  );
}

__global__ void
gkyl_dg_calc_pkpm_vars_grad_ln_n_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, const struct gkyl_array *rho,
  const struct gkyl_array *grad_rho, struct gkyl_array *grad_ln_n, int *num_null
)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);

    int null = up->pkpm_grad_ln_n(
      (const double *)gkyl_array_cfetch(rho, loc), (const double *)gkyl_array_cfetch(grad_rho, loc),
      (double *)gkyl_array_fetch(grad_ln_n, loc)
    );
    if (null > 0) {
      atomicAdd(num_null, null);
    }
  }
}

    int null = up->pkpm_grad_ln_n(
      (const double *)gkyl_array_cfetch(rho, loc), (const double *)gkyl_array_cfetch(grad_rho, loc),
      (double *)gkyl_array_fetch(grad_ln_n, loc)
    );
    if (null > 0) {
      atomicAdd(num_null, null);
    }
  }
}
// Host-side wrapper for grad(ln n)
int
gkyl_dg_calc_pkpm_vars_grad_ln_n_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *rho, const struct gkyl_array *grad_rho, struct gkyl_array *grad_ln_n
)
{
  gkyl_cu_memset(up->num_null_rho_dev, 0, sizeof(int));
  gkyl_dg_calc_pkpm_vars_grad_ln_n_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, rho->on_dev, grad_rho->on_dev, grad_ln_n->on_dev, up->num_null_rho_dev
  );
  int num_null = 0;
  gkyl_cu_memcpy(&num_null, up->num_null_rho_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  return num_null;
}

__global__ void
gkyl_dg_calc_pkpm_vars_normalize_dist_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array *n, const struct gkyl_array *f, struct gkyl_array *fbar, int *num_null
)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < phase_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    long pidx = gkyl_range_idx(&phase_range, idx);

    int null = up->pkpm_normalize(
      (const double *)gkyl_array_cfetch(n, cidx), (const double *)gkyl_array_cfetch(f, pidx),
      (double *)gkyl_array_fetch(fbar, pidx)
    );
    if (null > 0) {
      atomicAdd(num_null, null);
    }
  }
}

    int null = up->pkpm_normalize(
      (const double *)gkyl_array_cfetch(n, cidx), (const double *)gkyl_array_cfetch(f, pidx),
      (double *)gkyl_array_fetch(fbar, pidx)
    );
    if (null > 0) {
      atomicAdd(num_null, null);
    }
  }
}
// Host-side wrapper for the weak division of the distributions by n
int
gkyl_dg_calc_pkpm_vars_normalize_dist_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *n, const struct gkyl_array *f,
  struct gkyl_array *fbar
)
{
  gkyl_cu_memset(up->num_null_rho_dev, 0, sizeof(int));
  gkyl_dg_calc_pkpm_vars_normalize_dist_cu_kernel<<<phase_range->nblocks, phase_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, n->on_dev, f->on_dev, fbar->on_dev, up->num_null_rho_dev
  );
  int num_null = 0;
  gkyl_cu_memcpy(&num_null, up->num_null_rho_dev, sizeof(int), GKYL_CU_MEMCPY_D2H);
  return num_null;
}

__global__ void
gkyl_dg_calc_pkpm_vars_diag_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range,
  const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *euler_pkpm,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *theta, const struct gkyl_array *alpha_0,
  const struct gkyl_array *xi, const struct gkyl_array *c_h, struct gkyl_array *fluid_io,
  struct gkyl_array *pkpm_vars_io
)
{
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    up->pkpm_diag(
      up->conf_grid.dx, (const double *)gkyl_array_cfetch(vlasov_pkpm_moms, loc),
      (const double *)gkyl_array_cfetch(euler_pkpm, loc), (const double *)gkyl_array_cfetch(bvar, loc),
      (const double *)gkyl_array_cfetch(grad_b, loc),
      (const double *)gkyl_array_cfetch(u, loc), (const double *)gkyl_array_cfetch(grad_u, loc),
      (const double *)gkyl_array_cfetch(theta, loc), (const double *)gkyl_array_cfetch(alpha_0, loc),
      (const double *)gkyl_array_cfetch(xi, loc), (const double *)gkyl_array_cfetch(c_h, loc),
      (double *)gkyl_array_fetch(fluid_io, loc), (double *)gkyl_array_fetch(pkpm_vars_io, loc)
    );
  }
}

void
gkyl_dg_calc_pkpm_vars_diag_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *euler_pkpm,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *theta, const struct gkyl_array *alpha_0,
  const struct gkyl_array *xi, const struct gkyl_array *c_h, struct gkyl_array *fluid_io,
  struct gkyl_array *pkpm_vars_io
)
{
  gkyl_dg_calc_pkpm_vars_diag_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, vlasov_pkpm_moms->on_dev, euler_pkpm->on_dev, bvar->on_dev,
    grad_b->on_dev, u->on_dev, grad_u->on_dev, theta->on_dev, alpha_0->on_dev, xi->on_dev, c_h->on_dev,
    fluid_io->on_dev, pkpm_vars_io->on_dev
  );
}

// One thread per interior cell: the weak gradient of u in every direction from the
// cell and its two neighbors.
__global__ void
gkyl_dg_calc_pkpm_vars_grad_u_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, const struct gkyl_array *u,
  struct gkyl_array *grad_u
)
{
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long loc = gkyl_range_idx(&conf_range, idx);
    const double *u_c = (const double *)gkyl_array_cfetch(u, loc);
    double *gu = (double *)gkyl_array_fetch(grad_u, loc);
    for (int d = 0; d < cdim; ++d) {
      gkyl_copy_int_arr(cdim, idx, idx_l);
      gkyl_copy_int_arr(cdim, idx, idx_r);
      idx_l[d] = idx_l[d] - 1;
      idx_r[d] = idx_r[d] + 1;
      up->pkpm_grad_u[d](
        up->conf_grid.dx, (const double *)gkyl_array_cfetch(up->jacob_pos, loc),
        (const double *)gkyl_array_cfetch(u, gkyl_range_idx(&conf_range, idx_l)),
        u_c, (const double *)gkyl_array_cfetch(u, gkyl_range_idx(&conf_range, idx_r)), gu
      );
    }
  }
}

void
gkyl_dg_calc_pkpm_vars_grad_u_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *u, struct gkyl_array *grad_u
)
{
  gkyl_dg_calc_pkpm_vars_grad_u_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, u->on_dev, grad_u->on_dev
  );
}

__global__ void
gkyl_dg_calc_pkpm_vars_source_cu_kernel(
  struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, const struct gkyl_array *qmem,
  const struct gkyl_array *fluid, struct gkyl_array *rhs
)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);

    up->pkpm_source(
      (const double *)gkyl_array_cfetch(qmem, loc), (const double *)gkyl_array_cfetch(fluid, loc),
      (double *)gkyl_array_fetch(rhs, loc)
    );
  }
}

// Host-side wrapper for pkpm source term calculations
void
gkyl_dg_calc_pkpm_vars_source_cu(
  struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *qmem, const struct gkyl_array *fluid, struct gkyl_array *rhs
)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_pkpm_vars_source_cu_kernel<<<nblocks, nthreads>>>(
    up->on_dev, *conf_range, qmem->on_dev, fluid->on_dev, rhs->on_dev
  );
}

// CUDA kernel to set device pointers to pkpm vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void
dg_calc_pkpm_vars_set_cu_dev_ptrs(
  struct gkyl_dg_calc_pkpm_vars *up, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  up->pkpm_u_set = choose_pkpm_u_set_kern(b_type, cdim, poly_order);
  up->pkpm_u_copy = choose_pkpm_u_copy_kern(b_type, cdim, poly_order);
  for (int d = 0; d < cdim; ++d) {
    up->pkpm_grad_u[d] = choose_pkpm_grad_u_kern(b_type, cdim, poly_order, d);
  }
  up->pkpm_closure = choose_pkpm_closure_kern(b_type, cdim, poly_order);
  up->pkpm_prim = choose_pkpm_prim_kern(b_type, cdim, poly_order);
  up->pkpm_int = choose_pkpm_int_kern(b_type, cdim, poly_order);
  up->pkpm_grad_ln_n = choose_pkpm_grad_ln_n_kern(b_type, cdim, poly_order);
  up->pkpm_normalize = choose_pkpm_normalize_kern(b_type, cdim, poly_order);
  up->pkpm_source = choose_pkpm_source_kern(b_type, cdim, poly_order);
  up->pkpm_diag = choose_pkpm_diag_kern(b_type, cdim, poly_order);
}

gkyl_dg_calc_pkpm_vars *
gkyl_dg_calc_pkpm_vars_cu_dev_new(
  const struct gkyl_rect_grid *conf_grid, const struct gkyl_basis *cbasis,
  const struct gkyl_range *mem_range, const struct gkyl_vlasov_position_map *pos_map
)
{
  struct gkyl_dg_calc_pkpm_vars *up =
    (struct gkyl_dg_calc_pkpm_vars *)gkyl_malloc(sizeof(gkyl_dg_calc_pkpm_vars));

  up->conf_grid = *conf_grid;
  int nc = cbasis->num_basis;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->cdim = cdim;
  up->poly_order = poly_order;
  up->mem_range = *mem_range;
  assert(pos_map);
  up->pos_map = 0; // host pointer set after the device memcpy below
  up->jacob_pos = pos_map->jacob_pos->on_dev;

  // Linear system for the weak division of the momentum by the density (ux, uy, uz)
  up->As_u = gkyl_nmat_cu_dev_new(3 * mem_range->volume, nc, nc);
  up->xs_u = gkyl_nmat_cu_dev_new(3 * mem_range->volume, nc, 1);
  up->mem_u = gkyl_nmat_linsolve_lu_cu_dev_new(up->As_u->num, up->As_u->nr);

  up->num_null_rho = 0;
  up->num_null_rho_dev = (int *)gkyl_cu_malloc(sizeof(int));

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_pkpm_vars *up_cu =
    (struct gkyl_dg_calc_pkpm_vars *)gkyl_cu_malloc(sizeof(gkyl_dg_calc_pkpm_vars));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_pkpm_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_pkpm_vars_set_cu_dev_ptrs<<<1, 1>>>(up_cu, b_type, cdim, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;
  up->pos_map = gkyl_vlasov_position_map_acquire(pos_map);
  up->jacob_pos = pos_map->jacob_pos;

  return up;
}
