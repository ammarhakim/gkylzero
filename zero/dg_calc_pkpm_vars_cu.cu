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
#include <gkyl_wave_geom.h>
#include <gkyl_wv_eqn.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_dg_calc_pkpm_vars_set_cu_kernel(gkyl_dg_calc_pkpm_vars* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* p_ij, 
  const struct gkyl_array* pkpm_div_ppar, const struct gkyl_array* div_b, 
  struct gkyl_array* cell_avg_prim)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp_prim = 6 systems in each cell)
    long count = linc1*up->Ncomp_prim;

    const double *vlasov_pkpm_moms_d = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *p_ij_d = (const double*) gkyl_array_cfetch(p_ij, loc);
    const double *pkpm_div_ppar_d = (const double*) gkyl_array_cfetch(pkpm_div_ppar, loc);
    const double *div_b_d = (const double*) gkyl_array_cfetch(div_b, loc);

    int* cell_avg_prim_d = (int*) gkyl_array_fetch(cell_avg_prim, loc);
    // First index of cell_avg_prim is whether p = p_par + 2 p_perp < 0.0 at control points
    cell_avg_prim_d[1] = up->pkpm_set(count, As, xs, 
      vlasov_pkpm_moms_d, p_ij_d, pkpm_div_ppar_d, div_b_d);
  }
}

__global__ static void
gkyl_dg_calc_pkpm_vars_copy_cu_kernel(gkyl_dg_calc_pkpm_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range,
  struct gkyl_array* prim, struct gkyl_array* pkpm_accel)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp_prim = 6 systems in each cell)
    long count = linc1*up->Ncomp_prim;

    double* prim_d = (double*) gkyl_array_fetch(prim, loc);
    double* pkpm_accel_d = (double*) gkyl_array_fetch(pkpm_accel, loc);

    up->pkpm_copy(count, xs, prim_d, pkpm_accel_d);
  }
}

// Host-side wrapper for pkpm primitive variable calculation
void gkyl_dg_calc_pkpm_vars_advance_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* p_ij, 
  const struct gkyl_array* pkpm_div_ppar, const struct gkyl_array* div_b, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* prim, struct gkyl_array* pkpm_accel)
{
  struct gkyl_range conf_range = up->mem_range;
  
  gkyl_dg_calc_pkpm_vars_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range,
    vlasov_pkpm_moms->on_dev, p_ij->on_dev, pkpm_div_ppar->on_dev, div_b->on_dev, 
    cell_avg_prim->on_dev);

  bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
  assert(status);

  gkyl_dg_calc_pkpm_vars_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range, prim->on_dev, pkpm_accel->on_dev);
}

__global__ static void
gkyl_dg_calc_pkpm_vars_surf_set_cu_kernel(gkyl_dg_calc_pkpm_vars* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* p_ij)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp_surf systems in each cell)
    long count = linc1*up->Ncomp_surf;

    const double *vlasov_pkpm_moms_d = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *p_ij_d = (const double*) gkyl_array_cfetch(p_ij, loc);

    up->pkpm_surf_set(count, As, xs, 
      vlasov_pkpm_moms_d, p_ij_d);
  }
}

__global__ static void
gkyl_dg_calc_pkpm_vars_surf_copy_cu_kernel(gkyl_dg_calc_pkpm_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range,
  struct gkyl_array* prim_surf)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp_surf systems in each cell)
    long count = linc1*up->Ncomp_surf;

    double* prim_surf_d = (double*) gkyl_array_fetch(prim_surf, loc);

    up->pkpm_surf_copy(count, xs, prim_surf_d);
  }
}

// Host-side wrapper for pkpm surface primitive variable calculation
void gkyl_dg_calc_pkpm_vars_surf_advance_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* p_ij, 
  struct gkyl_array* prim_surf)
{
  struct gkyl_range conf_range = up->mem_range;
  
  gkyl_dg_calc_pkpm_vars_surf_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As_surf->on_dev, up->xs_surf->on_dev, conf_range,
    vlasov_pkpm_moms->on_dev, p_ij->on_dev);

  if (up->cdim > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem_surf, up->As_surf, up->xs_surf);
    assert(status);
  }

  gkyl_dg_calc_pkpm_vars_surf_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs_surf->on_dev, conf_range, prim_surf->on_dev);
}

__global__ static void
gkyl_dg_calc_pkpm_vars_u_set_cu_kernel(gkyl_dg_calc_pkpm_vars* up,
  struct gkyl_nmat *As_u, struct gkyl_nmat *xs_u, struct gkyl_range conf_range,
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* cell_avg_prim)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp_u = 3 systems in each cell)
    long count = linc1*up->Ncomp_u;

    const double *vlasov_pkpm_moms_d = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = (const double*) gkyl_array_cfetch(euler_pkpm, loc);

    int* cell_avg_prim_d = (int*) gkyl_array_fetch(cell_avg_prim, loc);
    // Zeroth index of cell_avg_prim is whether rho < 0.0 at control points
    cell_avg_prim_d[0] = up->pkpm_u_set(count, As_u, xs_u, 
      vlasov_pkpm_moms_d, euler_pkpm_d);
  }
}

__global__ static void
gkyl_dg_calc_pkpm_vars_u_copy_cu_kernel(gkyl_dg_calc_pkpm_vars* up, 
  struct gkyl_nmat *xs_u, struct gkyl_range conf_range,
  struct gkyl_array* pkpm_u)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);
    // fetch the correct count in the matrix (since we solve Ncomp_u = 3 systems in each cell)
    long count = linc1*up->Ncomp_u;

    double* pkpm_u_d = (double*) gkyl_array_fetch(pkpm_u, loc);

    up->pkpm_u_copy(count, xs_u, pkpm_u_d);
  }
}

// Host-side wrapper for pkpm flow velocity calculation
void gkyl_dg_calc_pkpm_vars_u_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* pkpm_u)
{
  struct gkyl_range conf_range = up->mem_range;
  
  gkyl_dg_calc_pkpm_vars_u_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As_u->on_dev, up->xs_u->on_dev, conf_range,
    vlasov_pkpm_moms->on_dev, euler_pkpm->on_dev, 
    cell_avg_prim->on_dev);

  bool status = gkyl_nmat_linsolve_lu_pa(up->mem_u, up->As_u, up->xs_u);
  assert(status);

  gkyl_dg_calc_pkpm_vars_u_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs_u->on_dev, conf_range, pkpm_u->on_dev);
}

__global__ void
gkyl_calc_pkpm_vars_u_surf_cu_kernel(struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* pkpm_u, struct gkyl_array* pkpm_u_surf)
{ 
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);

    const double* pkpm_u_d = (const double*) gkyl_array_cfetch(pkpm_u, loc);
    double* pkpm_u_surf_d = (double*) gkyl_array_fetch(pkpm_u_surf, loc);

    up->pkpm_u_surf(pkpm_u_d, pkpm_u_surf_d);
  }
}

// Host-side wrapper for pkpm surface flow velocity computation
void gkyl_dg_calc_pkpm_vars_u_surf_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* pkpm_u, struct gkyl_array* pkpm_u_surf)
{
  struct gkyl_range conf_range = up->mem_range;
  gkyl_calc_pkpm_vars_u_surf_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev, conf_range, 
    pkpm_u->on_dev, pkpm_u_surf->on_dev);
}

__global__ void
gkyl_calc_pkpm_vars_pressure_cu_kernel(struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* bb, const struct gkyl_array* vlasov_pkpm_moms, struct gkyl_array* p_ij)
{ 
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);

    const double *bb_d = (const double*) gkyl_array_cfetch(bb, loc);
    const double *vlasov_pkpm_moms_d = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, loc);

    double *p_ij_d = (double*) gkyl_array_fetch(p_ij, loc);
    up->pkpm_pressure(bb_d, vlasov_pkpm_moms_d, p_ij_d);
  }
}

// Host-side wrapper for pkpm pressure calculation
void gkyl_dg_calc_pkpm_vars_pressure_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bb, const struct gkyl_array* vlasov_pkpm_moms, struct gkyl_array* p_ij)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_calc_pkpm_vars_pressure_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    bb->on_dev, vlasov_pkpm_moms->on_dev, p_ij->on_dev);
}

__global__ void
gkyl_dg_calc_pkpm_vars_accel_cu_kernel(struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* pkpm_u_surf, const struct gkyl_array* pkpm_u, 
  const struct gkyl_array* prim, const struct gkyl_array* bb, 
  const struct gkyl_array* div_b, const struct gkyl_array* nu, 
  struct gkyl_array* pkpm_lax, struct gkyl_array* pkpm_accel)
{
  int cdim = up->cdim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idxc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&conf_range, idxc);

    const double *pkpm_u_surf_c = (const double*) gkyl_array_cfetch(pkpm_u_surf, linc);
    const double *prim_c = (const double*) gkyl_array_cfetch(prim, linc);

    const double *pkpm_u_d = (const double*) gkyl_array_cfetch(pkpm_u, linc);
    const double *bb_d = (const double*) gkyl_array_cfetch(bb, linc);
    const double *div_b_d = (const double*) gkyl_array_cfetch(div_b, linc);
    const double *nu_d = (const double*) gkyl_array_cfetch(nu, linc);

    double *pkpm_lax_d = (double*) gkyl_array_fetch(pkpm_lax, linc);
    double *pkpm_accel_d = (double*) gkyl_array_fetch(pkpm_accel, linc);

    // Compute T_perp/m div(b) and p_force
    up->pkpm_p_force(prim_c, div_b_d, pkpm_accel_d);

    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, idxc, idxl);
      gkyl_copy_int_arr(cdim, idxc, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(&conf_range, idxl); 
      long linr = gkyl_range_idx(&conf_range, idxr);

      const double *pkpm_u_surf_l = (const double*) gkyl_array_cfetch(pkpm_u_surf, linl);
      const double *pkpm_u_surf_r = (const double*) gkyl_array_cfetch(pkpm_u_surf, linr);
      const double *prim_l = (const double*) gkyl_array_cfetch(prim, linl);
      const double *prim_r = (const double*) gkyl_array_cfetch(prim, linr);
      
      up->pkpm_accel[dir](up->conf_grid.dx, 
        pkpm_u_surf_l, pkpm_u_surf_c, pkpm_u_surf_r, 
        prim_l, prim_c, prim_r, 
        pkpm_u_d, bb_d, nu_d,
        pkpm_lax_d, pkpm_accel_d);
    }
  }
}

// Host-side wrapper for pkpm acceleration variable calculations 
void
gkyl_dg_calc_pkpm_vars_accel_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* pkpm_u_surf, const struct gkyl_array* pkpm_u, 
  const struct gkyl_array* prim, const struct gkyl_array* bb, 
  const struct gkyl_array* div_b, const struct gkyl_array* nu, 
  struct gkyl_array* pkpm_lax, struct gkyl_array* pkpm_accel)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_pkpm_vars_accel_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    pkpm_u_surf->on_dev, pkpm_u->on_dev, prim->on_dev, 
    bb->on_dev, div_b->on_dev, nu->on_dev, 
    pkpm_lax->on_dev, pkpm_accel->on_dev);
}

__global__ void
gkyl_dg_calc_pkpm_integrated_vars_cu_kernel(struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, 
  struct gkyl_array* int_pkpm_vars)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);

    const double *vlasov_pkpm_moms_d = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *pkpm_u_d = (const double*) gkyl_array_cfetch(pkpm_u, loc);

    double *int_pkpm_vars_d = (double*) gkyl_array_fetch(int_pkpm_vars, loc);
    up->pkpm_int(vlasov_pkpm_moms_d, pkpm_u_d, int_pkpm_vars_d);
  }
}

// Host-side wrapper for pkpm integrated variables calculation
void
gkyl_dg_calc_pkpm_integrated_vars_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, 
  struct gkyl_array* int_pkpm_vars)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_pkpm_integrated_vars_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    vlasov_pkpm_moms->on_dev, pkpm_u->on_dev, 
    int_pkpm_vars->on_dev);
}

__global__ void
gkyl_dg_calc_pkpm_vars_explicit_source_cu_kernel(struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* qmem, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm,
  struct gkyl_array* rhs)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);

    const double *qmem_d = (const double*) gkyl_array_cfetch(qmem, loc);
    const double *vlasov_pkpm_moms_d = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = (const double*) gkyl_array_cfetch(euler_pkpm, loc);

    double *rhs_d = (double*) gkyl_array_fetch(rhs, loc);
    up->pkpm_explicit_source(qmem_d, vlasov_pkpm_moms_d, euler_pkpm_d, rhs_d);
  }
}

// Host-side wrapper for pkpm explicit source term calculations
void
gkyl_dg_calc_pkpm_vars_explicit_source_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array* qmem, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* rhs)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_pkpm_vars_explicit_source_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    qmem->on_dev, vlasov_pkpm_moms->on_dev, euler_pkpm->on_dev, 
    rhs->on_dev);
}

__global__ void
gkyl_dg_calc_pkpm_vars_io_cu_kernel(struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, const struct gkyl_array* p_ij, 
  const struct gkyl_array* prim, const struct gkyl_array* pkpm_accel, 
  struct gkyl_array* fluid_io, struct gkyl_array* pkpm_vars_io)
{
  int idx[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long loc = gkyl_range_idx(&conf_range, idx);

    const double *vlasov_pkpm_moms_d = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *pkpm_u_d = (const double*) gkyl_array_cfetch(pkpm_u, loc);
    const double *p_ij_d = (const double*) gkyl_array_cfetch(p_ij, loc);
    const double *prim_d = (const double*) gkyl_array_cfetch(prim, loc);
    const double *pkpm_accel_d = (const double*) gkyl_array_cfetch(pkpm_accel, loc);

    double *fluid_io_d = (double*) gkyl_array_fetch(fluid_io, loc);
    double *pkpm_vars_io_d = (double*) gkyl_array_fetch(pkpm_vars_io, loc);
    up->pkpm_io(vlasov_pkpm_moms_d, pkpm_u_d, p_ij_d, prim_d, pkpm_accel_d, 
      fluid_io_d, pkpm_vars_io_d);
  }
}

// Host-side wrapper for pkpm io. Computes conserved variables and copies primitive and acceleration variables to output array
void
gkyl_dg_calc_pkpm_vars_io_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, const struct gkyl_array* p_ij, 
  const struct gkyl_array* prim, const struct gkyl_array* pkpm_accel, 
  struct gkyl_array* fluid_io, struct gkyl_array* pkpm_vars_io)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_pkpm_vars_io_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    vlasov_pkpm_moms->on_dev, pkpm_u->on_dev, p_ij->on_dev, prim->on_dev, pkpm_accel->on_dev, 
    fluid_io->on_dev, pkpm_vars_io->on_dev);
}

__global__ void
gkyl_dg_calc_pkpm_vars_limiter_cu_kernel(struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, const struct gkyl_array* p_ij, 
  struct gkyl_array* euler_pkpm)
{
  int cdim = up->cdim;
  int idxl[GKYL_MAX_DIM], idxc[GKYL_MAX_DIM], idxr[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < conf_range.volume;
      linc1 += gridDim.x*blockDim.x)
  {
    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&conf_range, linc1, idxc);
    const struct gkyl_wave_cell_geom *geom = gkyl_wave_geom_get(up->geom, idxc);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&conf_range, idxc);

    const double *vlasov_pkpm_moms_c = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, linc);
    const double *pkpm_u_c = (const double*) gkyl_array_cfetch(pkpm_u, linc);
    const double *p_ij_c = (const double*) gkyl_array_cfetch(p_ij, linc);

    double *euler_pkpm_c = (double*) gkyl_array_fetch(euler_pkpm, linc);
    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, idxc, idxl);
      gkyl_copy_int_arr(cdim, idxc, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(&conf_range, idxl); 
      long linr = gkyl_range_idx(&conf_range, idxr);

      const double *vlasov_pkpm_moms_l = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, linl);
      const double *vlasov_pkpm_moms_r = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, linr);
      const double *pkpm_u_l = (const double*) gkyl_array_cfetch(pkpm_u, linl);
      const double *pkpm_u_r = (const double*) gkyl_array_cfetch(pkpm_u, linr);
      const double *p_ij_l = (const double*) gkyl_array_cfetch(p_ij, linl);
      const double *p_ij_r = (const double*) gkyl_array_cfetch(p_ij, linr);

      double *euler_pkpm_l = (double*) gkyl_array_fetch(euler_pkpm, linl);
      double *euler_pkpm_r = (double*) gkyl_array_fetch(euler_pkpm, linr);

      up->pkpm_limiter[dir](up->limiter_fac, up->wv_eqn, geom, 
        vlasov_pkpm_moms_l, vlasov_pkpm_moms_c, vlasov_pkpm_moms_r, 
        pkpm_u_l, pkpm_u_c, pkpm_u_r, 
        p_ij_l, p_ij_c, p_ij_r, 
        euler_pkpm_l, euler_pkpm_c, euler_pkpm_r); 
    }
  }
}

// Host-side wrapper for slope limiter of fluid variables
void
gkyl_dg_calc_pkpm_vars_limiter_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* pkpm_u, const struct gkyl_array* p_ij, 
  struct gkyl_array* euler_pkpm)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_pkpm_vars_limiter_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    vlasov_pkpm_moms->on_dev, pkpm_u->on_dev, p_ij->on_dev, 
    euler_pkpm->on_dev);
}

// CUDA kernel to set device pointers to pkpm vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_pkpm_vars_set_cu_dev_ptrs(struct gkyl_dg_calc_pkpm_vars *up, 
  int cdim, int poly_order, int poly_order_2p)
{
  // Set the function pointers for the order p quantities we solve for 
  up->pkpm_u_set = choose_pkpm_u_set_kern(cdim, poly_order);
  up->pkpm_u_copy = choose_pkpm_u_copy_kern(cdim, poly_order);  
  up->pkpm_u_surf = choose_pkpm_u_surf_kern(cdim, poly_order);
  up->pkpm_explicit_source = choose_pkpm_explicit_source_kern(cdim, poly_order);
  for (int d=0; d<cdim; ++d) {
    up->pkpm_limiter[d] = choose_pkpm_limiter_kern(d, cdim, poly_order);
  }
  // Set the function pointers for the order 2*p quantities we solve for
  up->pkpm_set = choose_pkpm_set_kern(cdim, poly_order_2p);
  up->pkpm_copy = choose_pkpm_copy_kern(cdim, poly_order_2p);
  up->pkpm_surf_set = choose_pkpm_surf_set_kern(cdim, poly_order_2p);
  up->pkpm_surf_copy = choose_pkpm_surf_copy_kern(cdim, poly_order_2p);
  up->pkpm_pressure = choose_pkpm_pressure_kern(cdim, poly_order_2p);
  up->pkpm_p_force = choose_pkpm_p_force_kern(cdim, poly_order_2p);
  up->pkpm_int = choose_pkpm_int_kern(cdim, poly_order_2p);
  up->pkpm_io = choose_pkpm_io_kern(cdim, poly_order_2p);
  // Fetch the kernels in each direction
  for (int d=0; d<cdim; ++d) {
    up->pkpm_accel[d] = choose_pkpm_accel_kern(d, cdim, poly_order_2p);
  }
}

gkyl_dg_calc_pkpm_vars*
gkyl_dg_calc_pkpm_vars_cu_dev_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_basis* cbasis_2p, 
  const struct gkyl_range *mem_range, 
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *wg, double limiter_fac)
{
  struct gkyl_dg_calc_pkpm_vars *up = (struct gkyl_dg_calc_pkpm_vars*) gkyl_malloc(sizeof(gkyl_dg_calc_pkpm_vars));

  up->conf_grid = *conf_grid;
  int cdim = cbasis->ndim;
  up->cdim = cdim;
  up->mem_range = *mem_range;

  // acquire pointer to wave equation object
  struct gkyl_wv_eqn *eqn = gkyl_wv_eqn_acquire(wv_eqn);
  up->wv_eqn = eqn->on_dev; // this is so the memcpy below has eqn on_dev

  // acquire pointer to wave equation object
  struct gkyl_wave_geom *geom = gkyl_wave_geom_acquire(wg);
  up->geom = geom->on_dev; // this is so the memcpy below has geom on_dev

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
  // Linear system for solving for ux, uy, uz
  int poly_order = cbasis->poly_order;
  up->Ncomp_u = 3;
  int nc = cbasis->num_basis;
  up->As_u = gkyl_nmat_cu_dev_new(up->Ncomp_u*mem_range->volume, nc, nc);
  up->xs_u = gkyl_nmat_cu_dev_new(up->Ncomp_u*mem_range->volume, nc, 1);
  up->mem_u = gkyl_nmat_linsolve_lu_cu_dev_new(up->As_u->num, up->As_u->nr);

  // 6 components: div(p_par b)/rho, p_perp/rho, rho/p_perp, 3*Txx/m, 3*Tyy/m, 3*Tzz/m
  int poly_order_2p = cbasis_2p->poly_order;
  up->Ncomp_prim = 7;
  int nc_2p = cbasis_2p->num_basis;
  up->As = gkyl_nmat_cu_dev_new(up->Ncomp_prim*mem_range->volume, nc_2p, nc_2p);
  up->xs = gkyl_nmat_cu_dev_new(up->Ncomp_prim*mem_range->volume, nc_2p, 1);
  up->mem = gkyl_nmat_linsolve_lu_cu_dev_new(up->As->num, up->As->nr);

  // 2*cdim components: 
  // 3*Txx/m at the left and right x surfaces 
  // 3*Tyy/m at the left and right y surfaces 
  // 3*Tzz/m at the left and right z surfaces
  up->Ncomp_surf = 2*cdim;
  int nc_surf_2p = cbasis_2p->num_basis/(poly_order_2p+1); 
  up->As_surf = gkyl_nmat_cu_dev_new(up->Ncomp_surf*mem_range->volume, nc_surf_2p, nc_surf_2p);
  up->xs_surf = gkyl_nmat_cu_dev_new(up->Ncomp_surf*mem_range->volume, nc_surf_2p, 1);
  up->mem_surf = gkyl_nmat_linsolve_lu_cu_dev_new(up->As_surf->num, up->As_surf->nr);

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_pkpm_vars *up_cu = (struct gkyl_dg_calc_pkpm_vars*) gkyl_cu_malloc(sizeof(gkyl_dg_calc_pkpm_vars));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_pkpm_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_pkpm_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, cdim, poly_order, poly_order_2p);

  // set parent on_dev pointer
  up->on_dev = up_cu;

  up->wv_eqn = eqn; // updater should store host pointer 
  up->geom = geom; 
  
  return up;
}
