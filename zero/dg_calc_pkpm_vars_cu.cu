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
gkyl_dg_calc_pkpm_vars_set_cu_kernel(gkyl_dg_calc_pkpm_vars* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* p_ij, const struct gkyl_array* pkpm_div_ppar, 
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
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1*up->Ncomp;

    const double *vlasov_pkpm_moms_d = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, loc);
    const double *euler_pkpm_d = (const double*) gkyl_array_cfetch(euler_pkpm, loc);
    const double *p_ij_d = (const double*) gkyl_array_cfetch(p_ij, loc);
    const double *pkpm_div_ppar_d = (const double*) gkyl_array_cfetch(pkpm_div_ppar, loc);

    int* cell_avg_prim_d = (int*) gkyl_array_fetch(cell_avg_prim, loc);

    cell_avg_prim_d[0] = up->pkpm_set(count, As, xs, 
      vlasov_pkpm_moms_d, euler_pkpm_d, p_ij_d, pkpm_div_ppar_d);
  }
}

__global__ static void
gkyl_dg_calc_pkpm_vars_copy_cu_kernel(gkyl_dg_calc_pkpm_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range,
  struct gkyl_array* prim, struct gkyl_array* prim_surf)
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
    // fetch the correct count in the matrix (since we solve Ncomp systems in each cell)
    long count = linc1*up->Ncomp;

    double* prim_d = (double*) gkyl_array_fetch(prim, loc);
    double* prim_surf_d = (double*) gkyl_array_fetch(prim_surf, loc);

    up->pkpm_copy(count, xs, prim_d, prim_surf_d);
  }
}

// Host-side wrapper for pkpm primitive variable calculation
void gkyl_dg_calc_pkpm_vars_advance_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* p_ij, const struct gkyl_array* pkpm_div_ppar, 
  struct gkyl_array* cell_avg_prim, struct gkyl_array* prim, struct gkyl_array* prim_surf)
{
  struct gkyl_range conf_range = up->mem_range;
  
  gkyl_dg_calc_pkpm_vars_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range,
    vlasov_pkpm_moms->on_dev, euler_pkpm->on_dev, 
    p_ij->on_dev, pkpm_div_ppar->on_dev, 
    cell_avg_prim->on_dev);

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_dg_calc_pkpm_vars_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range, prim->on_dev, prim_surf->on_dev);
}

__global__ void
gkyl_calc_pkpm_vars_pressure_cu_kernel(struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, struct gkyl_array* p_ij)
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

    const double *bvar_d = (const double*) gkyl_array_cfetch(bvar, loc);
    const double *vlasov_pkpm_moms_d = (const double*) gkyl_array_cfetch(vlasov_pkpm_moms, loc);

    double *p_ij_d = (double*) gkyl_array_fetch(p_ij, loc);
    up->pkpm_pressure(bvar_d, vlasov_pkpm_moms_d, p_ij_d);
  }
}

// Host-side wrapper for pkpm pressure calculation
void gkyl_dg_calc_pkpm_vars_pressure_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bvar, const struct gkyl_array* vlasov_pkpm_moms, struct gkyl_array* p_ij)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_calc_pkpm_vars_pressure_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    bvar->on_dev, vlasov_pkpm_moms->on_dev, p_ij->on_dev);
}

__global__ void
gkyl_dg_calc_pkpm_vars_accel_cu_kernel(struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* prim_surf, const struct gkyl_array* prim, 
  const struct gkyl_array* bvar, const struct gkyl_array* div_b, const struct gkyl_array* nu, 
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

    const double *prim_surf_c = (const double*) gkyl_array_cfetch(prim_surf, linc);

    const double *prim_d = (const double*) gkyl_array_cfetch(prim, linc);
    const double *bvar_d = (const double*) gkyl_array_cfetch(bvar, linc);
    const double *div_b_d = (const double*) gkyl_array_cfetch(div_b, linc);
    const double *nu_d = (const double*) gkyl_array_cfetch(nu, linc);

    double *pkpm_lax_d = (double*) gkyl_array_fetch(pkpm_lax, linc);
    double *pkpm_accel_d = (double*) gkyl_array_fetch(pkpm_accel, linc);

    // Compute T_perp/m div(b) and p_force
    up->pkpm_p_force(prim_d, div_b_d, pkpm_accel_d);

    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, idxc, idxl);
      gkyl_copy_int_arr(cdim, idxc, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(&conf_range, idxl); 
      long linr = gkyl_range_idx(&conf_range, idxr);

      const double *prim_surf_l = (const double*) gkyl_array_cfetch(prim_surf, linl);
      const double *prim_surf_r = (const double*) gkyl_array_cfetch(prim_surf, linr);
      
      up->pkpm_accel[dir](up->conf_grid.dx, 
        prim_surf_l, prim_surf_c, prim_surf_r, 
        prim_d, bvar_d, nu_d,
        pkpm_lax_d, pkpm_accel_d);
    }
  }
}

// Host-side wrapper for pkpm acceleration variable calculations with recovery or averaging
void
gkyl_dg_calc_pkpm_vars_accel_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* prim_surf, const struct gkyl_array* prim, 
  const struct gkyl_array* bvar, const struct gkyl_array* div_b, const struct gkyl_array* nu, 
  struct gkyl_array* pkpm_lax, struct gkyl_array* pkpm_accel)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_pkpm_vars_accel_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    prim_surf->on_dev, prim->on_dev, 
    bvar->on_dev, div_b->on_dev, nu->on_dev, 
    pkpm_lax->on_dev, pkpm_accel->on_dev);
}

__global__ void
gkyl_dg_calc_pkpm_integrated_vars_cu_kernel(struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* prim, struct gkyl_array* int_pkpm_vars)
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
    const double *euler_pkpm_d = (const double*) gkyl_array_cfetch(euler_pkpm, loc);
    const double *prim_d = (const double*) gkyl_array_cfetch(prim, loc);

    double *int_pkpm_vars_d = (double*) gkyl_array_fetch(int_pkpm_vars, loc);
    up->pkpm_int(vlasov_pkpm_moms_d, euler_pkpm_d, prim_d, int_pkpm_vars_d);
  }
}

// Host-side wrapper for pkpm integrated variables calculation
void
gkyl_dg_calc_pkpm_integrated_vars_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  const struct gkyl_array* prim, struct gkyl_array* int_pkpm_vars)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_pkpm_integrated_vars_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    vlasov_pkpm_moms->on_dev, euler_pkpm->on_dev, prim->on_dev, 
    int_pkpm_vars->on_dev);
}

__global__ void
gkyl_dg_calc_pkpm_vars_source_cu_kernel(struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, 
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
    up->pkpm_source(qmem_d, vlasov_pkpm_moms_d, euler_pkpm_d, rhs_d);
  }
}

// Host-side wrapper for pkpm source term calculations
void
gkyl_dg_calc_pkpm_vars_source_cu(struct gkyl_dg_calc_pkpm_vars *up, const struct gkyl_range *conf_range,
  const struct gkyl_array* qmem, const struct gkyl_array* vlasov_pkpm_moms, const struct gkyl_array* euler_pkpm, 
  struct gkyl_array* rhs)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_pkpm_vars_source_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    qmem->on_dev, vlasov_pkpm_moms->on_dev, euler_pkpm->on_dev, 
    rhs->on_dev);
}

__global__ void
gkyl_dg_calc_pkpm_vars_io_cu_kernel(struct gkyl_dg_calc_pkpm_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* vlasov_pkpm_moms, 
  const struct gkyl_array* euler_pkpm, const struct gkyl_array* p_ij, 
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
    const double *euler_pkpm_d = (const double*) gkyl_array_cfetch(euler_pkpm, loc);
    const double *p_ij_d = (const double*) gkyl_array_cfetch(p_ij, loc);
    const double *prim_d = (const double*) gkyl_array_cfetch(prim, loc);
    const double *pkpm_accel_d = (const double*) gkyl_array_cfetch(pkpm_accel, loc);

    double *fluid_io_d = (double*) gkyl_array_fetch(fluid_io, loc);
    double *pkpm_vars_io_d = (double*) gkyl_array_fetch(pkpm_vars_io, loc);
    up->pkpm_io(vlasov_pkpm_moms_d, euler_pkpm_d, p_ij_d, prim_d, pkpm_accel_d, 
      fluid_io_d, pkpm_vars_io_d);
  }
}

// Host-side wrapper for pkpm io. Computes conserved variables and copies primitive and acceleration variables to output array
void
gkyl_dg_calc_pkpm_vars_io_cu(struct gkyl_dg_calc_pkpm_vars *up, 
  const struct gkyl_range *conf_range, const struct gkyl_array* vlasov_pkpm_moms, 
  const struct gkyl_array* euler_pkpm, const struct gkyl_array* p_ij, 
  const struct gkyl_array* prim, const struct gkyl_array* pkpm_accel, 
  struct gkyl_array* fluid_io, struct gkyl_array* pkpm_vars_io)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_pkpm_vars_io_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    vlasov_pkpm_moms->on_dev, euler_pkpm->on_dev, p_ij->on_dev, prim->on_dev, pkpm_accel->on_dev, 
    fluid_io->on_dev, pkpm_vars_io->on_dev);
}

// CUDA kernel to set device pointers to pkpm vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_pkpm_vars_set_cu_dev_ptrs(struct gkyl_dg_calc_pkpm_vars *up, enum gkyl_basis_type b_type,
  int cdim, int poly_order)
{
  up->pkpm_set = choose_pkpm_set_kern(b_type, cdim, poly_order);
  up->pkpm_copy = choose_pkpm_copy_kern(b_type, cdim, poly_order);
  up->pkpm_pressure = choose_pkpm_pressure_kern(b_type, cdim, poly_order);
  up->pkpm_p_force = choose_pkpm_p_force_kern(b_type, cdim, poly_order);
  up->pkpm_source = choose_pkpm_source_kern(b_type, cdim, poly_order);
  up->pkpm_int = choose_pkpm_int_kern(b_type, cdim, poly_order);
  up->pkpm_io = choose_pkpm_io_kern(b_type, cdim, poly_order);
  // Fetch the kernels in each direction
  for (int d=0; d<cdim; ++d) 
    up->pkpm_accel[d] = choose_pkpm_accel_kern(d, b_type, cdim, poly_order);
}

gkyl_dg_calc_pkpm_vars*
gkyl_dg_calc_pkpm_vars_cu_dev_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis* cbasis, const struct gkyl_range *mem_range)
{
  struct gkyl_dg_calc_pkpm_vars *up = (struct gkyl_dg_calc_pkpm_vars*) gkyl_malloc(sizeof(gkyl_dg_calc_pkpm_vars));

  up->conf_grid = *conf_grid;
  int nc = cbasis->num_basis;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  enum gkyl_basis_type b_type = cbasis->b_type;
  up->cdim = cdim;
  up->poly_order = poly_order;
  up->Ncomp = 9;
  up->mem_range = *mem_range;

  // There are Ncomp*range->volume linear systems to be solved 
  // 6 components: ux, uy, uz, div(p_par b)/rho, p_perp/rho, rho/p_perp
  up->As = gkyl_nmat_cu_dev_new(up->Ncomp*mem_range->volume, nc, nc);
  up->xs = gkyl_nmat_cu_dev_new(up->Ncomp*mem_range->volume, nc, 1);
  up->mem = gkyl_nmat_linsolve_lu_cu_dev_new(up->As->num, up->As->nr);

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_pkpm_vars *up_cu = (struct gkyl_dg_calc_pkpm_vars*) gkyl_cu_malloc(sizeof(gkyl_dg_calc_pkpm_vars));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_pkpm_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_pkpm_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, b_type, cdim, poly_order);

  // set parent on_dev pointer
  up->on_dev = up_cu;
  
  return up;
}
