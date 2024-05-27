/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_calc_em_vars.h>
#include <gkyl_dg_calc_em_vars_priv.h>
#include <gkyl_util.h>
}

__global__ static void
gkyl_dg_calc_em_vars_set_cu_kernel(struct gkyl_dg_calc_em_vars* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* em, struct gkyl_array* temp_var)
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

    const double *em_d = (const double*) gkyl_array_cfetch(em, loc);

    up->em_calc_temp(em_d, (double*) gkyl_array_fetch(temp_var, loc));
    up->em_set(count, As, xs, (const double*) gkyl_array_cfetch(temp_var, loc));
  }
}

__global__ static void
gkyl_dg_calc_em_vars_copy_cu_kernel(struct gkyl_dg_calc_em_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* em, struct gkyl_array* cell_avg_bb, 
  struct gkyl_array* out, struct gkyl_array* out_surf)
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

    const double *em_d = (const double*) gkyl_array_cfetch(em, loc);
    int *cell_avg_bb_d = (int*) gkyl_array_fetch(cell_avg_bb, loc);
    double *out_d = (double*) gkyl_array_fetch(out, loc);
    double *out_surf_d = (double*) gkyl_array_fetch(out_surf, loc);

    up->em_copy(count, xs, em_d, cell_avg_bb_d, out_d, out_surf_d);
  }
}

void gkyl_dg_calc_em_vars_advance_cu(struct gkyl_dg_calc_em_vars *up, 
  const struct gkyl_array* em, struct gkyl_array* cell_avg_bb, 
  struct gkyl_array* out, struct gkyl_array* out_surf)
{
  gkyl_array_clear(up->temp_var, 0.0);
  struct gkyl_range conf_range = up->mem_range;
  
  gkyl_dg_calc_em_vars_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range,
    em->on_dev, up->temp_var->on_dev);

  bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
  assert(status);

  gkyl_dg_calc_em_vars_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range,
    em->on_dev, cell_avg_bb->on_dev, 
    out->on_dev, out_surf->on_dev);
}

__global__ static void
gkyl_dg_calc_em_vars_diag_set_cu_kernel(struct gkyl_dg_calc_em_vars* up,
  struct gkyl_nmat *As_diag, struct gkyl_nmat *xs_diag, struct gkyl_range conf_range,
  const struct gkyl_array* em, struct gkyl_array* temp_var)
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
    long count = linc1*up->Ncomp_diag;

    const double *em_d = (const double*) gkyl_array_cfetch(em, loc);

    up->em_calc_temp(em_d, (double*) gkyl_array_fetch(temp_var, loc));
    up->em_diag_set(count, As_diag, xs_diag, em_d, (const double*) gkyl_array_cfetch(temp_var, loc));
  }
}

__global__ static void
gkyl_dg_calc_em_vars_diag_copy_cu_kernel(struct gkyl_dg_calc_em_vars* up, 
  struct gkyl_nmat *xs_diag, struct gkyl_range conf_range,
  const struct gkyl_array* em, struct gkyl_array* cell_avg_bb, 
  struct gkyl_array* out)
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
    long count = linc1*up->Ncomp_diag;

    const double *em_d = (const double*) gkyl_array_cfetch(em, loc);
    int *cell_avg_bb_d = (int*) gkyl_array_fetch(cell_avg_bb, loc);
    double *out_d = (double*) gkyl_array_fetch(out, loc);

    up->em_diag_copy(count, xs_diag, em_d, cell_avg_bb_d, out_d);
  }
}

void gkyl_dg_calc_em_vars_diag_cu(struct gkyl_dg_calc_em_vars *up, 
  const struct gkyl_array* em, struct gkyl_array* cell_avg_bb, 
  struct gkyl_array* out)
{
  gkyl_array_clear(up->temp_var, 0.0);
  struct gkyl_range conf_range = up->mem_range;
  
  gkyl_dg_calc_em_vars_diag_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As_diag->on_dev, up->xs_diag->on_dev, conf_range,
    em->on_dev, up->temp_var->on_dev);

  bool status = gkyl_nmat_linsolve_lu_pa(up->mem_diag, up->As_diag, up->xs_diag);
  assert(status);

  gkyl_dg_calc_em_vars_diag_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs_diag->on_dev, conf_range,
    em->on_dev, cell_avg_bb->on_dev, 
    out->on_dev);
}

__global__ void
gkyl_dg_calc_em_vars_div_b_cu_kernel(struct gkyl_dg_calc_em_vars *up, struct gkyl_range conf_range, 
  const struct gkyl_array* bvar_surf, const struct gkyl_array* bvar, 
  struct gkyl_array* max_b, struct gkyl_array* div_b)
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

    const double *bvar_surf_c = (const double*) gkyl_array_cfetch(bvar_surf, linc);
    const double *bvar_d = (const double*) gkyl_array_cfetch(bvar, linc);

    double *max_b_d = (double*) gkyl_array_fetch(max_b, linc);
    double *div_b_d = (double*) gkyl_array_fetch(div_b, linc);

    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, idxc, idxl);
      gkyl_copy_int_arr(cdim, idxc, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(&conf_range, idxl); 
      long linr = gkyl_range_idx(&conf_range, idxr);

      const double *bvar_surf_l = (const double*) gkyl_array_cfetch(bvar_surf, linl);
      const double *bvar_surf_r = (const double*) gkyl_array_cfetch(bvar_surf, linr);
      
      up->em_div_b[dir](up->conf_grid.dx, 
        bvar_surf_l, bvar_surf_c, bvar_surf_r, 
        bvar_d, max_b_d, div_b_d);
    }
  }
}

// Host-side wrapper for div(b) and max(|b_i|) variable calculations 
void
gkyl_dg_calc_em_vars_div_b_cu(struct gkyl_dg_calc_em_vars *up, const struct gkyl_range *conf_range, 
  const struct gkyl_array* bvar_surf, const struct gkyl_array* bvar, 
  struct gkyl_array* max_b, struct gkyl_array* div_b)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_em_vars_div_b_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, 
    bvar_surf->on_dev, bvar->on_dev, 
    max_b->on_dev, div_b->on_dev);
}

__global__ void
gkyl_dg_calc_em_vars_limiter_cu_kernel(struct gkyl_dg_calc_em_vars *up, struct gkyl_range conf_range, 
  struct gkyl_array* em)
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

    double *em_c = (double*) gkyl_array_fetch(em, linc);

    for (int dir=0; dir<cdim; ++dir) {
      gkyl_copy_int_arr(cdim, idxc, idxl);
      gkyl_copy_int_arr(cdim, idxc, idxr);

      idxl[dir] = idxl[dir]-1; idxr[dir] = idxr[dir]+1;

      long linl = gkyl_range_idx(&conf_range, idxl); 
      long linr = gkyl_range_idx(&conf_range, idxr);

      double *em_l = (double*) gkyl_array_fetch(em, linl);
      double *em_r = (double*) gkyl_array_fetch(em, linr);
      
      up->em_limiter[dir](up->limiter_fac, up->wv_eqn, geom, em_l, em_c, em_r);
    }
  }
}

// Host-side wrapper for slope limiter of em variables
void
gkyl_dg_calc_em_vars_limiter_cu(struct gkyl_dg_calc_em_vars *up, const struct gkyl_range *conf_range, 
  struct gkyl_array* em)
{
  int nblocks = conf_range->nblocks;
  int nthreads = conf_range->nthreads;
  gkyl_dg_calc_em_vars_limiter_cu_kernel<<<nblocks, nthreads>>>(up->on_dev, *conf_range, em->on_dev);
}

// CUDA kernel to set device pointers to em vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_em_vars_set_cu_dev_ptrs(struct gkyl_dg_calc_em_vars *up, enum gkyl_basis_type b_type,
  int cdim, int poly_order, int poly_order_2p)
{
  up->em_calc_temp = choose_em_calc_BB_kern(cdim, poly_order_2p);
  up->em_set = choose_em_set_bvar_kern(cdim, poly_order_2p);
  up->em_copy = choose_em_copy_bvar_kern(cdim, poly_order_2p);  
  up->em_diag_set = choose_em_set_diag_kern(cdim, poly_order_2p);
  up->em_diag_copy = choose_em_copy_diag_kern(cdim, poly_order_2p);      
  // Fetch the kernels in each direction
  for (int d=0; d<cdim; ++d) {
    up->em_div_b[d] = choose_em_div_b_kern(d, cdim, poly_order_2p); 
    // Limiter is only applied to order p EM variables
    up->em_limiter[d] = choose_em_limiter_kern(d, cdim, poly_order);
  }
}

gkyl_dg_calc_em_vars*
gkyl_dg_calc_em_vars_cu_dev_new(const struct gkyl_rect_grid *conf_grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *cbasis_2p, 
  const struct gkyl_range *mem_range, 
  const struct gkyl_wv_eqn *wv_eqn, const struct gkyl_wave_geom *wg, 
  double limiter_fac)
{
  struct gkyl_dg_calc_em_vars *up = (struct gkyl_dg_calc_em_vars*) gkyl_malloc(sizeof(gkyl_dg_calc_em_vars));

  up->conf_grid = *conf_grid;
  int nc_2p = cbasis_2p->num_basis;
  int cdim = cbasis_2p->ndim;
  int poly_order_2p = cbasis_2p->poly_order;
  int poly_order = cbasis->poly_order;
  up->cdim = cdim;
  up->mem_range = *mem_range;

  // acquire pointer to wave equation object
  struct gkyl_wv_eqn *eqn = gkyl_wv_eqn_acquire(wv_eqn);
  up->wv_eqn = eqn->on_dev; // this is so the memcpy below has eqn on_dev

  // acquire pointer to wave equation object
  struct gkyl_wave_geom *geom = gkyl_wave_geom_acquire(wg);
  up->geom = geom->on_dev; // this is so the memcpy below has geom on_dev

  up->Ncomp = 6;
  up->Ncomp_diag = 9;

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

  // 6 component temporary variable for either storing B_i B_j (for computing bb and ExB/|B|^2) 
  up->temp_var = gkyl_array_cu_dev_new(GKYL_DOUBLE, 6*nc_2p, mem_range->volume);

  // There are 6 linear systems to be solved for bb = B_i B_j/|B|^2
  up->As = gkyl_nmat_cu_dev_new(up->Ncomp*mem_range->volume, nc_2p, nc_2p);
  up->xs = gkyl_nmat_cu_dev_new(up->Ncomp*mem_range->volume, nc_2p, 1);
  up->mem = gkyl_nmat_linsolve_lu_cu_dev_new(up->As->num, up->As->nr);

  // There are 9 linear systems to be solved for B_i B_j/|B|^2 and (E x B)_i/|B|^2
  up->As_diag = gkyl_nmat_cu_dev_new(up->Ncomp_diag*mem_range->volume, nc_2p, nc_2p);
  up->xs_diag = gkyl_nmat_cu_dev_new(up->Ncomp_diag*mem_range->volume, nc_2p, 1);
  up->mem_diag = gkyl_nmat_linsolve_lu_cu_dev_new(up->As_diag->num, up->As_diag->nr);

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_em_vars *up_cu = (struct gkyl_dg_calc_em_vars*) gkyl_cu_malloc(sizeof(gkyl_dg_calc_em_vars));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_em_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_em_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, cdim, poly_order, poly_order_2p);

  // set parent on_dev pointer
  up->on_dev = up_cu;

  up->wv_eqn = eqn; // updater should store host pointer 
  up->geom = geom;  
  
  return up;
}
