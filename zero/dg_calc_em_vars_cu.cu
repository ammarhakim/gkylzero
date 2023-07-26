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
gkyl_dg_calc_em_vars_set_cu_kernel(gkyl_dg_calc_em_vars* up,
  struct gkyl_nmat *As, struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2, struct gkyl_array* temp_var)
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
    printf("count = %ld\n", count);

    const double *em_d = (const double*) gkyl_array_cfetch(em, loc);
    int *cell_avg_magB2_d = (int*) gkyl_array_fetch(cell_avg_magB2, loc);

    up->em_calc_temp(em_d, (double*) gkyl_array_fetch(temp_var, loc));
    cell_avg_magB2_d[0] = up->em_set(count, As, xs, (double*) gkyl_array_fetch(temp_var, loc));
  }
}

__global__ static void
gkyl_dg_calc_em_vars_copy_cu_kernel(gkyl_dg_calc_em_vars* up, 
  struct gkyl_nmat *xs, struct gkyl_range conf_range,
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2, struct gkyl_array* out)
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
    int *cell_avg_magB2_d = (int*) gkyl_array_fetch(cell_avg_magB2, loc);
    double *out_d = (double*) gkyl_array_fetch(out, loc);

    up->em_copy(count, xs, em_d, cell_avg_magB2_d, out_d);
  }
}

void gkyl_dg_calc_em_vars_advance_cu(struct gkyl_dg_calc_em_vars *up, 
  const struct gkyl_array* em, struct gkyl_array* cell_avg_magB2, struct gkyl_array* out)
{
  gkyl_array_clear(up->temp_var, 0.0);
  struct gkyl_range conf_range = up->conf_range;
  
  gkyl_dg_calc_em_vars_set_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->As->on_dev, up->xs->on_dev, conf_range,
    em->on_dev, cell_avg_magB2->on_dev, up->temp_var->on_dev);

  if (up->poly_order > 1) {
    bool status = gkyl_nmat_linsolve_lu_pa(up->mem, up->As, up->xs);
    assert(status);
  }

  gkyl_dg_calc_em_vars_copy_cu_kernel<<<conf_range.nblocks, conf_range.nthreads>>>(up->on_dev,
    up->xs->on_dev, conf_range,
    em->on_dev, cell_avg_magB2->on_dev, out->on_dev);
}

// CUDA kernel to set device pointers to em vars kernel functions
// Doing function pointer stuff in here avoids troublesome cudaMemcpyFromSymbol
__global__ static void 
dg_calc_em_vars_set_cu_dev_ptrs(struct gkyl_dg_calc_em_vars *up, enum gkyl_basis_type b_type,
  int cdim,int poly_order, bool is_ExB)
{
  if (is_ExB) {
    up->em_calc_temp = choose_em_calc_num_ExB_kern(b_type, cdim, poly_order);
    up->em_set = choose_em_set_ExB_kern(b_type, cdim, poly_order);
    up->em_copy = choose_em_copy_ExB_kern(b_type, cdim, poly_order);
  }
  else {
    up->em_calc_temp = choose_em_calc_BB_kern(b_type, cdim, poly_order);
    up->em_set = choose_em_set_bvar_kern(b_type, cdim, poly_order);
    up->em_copy = choose_em_copy_bvar_kern(b_type, cdim, poly_order);    
  }
}

gkyl_dg_calc_em_vars*
gkyl_dg_calc_em_vars_cu_dev_new(const struct gkyl_basis* cbasis, 
  const struct gkyl_range *conf_range, 
  bool is_ExB)
{
  struct gkyl_dg_calc_em_vars *up = (struct gkyl_dg_calc_em_vars*) gkyl_malloc(sizeof(gkyl_dg_calc_em_vars));

  int nc = cbasis->num_basis;
  enum gkyl_basis_type b_type = cbasis->b_type;
  int cdim = cbasis->ndim;
  int poly_order = cbasis->poly_order;
  up->poly_order = poly_order;
  up->conf_range = *conf_range;

  if (is_ExB) 
    up->Ncomp = 3;
  else 
    up->Ncomp = 6;

  // There are Ncomp more linear systems to be solved 
  // 6 components of bb and 3 components of E x B
  up->As = gkyl_nmat_cu_dev_new(up->Ncomp*conf_range->volume, nc, nc);
  up->xs = gkyl_nmat_cu_dev_new(up->Ncomp*conf_range->volume, nc, 1);
  up->mem = gkyl_nmat_linsolve_lu_cu_dev_new(up->As->num, up->As->nr);
  // 6 component temporary variable for either storing B_i B_j (for computing bb) 
  // or (E x B)_i and B_i^2 (for computing E x B/|B|^2)
  up->temp_var = gkyl_array_cu_dev_new(GKYL_DOUBLE, 6*nc, conf_range->volume);

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_calc_em_vars *up_cu = (struct gkyl_dg_calc_em_vars*) gkyl_cu_malloc(sizeof(gkyl_dg_calc_em_vars));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_dg_calc_em_vars), GKYL_CU_MEMCPY_H2D);

  dg_calc_em_vars_set_cu_dev_ptrs<<<1,1>>>(up_cu, b_type, cdim, poly_order, is_ExB);

  // set parent on_dev pointer
  up->on_dev = up_cu;
  
  return up;
}
