/* -*- c++ -*- */
extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_mat.h>
#include <gkyl_prim_lbo_calc.h>
#include <gkyl_prim_lbo_calc_priv.h>
#include <gkyl_prim_lbo_vlasov_priv.h>
#include <gkyl_util.h>
}

static void
gkyl_get_array_range_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range range)
{
  int volume = range.volume;
  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  int ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;
  dimBlock->x = min(ac1, GKYL_DEFAULT_NUM_THREADS);
  dimGrid->x = gkyl_int_div_up(ac1, dimBlock->x);

  dimBlock->y = gkyl_int_div_up(GKYL_DEFAULT_NUM_THREADS, ac1);
  dimGrid->y = gkyl_int_div_up(volume, ac1*dimBlock->y);
}

__global__ static void
gkyl_prim_lbo_calc_set_cu_ker(gkyl_prim_lbo_calc* calc,
  struct gkyl_nmat *As, struct gkyl_nmat *xs,
  struct gkyl_basis cbasis, struct gkyl_range conf_rng,
  const struct gkyl_array* m0, const struct gkyl_array* m1, const struct gkyl_array* m2,
  const struct gkyl_array* cM, const struct gkyl_array* cE)
{
  int ndim = cbasis.ndim;
  int idx[GKYL_MAX_DIM];
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = conf_rng.iac[ndim-1] > 0 ? conf_rng.iac[ndim-1] : 1;

  // 2D thread grid
  // linc1 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
  // linc2 = idx2 + ac2*idx3 + ...
  for (unsigned long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
      linc2 < conf_rng.volume/ac1;
      linc2 += gridDim.y*blockDim.y)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc2.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    if (linc1 < ac1) {
      gkyl_sub_range_inv_idx(&conf_rng, linc1 + ac1*linc2, idx);
      long start = gkyl_range_idx(&conf_rng, idx);
    
      struct gkyl_mat lhs = gkyl_nmat_get(As, linc1 + ac1*linc2);
      struct gkyl_mat rhs = gkyl_nmat_get(xs, linc1 + ac1*linc2);
      const double *m0_d = (const double*) gkyl_array_cfetch(m0, start);
      const double *m1_d = (const double*) gkyl_array_cfetch(m1, start);
      const double *m2_d = (const double*) gkyl_array_cfetch(m2, start);
      const double *cM_d = (const double*) gkyl_array_cfetch(cM, start);
      const double *cE_d = (const double*) gkyl_array_cfetch(cE, start);

      gkyl_mat_clear(&lhs, 0.0); gkyl_mat_clear(&rhs, 0.0);

      calc->prim->self_prim(calc->prim, &lhs, &rhs, m0_d, m1_d, m2_d, cM_d, cE_d);
    }
  }
}

__global__ static void
gkyl_prim_lbo_copy_sol_cu_ker(struct gkyl_nmat *xs,
  struct gkyl_basis cbasis, struct gkyl_range conf_rng,
  int nc, int vdim,
  struct gkyl_array* uout, struct gkyl_array* vtSqout)
{
  int ndim = cbasis.ndim;
  int idx[GKYL_MAX_DIM];

  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = conf_rng.iac[ndim-1] > 0 ? conf_rng.iac[ndim-1] : 1;
  
  // 2D thread grid
  // linc1 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
  // linc2 = idx2 + ac2*idx3 + ...
  for (unsigned long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
      linc2 < conf_rng.volume/ac1;
      linc2 += gridDim.y*blockDim.y)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc2.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    if (linc1 < ac1) {
      gkyl_sub_range_inv_idx(&conf_rng, linc1 + ac1*linc2, idx);
      long start = gkyl_range_idx(&conf_rng, idx);

      struct gkyl_mat out_d = gkyl_nmat_get(xs, linc1 + ac1*linc2);
      double *u_d = (double*) gkyl_array_fetch(uout, start);
      double *vtSq_d = (double*) gkyl_array_fetch(vtSqout, start);
    
      prim_lbo_copy_sol(&out_d, nc, vdim, u_d, vtSq_d);
    }
  }
}

void
gkyl_prim_lbo_calc_advance_cu(gkyl_prim_lbo_calc* calc, struct gkyl_basis cbasis,
  struct gkyl_range conf_rng, 
  const struct gkyl_array* m0, const struct gkyl_array* m1,
  const struct gkyl_array* m2, const struct gkyl_array* cM, const struct gkyl_array* cE,
  struct gkyl_array* uout, struct gkyl_array* vtSqout)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, conf_rng);
  int nc = cbasis.num_basis;
  int vdim = calc->prim->pdim - calc->prim->cdim;
  int N = nc*(vdim + 1);      
  struct gkyl_nmat *A_d = gkyl_nmat_cu_dev_new(conf_rng.volume, N, N);
  struct gkyl_nmat *x_d = gkyl_nmat_cu_dev_new(conf_rng.volume, N, 1);

  gkyl_prim_lbo_calc_set_cu_ker<<<dimGrid, dimBlock>>>(calc->on_dev,
    A_d->on_dev, x_d->on_dev, cbasis, conf_rng,
    m0->on_dev, m1->on_dev, m2->on_dev,
    cM->on_dev, cE->on_dev);
  
  bool status = gkyl_nmat_linsolve_lu(A_d, x_d);

  gkyl_prim_lbo_copy_sol_cu_ker<<<dimGrid, dimBlock>>>(x_d->on_dev,
    cbasis, conf_rng, nc, vdim,
    uout->on_dev, vtSqout->on_dev);

  gkyl_nmat_release(A_d);
  gkyl_nmat_release(x_d);
}

gkyl_prim_lbo_calc*
gkyl_prim_lbo_calc_cu_dev_new(const struct gkyl_rect_grid *grid,
  struct gkyl_prim_lbo_type *prim)
{
  gkyl_prim_lbo_calc *up = (gkyl_prim_lbo_calc*) gkyl_malloc(sizeof(gkyl_prim_lbo_calc));
  up->grid = *grid;
  up->prim = prim;

  struct gkyl_prim_lbo_type *pt = gkyl_prim_lbo_type_acquire(prim);
  up->prim = pt->on_dev; // so memcpy below gets dev copy

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  gkyl_prim_lbo_calc *up_cu = (gkyl_prim_lbo_calc*) gkyl_cu_malloc(sizeof(gkyl_prim_lbo_calc));
  gkyl_cu_memcpy(up_cu, up, sizeof(gkyl_prim_lbo_calc), GKYL_CU_MEMCPY_H2D);

  up->prim = pt; // host portion of struct should have host copy
  up->on_dev = up_cu; // host pointer
  
  return up;
}
