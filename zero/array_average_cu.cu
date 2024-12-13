/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_average_priv.h>
#include <gkyl_array_average.h>
}

__global__ void
gkyl_array_average_set_ker_cu(struct gkyl_array_average *up)
{
  int ndim =  up->basis.ndim, poly_order = up->basis.poly_order;

  int op = -1; // -1 shifted to start with 0
  for (unsigned d = 0; d < ndim; d++)
    op += pow(2,d) * up->avg_dim[d];

  up->kernel = gkyl_array_average_ker_list[ndim-1].list[op].kernels[poly_order-1];

}

struct gkyl_array_average*
gkyl_array_average_cu_dev_new(struct gkyl_array_average *up)
{
  // Copy struct to device.
  struct gkyl_array_average *up_cu = (struct gkyl_array_average*) gkyl_cu_malloc(sizeof(struct gkyl_array_average));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_array_average), GKYL_CU_MEMCPY_H2D);

  // Set the kernel.
  gkyl_array_average_set_ker_cu<<<1,1>>>(up_cu);

  up->on_dev = up_cu;

  return up;
}

__global__ void
gkyl_array_average_advance_cu_ker(const struct gkyl_array_average *up, 
  const struct gkyl_array *fin, struct gkyl_array *avgout)
{
  int idx[GKYL_MAX_CDIM], idx_avg[GKYL_MAX_CDIM];

  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < up->local.volume; tid += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&up->local, tid, idx);

    // get the linear idx in the local range
    long lidx = gkyl_range_idx(&up->local, idx);

    // get avg-space linear index.
    for (unsigned int k = 0; k < up->local_avg.ndim; k++)
      idx_avg[k] = idx[k];
    long lidx_avg = gkyl_range_idx(&up->local_avg, idx_avg);

    // fetch the addresses where the weight and function are
    const double *fin_i = (const double*) gkyl_array_cfetch(fin, lidx);
    const double *win_i = up->isweighted? (const double*) gkyl_array_cfetch(up->weight, lidx) : 
        (const double*) gkyl_array_cfetch(up->weight, 0);

    // fetch the address where the avg is returned
    double *avg_i = (double*) gkyl_array_fetch(avgout, lidx_avg);
    
    // add (atomicAdd) the contribution of the local data to the avg
    up->kernel(up->subvol, win_i, fin_i, avg_i);
  }
}

void gkyl_array_average_advance_cu(const struct gkyl_array_average *up, 
  const struct gkyl_array *fin, struct gkyl_array *avgout)
{

  int nblocks = up->local.nblocks, nthreads = up->local.nthreads;
  // int nblocks_avg = up->local_avg.nblocks, nthreads_avg = up->local_avg.nthreads;

  gkyl_array_clear_range(avgout, 0.0, &up->local_avg);

  gkyl_array_average_advance_cu_ker<<<nblocks, nthreads>>>(up->on_dev,fin,avgout);

//   gkyl_cu_memset(out, 0, up->num_comp*sizeof(double));
//   const int nthreads = GKYL_DEFAULT_NUM_THREADS;
//   int nblocks = gkyl_int_div_up(range->volume, nthreads);
//   array_integrate_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(up->on_dev, fin->on_dev, *range, out);
  // device synchronize required because out may be host pinned memory
  cudaDeviceSynchronize();

  if (up->isweighted){
    gkyl_dg_div_op_range(up->div_mem, up->basis_avg,
      0, avgout, 0, avgout, 0, up->weight_avg, &up->local_avg);
  } else{
    // divide by the volume of the averaging domain
    gkyl_array_scale(avgout,up->vol_avg_inv);
  }

}
