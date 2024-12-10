/* -*- c++ -*- */

// CUB for reductions.
#include <cub/cub.cuh>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_array_average.h>
#include <gkyl_array_average_priv.h>
}

__global__ static void
gkyl_array_average_set_ker_cu(struct gkyl_array_average *up)
{
  int ndim =  up->tot_basis.ndim, poly_order = up->tot_basis.poly_order;

  int op = -1; // -1 shifted to start with 0
  for (unsigned d = 0; d < ndim; d++)
    op += pow(2,d) * up->avg_dim[d];

  up->kernel = gkyl_array_average_ker_list[ndim-1].list[op].kernels[poly_order-1];

}

struct gkyl_array_average*
gkyl_array_average_cu_dev_new(const gkyl_array_average *up)
{
  // Copy struct to device.
  struct gkyl_array_average *up_cu = (struct gkyl_array_average*) gkyl_cu_malloc(sizeof(struct gkyl_array_average));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_array_average), GKYL_CU_MEMCPY_H2D);

  // Set the kernel.
  gkyl_array_average_set_ker_cu<<<1,1>>>(up_cu);

  up->on_dev = up_cu;

  return up;
}

template <unsigned int BLOCKSIZE>
__global__ void
array_integrate_blockRedAtomic_cub(struct gkyl_array_average *up, 
  const struct gkyl_range *full_rng, const struct gkyl_range *sub_rng,
  const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT *avgout)
{
  int full_idx[GKYL_MAX_CDIM], sub_idx[GKYL_MAX_CDIM];

  /*
  The cuda code needs to loop over the full range and then
  converts the full dim index to the sub dim index
  */
  for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
      tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, tid, full_idx);

    long full_lidx = gkyl_range_idx(&phase_range, full_idx);
    const double* fptr = (const double*) gkyl_array_cfetch(fin, full_lidx);

    // get sub-space linear index using the full to sub dimension mapping
    for (unsigned int k = 0; k < sub_rng->ndim; k++)
      sub_idx[k] = full_idx[up->sub_dir[k]];
    long sub_lidx = gkyl_range_idx(sub_rng, sub_idx);

    double* avgptr = (double*) gkyl_array_fetch(avgout, sub_lidx);

    // reduce local f to local avg
    up->kernel(fptr, avgLocal, up->subvol);

    for (unsigned int k = 0; k < mout->ncomp; ++k) {
       if (tid < full_rng.volume)
         atomicAdd(&avgptr[k], avgLocal[k]);
    }
  }
  // Integrate in this cell
  

  //---- If we average over z, we need to perform a reduction

  // // Specialize BlockReduce for type double.
  // typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // // Allocate temporary storage in shared memory.
  // __shared__ typename BlockReduceT::TempStorage temp;

  // for (size_t k = 0; k < up->num_comp; ++k) {
  //   double f = 0;
  //   if (linc < range.volume) f = outLocal[k];
  //   double bResult = 0;
  //   bResult = BlockReduceT(temp).Reduce(f, cub::Sum());
  //   if (threadIdx.x == 0)
  //     atomicAdd(&out[k], bResult);
  // }
}

void gkyl_array_average_advance_cu(gkyl_array_average *up, 
  const struct gkyl_array *fin, struct gkyl_array *avgout)
{
  gkyl_cu_memset(out, 0, up->num_comp*sizeof(double));

  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  array_integrate_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(up->on_dev, fin->on_dev, *range, out);
  // device synchronize required because out may be host pinned memory
  cudaDeviceSynchronize();
}
