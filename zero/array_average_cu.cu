/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_array_average.h>
#include <gkyl_array_average_priv.h>
}

// __global__ static void
// gkyl_array_average_set_ker_cu(struct gkyl_array_average *up)
// {
//   int ndim =  up->basis.ndim, poly_order = up->basis.poly_order;

//   int op = -1; // -1 shifted to start with 0
//   for (unsigned d = 0; d < ndim; d++)
//     op += pow(2,d) * up->avg_dim[d];

//   up->kernel = gkyl_array_average_ker_list[ndim-1].list[op].kernels[poly_order-1];

// }

// struct gkyl_array_average*
// gkyl_array_average_cu_dev_new(struct gkyl_array_average *up)
// {
//   // // Copy struct to device.
//   // struct gkyl_array_average *up_cu = (struct gkyl_array_average*) gkyl_cu_malloc(sizeof(struct gkyl_array_average));
//   // gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_array_average), GKYL_CU_MEMCPY_H2D);

//   // // Set the kernel.
//   // gkyl_array_average_set_ker_cu<<<1,1>>>(up_cu);

//   // up->on_dev = up_cu;

//   return up;
// }

// // template <unsigned int BLOCKSIZE>
// // __global__ void
// // array_integrate_blockRedAtomic_cub(struct gkyl_array_average *up, 
// //   const struct gkyl_range *full_rng, const struct gkyl_range *local_avg,
// //   const struct gkyl_array *GKYL_RESTRICT fin, struct gkyl_array *GKYL_RESTRICT *avgout)
// // {
// //   int full_idx[GKYL_MAX_CDIM], sub_idx[GKYL_MAX_CDIM];

//   /*
//   The cuda code needs to loop over the full range and then
//   converts the full dim index to the sub dim index
//   */
//   // for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
//   //     tid < phase_range.volume; tid += blockDim.x*gridDim.x) {
//   //   gkyl_sub_range_inv_idx(&phase_range, tid, full_idx);

//   //   long full_lidx = gkyl_range_idx(&phase_range, full_idx);
//   //   const double* fptr = (const double*) gkyl_array_cfetch(fin, full_lidx);

//   //   // get sub-space linear index using the full to sub dimension mapping
//   //   for (unsigned int k = 0; k < local_avg->ndim; k++)
//   //     sub_idx[k] = full_idx[up->sub_dir[k]];
//   //   long sub_lidx = gkyl_range_idx(local_avg, sub_idx);

//   //   double* avgptr = (double*) gkyl_array_fetch(avgout, sub_lidx);

//   //   // reduce local f to local avg
//   //   up->kernel(fptr, avgLocal, up->subvol);

//   //   for (unsigned int k = 0; k < mout->ncomp; ++k) {
//   //      if (tid < full_rng.volume)
//   //        atomicAdd(&avgptr[k], avgLocal[k]);
//   //   }
//   // }
//   // Integrate in this cell
  

//   //---- If we average over z, we need to perform a reduction

//   // // Specialize BlockReduce for type double.
//   // typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

//   // // Allocate temporary storage in shared memory.
//   // __shared__ typename BlockReduceT::TempStorage temp;

//   // for (size_t k = 0; k < up->num_comp; ++k) {
//   //   double f = 0;
//   //   if (linc < range.volume) f = outLocal[k];
//   //   double bResult = 0;
//   //   bResult = BlockReduceT(temp).Reduce(f, cub::Sum());
//   //   if (threadIdx.x == 0)
//   //     atomicAdd(&out[k], bResult);
//   // }
// // }

// __global__ static void
// gkyl_array_average_advance_cu_ker(const struct gkyl_array_average *up, 
//   const struct gkyl_array *fin, struct gkyl_array *avgout)
// {
//   // int pidx[GKYL_MAX_DIM];

//   // const int num_cbasis = 20; // MF 2024/09/03: Hardcoded to p=2 3x ser for now.

//   // for(unsigned long tid = threadIdx.x + blockIdx.x*blockDim.x;
//   //     tid < up->tot_rng.volume; tid += blockDim.x*gridDim.x) {
//   //   gkyl_sub_range_inv_idx(&up->tot_rng, tid, pidx);

//   //   long clinidx = gkyl_range_idx(&up->sub_rng, pidx);
//   //   const int *shiftedf_c = (const int*) gkyl_array_cfetch(shiftedf, clinidx);

//   //   if (shiftedf_c[0]) {
//   //     const double *delta_m0_c = (const double*) gkyl_array_cfetch(delta_m0, clinidx);
//   //     if (kers->is_m0_positive(delta_m0_c)) {
//   //       // Rescale f so it has the same m0 at this conf-space cell.
//   //       const double *m0_c = (const double*) gkyl_array_cfetch(m0, clinidx);
//   //       double m0ratio_c[num_cbasis];
//   //       kers->conf_inv_op(m0_c, m0ratio_c);
//   //       kers->conf_mul_op(delta_m0_c, m0ratio_c, m0ratio_c);

//   //       long plinidx = gkyl_range_idx(&up->tot_rng, pidx);
//   //       double *fin_i = (double*) gkyl_array_fetch(fin, plinidx);
//   //       kers->conf_phase_mul_op(m0ratio_c, fin_i, distf_c);
//   //     }
//   //   }
//   // }
// }

void gkyl_array_average_advance_cu(const struct gkyl_array_average *up, 
  const struct gkyl_array *fin, struct gkyl_array *avgout)
{

//   // int nblocks_tot = up->local.nblocks, nthreads_tot = up->local.nthreads;
//   // int nblocks_sub = up->local_avg.nblocks, nthreads_sub = up->local_avg.nthreads;

//   // gkyl_array_clear_range(avgout, 0.0, &up->local_avg);

//   // gkyl_array_average_advance_cu_ker<<nblocks_tot, nthreads_tot>>(up,fin,avgout);

//   // gkyl_cu_memset(out, 0, up->num_comp*sizeof(double));
//   // const int nthreads = GKYL_DEFAULT_NUM_THREADS;
//   // int nblocks = gkyl_int_div_up(range->volume, nthreads);
//   // array_integrate_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(up->on_dev, fin->on_dev, *range, out);
//   // // device synchronize required because out may be host pinned memory
//   // cudaDeviceSynchronize();

}
