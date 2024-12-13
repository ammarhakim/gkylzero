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

  // Should we assign manually these variables or the gkyl_cu_memcpy is enough?
  // up_cu->local = up->local;
  // up_cu->local_avg = up->local_avg;
  // for (int d = 0; d < up->ndim; d++)
  //   up_cu->isdim_sub[d] = up->isdim_sub[d];
  // up_cu->basis = up->basis;
  // up_cu->isweighted = up->isweighted;

  // // allocate space for the weight
  // gkyl_array_cu_dev_new(GKYL_DOUBLE, up->weight->ncomp, up->weight->size);
  // gkyl_array_copy(up_cu->weight, up->weight);

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

    // build the idx_avg array with the subdim coordinates of the idx vector
    int cnter = 0;
    for (int i = 0; i < up->basis.ndim; i++){
      if (up->isdim_sub[i]) {
        idx_avg[i] = idx[cnter];
        cnter ++;
      }
    }

    long lidx_avg = gkyl_range_idx(&up->local_avg, idx_avg);

    // fetch the addresses where the weight and function are
    const double *fin_i = (const double*) gkyl_array_cfetch(fin, lidx);
    const double *win_i = up->isweighted? (const double*) gkyl_array_cfetch(up->weight, lidx) : 
        (const double*) gkyl_array_cfetch(up->weight, 0);
    // fetch the address where the avg is returned
    double *avg_i = (double*) gkyl_array_fetch(avgout, lidx_avg);
    
    // This line is creating a "CUDA error: an illegal memory access was encountered"
    up->kernel(up->subvol, win_i, fin_i, avg_i);
    // Even if we replace it by this line, the same error will rise
    // printf("win_i[%d] = %g\n", 0, win_i[0]);
    // It looks like the array are not well allocated. The win_i[0] should exist 
  }
}

void gkyl_array_average_advance_cu(const struct gkyl_array_average *up, 
  const struct gkyl_array *fin, struct gkyl_array *avgout)
{

  int nblocks = up->local.nblocks, nthreads = up->local.nthreads;

  gkyl_array_clear_range(avgout, 0.0, &up->local_avg);

  gkyl_array_average_advance_cu_ker<<<nblocks, nthreads>>>(up->on_dev, fin, avgout);

  if (up->isweighted){
    gkyl_dg_div_op_range(up->div_mem, up->basis_avg,
      0, avgout, 0, avgout, 0, up->weight_avg, &up->local_avg);
  } else{
    // divide by the volume of the averaging domain
    gkyl_array_scale(avgout,up->vol_avg_inv);
  }

}
