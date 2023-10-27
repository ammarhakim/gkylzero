/* -*- c++ -*- */

// CUB for reductions.
#include <cub/cub.cuh>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_array_integrate.h>
#include <gkyl_array_integrate_priv.h>
}

__global__ static void
gkyl_array_integrate_set_ker_cu(struct gkyl_array_integrate *up, enum gkyl_array_integrate_op op)
{
  up->kernel = gkyl_array_integrate_ker_list.kernels[op];
}

struct gkyl_array_integrate*
gkyl_array_integrate_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  int num_comp, enum gkyl_array_integrate_op op)
{
  // Allocate space for new updater.
  struct gkyl_array_integrate *up = (struct gkyl_array_integrate*) gkyl_malloc(sizeof(struct gkyl_array_integrate));

  up->num_basis = basis->num_basis;
  up->num_comp = num_comp;
  up->use_gpu = true;

  int ndim = basis->ndim;
  up->vol = 1.0;
  if (basis->poly_order > 0) {
    for (unsigned d=0; d<ndim; ++d)
      up->vol *= op == GKYL_ARRAY_INTEGRATE_OP_SQ? grid->dx[d]/2.0 : grid->dx[d]/sqrt(2.0);
  } else {
    // For polyOrder = 0 no normalization is applied.
    for (unsigned d=0; d<ndim; ++d)
      up->vol *= grid->dx[d];
  }

  // Copy struct to device.
  struct gkyl_array_integrate *up_cu = (struct gkyl_array_integrate*) gkyl_cu_malloc(sizeof(struct gkyl_array_integrate));
  gkyl_cu_memcpy(up_cu, up, sizeof(struct gkyl_array_integrate), GKYL_CU_MEMCPY_H2D);

  // Set the kernel.
  gkyl_array_integrate_set_ker_cu<<<1,1>>>(up_cu, op);

  up->on_dev = up_cu;

  return up;
}

template <unsigned int BLOCKSIZE>
__global__ void
array_integrate_blockRedAtomic_cub(const struct gkyl_array_integrate *up,
  const struct gkyl_array *inp, double weight, const struct gkyl_range range, double *out)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  int idx[GKYL_MAX_DIM];
  gkyl_sub_range_inv_idx(&range, linc, idx);
  long start = gkyl_range_idx(&range, idx);
  const double *fptr = (const double*) gkyl_array_cfetch(inp, start);

  double outLocal[10]; // Set to max of 10 (e.g. heat flux tensor).
  for (unsigned int k=0; k<up->num_comp; ++k)
    outLocal[k] = 0.0;

  // Integrate in this cell
  up->kernel(up->vol*weight, up->num_comp, up->num_basis, fptr, outLocal);

  for (size_t k = 0; k < up->num_comp; ++k) {
    double f = 0;
    if (linc < range.volume) f = outLocal[k];
    double bResult = 0;
    bResult = BlockReduceT(temp).Reduce(f, cub::Sum());
    if (threadIdx.x == 0)
      atomicAdd(&out[k], bResult);
  }
}

void gkyl_array_integrate_advance_cu(gkyl_array_integrate *up, const struct gkyl_array *arr,
  double weight, const struct gkyl_range *range, double *out)
{
  gkyl_cu_memset(out, 0, up->num_comp*sizeof(double));

  const int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks = gkyl_int_div_up(range->volume, nthreads);
  array_integrate_blockRedAtomic_cub<nthreads><<<nblocks, nthreads>>>(up->on_dev, arr->on_dev, weight, *range, out);
  // device synchronize required because out may be host pinned memory
  cudaDeviceSynchronize();
}
