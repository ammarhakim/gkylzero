/* -*- c++ -*- */

// CUB for reductions.
#include <cub/cub.cuh>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_array_reduce.h>
}

__device__ static double
atomicMax_double(double* address, double val)
{
  unsigned long long int* address_as_ull = (unsigned long long int*) address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
      __double_as_longlong(fmax(val, __longlong_as_double(assumed))));
  } while (assumed != old);
  return __longlong_as_double(old);
}

template <unsigned int BLOCKSIZE> __global__ void
arrayMax_blockRedAtomic_cub(const struct gkyl_array* inp, double* out)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  unsigned int nComp = inp->ncomp;
  unsigned long nCells = inp->size;

  const double *inp_d = (const double*) inp->data;

  for (unsigned int k = 0; k < nComp; ++k) {
    double f = inp_d[linc*nComp+k];
    double bResult;
    if (linc < nCells)
      bResult = BlockReduceT(temp).Reduce(f, cub::Max());
    if (threadIdx.x == 0)
      atomicMax_double(&out[k], bResult);
  };
}

template <unsigned int BLOCKSIZE>
__global__ void
arrayMax_range_blockRedAtomic_cub(const struct gkyl_array* inp, const struct gkyl_range range, double* out)
{
  unsigned long linc = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  long nCells = range.volume, nComp = inp->ncomp;

  int idx[GKYL_MAX_DIM];

  for (int k = 0; k < nComp; ++k) {
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    double* fptr = (double*) gkyl_array_cfetch(inp, start);
    double f = fptr[k];
    double bResult;
    if (linc < nCells)
      bResult = BlockReduceT(temp).Reduce(f, cub::Max());
    if (threadIdx.x == 0)
      atomicMax_double(&out[k], bResult);
  };
}

void
gkyl_array_reduce_max_cu(double *out_d, const struct gkyl_array* inp)
{
  int numCells = inp->size;
  const int blockSize = GKYL_DEFAULT_NUM_THREADS;
  arrayMax_blockRedAtomic_cub<blockSize><<<gkyl_int_div_up(numCells, blockSize), blockSize>>>(inp->on_device, out_d);
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_range_max_cu(double *out_d, const struct gkyl_array* inp, struct gkyl_range range)
{
  int numCells = range.volume;
  const int blockSize = GKYL_DEFAULT_NUM_THREADS;
  arrayMax_range_blockRedAtomic_cub<blockSize><<<gkyl_int_div_up(numCells, blockSize), blockSize>>>(
    inp->on_device, range, out_d);
}
