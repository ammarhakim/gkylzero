// CUB for reductions.
#include <cub/cub.cuh>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_util.h>
#include <gkyl_array_reduce.h>
}

int iDivUp(int a, int b) {
  // Round a / b to nearest higher integer value
  return (a % b != 0) ? (a / b + 1) : (a / b);
}

__device__ static float atomicMax_float(float* address, float val)
{
  int* address_as_i = (int*) address;
  int old = *address_as_i, assumed;
  do {
    assumed = old;
    old = ::atomicCAS(address_as_i, assumed,
      __float_as_int(::fmaxf(val, __int_as_float(assumed))));
  } while (assumed != old);
  return __int_as_float(old);
}

__device__ double atomicMax_double(double* address, double val)
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

template <unsigned int BLOCKSIZE>
__global__ void arrayMax_blockRedAtomic_cub(const struct gkyl_array* inp, double* out) {

  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

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
    if (linc < nCells) bResult = BlockReduceT(temp).Reduce(f, cub::Max());
    if (threadIdx.x == 0) {
      atomicMax_double(&out[k], bResult);
    }
  };
}

template <unsigned int BLOCKSIZE>
__global__ void arrayMax_range_blockRedAtomic_cub(const struct gkyl_array* inp, const struct gkyl_range range, double* out) {

  unsigned long linc = blockIdx.x * blockDim.x + threadIdx.x;

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
    if (linc < nCells) bResult = BlockReduceT(temp).Reduce(f, cub::Max());
    if (threadIdx.x == 0) atomicMax_double(&out[k], bResult);
  };
}

void
gkyl_array_reduce_init_cu(struct gkyl_array* arrIn, struct gkyl_array_reduce_util* redutil)
{
  // Current implementation (CUB block reduce + atomic reduce) does not need to preallocate aything.
}

void
gkyl_array_reduce_free_cu(struct gkyl_array_reduce_util* redutil)
{
  // Current implementation (CUB block reduce + atomic reduce) does not allocate additional memory.
}

void
gkyl_array_reduce_max_cu(double *out_d, const struct gkyl_array* inp)
{
  // Reduce a gkyl_array component-wise. The output must be a device array
  // with as many elements as there are components in the input array.

  int numCells = inp->size;
  const int blockSize = GKYL_DEFAULT_NUM_THREADS;

  arrayMax_blockRedAtomic_cub<blockSize><<<iDivUp(numCells, blockSize), blockSize>>>(inp->on_device, out_d);
  cudaDeviceSynchronize();
}

void
gkyl_array_reduce_range_max_cu(double *out_d, const struct gkyl_array* inp, const struct gkyl_range range)
{
  // Reduce a gkyl_array component-wise within a specified range. The output must be
  // a device array with as many elements as there are components in the input array.

  int numCells = range.volume;
  const int blockSize = GKYL_DEFAULT_NUM_THREADS;

  arrayMax_range_blockRedAtomic_cub<blockSize><<<iDivUp(numCells, blockSize), blockSize>>>(inp->on_device, range, out_d);
}
