#include <cuda_runtime.h>
#include <cub/cub.cuh>

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
cub_blockReduce_max(const double* inp, int size, int offset, double* out)
{
  unsigned long idx = blockIdx.x*blockDim.x + threadIdx.x;

  // Specialize BlockReduce for type double.
  typedef cub::BlockReduce<double, BLOCKSIZE> BlockReduceT;

  // Allocate temporary storage in shared memory.
  __shared__ typename BlockReduceT::TempStorage temp;

  double bResult = 0.;
  double data = -DBL_MAX;
  if (idx < size) data = inp[idx+offset];
  bResult = BlockReduceT(temp).Reduce(data, cub::Max());
  if (threadIdx.x == 0) {
    atomicMax_double(out, bResult);
  }
}

int main()
{
  // array sizes
  int a1_size = 100;
  int b1_size = 100;
 
  // allocate device arrays
  double* a1;
  cudaMalloc(&a1, sizeof(double)*a1_size);
  
  double* b1;
  cudaMalloc(&b1, sizeof(double)*b1_size);

  // allocate host arrays
  double* a1_h;
  cudaMallocHost(&a1_h, sizeof(double)*a1_size);
  
  double* b1_h;
  cudaMallocHost(&b1_h, sizeof(double)*b1_size);

  // allocate host-pinned memory for result
  double *a1_max, *b1_max;
  cudaMallocHost(&a1_max, sizeof(double));
  cudaMallocHost(&b1_max, sizeof(double));

  // initialize arrays on host
  for (unsigned i=0; i<a1_size; ++i) {
    a1_h[i] = (double) i;
  }

  for (unsigned i=0; i<b1_size; ++i) {
    b1_h[i] = (double) i+1000;
  }

  // copy host arrays to device
  cudaMemcpy(a1, a1_h, sizeof(double)*a1_size, cudaMemcpyHostToDevice);
  cudaMemcpy(b1, b1_h, sizeof(double)*b1_size, cudaMemcpyHostToDevice);

  // do reduction on b1
  const int nthreads = 64;
  int nblocks = b1_size/nthreads+1;
  cub_blockReduce_max<nthreads><<<nblocks, nthreads>>>(b1, b1_size, 0, b1_max);
  cudaDeviceSynchronize();

  // check result
  printf("b1_max = %f, correct = %f\n", b1_max[0], (double) b1_size-1+1000);

  // do reduction on subset of a1
  const int nthreads2 = 64;
  int offset = 10; // e.g. for halo
  int size = (a1_size-2*offset);
  nblocks = size/nthreads2+1;
  cub_blockReduce_max<nthreads2><<<nblocks, nthreads2>>>(a1, size, offset, a1_max);
  cudaDeviceSynchronize();

  // check result
  printf("a1_max = %f, correct = %f\n", a1_max[0], (double) size+offset-1);

  // clean up
  cudaFree(a1);
  cudaFree(b1);
  cudaFreeHost(a1_h);
  cudaFreeHost(b1_h);
  cudaFreeHost(a1_max);
  cudaFreeHost(b1_max);
}
