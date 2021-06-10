#include <acutest.h>
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_reduce.h>
#include <gkyl_util.h>

// CUDA specific tests
#ifdef GKYL_HAVE_CUDA

void test_cu_array_reduce_max()
{

  unsigned long numComp = 3, numCells = 10;
  // create host array and device copy
  struct gkyl_array *a1 = gkyl_array_new(GKYL_DOUBLE, numComp, numCells);
  struct gkyl_array *a1_cu = gkyl_array_cu_dev_new(GKYL_DOUBLE, numComp, numCells);
  // initialize data
  double *a1_d = a1->data;
  for (unsigned i=0; i<numCells; ++i) {
    for (unsigned k=0; k<numComp; ++k) {
      a1_d[i*numComp+k] = (double)i+(double)k*0.10;
  }}
  gkyl_array_copy(a1_cu, a1);
  // create a device and host arrays to store reduction
  double* a1max = gkyl_malloc(numComp*sizeof(double));
  double* a1max_cu = (double*) gkyl_cu_malloc(numComp*sizeof(double));

  // component-wise reduce array and copy the result back to host 
  gkyl_array_reduce(a1max_cu, a1_cu, GKYL_MAX);

  // copy to host and check values.
  gkyl_cu_memcpy(a1max, a1max_cu, numComp*sizeof(double), GKYL_CU_MEMCPY_D2H);
  for (unsigned k=0; k<numComp; ++k) {
    double amax;
    TEST_CHECK( gkyl_compare(a1max[k], (double)(numCells-1)+(double)k*0.10, 1e-14) );
  }

  gkyl_free(a1max);
  gkyl_cu_free(a1max_cu);
  gkyl_array_release(a1);
  gkyl_array_release(a1_cu);

}

#endif

TEST_LIST = {
#ifdef GKYL_HAVE_CUDA
  { "cu_array_reduce_max", test_cu_array_reduce_max },
#endif
  { NULL, NULL },
};
