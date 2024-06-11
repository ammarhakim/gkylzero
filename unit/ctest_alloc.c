#include <acutest.h>
#include <gkyl_alloc.h>

void
test_aligned_alloc()
{
  int *d1 = gkyl_aligned_alloc(8, 100*sizeof(int));
  TEST_CHECK( (ptrdiff_t) d1 % 8 == 0 );
  gkyl_aligned_free(d1);

  int *d2 = gkyl_aligned_alloc(16, 100*sizeof(int));
  TEST_CHECK( (ptrdiff_t) d2 % 16 == 0 );
  gkyl_aligned_free(d2);

  int *d3 = gkyl_aligned_alloc(32, 100*sizeof(int));
  TEST_CHECK( (ptrdiff_t) d3 % 32 == 0 );
  gkyl_aligned_free(d3);

  int *d4 = gkyl_aligned_alloc(64, 100*sizeof(int));
  TEST_CHECK( (ptrdiff_t) d4 % 64 == 0 );
  gkyl_aligned_free(d4);
}

void
test_aligned_realloc()
{
  int n = 10;
  int *d = gkyl_aligned_alloc(16, n*sizeof(int));

  for (int i=0; i<n; ++i)
    d[i] = 2*i;

  int *rd = gkyl_aligned_realloc(d, 16, n*sizeof(int), 2*n*sizeof(int));

  TEST_CHECK( (ptrdiff_t) rd % 16 == 0 );

  for (int i=0; i<n; ++i)
    TEST_CHECK( rd[i] == 2*i );

  gkyl_aligned_free(rd);
}

#ifdef GKYL_HAVE_CUDA
void
test_cu_malloc()
{
  // Test a simple allocation on the GPU.
  int nelem = 6;

  double *arr = gkyl_cu_malloc(nelem*sizeof(double));
  double *arr_ho = gkyl_malloc(nelem*sizeof(double));

  for (int i=0; i<nelem; i++)
    arr_ho[i] = 11.*i;

  gkyl_cu_memcpy(arr, arr_ho, nelem*sizeof(double), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(arr_ho, arr, nelem*sizeof(double), GKYL_CU_MEMCPY_D2H);

  for (int i=0; i<nelem; i++)
    TEST_CHECK( arr_ho[i] == 11.*i);

  gkyl_free(arr_ho);
  gkyl_cu_free(arr);
}

int dev_cu_malloc_array(double **arr, int narr, int nelem);

void
test_cu_malloc_array()
{
  // Test allocation of arrays of arrays on the GPU.
  int narr = 2;
  int nelem = 6;

  double **arr     = gkyl_cu_malloc(narr*sizeof(double*));
  double **arr_mem = gkyl_malloc(narr*sizeof(double*));
  double **arr_ho  = gkyl_malloc(narr*sizeof(double*));

  for (int k=0; k<narr; k++) {
    arr_mem[k] = gkyl_cu_malloc(nelem*sizeof(double));
    arr_ho[k] = gkyl_malloc(nelem*sizeof(double));

    double *arr_d = arr_ho[k];
    for (int i=0; i<nelem; i++)
      arr_d[i] = (double)(nelem*k+i);

    gkyl_cu_memcpy(arr_mem[k], arr_ho[k], nelem*sizeof(double), GKYL_CU_MEMCPY_H2D);

  }
  gkyl_cu_memcpy(arr, arr_mem, narr*sizeof(double*), GKYL_CU_MEMCPY_H2D);

  // Check arrays in device kernel
  int nfail = dev_cu_malloc_array(arr, narr, nelem);
  TEST_CHECK(nfail == 0);

  for (int k=0; k<narr; k++) {
    gkyl_free(arr_ho[k]);
    gkyl_cu_free(arr_mem[k]);
  }
  gkyl_free(arr_ho);
  gkyl_free(arr_mem);
  gkyl_cu_free(arr);
}
#endif

TEST_LIST = {
  { "aligned_alloc", test_aligned_alloc },
  { "aligned_realloc", test_aligned_realloc },
#ifdef GKYL_HAVE_CUDA
  { "cu_malloc", test_cu_malloc },
  { "cu_malloc_array", test_cu_malloc_array },
#endif
  { NULL, NULL },
};
