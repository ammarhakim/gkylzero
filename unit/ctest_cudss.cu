#include <cuDSS.h>

extern "C" {
#include <gkyl_alloc.h>
#include <assert.h>
}

#define TEST_NO_MAIN
#include <acutest.h>

extern "C" {
void test_cudss_simple();
}

#define checkCUDSS(call, status, msg) \
do { \
    status = call; \
    if (status != CUDSS_STATUS_SUCCESS) { \
        printf("Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: " #msg "\n", status); \
        exit(EXIT_FAILURE); \
    } \
} while(0);

void test_cudss_simple()
{
  // This is meant to replicate the "simple" example in the cuDSS folder of the CUDA samples repo:
  //   https://github.com/NVIDIA/CUDALibrarySamples/blob/master/cuDSS/simple/simple.cpp

  // ---------------------------------------------------------
  // cuDSS example: solving a real linear 5x5 system
  // with a symmetric positive-definite matrix
  // ---------------------------------------------------------
  cudssStatus_t status = CUDSS_STATUS_SUCCESS;

  int n = 5;
  int nnz = 8;
  int nrhs = 1;

  int *csr_offsets_h = NULL;
  int *csr_columns_h = NULL;
  double *csr_values_h = NULL;
  double *x_values_h = NULL, *b_values_h = NULL;

  int *csr_offsets_d = NULL;
  int *csr_columns_d = NULL;
  double *csr_values_d = NULL;
  double *x_values_d = NULL, *b_values_d = NULL;

  csr_offsets_h = (int*)malloc((n + 1) * sizeof(int));
  csr_columns_h = (int*)malloc(nnz * sizeof(int));
  csr_values_h = (double*)malloc(nnz * sizeof(double));
  x_values_h = (double*)malloc(nrhs * n * sizeof(double));
  b_values_h = (double*)malloc(nrhs * n * sizeof(double));

  if (!csr_offsets_h || ! csr_columns_h || !csr_values_h ||
      !x_values_h || !b_values_h) {
      printf("Error: host memory allocation failed\n");
      assert(false);
  }

  /* Initialize host memory for A and b */
  int i = 0;
  csr_offsets_h[i++] = 0;
  csr_offsets_h[i++] = 2;
  csr_offsets_h[i++] = 4;
  csr_offsets_h[i++] = 6;
  csr_offsets_h[i++] = 7;
  csr_offsets_h[i++] = 8;

  i = 0;
  csr_columns_h[i++] = 0; csr_columns_h[i++] = 2;
  csr_columns_h[i++] = 1; csr_columns_h[i++] = 2;
  csr_columns_h[i++] = 2; csr_columns_h[i++] = 4;
  csr_columns_h[i++] = 3;
  csr_columns_h[i++] = 4;

  i = 0;
  csr_values_h[i++] = 4.0; csr_values_h[i++] = 1.0;
  csr_values_h[i++] = 3.0; csr_values_h[i++] = 2.0;
  csr_values_h[i++] = 5.0; csr_values_h[i++] = 1.0;
  csr_values_h[i++] = 1.0;
  csr_values_h[i++] = 2.0;

  /* Note: Right-hand side b is initialized with values which correspond
     to the exact solution vector {1, 2, 3, 4, 5} */
  i = 0;
  b_values_h[i++] = 7.0;
  b_values_h[i++] = 12.0;
  b_values_h[i++] = 25.0;
  b_values_h[i++] = 4.0;
  b_values_h[i++] = 13.0;

  /* Allocate device memory for A, x and b */
  checkCuda(cudaMalloc(&csr_offsets_d, (n + 1) * sizeof(int)));
  checkCuda(cudaMalloc(&csr_columns_d, nnz * sizeof(int)));
  checkCuda(cudaMalloc(&csr_values_d, nnz * sizeof(double)));
  checkCuda(cudaMalloc(&b_values_d, nrhs * n * sizeof(double)));
  checkCuda(cudaMalloc(&x_values_d, nrhs * n * sizeof(double)));

  /* Copy host memory to device for A and b */
  checkCuda(cudaMemcpy(csr_offsets_d, csr_offsets_h, (n + 1) * sizeof(int),
                      cudaMemcpyHostToDevice));
  checkCuda(cudaMemcpy(csr_columns_d, csr_columns_h, nnz * sizeof(int),
                      cudaMemcpyHostToDevice));
  checkCuda(cudaMemcpy(csr_values_d, csr_values_h, nnz * sizeof(double),
                      cudaMemcpyHostToDevice));
  checkCuda(cudaMemcpy(b_values_d, b_values_h, nrhs * n * sizeof(double),
                      cudaMemcpyHostToDevice));

  /* Create a CUDA stream */
  cudaStream_t stream = NULL;
  checkCuda(cudaStreamCreate(&stream));

  /* Creating the cuDSS library handle */
  cudssHandle_t handle;

  checkCUDSS(cudssCreate(&handle), status, "cudssCreate");

  /* (optional) Setting the custom stream for the library handle */
  checkCUDSS(cudssSetStream(handle, stream), status, "cudssSetStream");

  /* Creating cuDSS solver configuration and data objects */
  cudssConfig_t solverConfig;
  cudssData_t solverData;

  checkCUDSS(cudssConfigCreate(&solverConfig), status, "cudssConfigCreate");
  checkCUDSS(cudssDataCreate(handle, &solverData), status, "cudssDataCreate");

  /* Create matrix objects for the right-hand side b and solution x (as dense matrices). */
  cudssMatrix_t x, b;

  int64_t nrows = n, ncols = n;
  int ldb = ncols, ldx = nrows;
  checkCUDSS(cudssMatrixCreateDn(&b, ncols, nrhs, ldb, b_values_d, CUDA_R_64F,
                       CUDSS_LAYOUT_COL_MAJOR), status, "cudssMatrixCreateDn for b");
  checkCUDSS(cudssMatrixCreateDn(&x, nrows, nrhs, ldx, x_values_d, CUDA_R_64F,
                       CUDSS_LAYOUT_COL_MAJOR), status, "cudssMatrixCreateDn for x");

  /* Create a matrix object for the sparse input matrix. */
  cudssMatrix_t A;
  cudssMatrixType_t mtype     = CUDSS_MTYPE_SPD;
  cudssMatrixViewType_t mview = CUDSS_MVIEW_UPPER;
  cudssIndexBase_t base       = CUDSS_BASE_ZERO;
  checkCUDSS(cudssMatrixCreateCsr(&A, nrows, ncols, nnz, csr_offsets_d, NULL,
                       csr_columns_d, csr_values_d, CUDA_R_32I, CUDA_R_64F, mtype, mview,
                       base), status, "cudssMatrixCreateCsr");

  /* Symbolic factorization */
  checkCUDSS(cudssExecute(handle, CUDSS_PHASE_ANALYSIS, solverConfig, solverData,
                       A, x, b), status, "cudssExecute for analysis");

  /* Factorization */
  checkCUDSS(cudssExecute(handle, CUDSS_PHASE_FACTORIZATION, solverConfig,
                       solverData, A, x, b), status, "cudssExecute for factor");

  /* Solving */
  checkCUDSS(cudssExecute(handle, CUDSS_PHASE_SOLVE, solverConfig, solverData,
                       A, x, b), status, "cudssExecute for solve");

  /* Destroying opaque objects, matrix wrappers and the cuDSS library handle */
  checkCUDSS(cudssMatrixDestroy(A), status, "cudssMatrixDestroy for A");
  checkCUDSS(cudssMatrixDestroy(b), status, "cudssMatrixDestroy for b");
  checkCUDSS(cudssMatrixDestroy(x), status, "cudssMatrixDestroy for x");
  checkCUDSS(cudssDataDestroy(handle, solverData), status, "cudssDataDestroy");
  checkCUDSS(cudssConfigDestroy(solverConfig), status, "cudssConfigDestroy");
  checkCUDSS(cudssDestroy(handle), status, "cudssHandleDestroy");

  checkCuda(cudaStreamSynchronize(stream));

  /* Print the solution and compare against the exact solution */
  checkCuda(cudaMemcpy(x_values_h, x_values_d, nrhs * n * sizeof(double),
                      cudaMemcpyDeviceToHost));

  int passed = 1;
  for (int i = 0; i < n; i++) {
    printf("x[%d] = %1.4f expected %1.4f\n", i, x_values_h[i], double(i+1));
    if (fabs(x_values_h[i] - (i + 1)) > 2.e-15)
      passed = 0;
  }

  /* Release the data allocated on the user side */
  free(csr_offsets_h);
  free(csr_columns_h);
  free(csr_values_h);
  free(x_values_h);
  free(b_values_h);
  cudaFree(csr_offsets_d);
  cudaFree(csr_columns_d);
  cudaFree(csr_values_d);
  cudaFree(x_values_d);
  cudaFree(b_values_d);

  if (status == CUDSS_STATUS_SUCCESS && passed)
    printf("Example PASSED\n");
  else
    printf("Example FAILED\n");
}
