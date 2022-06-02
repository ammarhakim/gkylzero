#include <gkyl_util.h>
#include <gkylzero.h>
#include <cusparse.h>
#include <cusolverSp.h>
#include <stdbool.h>

#define TEST_NO_MAIN
#include <acutest.h>

extern "C" {
void test_cusolver_qr();
}

void test_cusolver_qr()
{
/*  
 * This is the small 5x5 example used in the Sections 2 and 3 of the 
 * Users' Guide to illustrate how to call a SuperLU routine, and the
 * matrix data structures used by SuperLU.
 *
 */
  double   s, u, p, e, r, l;
  int      nrhs, m, n, nnz;
  cusolverSpHandle_t cusolverH = NULL;
  cusparseHandle_t cusparseH = NULL;
  csrqrInfo_t info = NULL;
  cusparseMatDescr_t A = NULL;
  cudaStream_t stream = NULL;
  size_t size_qr = 0;
  size_t size_internal = 0;
  void *buffer_qr = nullptr; // working space for numerical factorization
  

  /* Initialize matrix A. */
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  m = n = 5;
  nnz = 12;
  nrhs = 1;

  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
  int cooRowInd[] = {0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4};
  int colInd[] = {0, 2, 3, 0, 1, 1, 2, 3, 4, 0, 1, 4};
  double Aval[] = {s, u, u, l, u, l, p, e, u, l, l, r};

  cusolverSpCreate(&cusolverH);
  cusparseCreate(&cusparseH);

  cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
  cusolverSpSetStream(cusolverH, stream);

  // allocate vectors on device, and copy to device
  int *cooRowInd_cu, *colInd_cu;
  double *Aval_cu;
  cooRowInd_cu = (int*) gkyl_cu_malloc(sizeof(int)*nnz);
  colInd_cu = (int*) gkyl_cu_malloc(sizeof(int)*nnz);
  Aval_cu = (double*) gkyl_cu_malloc(sizeof(double)*nnz);

  gkyl_cu_memcpy(cooRowInd_cu, cooRowInd, sizeof(int)*nnz, GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(colInd_cu, colInd, sizeof(int)*nnz, GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(Aval_cu, Aval, sizeof(double)*nnz, GKYL_CU_MEMCPY_H2D);

  // convert to csr format
  int *csrRowPtr_cu;
  csrRowPtr_cu = (int*) gkyl_cu_malloc(sizeof(int)*(m+1));
  cusparseXcoo2csr(cusparseH, cooRowInd_cu, nnz, m, csrRowPtr_cu, CUSPARSE_INDEX_BASE_ZERO);

  cusparseCreateMatDescr(&A);

  cusparseSetMatType(A, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(A, CUSPARSE_INDEX_BASE_ZERO); 

  cusolverSpCreateCsrqrInfo(&info);
  
  // allocate rhs vector
  double b[] = {1.0, 1.0, 1.0, 1.0, 1.0};
  double *b_cu;
  b_cu = (double*) gkyl_cu_malloc(sizeof(double)*m);
  gkyl_cu_memcpy(b_cu, b, sizeof(double)*m, GKYL_CU_MEMCPY_H2D);

  // allocate solution vector
  double *x = (double*) gkyl_malloc(sizeof(double)*n);
  double *x_cu = (double*) gkyl_cu_malloc(sizeof(double)*n);

  // symbolic analysis
  cusolverSpXcsrqrAnalysisBatched(cusolverH, m, n, nnz, A, csrRowPtr_cu, colInd_cu, info);

  // prepare working space
  cusolverSpDcsrqrBufferInfoBatched(cusolverH, m, n, nnz, A, Aval_cu,
                                                     csrRowPtr_cu, colInd_cu, nrhs, info,
                                                     &size_internal, &size_qr);
  
  // numerical factorization
  buffer_qr = (void*) gkyl_cu_malloc(size_qr);
  cusolverSpDcsrqrsvBatched(cusolverH, m, n, nnz, A, Aval_cu, csrRowPtr_cu,
                                             colInd_cu, b_cu, x_cu, nrhs, info, buffer_qr);

  cudaStreamSynchronize(stream);

  gkyl_cu_memcpy(x, x_cu, sizeof(double)*n, GKYL_CU_MEMCPY_D2H);

  /* Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336] */
  TEST_CHECK( gkyl_compare_double(-1.0/32.0, x[0], 1e-14) );
  TEST_CHECK( gkyl_compare_double( 11.0/168.0, x[1], 1e-14) );
  TEST_CHECK( gkyl_compare_double( 3.0/224.0, x[2], 1e-14) );
  TEST_CHECK( gkyl_compare_double( 1.0/16.0, x[3], 1e-14) );
  TEST_CHECK( gkyl_compare_double( 11.0/336.0, x[4], 1e-14) );

  gkyl_cu_free(x_cu);
  gkyl_cu_free(b_cu);
  gkyl_cu_free(buffer_qr);
  gkyl_cu_free(cooRowInd_cu);
  gkyl_cu_free(colInd_cu);
  gkyl_cu_free(Aval_cu);
  gkyl_cu_free(csrRowPtr_cu);

  cusolverSpDestroy(cusolverH);
  cudaStreamDestroy(stream);
}

