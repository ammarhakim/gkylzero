#include <gkyl_cusolver_ops.h>

#define TEST_NO_MAIN
#include <acutest.h>

extern "C" {
void test_cusolver_qr();
void test_cusolver_rf();
void test_cusolver_ops();
}

void test_cusolver_qr()
{
/*  
 * This is the small 5x5 example used in the Sections 2 and 3 of the 
 * SuperLU Users' Guide to illustrate how to solve a linear problem
 * with cuSolver.
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

  int *csrRowPtr = (int*) gkyl_malloc(sizeof(int)*(m+1));
  gkyl_cu_memcpy(csrRowPtr, csrRowPtr_cu, sizeof(int)*(m+1), GKYL_CU_MEMCPY_D2H);

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

void test_cusolver_ops()
{
  double s, u, p, e, r, l;
  int    nrhs, m, n, nnz;

  /* Initialize matrix A. */
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  m = n = 5;
  nrhs = 1;
  
  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  gkyl_mat_triples *tri = gkyl_mat_triples_new(m, n);
  gkyl_mat_triples_set_rowmaj_order(tri);
  // row 0
  gkyl_mat_triples_insert(tri, 0, 0, s);
  gkyl_mat_triples_insert(tri, 0, 2, u);
  gkyl_mat_triples_insert(tri, 0, 3, u);
  // row 1
  gkyl_mat_triples_insert(tri, 1, 0, l);
  gkyl_mat_triples_insert(tri, 1, 1, u);
  // row 2
  gkyl_mat_triples_insert(tri, 2, 1, l);
  gkyl_mat_triples_insert(tri, 2, 2, p);
  // row 3
  gkyl_mat_triples_insert(tri, 3, 3, e);
  gkyl_mat_triples_insert(tri, 3, 4, u);
  // row 4
  gkyl_mat_triples_insert(tri, 4, 0, l);
  gkyl_mat_triples_insert(tri, 4, 1, l);
  gkyl_mat_triples_insert(tri, 4, 4, r);

  // Create the cuSolver linear problem setup.
  gkyl_cusolver_prob *prob = gkyl_cusolver_prob_new(m, n, nrhs);

  // Allocate the A matrix from triples.
  gkyl_cusolver_amat_from_triples(prob, tri);
  gkyl_mat_triples_release(tri);

  // Create right-hand side matrix B = transpose([1,1,1,1,1]).
  gkyl_mat_triples *triRHS = gkyl_mat_triples_new(m, nrhs);
  gkyl_mat_triples_insert(triRHS, 0, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 1, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 2, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 3, 0, 1.0);
  gkyl_mat_triples_insert(triRHS, 4, 0, 1.0);
  gkyl_cusolver_brhs_from_triples(prob, triRHS);
  gkyl_mat_triples_release(triRHS);


  gkyl_cusolver_solve(prob);
  gkyl_cusolver_finish_host(prob);

  // Solution is: [-1/32, 11/168, 3/224, 1/16, 11/336].
  TEST_CHECK( gkyl_compare_double(-1.0/32.0,   gkyl_cusolver_get_sol_lin(prob,0), 1e-14) );
  TEST_CHECK( gkyl_compare_double( 11.0/168.0, gkyl_cusolver_get_sol_lin(prob,1), 1e-14) );
  TEST_CHECK( gkyl_compare_double( 3.0/224.0,  gkyl_cusolver_get_sol_lin(prob,2), 1e-14) );
  TEST_CHECK( gkyl_compare_double( 1.0/16.0,   gkyl_cusolver_get_sol_lin(prob,3), 1e-14) );
  TEST_CHECK( gkyl_compare_double( 11.0/336.0, gkyl_cusolver_get_sol_lin(prob,4), 1e-14) );

  gkyl_cusolver_prob_release(prob);
}


void test_cusolver_rf()
{
/*  
 * This is the small 5x5 example used in the Sections 2 and 3 of the 
 * SuperLU Users' Guide to illustrate how to solve a linear problem
 * with cuSolverRF, which allows us to update the coefficients of the
 * A matrix and reuse some of the information from the first LU
 * decomposition.
 * See cuda-samples/7_CUDALibraries/cuSolverRf/cuSolverRf.cpp
 *
 */
  double   s, u, p, e, r, l;

  /* Initialize matrix A. */
  /*  A : matrix([s,0,u,u,0],[l,u,0,0,0],[0,l,p,0,0],[0,0,0,e,u],[l,l,0,0,r]); */
  s = 19.0; u = 21.0; p = 16.0; e = 5.0; r = 18.0; l = 12.0;
  double csrValA[] = {s, u, u, l, u, l, p, e, u, l, l, r};
  int csrRowIndA[] = {0, 3, 5, 7, 9};
  int csrColIndA[] = {0, 2, 3, 0, 1, 1, 2, 3, 4, 0, 1, 4};

  int rowsA = sizeof(csrRowIndA)/sizeof(csrRowIndA[0]); // number of rows of A
  int colsA = rowsA; // number of columns of A
  int nnzA  = sizeof(csrValA)/sizeof(csrValA[0]); // number of nonzeros of A


  // Constants used in cusolverRf:
  //   nzero is the value below which zero pivot is flagged.
  //   nboost is the value which is substitured for zero pivot.
  double nzero = 0.0;
  double nboost= 0.0;
  // Constants used in cusolverSp:
  //   singularity is -1 if A is invertible under tol
  //   tol determines the condition of singularity
  //   pivot_threshold decides pivoting strategy
  int singularity = 0;
  const double tol = 1.e-14;
  const double pivot_threshold = 1.0;
  // Constants used in cusolverRf:
  const cusolverRfFactorization_t fact_alg = CUSOLVERRF_FACTORIZATION_ALG0; // default
  const cusolverRfTriangularSolve_t solve_alg = CUSOLVERRF_TRIANGULAR_SOLVE_ALG1; // default
