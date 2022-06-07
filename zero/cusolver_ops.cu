#include <cusparse.h>
#include <cusolverSp.h>
#include <cusolverRf.h>
#include <cusolverSp_LOWLEVEL_PREVIEW.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_mat_triples.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_cusolver_ops.h>
}


#ifdef GKYL_HAVE_CUDA
#include <cusparse.h>
#include <cusolverSp.h>
#include <cusolverRf.h>
#include <cusolverSp_LOWLEVEL_PREVIEW.h>
#endif

struct gkyl_cusolver_prob {
  double *rhs, *rhs_cu; // right-hand side entries. 
  double *x, *x_cu;
  double *nzval_cu;
  int *rowptr_cu, *colind_cu;
  int mrow, ncol; // A is a mrow x ncol matrix.
  int nnz; // number of non-zero entries in A.
  int nrhs; // number of problems to solve (B is an mrow x nrhs matrix).
  cusolverSpHandle_t cusolverH;
  csrqrInfo_t info;
  cusparseMatDescr_t A;
  cudaStream_t stream;
  size_t size_qr;
  size_t size_internal;
  void *buffer_qr; // working space for numerical factorization
};

gkyl_cusolver_prob*
gkyl_cusolver_prob_new(const int mrow, const int ncol, const int nprob)
{
  struct gkyl_cusolver_prob *prob = (struct gkyl_cusolver_prob*) gkyl_malloc(sizeof(*prob));

  prob->mrow = mrow;
  prob->ncol = ncol;
  prob->nrhs = nprob;

  prob->rhs = (double*) gkyl_malloc(sizeof(double)*mrow*nprob);
  prob->rhs_cu = (double*) gkyl_cu_malloc(sizeof(double)*mrow*nprob);

  prob->x = (double*) gkyl_malloc(sizeof(double)*mrow*nprob);
  prob->x_cu = (double*) gkyl_cu_malloc(sizeof(double)*mrow*nprob);

  // create solver
  cusolverSpCreate(&prob->cusolverH);

  cudaStreamCreateWithFlags(&prob->stream, cudaStreamNonBlocking);
  cusolverSpSetStream(prob->cusolverH, prob->stream);

  // create matrix A
  cusparseCreateMatDescr(&prob->A);

  cusparseSetMatType(prob->A, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(prob->A, CUSPARSE_INDEX_BASE_ZERO); 

  cusolverSpCreateCsrqrInfo(&prob->info);

  return prob;
}

void
gkyl_cusolver_amat_from_triples(gkyl_cusolver_prob *prob, gkyl_mat_triples *tri)
{
  // triples must be in rowmaj order for cusolver
  assert(gkyl_mat_triples_is_rowmaj(tri));
  prob->nnz = gkyl_mat_triples_size(tri);

  // Use CSR format
  double *nzval = (double*) gkyl_malloc(sizeof(double)*(prob->nnz)); // non-zero matrix elements.
  int *colind = (int*) gkyl_malloc(sizeof(int)*prob->nnz); // col index of entries in nzval.
  int *rowptr = (int*) gkyl_malloc(sizeof(int)*(prob->mrow+1)); // 1st entry of each row as index in nzval.

  bool *rowptr_assigned = (bool*) gkyl_malloc(sizeof(bool)*prob->mrow);
  for (size_t i=0; i<prob->mrow; i++)
    rowptr_assigned[i] = false;

  // Sorted (row-major order) keys (linear indices to flattened matrix).
  gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri);
  for (size_t i=0; i<prob->nnz; ++i) {
    gkyl_mat_triples_iter_next(iter); // bump iterator.
    struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);
    size_t idx[2] = { mt.row, mt.col };
    
    nzval[i] = mt.val;
    colind[i] = idx[1];
    if (!rowptr_assigned[idx[0]]) {
      rowptr[idx[0]] = i;
      rowptr_assigned[idx[0]] = true;
    }
  }
  rowptr[prob->mrow] = prob->nnz;

  gkyl_mat_triples_iter_release(iter);
  gkyl_free(rowptr_assigned);

  // copy arrays to device
  prob->nzval_cu = (double*) gkyl_cu_malloc(sizeof(double)*(prob->nnz)); // non-zero matrix elements.
  prob->colind_cu = (int*) gkyl_cu_malloc(sizeof(int)*prob->nnz); // col index of entries in nzval.
  prob->rowptr_cu = (int*) gkyl_cu_malloc(sizeof(int)*(prob->mrow+1)); // 1st entry of each row as index in nzval.
  gkyl_cu_memcpy(prob->nzval_cu, nzval, sizeof(double)*prob->nnz, GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(prob->colind_cu, colind, sizeof(int)*prob->nnz, GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(prob->rowptr_cu, rowptr, sizeof(int)*(prob->mrow+1), GKYL_CU_MEMCPY_H2D);

  // symbolic analysis
  cusolverSpXcsrqrAnalysisBatched(prob->cusolverH, prob->mrow, prob->ncol, prob->nnz, prob->A, prob->rowptr_cu, prob->colind_cu, prob->info);

  // prepare working space
  cusolverSpDcsrqrBufferInfoBatched(prob->cusolverH, prob->mrow, prob->ncol, prob->nnz, prob->A, prob->nzval_cu,
                                    prob->rowptr_cu, prob->colind_cu, prob->nrhs, prob->info,
                                    &prob->size_internal, &prob->size_qr);

  
  prob->buffer_qr = (void*) gkyl_cu_malloc(prob->size_qr);

  gkyl_free(colind);
  gkyl_free(nzval);
  gkyl_free(rowptr);
}

void
gkyl_cusolver_brhs_from_triples(gkyl_cusolver_prob *prob, gkyl_mat_triples *tri)
{
  long nnz_rhs = gkyl_mat_triples_size(tri);  // number of non-zero entries in RHS matrix B
  
  // sorted (column-major order) keys (linear indices to flattened matrix)
  gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri);
  for (size_t i=0; i<nnz_rhs; i++) {
    gkyl_mat_triples_iter_next(iter); // bump iterator
    struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);    
    prob->rhs[i] = mt.val;
  }
  gkyl_mat_triples_iter_release(iter);
  
  gkyl_cu_memcpy(prob->rhs_cu, prob->rhs, sizeof(double)*prob->mrow*prob->nrhs, GKYL_CU_MEMCPY_H2D);
}

void
gkyl_cusolver_solve(gkyl_cusolver_prob *prob)
{
  cusolverSpDcsrqrsvBatched(prob->cusolverH, prob->mrow, prob->ncol, prob->nnz, prob->A, prob->nzval_cu, prob->rowptr_cu,
                            prob->colind_cu, prob->rhs_cu, prob->x_cu, prob->nrhs, prob->info, prob->buffer_qr);
}

void
gkyl_cusolver_finish_host(gkyl_cusolver_prob *prob)
{
  cudaStreamSynchronize(prob->stream);
  gkyl_cu_memcpy(prob->x, prob->x_cu, sizeof(double)*prob->mrow*prob->nrhs, GKYL_CU_MEMCPY_D2H);
}

double*
gkyl_cusolver_get_rhs_ptr(gkyl_cusolver_prob *prob, const long loc)
{
  return &prob->rhs_cu[loc];
}

double*
gkyl_cusolver_get_sol_ptr(gkyl_cusolver_prob *prob, const long loc)
{
  return &prob->x_cu[loc];
}

double
gkyl_cusolver_get_sol_ij(gkyl_cusolver_prob *prob, const long ielement, const long jprob)
{
  return prob->x[jprob*prob->mrow+ielement];
}


double
gkyl_cusolver_get_sol_lin(gkyl_cusolver_prob *prob, const long loc)
{
  return prob->x[loc];
}

void
gkyl_cusolver_prob_release(gkyl_cusolver_prob *prob)
{
  gkyl_cu_free(prob->x_cu);
  gkyl_cu_free(prob->rhs_cu);
  gkyl_cu_free(prob->buffer_qr);
  gkyl_cu_free(prob->colind_cu);
  gkyl_cu_free(prob->nzval_cu);
  gkyl_cu_free(prob->rowptr_cu);

  gkyl_free(prob->x);
  gkyl_free(prob->rhs);

  cusolverSpDestroy(prob->cusolverH);
  cudaStreamDestroy(prob->stream);
  gkyl_free(prob);
}
