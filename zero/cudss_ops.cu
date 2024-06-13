#include <cuDSS.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_mat_triples.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_cudss_ops.h>
}

#define checkCUDSS(call, status, msg) \
do { \
    status = call; \
    if (status != CUDSS_STATUS_SUCCESS) { \
        printf("Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: " #msg "\n", status); \
        exit(EXIT_FAILURE); \
    } \
} while(0);

struct gkyl_cudss_prob {
  double *rhs, *rhs_cu; // right-hand side vector (reused to store the answer x).
  int nprob; // number of problems to solve.
  int mrow, ncol; // A is a mrow x ncol matrix.
  int nnz; // number of non-zero entries in A.
  int nrhs; // number of columns in B (B is an mrow x nrhs matrix).

  cudaStream_t stream; // CUDA stream cuDSS runs on.
  cudssHandle_t handle; // cuDSS handle.
  cudssConfig_t solverConfig;
  cudssData_t solverData;

  cudssMatrix_t A; // cuDSS object holding the LHS matrix.
  cudssMatrix_t x, b; // cuDSS objects holding the unknowns vector and RHS vector.

  // Arrays used to populate the LHS matrix in CSR format.
  double *csr_val_cu;
  int *csr_rowptr_cu, *csr_colind_cu;
};

gkyl_cudss_prob*
gkyl_cudss_prob_new(int nprob, int mrow, int ncol, int nrhs)
{
  struct gkyl_cudss_prob *prob = (struct gkyl_cudss_prob*) gkyl_malloc(sizeof(*prob));

  prob->nprob = nprob;
  prob->mrow = mrow;
  prob->ncol = ncol;
  prob->nrhs = nrhs;

  prob->rhs = (double*) gkyl_malloc(nrhs * mrow * sizeof(double));
  prob->rhs_cu = (double*) gkyl_cu_malloc(nrhs * mrow * sizeof(double));

  cudssStatus_t status = CUDSS_STATUS_SUCCESS;

  /* Create a CUDA stream */
  prob->stream = NULL;
  checkCuda(cudaStreamCreate(&prob->stream));

  /* Creating the cuDSS library handle */
  checkCUDSS(cudssCreate(&prob->handle), status, "cudssCreate");

  /* (optional) Setting the custom stream for the library handle */
  checkCUDSS(cudssSetStream(prob->handle, prob->stream), status, "cudssSetStream");

  /* Creating cuDSS solver configuration and data objects */
  checkCUDSS(cudssConfigCreate(&prob->solverConfig), status, "cudssConfigCreate");
  checkCUDSS(cudssDataCreate(prob->handle, &prob->solverData), status, "cudssDataCreate");

  /* Create matrix objects for the right-hand side b and solution x (as dense matrices). */
  int64_t mrow_64 = mrow, ncol_64 = ncol;
  int ldb = ncol_64, ldx = mrow_64;
  // Temporary vector. True values set later.
  double *tmp_vec_h = (double*) gkyl_malloc(nrhs * mrow * sizeof(double));
  double *tmp_vec_d = (double*) gkyl_cu_malloc(nrhs * mrow * sizeof(double));
  for (int i=0; i<nrhs * mrow; i++)
    tmp_vec_h[i] = 1.0;
  gkyl_cu_memcpy(tmp_vec_d, tmp_vec_h, nrhs * mrow * sizeof(double), GKYL_CU_MEMCPY_H2D);

  checkCUDSS(cudssMatrixCreateDn(&prob->b, ncol_64, nrhs, ldb, tmp_vec_d, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR),
    status, "cudssMatrixCreateDn for b");
  checkCUDSS(cudssMatrixCreateDn(&prob->x, mrow_64, nrhs, ldx, tmp_vec_d, CUDA_R_64F, CUDSS_LAYOUT_COL_MAJOR),
    status, "cudssMatrixCreateDn for x");

  gkyl_cu_free(tmp_vec_d);
  gkyl_free(tmp_vec_h);

  return prob;
}

void
gkyl_cudss_amat_from_triples(struct gkyl_cudss_prob *prob, struct gkyl_mat_triples **tri)
{
  prob->nnz = gkyl_mat_triples_size(tri[0]);
  for (size_t k=0; k<prob->nprob; k++) {
    assert(gkyl_mat_triples_size(tri[k]) == prob->nnz);  // No. of nonzeros must be the same for every problem.
    assert(gkyl_mat_triples_is_rowmaj(tri[k]));  // Triples must be in rowmaj order for cusolver.
  }

  // Convert triples to CSR arrays on device.
  // Use CSR format
  double *csr_val = (double*) gkyl_malloc(prob->nprob*prob->nnz*sizeof(double)); // non-zero matrix elements.
  int *csr_colind = (int*) gkyl_malloc(sizeof(int)*prob->nnz); // col index of entries in csrvalA.
  int *csr_rowptr = (int*) gkyl_malloc(sizeof(int)*(prob->mrow+1)); // 1st entry of each row as index in csrvalA.

  bool *csr_rowptr_assigned = (bool*) gkyl_malloc(sizeof(bool)*prob->mrow);
  for (size_t i=0; i<prob->mrow; i++) csr_rowptr_assigned[i] = false;

  // Sorted (row-major order) keys (linear indices to flattened matrix).
  for (size_t k=0; k<prob->nprob; k++) {
    gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri[k]);
    for (size_t i=0; i<prob->nnz; ++i) {
      gkyl_mat_triples_iter_next(iter); // bump iterator.
      struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);
      size_t idx[2] = { mt.row, mt.col };

      csr_val[k*prob->nnz+i] = mt.val;
      if (k==0) {
        csr_colind[i] = idx[1];
        if (!csr_rowptr_assigned[idx[0]]) {
          csr_rowptr[idx[0]] = i;
          csr_rowptr_assigned[idx[0]] = true;
        }
      }
    }
    gkyl_mat_triples_iter_release(iter);
  }
  csr_rowptr[prob->mrow] = prob->nnz;
  gkyl_free(csr_rowptr_assigned);

  // Copy arrays to device.
  prob->csr_val_cu = (double*) gkyl_cu_malloc(prob->nprob*prob->nnz*sizeof(double)); // Non-zero matrix elements.
  prob->csr_colind_cu = (int*) gkyl_cu_malloc(sizeof(int)*prob->nnz); // Col index of entries in csrvalA.
  prob->csr_rowptr_cu = (int*) gkyl_cu_malloc(sizeof(int)*(prob->mrow+1)); // 1st entry of each row as index in csrvalA.
  gkyl_cu_memcpy(prob->csr_val_cu, csr_val, prob->nprob*prob->nnz*sizeof(double), GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(prob->csr_colind_cu, csr_colind, sizeof(int)*prob->nnz, GKYL_CU_MEMCPY_H2D);
  gkyl_cu_memcpy(prob->csr_rowptr_cu, csr_rowptr, sizeof(int)*(prob->mrow+1), GKYL_CU_MEMCPY_H2D);

  // Create a matrix object for the sparse input matrix.
  cudssStatus_t status = CUDSS_STATUS_SUCCESS;
//  cudssMatrixType_t mtype     = CUDSS_MTYPE_SPD;
  cudssMatrixType_t mtype     = CUDSS_MTYPE_GENERAL;
  cudssMatrixViewType_t mview = CUDSS_MVIEW_UPPER;
  cudssIndexBase_t base       = CUDSS_BASE_ZERO;
  checkCUDSS(cudssMatrixCreateCsr(&prob->A, prob->mrow, prob->ncol, prob->nnz, prob->csr_rowptr_cu, NULL,
    prob->csr_colind_cu, prob->csr_val_cu, CUDA_R_32I, CUDA_R_64F, mtype, mview,
    base), status, "cudssMatrixCreateCsr");

  // Symbolic factorization.
  checkCUDSS(cudssExecute(prob->handle, CUDSS_PHASE_ANALYSIS, prob->solverConfig, prob->solverData,
    prob->A, prob->x, prob->b), status, "cudssExecute for analysis");

  // Factorization.
  checkCUDSS(cudssExecute(prob->handle, CUDSS_PHASE_FACTORIZATION, prob->solverConfig,
    prob->solverData, prob->A, prob->x, prob->b), status, "cudssExecute for factor");

  gkyl_free(csr_val);
  gkyl_free(csr_colind);
  gkyl_free(csr_rowptr);
}

void
gkyl_cudss_brhs_from_triples(struct gkyl_cudss_prob *prob, gkyl_mat_triples *tri)
{
  long nnz_rhs = gkyl_mat_triples_size(tri);  // Number of non-zero entries in RHS matrix B.

  // Sorted (column-major order) keys (linear indices to flattened matrix).
  gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri);
  for (size_t i=0; i<nnz_rhs; i++) {
    gkyl_mat_triples_iter_next(iter); // bump iterator
    struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);
    prob->rhs[i] = mt.val;
  }
  gkyl_mat_triples_iter_release(iter);

  gkyl_cu_memcpy(prob->rhs_cu, prob->rhs, sizeof(double)*prob->mrow*prob->nrhs, GKYL_CU_MEMCPY_H2D);

  cudssStatus_t status = CUDSS_STATUS_SUCCESS;
  checkCUDSS(cudssMatrixSetValues(prob->b, prob->rhs_cu), status, "cudssMatrixSetValues for setting brhs_from_triples");
}

void
gkyl_cudss_brhs_from_vec(struct gkyl_cudss_prob *prob, double *rhs)
{
  cudssStatus_t status = CUDSS_STATUS_SUCCESS;
  checkCUDSS(cudssMatrixSetValues(prob->b, rhs), status, "cudssMatrixSetValues for setting brhs_from_vec");
}

void
gkyl_cudss_solve(struct gkyl_cudss_prob *prob)
{
  cudssStatus_t status = CUDSS_STATUS_SUCCESS;

  checkCUDSS(cudssExecute(prob->handle, CUDSS_PHASE_SOLVE, prob->solverConfig, prob->solverData,
    prob->A, prob->x, prob->b), status, "cudssExecute for solve");
}

void
gkyl_cudss_sync(struct gkyl_cudss_prob *prob)
{
  cudaStreamSynchronize(prob->stream);
}

void
gkyl_cudss_finish_host(struct gkyl_cudss_prob *prob)
{
  //cudaStreamSynchronize(prob->stream); // not needed when using blocking stream
  cudssStatus_t status = CUDSS_STATUS_SUCCESS;
  checkCUDSS(cudssMatrixGetDn(prob->x, NULL, NULL, NULL, (void**) &(prob->rhs_cu), NULL, NULL), status, "cudssMatrixGetDn in finish_host");
  gkyl_cu_memcpy(prob->rhs, prob->rhs_cu, sizeof(double)*prob->mrow*prob->nrhs, GKYL_CU_MEMCPY_D2H);
}

void
gkyl_cudss_clear_rhs(struct gkyl_cudss_prob *prob, double val)
{
  gkyl_cu_memset(prob->rhs_cu, val, prob->mrow*prob->nrhs*sizeof(double));
}

double*
gkyl_cudss_get_rhs_ptr(struct gkyl_cudss_prob *prob, long loc)
{
  return prob->rhs_cu+loc;
}

double*
gkyl_cudss_get_sol_ptr(struct gkyl_cudss_prob *prob, long loc)
{
  return prob->rhs_cu+loc;
}

double
gkyl_cudss_get_sol_lin(struct gkyl_cudss_prob *prob, long loc)
{
  return prob->rhs[loc];
}

void
gkyl_cudss_prob_release(struct gkyl_cudss_prob *prob)
{
  cudssStatus_t status = CUDSS_STATUS_SUCCESS;

  checkCUDSS(cudssMatrixDestroy(prob->A), status, "cudssMatrixDestroy for A");
  checkCUDSS(cudssMatrixDestroy(prob->b), status, "cudssMatrixDestroy for b");
  checkCUDSS(cudssMatrixDestroy(prob->x), status, "cudssMatrixDestroy for x");
  checkCUDSS(cudssDataDestroy(prob->handle, prob->solverData), status, "cudssDataDestroy");
  checkCUDSS(cudssConfigDestroy(prob->solverConfig), status, "cudssConfigDestroy");
  checkCUDSS(cudssDestroy(prob->handle), status, "cudssHandleDestroy");

  gkyl_cu_free(prob->csr_colind_cu);
  gkyl_cu_free(prob->csr_rowptr_cu);
  gkyl_cu_free(prob->csr_val_cu);

  checkCuda(cudaStreamSynchronize(prob->stream));
  cudaStreamDestroy(prob->stream);

  gkyl_free(prob->rhs);
  gkyl_cu_free(prob->rhs_cu);

  gkyl_free(prob);
}
