#include <gkyl_alloc.h>
#include <gkyl_mat_triples.h>
#include <gkyl_superlu_ops.h>

#include <stdbool.h>

struct gkyl_superlu_prob {
  SuperMatrix A, B; // matrices in Ax=B problem.
  SuperMatrix L, U; // L and U factors in LU decomposition.
  double *rhs; // right-hand side entries. 
  int *perm_r; // row permutations from partial pivoting.
  int *perm_c; // column permutation vector.
  int mrow, ncol; // A is a mrow x ncol matrix.
  int nnz; // number of non-zero entries in A.
  int nrhs; // number of problems to solve (B is an mrow x nrhs matrix).
  int info, permc_spec;
  superlu_options_t options;
  SuperLUStat_t stat;
};

gkyl_superlu_prob*
gkyl_superlu_prob_new(const int mrow, const int ncol, const int nprob)
{
  struct gkyl_superlu_prob *prob = gkyl_malloc(sizeof(*prob));

  prob->mrow = mrow;
  prob->ncol = ncol;
  prob->nrhs = nprob;

  prob->rhs = doubleMalloc(mrow*nprob);
  prob->perm_r = intMalloc(mrow);
  prob->perm_c = intMalloc(ncol);

  // Set the default input options.
  set_default_options(&prob->options);
  prob->options.ColPerm = NATURAL;

  // Initialize the statistics variables.
  StatInit(&prob->stat);

  return prob;
}

void
gkyl_superlu_amat_from_triples(gkyl_superlu_prob *prob, gkyl_mat_triples *tri)
{
  prob->nnz = gkyl_mat_triples_size(tri);

  // allocate some memory needed in superlu. NOTE: this memory is
  // deleted when Destroy_CompCol_Matrix is called, and so we do not
  // need to do it ourselves.
  double *nzval = doubleMalloc(prob->nnz);
  int *rowind = intMalloc(prob->nnz);
  int *colptr = intMalloc(prob->ncol+1);

  bool *colptr_assigned = gkyl_malloc(sizeof(bool[prob->ncol]));
  for (size_t i=0; i<prob->ncol; i++)
    colptr_assigned[i] = false;

  // sorted (column-major order) keys (linear indices to flattened matrix).
  gkyl_mat_triples_iter *iter = gkyl_mat_triples_iter_new(tri);
  for (size_t i=0; i<prob->nnz; ++i) {
    gkyl_mat_triples_iter_next(iter); // bump iterator
    struct gkyl_mtriple mt = gkyl_mat_triples_iter_at(iter);
    size_t idx[2] = { mt.row, mt.col };
    
    nzval[i] = mt.val;
    rowind[i] = idx[0];
    if (!colptr_assigned[idx[1]]) {
      colptr[idx[1]] = i;
      colptr_assigned[idx[1]] = true;
    }
  }
  colptr[prob->ncol] = prob->nnz;
  
  gkyl_mat_triples_iter_release(iter);
  gkyl_free(colptr_assigned);

  // Create matrix A. See SuperLU manual for definitions.
  dCreate_CompCol_Matrix(&prob->A, prob->mrow, prob->ncol, prob->nnz,
    nzval, rowind, colptr, SLU_NC, SLU_D, SLU_GE);
}

void gkyl_superlu_print_amat(gkyl_superlu_prob *prob)
{
  dPrint_CompCol_Matrix("A", &prob->A);
}

void
gkyl_superlu_brhs_from_triples(gkyl_superlu_prob *prob, gkyl_mat_triples *tri)
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
  
  // Create RHS matrix B. See SuperLU manual for definitions.
  dCreate_Dense_Matrix(&prob->B, prob->mrow, prob->nrhs, prob->rhs, prob->mrow,
    SLU_DN, SLU_D, SLU_GE);
}

void
gkyl_superlu_solve(gkyl_superlu_prob *prob)
{
  dgssv(&prob->options, &prob->A, prob->perm_c, prob->perm_r, &prob->L, &prob->U,
    &prob->B, &prob->stat, &prob->info);
}

double
gkyl_superlu_get_rhs(gkyl_superlu_prob *prob, const long loc)
{
  return prob->rhs[loc];
}

void
gkyl_superlu_prob_release(gkyl_superlu_prob *prob)
{
  SUPERLU_FREE (prob->rhs);
  SUPERLU_FREE (prob->perm_r);
  SUPERLU_FREE (prob->perm_c);
  Destroy_CompCol_Matrix(&prob->A);
  Destroy_SuperMatrix_Store(&prob->B);
  Destroy_SuperNode_Matrix(&prob->L);
  Destroy_CompCol_Matrix(&prob->U);
  StatFree(&prob->stat);
  gkyl_free(prob);
}
