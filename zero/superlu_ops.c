#include <gkyl_superlu_ops.h>
#include <gkyl_alloc.h>
#include <stdbool.h>

struct gkyl_superlu_prob {
  SuperMatrix A, B; // matrices in Ax=B problem.
  SuperMatrix L, U; // L and U factors in LU decomposition.
  double   *a; // non-zero elements in A.
  int      *asub, *xa; // row index of entries in a, & element in xa corresponding to 1st entry of each column.
  double   *rhs; // right-hand side entries. 
  int      *perm_r; // row permutations from partial pivoting.
  int      *perm_c; // column permutation vector.
  int      mrow, ncol; // A is a mrow x ncol matrix.
  int      nnz; // number of non-zero entries in A.
  int      nrhs; // number of problems to solve (B is an mrow x nrhs matrix).
  int      info, permc_spec;
  superlu_options_t options;
  SuperLUStat_t stat;
};

gkyl_superlu_prob*
gkyl_superlu_prob_new(const int mrow, const int ncol, const int nprob)
{
  struct gkyl_superlu_prob *prob = gkyl_malloc(sizeof(struct gkyl_superlu_prob));

  prob->mrow = mrow;
  prob->ncol = ncol;
  prob->nrhs = nprob;

  if ( !(prob->rhs = doubleMalloc(mrow * nprob)) ) ABORT("superlu_ops: Malloc fails for rhs[].");

  if ( !(prob->perm_r = intMalloc(mrow)) ) ABORT("superlu_ops: Malloc fails for perm_r[].");
  if ( !(prob->perm_c = intMalloc(ncol)) ) ABORT("superlu_ops: Malloc fails for perm_c[].");

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

  if ( !(prob->a = doubleMalloc(prob->nnz)) ) ABORT("superlu_ops: Malloc fails for a[].");
  if ( !(prob->asub = intMalloc(prob->nnz)) ) ABORT("superlu_ops: Malloc fails for asub[].");
  if ( !(prob->xa = intMalloc(prob->ncol+1)) ) ABORT("superlu_ops: Malloc fails for xa[].");

  bool *colptr_assigned = gkyl_malloc(prob->ncol*sizeof(bool));
  for (size_t i=0; i<prob->ncol; i++) {
    colptr_assigned[i] = false;
  }

  long *skeys = gkyl_mat_triples_keys_colmo(tri); // sorted (column-major order) keys (linear indices to flattened matrix).
  for (size_t i=0; i<prob->nnz; i++) {
    int idx[2];
    gkyl_mat_triples_key_to_idx(tri, skeys[i], idx);

    prob->a[i] = gkyl_mat_triples_get(tri, idx[0], idx[1]);
    prob->asub[i] = idx[0];
    if (!colptr_assigned[idx[1]]) {
      prob->xa[idx[1]] = i;
      colptr_assigned[idx[1]] = true;
    }
  }
  prob->xa[prob->ncol] = prob->nnz;
  gkyl_free(skeys);
  gkyl_free(colptr_assigned);

  // Create matrix A. See SuperLU manual for definitions.
  dCreate_CompCol_Matrix(&prob->A, prob->mrow, prob->ncol, prob->nnz,
    prob->a, prob->asub, prob->xa, SLU_NC, SLU_D, SLU_GE);

}

void gkyl_superlu_print_amat(gkyl_superlu_prob *prob)
{
  dPrint_CompCol_Matrix("A", &prob->A);
}

void
gkyl_superlu_brhs_from_triples(gkyl_superlu_prob *prob, gkyl_mat_triples *tri)
{

  if ( !(prob->rhs = doubleMalloc(prob->mrow*prob->nrhs)) ) ABORT("superlu_ops: Malloc fails for rhs[].");

  long nnz_rhs = gkyl_mat_triples_size(tri);  // number of non-zero entries in RHS matrix B.
  long *skeys = gkyl_mat_triples_keys_colmo(tri); // sorted (column-major order) keys (linear indices to flattened matrix).
  for (size_t i=0; i<nnz_rhs; i++) {
    prob->rhs[i] = gkyl_mat_triples_val_at_key(tri, skeys[i]);
  }
  gkyl_free(skeys);

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
}
