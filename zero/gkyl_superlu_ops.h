#pragma once

#include <slu_ddefs.h>
#include <gkyl_mat_triples.h>

// Object type
typedef struct gkyl_superlu_prob gkyl_superlu_prob;

/**
 * Create a new SuperLU Ax=B problem.
 */
struct gkyl_superlu_prob* gkyl_superlu_prob_new(const int mrow, const int ncol, const int nprob);

/**
 * Initialize SuperLU matrix A in Ax=B problem from a list of triples.
 *
 * @param prob SuperLu struct holding arrays used in problem.
 * @param tri coordinates & values of non-zero entries in A matrix (triplets).
 */
void gkyl_superlu_amat_from_triples(gkyl_superlu_prob *prob, gkyl_mat_triples *tri);

/**
 * Initialize right-hand-side SuperLU matrix B in Ax=B problem from a list of
 * triples.
 *
 * @param prob SuperLu struct holding arrays used in problem.
 * @param tri coordinates & values of non-zero entries in B matrix (triplets).
 */
void gkyl_superlu_brhs_from_triples(gkyl_superlu_prob *prob, gkyl_mat_triples *tri);

/**
 * Solve Ax=B problem.
 *
 * @param prob SuperLu struct holding arrays used in problem.
 */
void gkyl_superlu_solve(gkyl_superlu_prob *prob);

/**
 * Obtain the RHS value at location loc.
 *
 * @param linear index into the RHS flattened array for the desired value.
 * @return RHS value.
 */
double gkyl_superlu_get_rhs(gkyl_superlu_prob *prob, const long loc);

/**
 * Release SuperLU problem
 *
 * @param prob Pointer to SuperLU problem to release.
 */
void gkyl_superlu_prob_release(gkyl_superlu_prob *prob);
