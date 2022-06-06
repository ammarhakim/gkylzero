#pragma once

#include <gkyl_mat_triples.h>
#include <gkyl_superlu.h>

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
 * Method to print the matrix A to screen.
 */
void gkyl_superlu_print_amat(gkyl_superlu_prob *prob);

/**
 * Perform the LU decomposition of the A matrix.
 * The _solve method will use these if they are pre-computed.
 */
void gkyl_superlu_ludecomp(gkyl_superlu_prob *prob);

/**
 * Initialize right-hand-side SuperLU matrix B in Ax=B problem from a list of
 * triples.
 *
 * @param prob SuperLu struct holding arrays used in problem.
 * @param tri coordinates & values of non-zero entries in B matrix (triplets).
 */
void gkyl_superlu_brhs_from_triples(gkyl_superlu_prob *prob, gkyl_mat_triples *tri);

/**
 * Initialize right-hand-side SuperLU matrix B in Ax=B problem from a simple
 * array. This is enough especially when solving a single Ax=B problem (with B
 * a vector) and makes the GPU interface implementation easier.
 *
 * @param prob SuperLu struct holding arrays used in problem.
 * @param bin column array for the right side source.
 */
void gkyl_superlu_brhs_from_array(gkyl_superlu_prob *prob, const double *bin);

/**
 * Solve Ax=B problem.
 *
 * @param prob SuperLu struct holding arrays used in problem.
 */
void gkyl_superlu_solve(gkyl_superlu_prob *prob);

/**
 * Obtain the RHS ielement-th value of the jprob-th linear problem.
 * Recall the RHS is a mrow x nrhs matrix.
 *
 * @param linear index into the RHS flattened array for the desired value.
 * @return RHS value.
 */
double gkyl_superlu_get_rhs_ij(gkyl_superlu_prob *prob, const long ielement, const long jprob);

/**
 * Obtain the RHS value at location loc (a linear index into the RHS matrix).
 *
 * @param linear index into the RHS flattened array for the desired value.
 * @return RHS value.
 */
double gkyl_superlu_get_rhs_lin(gkyl_superlu_prob *prob, const long loc);

/**
 * Release SuperLU problem
 *
 * @param prob Pointer to SuperLU problem to release.
 */
void gkyl_superlu_prob_release(gkyl_superlu_prob *prob);
