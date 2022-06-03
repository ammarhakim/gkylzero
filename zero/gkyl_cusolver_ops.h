#pragma once

#include <assert.h>
#include <stdbool.h>
#include "gkyl_alloc_flags_priv.h"

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_mat_triples.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
}
#include <cusparse.h>
#include <cusolverSp.h>

// Object type
typedef struct gkyl_cusolver_prob gkyl_cusolver_prob;

/**
 * Create a new cuSolver Ax=B problem.
 */
struct gkyl_cusolver_prob* gkyl_cusolver_prob_new(const int mrow, const int ncol, const int nprob);

/**
 * Initialize cuSolver matrix A in Ax=B problem from a list of triples.
 *
 * @param prob cuSolver struct holding arrays used in problem.
 * @param tri coordinates & values of non-zero entries in A matrix (triplets).
 */
void gkyl_cusolver_amat_from_triples(gkyl_cusolver_prob *prob, gkyl_mat_triples *tri);

/**
 * Method to print the matrix A to screen.
 */
void gkyl_cusolver_print_amat(gkyl_cusolver_prob *prob);

/**
 * Perform the LU decomposition of the A matrix.
 * The _solve method will use these if they are pre-computed.
 */
void gkyl_cusolver_ludecomp(gkyl_cusolver_prob *prob);

/**
 * Initialize right-hand-side cuSolver matrix B in Ax=B problem from a list of
 * triples.
 *
 * @param prob cuSolver struct holding arrays used in problem.
 * @param tri coordinates & values of non-zero entries in B matrix (triplets).
 */
void gkyl_cusolver_brhs_from_triples(gkyl_cusolver_prob *prob, gkyl_mat_triples *tri);

/**
 * Solve Ax=B problem.
 *
 * @param prob cuSolver struct holding arrays used in problem.
 */
void gkyl_cusolver_solve(gkyl_cusolver_prob *prob);

/**
 * Copy solution back to host
 *
 * @param prob cuSolver struct holding arrays used in problem.
 */
void gkyl_cusolver_finish_host(gkyl_cusolver_prob *prob);

/**
 * Obtain the RHS ielement-th value of the jprob-th linear problem.
 * Recall the RHS is a mrow x nrhs matrix.
 *
 * @param linear index into the RHS flattened array for the desired value.
 * @return RHS value.
 */
double gkyl_cusolver_get_sol_ij(gkyl_cusolver_prob *prob, const long ielement, const long jprob);

/**
 * Obtain the RHS value at location loc (a linear index into the RHS matrix).
 *
 * @param linear index into the RHS flattened array for the desired value.
 * @return RHS value.
 */
double gkyl_cusolver_get_sol_lin(gkyl_cusolver_prob *prob, const long loc);

/**
 * Release cuSolver problem
 *
 * @param prob Pointer to cuSolver problem to release.
 */
void gkyl_cusolver_prob_release(gkyl_cusolver_prob *prob);
