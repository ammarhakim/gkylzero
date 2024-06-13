#pragma once

#include <assert.h>
#include <stdbool.h>

#include <gkyl_alloc.h>
#include <gkyl_mat_triples.h>
#include <gkyl_array_ops.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Object type
typedef struct gkyl_cudss_prob gkyl_cudss_prob;

/**
 * Create a new cuDSS A_i x_i = B_i where i \in {0,1,...,nprob-1},
 * for solving nprob linear problems. This module can be used in two ways:
 *   a) Using nprob=1 and nrhs>=1, so that one can solve nrhs linear problems
 *      with the same left-side matrix A_i for each problem, i.e. B is a matrix
 *      of nprob columns.
 *   b) Using nprob=>1 and nrhs=1, so we are solving nprob separate linear
 *      problems (each with a different left-side Matrix A_i, but with the same
 *      dimensions and sparsity pattern) where each problem only has a
 *      right-side vector with a single column.
 */
struct gkyl_cudss_prob* gkyl_cudss_prob_new(int nprob, int mrow, int ncol, int nrhs);

/**
 * Initialize cuDSS matrix A in Ax=B problem from a list of triples.
 *
 * @param prob cuDSS struct holding arrays used in problem.
 * @param tri (array of) coordinates & values of non-zero entries in A matrix (triplets).
 */
void gkyl_cudss_amat_from_triples(struct gkyl_cudss_prob *prob, struct gkyl_mat_triples **tri);

/**
 * Initialize right-hand-side cuDSS matrix B in Ax=B problem from a list of
 * triples.
 *
 * @param prob cuDSS struct holding arrays used in problem.
 * @param tri coordinates & values of non-zero entries in B matrix (triplets).
 */
void gkyl_cudss_brhs_from_triples(struct gkyl_cudss_prob *prob, gkyl_mat_triples *tri);

/**
 * Solve Ax=B problem.
 *
 * @param prob cuDSS struct holding arrays used in problem.
 */
void gkyl_cudss_solve(struct gkyl_cudss_prob *prob);

/**
 * Copy solution back to host
 *
 * @param prob cuDSS struct holding arrays used in problem.
 */
void gkyl_cudss_finish_host(struct gkyl_cudss_prob *prob);

/**
 * Synchronize the stream cuDSS ran on.
 *
 * @param prob cuDSS struct holding arrays used in problem.
 */
void gkyl_cudss_sync(struct gkyl_cudss_prob *prob);

/**
 * Obtain the RHS value at location loc (a linear index into the RHS matrix).
 *
 * @param linear index into the RHS flattened array for the desired value.
 * @return RHS value.
 */
double gkyl_cudss_get_sol_lin(struct gkyl_cudss_prob *prob, long loc);

/**
 * Release cuDSS problem
 *
 * @param prob Pointer to cuDSS problem to release.
 */
void gkyl_cudss_prob_release(struct gkyl_cudss_prob *prob);
