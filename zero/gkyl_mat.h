#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

/** Flags for indicating (conjugate) transpose */
enum gkyl_mat_trans { GKYL_NO_TRANS, GKYL_TRANS, GKYL_CONJ_TRANS };

/**
 * Matrix object. Stored in column major order.
 */
struct gkyl_mat {
  size_t nr, nc; // Number of rows, columns
  double *data; // Pointer to data
};

/**
 * Multiple matrices, each stored in column major order
 */
struct gkyl_nmat {
  size_t num; // Number of matrices
  size_t nr, nc; // Number of rows, columns
  double data[]; // Pointer to data
};

/**
 * Construct new matrix with all elements initialized to @a
 * val. Delete using gkyl_mat_release method.
 *
 * @param nr Number of rows
 * @param nc Number of cols
 * @param val Initial value
 * @return Pointer to new matrix.
 */
struct gkyl_mat* gkyl_mat_new(size_t nr, size_t nc, double val);

/**
 * Construct new multi-matrix (batch of matrices) with all elements in
 * each matric initialized to @a val. Delete using gkyl_nmat_release
 * method. Each matrix has the same shape.
 *
 * @param num Number of matrices
 * @param nr Number of rows
 * @param nc Number of cols
 * @param val Initial value
 * @return Pointer to new multi-matrix.
 */
struct gkyl_nmat* gkyl_nmat_new(size_t num, size_t nr, size_t nc, double val);

/**
 * Clone matrix.
 */
struct gkyl_mat* gkyl_mat_clone(const struct gkyl_mat *in);

/**
 * Get a matrix from multi-matrix. DO NOT free the returned matrix!
 *
 * @param n Matrix to fetch
 * @return Matrix (DO NOT free/release this)
 */
static inline struct gkyl_mat
gkyl_nmat_get(struct gkyl_nmat *mat, size_t num)
{
  return (struct gkyl_mat) {
    .nr = mat->nr,
    .nc = mat->nc,
    .data = mat->data+num*mat->nr*mat->nc
  };
}

/**
 * Set value in matrix.
 */
static inline void
gkyl_mat_set(struct gkyl_mat *mat, size_t r, size_t c, double val)
{
  mat->data[c*mat->nr+r] = val;
}

/**
 * Get value from matrix.
 */
static inline double
gkyl_mat_get(const struct gkyl_mat *mat, size_t r, size_t c)
{
  return mat->data[c*mat->nr+r];
}

/**
 * Get column of matrix as const pointer.
 */
static inline const double*
gkyl_mat_get_ccol(const struct gkyl_mat *mat, size_t c)
{
  return mat->data+c*mat->nr;
}

/**
 * Get column of matrix as pointer.
 */
static inline double*
gkyl_mat_get_col(struct gkyl_mat *mat, size_t c)
{
  return mat->data+c*mat->nr;
}

/**
 * Set all elements of matrix to specified value. Returns pointer to @a mat.
 */ 
struct gkyl_mat* gkyl_mat_clear(struct gkyl_mat *mat, double val);

/**
 * Set all elements on diagonal to specified value. All other elements
 * are set to 0.0. Returns pointer to @a mat.
 */ 
struct gkyl_mat* gkyl_mat_diag(struct gkyl_mat *mat, double val);

/**
 * Write matrix to file. Output is in Maxima matrix format
 */
void gkyl_mat_show(const char *name, FILE *fp, const struct gkyl_mat *mat);

/**
 * Computes matrix-matrix product:
 *
 * C = alpha*OP(A)*OP(B) + beta*C
 *
 * where OP(A) indicates transpose/no-transpose based on the
 * transa/transb flags.
 *
 * C is returned
 */
struct gkyl_mat* gkyl_mat_mm(double alpha, double beta,
  enum gkyl_mat_trans transa, const struct gkyl_mat *A,
  enum gkyl_mat_trans transb, const struct gkyl_mat *B, struct gkyl_mat *C);

/**
 * Solve system of linear equations using LU decomposition. On input
 * the RHS must be in the "x" matrix (each column represents a RHS
 * vector) and on output "x" is replaced with the solution(s). Returns
 * true on success, false otherwise. Note that on output A is replaced
 * by its LU factors.
 *
 * The ipiv input is an chunk of memory that is sizeof(lapack_int[N]),
 * where N is the number of equations. It is safest to assume
 * lapack_int is long (it may be smaller). You must
 * allocate/deallocate ipiv yourself! Use:
 *
 * ipiv = gkyl_mem_buff_new(sizeof(long[N]));
 * gkyl_mat_linsolve_lu(A, x, gkyl_mem_buff_data(ipiv));
 * gkyl_mem_buff_release(ipiv);
 *
 * The reason for passing ipiv to this function is that it avoids
 * allocations inside this function.
 */
bool gkyl_mat_linsolve_lu(struct gkyl_mat *A, struct gkyl_mat *x, void* ipiv);

/**
 * Release matrix
 *
 * @param mat Pointer to matrix to release
 */
void gkyl_mat_release(struct gkyl_mat *mat);

/**
 * Release multi-matrix
 *
 * @param mat Pointer to multi-matrix to release
 */
void gkyl_nmat_release(struct gkyl_nmat *nmat);
