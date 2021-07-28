#pragma once

/** Flags for indicating (conjugate) transpose */
enum gkyl_mat_trans { GKYL_NO_TRANS, GKYL_TRANS, GKYL_CONJ_TRANS };

/**
 * Matrix object. Stored in column major order.
 */
struct gkyl_mat {
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
 * Clone matrix.
 */
struct gkyl_mat* gkyl_mat_clone(const struct gkyl_mat *in);

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
 * Release matrix
 *
 * @param mat Pointer to matrix to release
 */
void gkyl_mat_release(struct gkyl_mat *mat);
