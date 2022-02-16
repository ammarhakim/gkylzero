#pragma once

#include <stddef.h>

// A ,triple stores 'val' at (row, col)
struct gkyl_mtriple {
  size_t row, col;
  double val;  
};

/** Triples stores list of (i,j,val) for use in contructing sparse matrices */
typedef struct gkyl_mat_triples gkyl_mat_triples;

/**
 * Create new triples list in a matrix of shape (nr, nc).
 *
 * @param nr Number of rows
 * @param nc Number of cols
 * @return Pointer to new empty triples object.
 */
gkyl_mat_triples* gkyl_mat_triples_new(size_t nr, size_t nc);

/**
 * Insert value 'val' in triples list at location (i,j)
 *
 * @return value inserted (val)
 */
double gkyl_mat_triples_insert(gkyl_mat_triples *tri, size_t i, size_t j, double val);

/**
 * Accumulate value 'val' in triples list at location (i,j). If an
 * element at this location exists it is incremented by 'val'. New
 * value at the location is returned.
 */
double gkyl_mat_triples_accum(gkyl_mat_triples *tri, size_t i, size_t j, double val);

/**
 * Returns value in triples list at location (i,j)
 *
 * @return value at (i,j)
 */
double gkyl_mat_triples_get(const gkyl_mat_triples *tri, size_t i, size_t j);

/**
 * Return number of elements inserted.
 */
size_t gkyl_mat_triples_size(const gkyl_mat_triples *tri);

/**
 * Return an array with the keys (location in flattened matrix) of each triple.
 * 
 * @return array with keys.
 */
long *gkyl_mat_triples_keys(const gkyl_mat_triples *tri);

/**
 * Return an array with the keys (location in flattened matrix) of each triple sorted in column-major order.
 * 
 * @return array with keys.
 */
long* gkyl_mat_triples_keys_colmo(const gkyl_mat_triples *tri);

/**
 * Returns value in triples list given its key (location in flattened matrix).
 *
 * @return value with key loc.
 */
double gkyl_mat_triples_val_at_key(const gkyl_mat_triples *tri, long loc);

/**
 * Returns an array with row and column indices corresponding to key loc (location in flattened matrix).
 *
 * @return int[2] with row/column indexes.
 */
void gkyl_mat_triples_key_to_idx(const gkyl_mat_triples *tri, long loc, int idx[2]);

/**
 * Release triples
 *
 * @param tri Pointer to triples to release.
 */
void gkyl_mat_triples_release(gkyl_mat_triples *tri);
