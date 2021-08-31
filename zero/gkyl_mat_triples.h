#pragma once

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
 * Release triples
 *
 * @param tri Pointer to triples to release.
 */
void gkyl_mat_triples_release(gkyl_mat_triples *tri);
