#pragma once

/** Triples stores list of (i,j,val) for use in contructing sparse matrices */
typedef struct gkyl_mat_triples gkyl_mat_triples;

/**
 * Create new triples list.
 */
gkyl_mat_triples* gkyl_mat_triples_new();

/**
 * Insert value 'val' in triples list at location (i,j)
 */
void gkyl_mat_triples_insert(int i, int j, double val);

/**
 * Accumulate value 'val' in triples list at location (i,j). If an
 * element at this location exists it is incremented by 'val'. New
 * value at the location is returned.
 */
double gkyl_mat_triples_accum(int i, int j, double val);

/**
 * Release triples
 *
 * @param tri Pointer to triples to release.
 */
void gkyl_mat_triples_release(gkyl_mat_triples *tri);
