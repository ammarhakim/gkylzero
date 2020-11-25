#ifndef GKYL_RANGE_H
#define GKYL_RANGE_H

#include <gkyl_vargm.h>

// Maximum dimension of grids
#ifndef GKYL_MAX_DIM
# define GKYL_MAX_DIM 6
#endif

/**
 * Series of indexing "functions" to compute linear index into range
 */
#define gkyl_ridx1(r, i1)                       \
    ((r).ac[0]+(i1)*(r).ac[1])
#define gkyl_ridx2(r, i1, i2)                   \
    ((r).ac[0]+((i1)*(r).ac[1]+(i2)*(r).ac[2]))
#define gkyl_ridx3(r, i1, i2, i3)                               \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3]))
#define gkyl_ridx4(r, i1, i2, i3, i4)                                   \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3]+(i4)*(r).ac[4]))
#define gkyl_ridx5(r, i1, i2, i3, i4, i5)                               \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3])+((i4)*(r).ac[4]+(i5)*(r).ac[5]))
#define gkyl_ridx6(r, i1, i2, i3, i4, i5, i6)                           \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3])+((i4)*(r).ac[4]+(i5)*(r).ac[5]+(i6)*(r).ac[6]))

/** Generic indexing macro : works for 1D-6D */
#define gkyl_ridx(r, ...) VFUNC1(gkyl_ridx, r, __VA_ARGS__)

/** Indexing macros taking index defined as array of int */
#define gkyl_ridxn(r, idx) gkyl_range_idx(&(r), idx)

/**
 * Range object, representing an N-dimensional integer index
 * set. Lower and upper limits are inclusive.
 */
struct gkyl_range {
    int ndim; // number of dimension
    int lower[GKYL_MAX_DIM]; // lower bound
    int upper[GKYL_MAX_DIM]; // upper bound (inclusive)
    long volume; // total volume of range

    /* do not access directly */
    int ac[GKYL_MAX_DIM+1]; // coefficients for indexing
    int linIdxZero; // linear index of {0,0,...}
};

/**
 * Iterator object into the range. You can read the 'idx' pointer but
 * must not modify it or any other members of this struct.
 */
struct gkyl_range_iter {
    int idx[GKYL_MAX_DIM]; // current index

    /* do not access directly */
    int isFirst;
    int isEmpty;
    int startIdx[GKYL_MAX_DIM];
    int num_bumps;
    int max_bumps; 
    struct gkyl_range range;
};

/**
 * Initialize new range object.
 *
 * @param rng Range object to initialize
 * @param ndim Dimension of range to create.
 * @param lower Lower indices of range
 * @param upper Upper indices of range
 */
void gkyl_range_init(struct gkyl_range *rng, int ndim,
  const int *lower, const int *upper);

/**
 * Create new range object from specified shape. This sets the lower
 * indices to [0,...] and upper indices to [shape[0]-1, ...].
 *
 * @param rng Range object to initialize
 * @param ndim Dimensiom of range to create.
 * @param shape Shape of region
 */
void gkyl_range_init_from_shape(struct gkyl_range *rng, int ndim,
  const int *shape);

/**
 * Shape in direction dir
 *
 * @param rng Range object
 * @param dir Direction to compute shape
 * @return Shape in direction dit
 */
int gkyl_range_shape(const struct gkyl_range *rng, int dir);

/**
 * Return range which has 'dir' direction shortened to length
 * 'len'. The shortened range has the same dimensions and the same
 * start index in 'dir'.
 *
 * @param srng Shortned range.
 * @param rng Range object to shorten
 * @param dir Direction to shorten
 * @param len Length of shortened direction
 */
void gkyl_range_shorten(struct gkyl_range* srng,
  const struct gkyl_range* rng, int dir, int len);

/**
 * Return range in direction 'dir' which corresponds to the "lower
 * skin" cells.  Lower skin cells refer to the second inner-most layer
 * of cells on the lower end of the range. Location of skin cells is
 * with respect to the number of ghost cells in the range.
 *
 * @param srng Skin range
 * @param range Range object to find lower skin cells of
 * @param dir Direction to find lower skin cells in
 * @param nghost Number of ghost cells to determine location of skin cell region
 */
void gkyl_range_lower_skin(struct gkyl_range* srng,
  const struct gkyl_range* range, int dir, int nghost);

/**
 * Return range in direction 'dir' which corresponds to the "upper
 * skin" cells.  Upper skin cells refer to the second inner-most layer
 * of cells on the upper end of the range. Location of skin cells is
 * with respect to the number of ghost cells in the range.
 *
 * @param srng Skin range
 * @param range Range object to find upper skin cells of
 * @param dir Direction to find upper skin cells in
 * @param nghost Number of ghost cells to determine location of skin cell region
 */
void gkyl_range_upper_skin(struct gkyl_range* srng,
  const struct gkyl_range* range, int dir, int nghost);

/**
 * Return range in direction 'dir' which corresponds to the "lower
 * ghost" cells.  Lower ghost cells refer to the outer-most layer of
 * cells on the lower end of the range. Note that this outermost layer
 * is of size nghost in direction dir.
 *
 * @param srng Ghost range
 * @param range Range object to find lower ghost cells of
 * @param dir Direction to find lower ghost cells in
 * @param nghost Number of ghost cells to determine size of ghost cell region
 */
void gkyl_range_lower_ghost(struct gkyl_range* srng,
  const struct gkyl_range* range, int dir, int nghost);

/**
 * Return range in direction 'dir' which corresponds to the "upper
 * ghost" cells.  Upper ghost cells refer to the outer-most layer of
 * cells on the upper end of the range. Note that this outermost layer
 * is of size nghost in direction dir.
 *
 * @param srng Ghost range
 * @param range Range object to find upper ghost cells of
 * @param dir Direction to find upper ghost cells in
 * @param nghost Number of ghost cells to determine size of ghost cell region
 * @return Pointer to newly created range. Call release() to free.
 */
void gkyl_range_upper_ghost(struct gkyl_range* srng,
  const struct gkyl_range* range, int dir, int nghost);

/**
 * Compute offset given relative index. So for a 2D range, idx[2] = {1,
 * 0} will compute the relative offset from {i, j} to {i+1, j} in the
 * mapping from indices to a linear integer space.
 *
 * @param range Range to find offset in.
 * @param idx Relative index for offset calculation
 * @return Relatice offset to idx.
 */
int gkyl_range_offset(const struct gkyl_range* range, const int *idx);

/**
 * General indexing function. Returns linear index into the index
 * range mapped by 'range'.
 *
 * @param range Range object to index
 * @param idx Index for which to compute linear index
 */
int gkyl_range_idx(const struct gkyl_range* range, const int *idx);

/**
 * Create iterator. The returned iterator can be used in a 'while'
 * loop using gkyl_range_iter_next method to loop over the index set
 * spanned by the region.
 *
 * @param range Range object.
 * @return New iterator object for 'range'
 */
void gkyl_range_new_iter(struct gkyl_range_iter *iter,
  const struct gkyl_range* range);


#endif // GKYL_RANGE_H
