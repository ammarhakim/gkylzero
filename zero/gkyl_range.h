#pragma once

#include <gkyl_util.h>
#include <gkyl_vargm.h>

/**
 * Series of indexing "functions" to compute linear index into range
 */
#define gkyl_ridx1(r, i1)                       \
    ((r).ac[0]+(i1)*(r).ac[1])
#define gkyl_ridx2(r, i1, i2)                   \
    ((r).ac[0]+((i1)*(r).ac[1]+(i2)*(r).ac[2]))
#define gkyl_ridx3(r, i1, i2, i3)                                       \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3]))
#define gkyl_ridx4(r, i1, i2, i3, i4)                                   \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3]+(i4)*(r).ac[4]))
#define gkyl_ridx5(r, i1, i2, i3, i4, i5)                               \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3])+((i4)*(r).ac[4]+(i5)*(r).ac[5]))
#define gkyl_ridx6(r, i1, i2, i3, i4, i5, i6)                           \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3])+((i4)*(r).ac[4]+(i5)*(r).ac[5]+(i6)*(r).ac[6]))
#define gkyl_ridx7(r, i1, i2, i3, i4, i5, i6, i7)                       \
    (((r).ac[0]+(i1)*(r).ac[1])+((i2)*(r).ac[2]+(i3)*(r).ac[3])+((i4)*(r).ac[4]+(i5)*(r).ac[5]+(i6)*(r).ac[6]) + (i7)*(r).ac[7])

/** Generic indexing: works for 1D-7D (VFUNC1 is defined-ed in
 * gkyl_vargm.h) */
#define gkyl_ridx(r, ...) VFUNC1(gkyl_ridx, r, __VA_ARGS__)

/** Indexing macro taking index defined as array of int */
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
    // do not access directly
    int ilo[GKYL_MAX_DIM]; // for use in inverse indexer
    long ac[GKYL_MAX_DIM+1]; // coefficients for indexing
    long linIdxZero; // linear index of {0,0,...}
};

/**
 * Iterator object into the range. You can read the 'idx' pointer but
 * must not modify it or any other members of this struct.
 */
struct gkyl_range_iter {
    int idx[GKYL_MAX_DIM]; // current index
    // do not access directly
    int is_first, ndim;
    int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
};

/**
 * Skip-list based iterator object.
 */
struct gkyl_range_skip_iter {
    long delta; // number of contiguous elements
    struct gkyl_range range; // outer range for iteration
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
 * Create a sub-range from a given range. The sub-range is completely
 * contained inside the parent range.
 *
 * @param rng New range object to initialize
 * @param bigrng Parent range object 
 * @param sublower Lower indices of sub-range
 * @param subupper Upper indices of sub-range
 */
void gkyl_sub_range_init(struct gkyl_range *rng,
  const struct gkyl_range *bigrng, const int *sublower, const int *subupper);

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
 * Return range which has some directions removed by setting the index
 * in those directions to fixed values. The "deflated" range has lower
 * dimension than the parent 'rng' object. The indexing into the
 * returned lower dimensional range gives the same index as the
 * corresponding location in the parent range (with the missing
 * indices set to 'locDir[dir]'). 
 *
 * @param srng Deflated range.
 * @param rng Range object to deflate
 * @param remDir 'ndim' Array of flags: 0 to keep direction, 1 to remove
 * @param loc Index to set removed direction.
 */
void gkyl_range_deflate(struct gkyl_range* srng,
  const struct gkyl_range* rng, const int *remDir, const int *locDir);

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
long gkyl_range_offset(const struct gkyl_range* range, const int *idx);

/**
 * General indexing function. Returns linear index into the index
 * range mapped by 'range'.
 *
 * @param range Range object to index
 * @param idx Index for which to compute linear index
 */
long gkyl_range_idx(const struct gkyl_range* range, const int *idx);

/**
 * Inverse indexer, mapping a linear index to an N-dimension index
 * into 'range' object.
 *
 * @param range Range object to map into
 * @param loc Linear index in [0, range->volume)
 * @param idx On output, the N-dimensional index into 'range'
 */
void gkyl_range_inv_idx(const struct gkyl_range *range, long loc, int *idx);

/**
 * Create iterator. The returned iterator can be used in a 'while'
 * loop using gkyl_range_iter_next method to loop over the index set
 * spanned by the region.
 *
 * @param range Range object.
 * @return New iterator object for 'range'
 */
void gkyl_range_iter_init(struct gkyl_range_iter *iter,
  const struct gkyl_range* range);

/**
 * Reset iterator back to start index. Call this method if you want to
 * reuse an iterator that was already used in a while() loop.
 * 
 * @param iter Iterator to reset.
 */
void gkyl_range_iter_reset(struct gkyl_range_iter *iter);

/**
 * Get next index into range. The iter->idx array holds the next
 * index. This should not be modified by the user!
 *
 * @param iter Iterator object. On exit, iter->idx has the next index
 * @return 1 if there are more indices remaining, 0 if done.
 */
int gkyl_range_iter_next(struct gkyl_range_iter *iter);

/**
 * Create skip-iterator. The returned iterator can be used in a nested
 * loop structure: a for loop inside a 'while' loop.
 *
 * @param range Range object.
 * @return New iterator object for 'range'
 */
void gkyl_range_skip_iter_init(struct gkyl_range_skip_iter *iter,
  const struct gkyl_range* range);
