#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_range.h>
#include <gkyl_util.h>

// flags and corresponding bit-masks
enum range_flags { R_IS_SUB_RANGE };
static const uint32_t masks[] =
{ 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

// sub-range flags
#define SET_SUB_RANGE(flags) (flags) |= masks[R_IS_SUB_RANGE]
#define CLEAR_SUB_RANGE(flags) (flags) &= ~masks[R_IS_SUB_RANGE]
#define IS_SUB_RANGE(flags) (((flags) & masks[R_IS_SUB_RANGE]) != 0)

// Computes coefficients for mapping indices in row-major order
static void
calc_rowmajor_ac(struct gkyl_range* range, long ac[])
{
  int ndim = range->ndim;
  ac[ndim] = 1L;
  for (int i=ndim-1; i>=1; --i)
    ac[i] = ac[i+1]*gkyl_range_shape(range, i);
  long start = 0L;
  for (int i=0; i<ndim; ++i)
    start += ac[i+1]*range->lower[i];
  ac[0] = -start;
}

// Computes stuff needed for "skip iterator"
static long
calc_skip_iter(const struct gkyl_range *rng, int *remDir)
{
  int up[GKYL_MAX_DIM];
  for (int d=0; d<rng->ndim; ++d) {
    remDir[d] = 1;
    up[d] = rng->upper[d];
  }
  int d = 0;
  long vol = rng->volume;
  long loidx = gkyl_range_idx(rng, rng->lower);
  long del = gkyl_range_idx(rng,up)-loidx+1;
  while (del != vol) {
    up[d] = rng->lower[d];
    vol /= gkyl_range_shape(rng, d);
    del = gkyl_range_idx(rng,up)-loidx+1;
    remDir[d] = 0;
    d += 1;
  }
  return del;
}

// compute volume, safely (for malformed ranges
static long
calc_volume_safely(int ndim, const int *lower, const int *upper)
{
  int is_zero_vol = 0;
  long vol = 1L;
  for (int i=0; i<ndim; ++i) {
    vol *= upper[i]-lower[i]+1;
    is_zero_vol = GKYL_MAX2(is_zero_vol, upper[i]<lower[i] ? 1 : 0);
  }
  if (is_zero_vol) vol = 0;
  return vol;
}

void
gkyl_range_init(struct gkyl_range *rng, int ndim,
  const int *lower, const int *upper)
{
//  // MF 2023/07/07: commenting this out because it causes seg faults in g2.
//  *rng = (struct gkyl_range) { };
  
  int is_zero_vol = 0;
  rng->ndim = ndim;
  rng->volume = 1L;
  for (int i=0; i<ndim; ++i) {
    rng->ilo[i] = rng->lower[i] = lower[i];
    rng->upper[i] = upper[i];
    rng->volume *= upper[i]-lower[i]+1;
    // need to handle case when upper[i]<lower[i]
    is_zero_vol = GKYL_MAX2(is_zero_vol, upper[i]<lower[i] ? 1 : 0);
  }
  // reset volume if any lower[d] <= upper[d]
  if (is_zero_vol) rng->volume = 0;
  
  calc_rowmajor_ac(rng, rng->ac);
  gkyl_copy_long_arr(GKYL_MAX_DIM+1, rng->ac, rng->iac);

  int idxZero[GKYL_MAX_DIM];
  for (int i=0; i<ndim; ++i) idxZero[i] = 0;
  rng->linIdxZero = gkyl_range_idx(rng, idxZero);

  rng->nsplit = 1;
  rng->tid = 0;

  rng->flags = 0;

  // for CUDA ops
  rng->nthreads = GKYL_DEFAULT_NUM_THREADS;
  rng->nblocks = rng->volume/rng->nthreads + 1;
}

void
gkyl_range_init_from_shape(struct gkyl_range *rng, int ndim, const int *shape)
{
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  for (int i=0; i<ndim; ++i) {
    lo[i] = 0; // lower-left corner has index (0,0,...)
    up[i] = shape[i]-1;
  }
  gkyl_range_init(rng, ndim, lo, up);
}

void
gkyl_range_ten_prod(struct gkyl_range *rng, const struct gkyl_range *a, const struct gkyl_range *b)
{
  int adim = a->ndim, bdim = b->ndim;
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  for (int d=0; d<adim; ++d) {
    lower[d] = a->lower[d];
    upper[d] = a->upper[d];
  }
  for (int d=0; d<bdim; ++d) {
    lower[adim+d] = b->lower[d];
    upper[adim+d] = b->upper[d];
  }
  gkyl_range_init(rng, adim+bdim, lower, upper);
}

void
gkyl_range_shift(struct gkyl_range *rng, const struct gkyl_range *inp,
  const int *delta)
{
  int lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];

  for (int d=0; d<inp->ndim; ++d) {
    lower[d] = inp->lower[d] + delta[d];
    upper[d] = inp->upper[d] + delta[d];
  }
  gkyl_range_init(rng, inp->ndim, lower, upper);
}

int
gkyl_range_is_sub_range(const struct gkyl_range *rng)
{
  return IS_SUB_RANGE(rng->flags);
}

int
gkyl_range_contains_idx(const struct gkyl_range *rng, const int *idx)
{
  for (int i=0; i<rng->ndim; ++i) {
    if ( (idx[i] < rng->lower[i]) || (idx[i] > rng->upper[i]) )
      return 0;
  }
  return 1;
}

void
gkyl_sub_range_init(struct gkyl_range *rng,
  const struct gkyl_range *bigrng, const int *sublower, const int *subupper)
{
  rng->ndim = bigrng->ndim;
  rng->volume = 1L;
  for (int i=0; i<rng->ndim; ++i) {
    rng->lower[i] = sublower[i] >= bigrng->lower[i] ? sublower[i] : bigrng->lower[i];
    rng->upper[i] = subupper[i] <= bigrng->upper[i] ? subupper[i] : bigrng->upper[i];
    rng->ilo[i] = bigrng->ilo[i]; // so inv indexer works correctly
    rng->volume *= rng->upper[i]-rng->lower[i]+1;
  }
  for (int i=0; i<rng->ndim+1; ++i)
    rng->ac[i] = bigrng->ac[i];
  rng->linIdxZero = bigrng->linIdxZero;

  rng->nsplit = bigrng->nsplit;
  rng->tid = bigrng->tid;
  
  rng->flags = bigrng->flags;
  SET_SUB_RANGE(rng->flags);

  // we need to construct iac such that sub_range_inv_idx works
  // properly
  struct gkyl_range sub_range;
  gkyl_range_init(&sub_range, rng->ndim, rng->lower, rng->upper);
  gkyl_copy_long_arr(GKYL_MAX_DIM+1, sub_range.ac, rng->iac);

  // for CUDA ops
  rng->nthreads = GKYL_DEFAULT_NUM_THREADS;
  rng->nblocks = rng->volume/rng->nthreads + 1;
}

struct gkyl_range
gkyl_range_split(struct gkyl_range *rng, int nsplit, int tid)
{
  struct gkyl_range r = *rng;
  r.nsplit = nsplit;
  r.tid = tid;
  return r;
}

// Computes split and returns number of elements handled locally and
// the initial index into the range. Number of elements is returned
// and start index set in 'lower'
static long
range_calc_split(const struct gkyl_range *rng, int *lower)
{
  const int nsplit = rng->nsplit, tid = rng->tid;
  long quot = rng->volume/nsplit, rem = rng->volume % nsplit;
  
  long len = gkyl_range_split_len(rng);
  long start = tid < rem ? tid*(quot+1) : rem*(quot+1) + (tid-rem)*quot;

  if (IS_SUB_RANGE(rng->flags)) {
    // as 'start' in sub-range we need to use an additional
    // indirection to compute the 'lower' bounds
    struct gkyl_range subrange;
    gkyl_range_init(&subrange, rng->ndim, rng->lower, rng->upper);
    gkyl_range_inv_idx(&subrange, start, lower);
  }
  else {
    gkyl_range_inv_idx(rng, start, lower);
  }

  return len;
}

long
gkyl_range_split_len(const struct gkyl_range *rng)
{
  const long quot = rng->volume/rng->nsplit, rem = rng->volume % rng->nsplit;
  return rng->tid < rem ? quot+1 : quot;
}

void
gkyl_range_deflate(struct gkyl_range* srng,
  const struct gkyl_range* rng, const int *remDir, const int *locDir)
{
  srng->linIdxZero = rng->linIdxZero;
  srng->ndim = 0;
  srng->volume = 1;  
  for (int i=0, j=0; i<rng->ndim; ++i) {
    if (!remDir[i]) {
      srng->lower[j] = rng->lower[i];
      srng->upper[j] = rng->upper[i];
      srng->ilo[j] = rng->ilo[i];
      srng->ac[j+1] = rng->ac[i+1];      
      srng->ndim += 1;
      srng->volume *= gkyl_range_shape(rng, i);
      j += 1;
    }
  }
  long adel = 0; // need to adjust ac[0]
  for (int i=0; i<rng->ndim; ++i)
    if (remDir[i])
      adel += locDir[i]*rng->ac[i+1];
  srng->ac[0] = rng->ac[0] + adel;

  srng->nsplit = rng->nsplit;
  srng->tid = rng->tid;

  srng->flags = rng->flags;
  SET_SUB_RANGE(srng->flags);

  // for CUDA ops
  srng->nthreads = GKYL_DEFAULT_NUM_THREADS;
  srng->nblocks = srng->volume/srng->nthreads + 1;
}

void
gkyl_range_shorten_from_above(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int len)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM] = {0}, up[GKYL_MAX_DIM] = {0};
  
  for (int i=0; i<ndim; ++i) {
    lo[i] = range->lower[i];
    up[i] = range->upper[i];
  }
  up[dir] = lo[dir]+len-1;
  gkyl_sub_range_init(rng, range, lo, up);
}

void
gkyl_range_shorten_from_below(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int len)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM] = {0}, up[GKYL_MAX_DIM] = {0};
  
  for (int i=0; i<ndim; ++i) {
    lo[i] = range->lower[i];
    up[i] = range->upper[i];
  }
  lo[dir] = up[dir]-len+1;
  gkyl_sub_range_init(rng, range, lo, up);
}

void
gkyl_range_extend(struct gkyl_range *erng,
  const struct gkyl_range* range, const int *elo, const int *eup)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM] = {0}, up[GKYL_MAX_DIM] = {0};

  for (int i=0; i<ndim; ++i) {
    lo[i] = range->lower[i]-elo[i];
    up[i] = range->upper[i]+eup[i];
  }
  gkyl_range_init(erng, ndim, lo, up);
}

void
gkyl_range_lower_skin(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int nskin)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  
  for (int i=0; i<ndim; ++i) {
    lo[i] = range->lower[i];
    up[i] = range->upper[i];
  }
  up[dir] = range->lower[dir]+nskin-1;
  gkyl_sub_range_init(rng, range, lo, up);
}

void
gkyl_range_upper_skin(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int nskin)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  
  for (int i=0; i<ndim; ++i) {
    lo[i] = range->lower[i];
    up[i] = range->upper[i];
  }
  lo[dir] = range->upper[dir]-nskin+1;
  gkyl_sub_range_init(rng, range, lo, up);
}

// Increment an int vector by fact*del[d] in each direction d.
static inline void
incr_int_array(int ndim, int fact, const int * GKYL_RESTRICT del,
  const int * GKYL_RESTRICT inp, int *GKYL_RESTRICT out)
{
  for (int i=0; i<ndim; ++i)
    out[i] = inp[i] + fact*del[i];
}

/**
 * Create ghost and skin sub-ranges given parent (extended
 * range). This code is somewhat convoluted as the skin and ghost
 * ranges need to be sub-ranges of the extended range on the grid and
 * not include corners. I am not sure how to handle corners on
 * physical boundaries. Also, perhaps this code could be simplified.
 */
void
gkyl_skin_ghost_ranges(struct gkyl_range *skin, struct gkyl_range *ghost,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_range *parent, const int *nghost)
{
  int ndim = parent->ndim;
  int lo[GKYL_MAX_DIM] = {0}, up[GKYL_MAX_DIM] = {0};

  if (edge == GKYL_LOWER_EDGE) {

    incr_int_array(ndim, 1, nghost, parent->lower, lo);
    incr_int_array(ndim, -1, nghost, parent->upper, up);
    
    up[dir] = lo[dir]+nghost[dir]-1;
    gkyl_sub_range_init(skin, parent, lo, up);

    incr_int_array(ndim, 1, nghost, parent->lower, lo);
    incr_int_array(ndim, -1, nghost, parent->upper, up);
    
    lo[dir] = lo[dir]-nghost[dir];
    up[dir] = lo[dir]+nghost[dir]-1;
    gkyl_sub_range_init(ghost, parent, lo, up);
  }
  else {

    incr_int_array(ndim, 1, nghost, parent->lower, lo);
    incr_int_array(ndim, -1, nghost, parent->upper, up);
    
    lo[dir] = up[dir]-nghost[dir]+1;
    gkyl_sub_range_init(skin, parent, lo, up);

    incr_int_array(ndim, 1, nghost, parent->lower, lo);
    incr_int_array(ndim, -1, nghost, parent->upper, up);
    
    up[dir] = up[dir]+nghost[dir]+1;
    lo[dir] = up[dir]-nghost[dir];
    gkyl_sub_range_init(ghost, parent, lo, up);
  }
}

int
gkyl_range_intersect(struct gkyl_range* irng,
  const struct gkyl_range *r1, const struct gkyl_range *r2)
{
  int ndim = r1->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) {
    lo[d] = r1->lower[d] > r2->lower[d] ? r1->lower[d] : r2->lower[d];
    up[d] = r1->upper[d] < r2->upper[d] ? r1->upper[d] : r2->upper[d];
  }
  gkyl_range_init(irng, ndim, lo, up);
  return irng->volume > 0 ? 1 : 0;
}

int
gkyl_sub_range_intersect(struct gkyl_range* irng,
  const struct gkyl_range *r1, const struct gkyl_range *r2)
{
  int ndim = r1->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) {
    lo[d] = r1->lower[d] > r2->lower[d] ? r1->lower[d] : r2->lower[d];
    up[d] = r1->upper[d] < r2->upper[d] ? r1->upper[d] : r2->upper[d];
  }
  
  long vol = irng->volume = calc_volume_safely(ndim, lo, up);
  if (vol > 0)
    gkyl_sub_range_init(irng, r1, lo, up);
  else
    gkyl_range_init(irng, ndim, lo, up);
  return irng->volume > 0 ? 1 : 0;
}

bool
gkyl_range_is_on_lower_edge(int dir, const struct gkyl_range *range,
  const struct gkyl_range *parent)
{
  if (range->lower[dir] == parent->lower[dir])
    return true;
  return false;
  
}

bool
gkyl_range_is_on_upper_edge(int dir, const struct gkyl_range *range,
  const struct gkyl_range *parent)
{
  if (range->upper[dir] == parent->upper[dir])
    return true;
  return false;  
}
    
void
gkyl_range_iter_init(struct gkyl_range_iter *iter,
  const struct gkyl_range* range)
{
  iter->is_first = 1;
  iter->ndim = range->ndim;
  iter->bumps_left = range->volume > 0? range_calc_split(range, iter->idx) : 0;
  
  for (int i=0; i<range->ndim; ++i) {
    iter->lower[i] = range->lower[i];
    iter->upper[i] = range->upper[i];
  }
}

void
gkyl_range_iter_no_split_init(struct gkyl_range_iter *iter,
  const struct gkyl_range* range)
{
  iter->is_first = 1;
  iter->ndim = range->ndim;
  iter->bumps_left = range->volume > 0? range_calc_split(range, iter->idx) : 0;
  
  for (int i=0; i<range->ndim; ++i) {
    iter->idx[i] = iter->lower[i] = range->lower[i];
    iter->upper[i] = range->upper[i];
  }  
}

int
gkyl_range_iter_next(struct gkyl_range_iter *iter)
{
  if (iter->bumps_left-- < 1) return 0;
  
  if (iter->is_first) {
    iter->is_first = 0;
    return 1;
  }
  for (int dir=iter->ndim-1; dir>=0; --dir) {
    iter->idx[dir] += 1;
    if (iter->idx[dir] > iter->upper[dir])
      iter->idx[dir] = iter->lower[dir];
    else
      return 1;
  }
  return 0;
}

void
gkyl_range_skip_iter_init(struct gkyl_range_skip_iter *iter,
  const struct gkyl_range* range)
{
  int remDir[GKYL_MAX_DIM];
  iter->delta = calc_skip_iter(range, remDir);
  gkyl_range_deflate(&iter->range, range, remDir, range->lower);
}

void
gkyl_print_range(const struct gkyl_range* range, const char *nm, FILE *fp)
{
  fprintf(fp, "%s = { ", nm);

  fprintf(fp, " lower = { ");
  for (int d=0; d<range->ndim; ++d)
    fprintf(fp, "%d%c ", range->lower[d], d==range->ndim-1 ? ' ' : ',');
  fprintf(fp, "}, ");

  fprintf(fp, "upper = { ");
  for (int d=0; d<range->ndim; ++d)
    fprintf(fp, "%d%c ", range->upper[d] , d==range->ndim-1 ? ' ' : ',');
  fprintf(fp, "}, ");

  fprintf(fp, " volume = %ld, ", range->volume );
  fprintf(fp, " is_sub_range = %d", gkyl_range_is_sub_range(range) );
  
  fprintf(fp, " }\n");
  fflush(fp);
}
