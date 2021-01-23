#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_range.h>
#include <gkyl_util.h>

// flags and corresponding bit-masks
enum range_flags { R_IS_SUB_RANGE, R_IS_THREADED };
static uint32_t masks[] =
{ 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

// sub-range flags
#define SET_SUB_RANGE(flags) (flags) |= masks[R_IS_SUB_RANGE]
#define CLEAR_SUB_RANGE(flags) (flags) &= ~masks[R_IS_SUB_RANGE]
#define IS_SUB_RANGE(flags) (((flags) & masks[R_IS_SUB_RANGE]) != 0)

// threading flags
#define SET_THREADED(flags) (flags) |= masks[R_IS_THREADED]
#define CLEAR_THREADED(flags) (flags) &= ~masks[R_IS_THREADED]
#define IS_THREADED(flags) (((flags) & masks[R_IS_THREADED]) != 0)

// Computes coefficients for mapping indices in row-major order
static void
calc_rowmajor_ac(struct gkyl_range* range)
{
  int ndim = range->ndim;
  range->ac[ndim] = 1L;
  for (int i=ndim-1; i>=1; --i)
    range->ac[i] = range->ac[i+1]*gkyl_range_shape(range, i);
  long start = 0L;
  for (int i=0; i<ndim; ++i)
    start += range->ac[i+1]*range->lower[i];
  range->ac[0] = -start;
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

void
gkyl_range_init(struct gkyl_range *rng, int ndim,
  const int *lower, const int *upper)
{
  rng->ndim = ndim;
  rng->volume = 1L;
  for (int i=0; i<ndim; ++i) {
    rng->ilo[i] = rng->lower[i] = lower[i];
    rng->upper[i] = upper[i];
    rng->volume *= upper[i]-lower[i]+1;
  }
  calc_rowmajor_ac(rng);

  int idxZero[GKYL_MAX_DIM];
  for (int i=0; i<ndim; ++i) idxZero[i] = 0;
  rng->linIdxZero = gkyl_range_idx(rng, idxZero);

  rng->th_start = 0;
  rng->th_len = rng->volume;

  rng->flags = 0;
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

int
gkyl_range_shape(const struct gkyl_range *range, int dir)
{
  return range->upper[dir]-range->lower[dir]+1;
}

int gkyl_range_is_sub_range(const struct gkyl_range *rng)
{
  return IS_SUB_RANGE(rng->flags);
}

int gkyl_range_is_threaded(const struct gkyl_range *rng)
{
  return IS_THREADED(rng->flags);
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

  rng->th_start = gkyl_range_idx(bigrng, rng->lower);
  rng->th_len = rng->volume;  

  rng->flags = bigrng->flags;
  SET_SUB_RANGE(rng->flags);
}


void
gkyl_range_thread(struct gkyl_range *rng, int nthreads, int tid)
{
  long offset = gkyl_range_idx(rng, rng->lower);
  long quot = rng->volume/nthreads, rem = rng->volume % nthreads;
  if (tid < rem) {
    rng->th_len = quot+1;
    rng->th_start = offset + tid*(quot+1);
  }
  else {
    rng->th_len = quot;
    rng->th_start = offset + rem*(quot+1) + (tid-rem)*quot;
  }

  // for sub-ranges we need to adjust things
  if (IS_SUB_RANGE(rng->flags)) {
    struct gkyl_range_iter iter;
    gkyl_range_iter_init_ignore_threading(&iter, rng);
    // we need to ignore previous threading information

    long nbumps = rng->th_start-offset+1;
    for (int i=0; i<nbumps; ++i) gkyl_range_iter_next(&iter);
    rng->th_start = gkyl_range_idx(rng, iter.idx);
  }
  SET_THREADED(rng->flags);
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

  srng->th_start = gkyl_range_idx(srng, srng->lower);
  srng->th_len = srng->volume;

  srng->flags = rng->flags;
  SET_SUB_RANGE(srng->flags);
}

void
gkyl_range_shorten(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int len)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  
  for (int i=0; i<ndim; ++i) {
    lo[i] = range->lower[i];
    up[i] = range->upper[i];
  }
  up[dir] = lo[dir]+len-1;
  gkyl_sub_range_init(rng, range, lo, up);
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

long
gkyl_range_offset(const struct gkyl_range* range, const int *idx)
{
  return gkyl_range_idx(range, idx) - range->linIdxZero;
}

long
gkyl_range_idx(const struct gkyl_range* range, const int *idx)
{
#define RI(...) gkyl_ridx(*range, __VA_ARGS__)
  switch (range->ndim) {
    case 0:
        return range->ac[0];
        break;    
    case 1:
        return RI(idx[0]); 
        break;
    case 2:
        return RI(idx[0], idx[1]);
        break;
    case 3:
        return RI(idx[0], idx[1], idx[2]);
        break;
    case 4:
        return RI(idx[0], idx[1], idx[2], idx[3]);
        break;
    case 5:
        return RI(idx[0], idx[1], idx[2], idx[3], idx[4]);
        break;
    case 6:
        return RI(idx[0], idx[1], idx[2], idx[3], idx[4], idx[5]);
        break;
    case 7:
        return RI(idx[0], idx[1], idx[2], idx[3], idx[4], idx[5], idx[6]);
        break;
  }
  return 0;
#undef RI
}

void
gkyl_range_inv_idx(const struct gkyl_range *range, long loc, int *idx)
{
  long n = loc;
  for (int i=1; i<=range->ndim; ++i) {
    long quot = n/range->ac[i];
    long rem = n % range->ac[i];
    idx[i-1] = quot + range->ilo[i-1];
    n = rem;
  }
}

void
gkyl_range_iter_init(struct gkyl_range_iter *iter,
  const struct gkyl_range* range)
{
  iter->is_first = 1;
  iter->ndim = range->ndim;
  iter->bumps_left = range->th_len;
  
  gkyl_range_inv_idx(range, range->th_start, iter->idx);
  for (int i=0; i<range->ndim; ++i) {
    iter->lower[i] = range->lower[i];
    iter->upper[i] = range->upper[i];
  }
}

void
gkyl_range_iter_init_ignore_threading(struct gkyl_range_iter *iter,
  const struct gkyl_range* range)
{
  iter->is_first = 1;
  iter->ndim = range->ndim;
  iter->bumps_left = range->volume;
  
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
