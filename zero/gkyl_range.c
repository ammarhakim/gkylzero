#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_range.h>
#include <gkyl_util.h>

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

// Computes stuff needd for "skip iterator"
long
calc_skip_iter(const struct gkyl_range *rng, int *remDir)
{
  int up[GKYL_MAX_DIM];
  for (unsigned d=0; d<rng->ndim; ++d) up[d] = rng->upper[d];
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
  for (unsigned i=0; i<ndim; ++i) {
    rng->ilo[i] = rng->lower[i] = lower[i];
    rng->upper[i] = upper[i];
    rng->volume *= upper[i]-lower[i]+1;
  }
  calc_rowmajor_ac(rng);

  int idxZero[GKYL_MAX_DIM];
  for (unsigned i=0; i<ndim; ++i) idxZero[i] = 0;
  rng->linIdxZero = gkyl_range_idx(rng, idxZero);
}

void
gkyl_range_init_from_shape(struct gkyl_range *rng, int ndim, const int *shape)
{
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  for (unsigned i=0; i<ndim; ++i) {
    lo[i] = 0; // lower-left corner has index (0,0,...)
    up[i] = shape[i]-1;
  }
  gkyl_range_init(rng, ndim, lo, up);
}

void
gkyl_sub_range_init(struct gkyl_range *rng,
  const struct gkyl_range *bigrng, const int *sublower, const int *subupper)
{
  rng->ndim = bigrng->ndim;
  rng->volume = 1L;
  for (unsigned i=0; i<rng->ndim; ++i) {
    rng->lower[i] = sublower[i] >= bigrng->lower[i] ? sublower[i] : bigrng->lower[i];
    rng->upper[i] = subupper[i] <= bigrng->upper[i] ? subupper[i] : bigrng->upper[i];
    rng->ilo[i] = bigrng->ilo[i]; // so inv indexer works correctly
    rng->volume *= rng->upper[i]-rng->lower[i]+1;
  }
  for (unsigned i=0; i<rng->ndim+1; ++i)
    rng->ac[i] = bigrng->ac[i];
  rng->linIdxZero = bigrng->linIdxZero;
}

int
gkyl_range_shape(const struct gkyl_range *range, int dir)
{
  return range->upper[dir]-range->lower[dir]+1;
}

void
gkyl_range_shorten(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int len)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  
  for (unsigned i=0; i<ndim; ++i) {
    lo[i] = range->lower[i];
    up[i] = range->upper[i];
  }
  up[dir] = lo[dir]+len-1;
  gkyl_range_init(rng, ndim, lo, up);
}

void
gkyl_range_deflate(struct gkyl_range* srng,
  const struct gkyl_range* rng, const int *remDir, const int *locDir)
{
  srng->linIdxZero = rng->linIdxZero;

  srng->ndim = 0;
  srng->volume = 1;  
  for (unsigned i=0, j=0; i<rng->ndim; ++i) {
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
  for (unsigned i=0; i<rng->ndim; ++i) {
    if (remDir[i])
      adel += locDir[i]*rng->ac[i+1];
  }
  srng->ac[0] = rng->ac[0] + adel;
}

void
gkyl_range_lower_skin(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int nghost)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  
  for (unsigned i=0; i<ndim; ++i) {
    lo[i] = range->lower[i];
    up[i] = range->upper[i];
  }
  up[dir] = range->lower[dir]+nghost-1;
  gkyl_range_init(rng, ndim, lo, up);
}

void
gkyl_range_upper_skin(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int nghost)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  
  for (unsigned i=0; i<ndim; ++i) {
    lo[i] = range->lower[i];
    up[i] = range->upper[i];
  }
  lo[dir] = range->upper[dir]-nghost+1;
  gkyl_range_init(rng, ndim, lo, up);
}

void
gkyl_range_lower_ghost(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int nghost)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  
  for (unsigned i=0; i<ndim; ++i) {
    lo[i] = range->lower[i];
    up[i] = range->upper[i];
  }
  lo[dir] = range->lower[dir]-nghost;
  up[dir] = range->lower[dir]-1;
  gkyl_range_init(rng, ndim, lo, up);
}

void
gkyl_range_upper_ghost(struct gkyl_range *rng,
  const struct gkyl_range* range, int dir, int nghost)
{
  int ndim = range->ndim;
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  
  for (unsigned i=0; i<ndim; ++i) {
    lo[i] = range->lower[i];
    up[i] = range->upper[i];
  }
  lo[dir] = range->upper[dir]+1;
  up[dir] = range->upper[dir]+nghost;
  gkyl_range_init(rng, ndim, lo, up);
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
  for (unsigned i=0; i<range->ndim; ++i) {
    iter->idx[i] = iter->lower[i] = range->lower[i];
    iter->upper[i] = range->upper[i];
  }
}

void gkyl_range_iter_reset(struct gkyl_range_iter *iter)
{
  iter->is_first = 1;
  for (unsigned i=0; i<iter->ndim; ++i)
    iter->idx[i] = iter->lower[i];
}

int
gkyl_range_iter_next(struct gkyl_range_iter *iter)
{
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
  for (unsigned d=0; d<range->ndim; ++d) remDir[d] = 1;
  iter->delta = calc_skip_iter(range, remDir);
  gkyl_range_deflate(&iter->range, range, remDir, range->lower);
}
