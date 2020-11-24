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
  range->ac[ndim] = 1;
  for (int i=ndim-1; i>=1; --i)
    range->ac[i] = range->ac[i+1]*gkyl_range_shape(range, i);
  int start = 0;
  for (int i=0; i<ndim; ++i)
    start += range->ac[i+1]*range->lower[i];
  range->ac[0] = -start;
}

void
gkyl_range_init(struct gkyl_range *rng, int ndim, int *lower, int *upper)
{
  rng->ndim = ndim;
  rng->volume = 1;
  for (unsigned i=0; i<ndim; ++i) {
    rng->lower[i] = lower[i];
    rng->upper[i] = upper[i];
    rng->volume *= (upper[i]-lower[i]+1);
  }
  calc_rowmajor_ac(rng);

  int idxZero[GKYL_MAX_DIM];
  for (unsigned i=0; i<ndim; ++i) idxZero[i] = 0;
  // TODO: Compute the index of the zero element
}

void
gkyl_range_init_from_shape(struct gkyl_range *rng, int ndim, int *shape)
{
  int lo[GKYL_MAX_DIM], up[GKYL_MAX_DIM];
  for (unsigned i=0; i<ndim; ++i) {
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

int
gkyl_range_offset(const struct gkyl_range* range, const int *idx)
{
  return 0;
  // TODO: gkyl_range_index(range, idx) - range->linIdxZero;
}

void
gkyl_range_new_iter(struct gkyl_range_iter *iter,
  const struct gkyl_range* range)
{
  iter->num_bumps = 0;
  iter->max_bumps = range->volume > 0 ? range->volume : 0;

  for (unsigned i=0; i<range->ndim; ++i) {
    iter->startIdx[i] = range->lower[i];
    iter->idx[i] = range->lower[i];
  }
  memcpy(&iter->range, range, sizeof(struct gkyl_range));
}

