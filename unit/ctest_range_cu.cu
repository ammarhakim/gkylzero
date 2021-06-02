/* -*- c -*- */

#include <gkyl_alloc.h>
extern "C" {
#include <gkyl_range.h>
}

int
my_gkyl_range_shape(const struct gkyl_range *range, int dir)
{
  return range->upper[dir]-range->lower[dir]+1;
}

#include <stdio.h>

extern "C" {
    void cu_range_test(const struct gkyl_range *rng);
}

__global__
void ker_cu_range_test(const struct gkyl_range *rng)
{
  printf("%d. (%d %d)\n", rng->ndim, rng->lower[0], rng->upper[0]);
  int s0 = gkyl_range_shape(rng, 0);
}

void cu_range_test(const struct gkyl_range *rng)
{
  ker_cu_range_test<<<1,1>>>(rng);
}


