#include <assert.h>
#include <stdint.h>
#include <gkyl_array_rio.h>

void
gkyl_array_write(const struct gkyl_array *arr, FILE *fp)
{
  uint64_t rank = arr->rank;
  fwrite(&rank, sizeof(uint64_t), 1, fp);
  fwrite(arr->shape, sizeof(size_t), rank, fp);
  fwrite(arr->data, arr->elemSz*arr->size, 1, fp);
}

void
gkyl_sub_array_write(const struct gkyl_range *range,
  const struct gkyl_array *arr, FILE *fp)
{
  assert(range->volume <= arr->size);
#define _F(loc) gkyl_array_fetch(arr, loc)

  uint64_t rank = arr->rank;
  fwrite(&rank, sizeof(uint64_t), 1, fp);
  fwrite(arr->shape, sizeof(size_t), rank, fp);  
  
  // construct a skip iterator to allow writing (potentially) in
  // chunks rather than element by element or requiring a copy of the
  // data to be made
  struct gkyl_range_skip_iter skip;
  gkyl_range_skip_iter_init(&skip, range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skip.range);
  
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&skip.range, iter.idx);
    fwrite(_F(start), arr->elemSz*skip.delta, 1, fp);
  }
#undef _F
}
