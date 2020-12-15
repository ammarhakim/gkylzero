#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>

#include <gkyl_array_rio.h>

void
gkyl_array_write(const struct gkyl_array *arr, FILE *fp)
{
  uint64_t elemSz = arr->elemSz, size = arr->size;
  fwrite(&elemSz, sizeof(uint64_t), 1, fp);
  fwrite(&size, sizeof(uint64_t), 1, fp);
  fwrite(arr->data, arr->elemSz*arr->size, 1, fp);
}

void
gkyl_sub_array_write(const struct gkyl_range *range,
  const struct gkyl_array *arr, FILE *fp)
{
  assert(range->volume <= arr->size);
#define _F(loc) gkyl_array_fetch(arr, loc)
  
  uint64_t elemSz = arr->elemSz, size = arr->size;
  fwrite(&elemSz, sizeof(uint64_t), 1, fp);
  fwrite(&size, sizeof(uint64_t), 1, fp);
  
  // construct skip iterator to allow writing (potentially) in chunks
  // rather than element by element or requiring a copy of data
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

int
gkyl_grid_array_write(const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  const struct gkyl_array *arr, const char *fname)
{
  FILE *fp = fopen(fname, "wb"); if (!fp) return errno;
  gkyl_rect_grid_write(grid, fp);
  gkyl_sub_array_write(range, arr, fp);
  fclose(fp);
  return 0;
}
