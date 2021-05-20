#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>

#include <gkyl_array_rio.h>

// code for array datatype for use in IO
static const uint64_t array_data_type[] = {
  [GKYL_INT] = 0,
  [GKYL_FLOAT] = 1,
  [GKYL_DOUBLE] = 2,
  [GKYL_USER] = 32,
};

void
gkyl_array_write(const struct gkyl_array *arr, FILE *fp)
{
  uint64_t esznc = arr->esznc, size = arr->size;
  fwrite(&esznc, sizeof(uint64_t), 1, fp);
  fwrite(&size, sizeof(uint64_t), 1, fp);
  fwrite(arr->data, arr->esznc*arr->size, 1, fp);
}

void
gkyl_sub_array_write(const struct gkyl_range *range,
  const struct gkyl_array *arr, FILE *fp)
{
#define _F(loc) gkyl_array_cfetch(arr, loc)
  
  uint64_t esznc = arr->esznc, size = arr->size;
  fwrite(&esznc, sizeof(uint64_t), 1, fp);
  fwrite(&size, sizeof(uint64_t), 1, fp);
  
  // construct skip iterator to allow writing (potentially) in chunks
  // rather than element by element or requiring a copy of data
  struct gkyl_range_skip_iter skip;
  gkyl_range_skip_iter_init(&skip, range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skip.range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&skip.range, iter.idx);
    fwrite(_F(start), arr->esznc*skip.delta, 1, fp);
  }
#undef _F
}

int
gkyl_grid_array_write(const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  const struct gkyl_array *arr, const char *fname)
{
  FILE *fp = fopen(fname, "wb"); if (!fp) return errno;

  // data double
  uint64_t real_type = array_data_type[arr->type];
  fwrite(&real_type, sizeof(uint64_t), 1, fp);
  
  gkyl_rect_grid_write(grid, fp);
  gkyl_sub_array_write(range, arr, fp);
  fclose(fp);
  return 0;
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
  
  fprintf(fp, " }\n ");
}
