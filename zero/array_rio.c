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

// size in bytes for various data-types
static const size_t array_elem_size[] = {
  [GKYL_INT] = sizeof(int),
  [GKYL_FLOAT] = sizeof(float),
  [GKYL_DOUBLE] = sizeof(double),
  [GKYL_USER] = 1,
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
  
  uint64_t esznc = arr->esznc, size = range->volume;
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
gkyl_grid_sub_array_write(const struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  const struct gkyl_array *arr, const char *fname)
{
  FILE *fp = 0;
  with_file (fp, fname, "w") {
    uint64_t real_type = array_data_type[arr->type];
    fwrite(&real_type, sizeof(uint64_t), 1, fp);
    gkyl_rect_grid_write(grid, fp);
    gkyl_sub_array_write(range, arr, fp);
  }
  return errno;
}

struct gkyl_array*
gkyl_array_new_from_file(enum gkyl_elem_type type, FILE *fp)
{
  struct gkyl_array* arr = 0;
  
  uint64_t esznc, size;
  if (1 != fread(&esznc, sizeof(uint64_t), 1, fp))
    return 0;
  if (1 != fread(&size, sizeof(uint64_t), 1, fp))
    return 0;

  int ncomp = esznc/array_elem_size[type];
  arr = gkyl_array_new(type, ncomp, size);
  if (1 != fread(arr->data, arr->esznc*arr->size, 1, fp)) {
    gkyl_array_release(arr);
    arr = 0;
  }
  return arr;
}

bool
gkyl_sub_array_read(const struct gkyl_range *range, struct gkyl_array *arr, FILE *fp)
{
#define _F(loc) gkyl_array_fetch(arr, loc)
  
  uint64_t esznc, size;
  if (1 != fread(&esznc, sizeof(uint64_t), 1, fp))
    return false;
  if (1 != fread(&size, sizeof(uint64_t), 1, fp))
    return false;

  if ((size != range->volume) || (size > arr->size))
    return false;
  
  // construct skip iterator to allow reading (potentially) in chunks
  // rather than element by element or requiring a copy of data
  struct gkyl_range_skip_iter skip;
  gkyl_range_skip_iter_init(&skip, range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skip.range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&skip.range, iter.idx);
    if (1 != fread(_F(start), arr->esznc*skip.delta, 1, fp))
      return false;
  }

  return true;
#undef _F
}

int
gkyl_grid_sub_array_read(struct gkyl_rect_grid *grid, const struct gkyl_range *range,
  struct gkyl_array *arr, const char* fname)
{
  FILE *fp = 0;
  with_file (fp, fname, "r") {
    uint64_t real_type = 0;
    if (1 != fread(&real_type, sizeof(uint64_t), 1, fp))
      break;
    
    gkyl_rect_grid_read(grid, fp);
    gkyl_sub_array_read(range, arr, fp);
  }
  return errno;  
}

struct gkyl_array*
gkyl_grid_array_new_from_file(struct gkyl_rect_grid *grid, const char* fname)
{
  struct gkyl_array *arr = 0;
  FILE *fp = 0;
  with_file (fp, fname, "r") {
    uint64_t real_type = 0;
    if (1 != fread(&real_type, sizeof(uint64_t), 1, fp))
      break;
    
    gkyl_rect_grid_read(grid, fp);
    arr = gkyl_array_new_from_file(real_type, fp);
  }
  return arr;
}
