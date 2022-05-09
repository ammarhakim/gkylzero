#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>

/*
  IF THIS FORMAT IF MODIFIED, PLEASE COPY AND THEN CHANGE THE
  DESCRIPTION SO WE HAVE THE OLDER VERSIONS DOCUMENTED HERE. UPDATE
  VERSION BY 1 EACH TIME YOU CHANGE THE FORMAT.

  The format of the gkyl binary output is as follows.

  ## Version 0: Jan 2021. Created by A.H. Note Version 0 has no header
     information

  Data      Type and meaning
  --------------------------
  ndim      uint64_t Dimension of field
  cells     uint64_t[ndim] number of cells in each direction
  lower     float64[ndim] Lower bounds of grid
  upper     float64[ndim] Upper bounds of grid
  esznc     uint64_t Element-size * number of components in field
  size      uint64_t Total number of cells in field
  DATA      size*esznc bytes of data  
  
  ## Version 1: May 9th 2022. Created by A.H

  Data      Type and meaning
  --------------------------
  gkyl0     5 bytes
  version   uint64_t 
  file_type uint64_t (0: field data, 1: diagnostic data)
  meta_size uint64_t Number of bytes of meta-data
  DATA      meta_size bytes of data. This is in msgpack format

  For file_type = 0 (field) the above header is followed by

  real_type uint64_t. Indicates real type of data 
  ndim      uint64_t Dimension of field
  cells     uint64_t[ndim] number of cells in each direction
  lower     float64[ndim] Lower bounds of grid
  upper     float64[ndim] Upper bounds of grid
  esznc     uint64_t Element-size * number of components in field
  size      uint64_t Total number of cells in field
  DATA      size*esznc bytes of data
  
 */

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

// type-id for field data
static const uint64_t field_file_type = 1;

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
    const char g0[5] = "gkyl0";

    // Version 1 header
    fwrite(g0, sizeof(char[5]), 1, fp);
    uint64_t version = 1;
    fwrite(&version, sizeof(uint64_t), 1, fp);
    fwrite(&field_file_type, sizeof(uint64_t), 1, fp);
    uint64_t meta_size = 0; // THIS WILL CHANGE ONCE METADATA IS EMBEDDED
    fwrite(&meta_size, sizeof(uint64_t), 1, fp);

    // Version 0 format is used for rest of the file
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
