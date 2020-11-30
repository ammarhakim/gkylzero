// Gkeyll includes
#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_util.h>

// std includes
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static void
array_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_array *arr = container_of(ref, struct gkyl_array, ref_count);
  free(arr->data);
  free(arr);
}

struct gkyl_array*
gkyl_array_new(int rank, size_t elemSz, const size_t *shape)
{
  struct gkyl_array* arr = gkyl_malloc(sizeof(struct gkyl_array));
  arr->rank = rank; arr->elemSz = elemSz;
  arr->size = 1;
  for (unsigned i=0; i<rank; ++i) {
    arr->size *= shape[i];
    arr->shape[i] = shape[i];
  }
  arr->ref_count = (struct gkyl_ref_count) { array_free, 1 };  
  arr->data = gkyl_calloc(arr->size, elemSz);
  return arr;
}

struct gkyl_array*
gkyl_array_reshape(struct gkyl_array* arr, int rank, const size_t *shape)
{
  arr->rank = rank;

  size_t size = 1;
  for (unsigned i=0; i<rank; ++i) {
    size *= shape[i];
    arr->shape[i] = shape[i];
  }
  size_t origSize = arr->size;
  if (size != origSize) {
    arr->size = size;
    arr->data = realloc(arr->data, size*arr->elemSz);
  }
  
  if (origSize > 0 && size > origSize) {
    // replicate data to fill newly allocated memory
    size_t i, elemSz = arr->elemSz;
    for (i=0; i<size/origSize-1; ++i)
      memcpy((char*) arr->data+(i+1)*origSize*elemSz, arr->data, origSize*elemSz);
    if ((size % origSize) > 0)
      memcpy((char*) arr->data+(i+1)*origSize*elemSz, arr->data, (size % origSize)*elemSz);
  }
  return arr;
}

struct gkyl_array*
gkyl_array_copy(struct gkyl_array* dest, const struct gkyl_array* src)
{
  size_t ncopy = src->size < dest->size ? src->size : dest->size;
  memcpy(dest->data, src->data, ncopy*src->elemSz);
  return dest;
}

struct gkyl_array*
gkyl_array_clone(const struct gkyl_array* src)
{
  struct gkyl_array* arr = gkyl_malloc(sizeof(struct gkyl_array));
  arr->rank = src->rank; arr->elemSz = src->elemSz;
  arr->size = src->size;
  for (unsigned i=0; i<arr->rank; ++i)
    arr->shape[i] = src->shape[i];
  
  arr->data = gkyl_calloc(arr->size, arr->elemSz);
  memcpy(arr->data, src->data, arr->size*arr->elemSz);
  arr->ref_count = (struct gkyl_ref_count) { array_free, 1 };

  return arr;
}

void
gkyl_array_write(const struct gkyl_array *arr, FILE *fp)
{
  uint64_t rank = arr->rank;
  fwrite(&rank, sizeof(uint64_t), 1, fp); // write as 64-bit integer
  fwrite(arr->shape, sizeof(size_t), rank, fp);
  fwrite(arr->data, arr->elemSz*arr->size, 1, fp);
}

struct gkyl_array*
gkyl_array_aquire(const struct gkyl_array* arr)
{
  gkyl_ref_count_inc(&arr->ref_count);
  return (struct gkyl_array*) arr;
}

void
gkyl_array_release(const struct gkyl_array* arr)
{
  gkyl_ref_count_dec(&arr->ref_count);
}
