#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_util.h>

static void
array_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_array *arr = container_of(ref, struct gkyl_array, ref_count);
  gkyl_free(arr->data);
  gkyl_free(arr);
}

struct gkyl_array*
gkyl_array_new(size_t elemSz, size_t size)
{
  struct gkyl_array* arr = gkyl_malloc(sizeof(struct gkyl_array));
  arr->elemSz = elemSz;
  arr->size = size;
  arr->ref_count = (struct gkyl_ref_count) { array_free, 1 };  
  arr->data = gkyl_calloc(arr->size, elemSz);
  return arr;
}

struct gkyl_array*
gkyl_array_copy(struct gkyl_array* dest, const struct gkyl_array* src)
{
  assert(dest->elemSz == src->elemSz);
  
  long ncopy = src->size < dest->size ? src->size : dest->size;
  memcpy(dest->data, src->data, ncopy*src->elemSz);
  return dest;
}

struct gkyl_array*
gkyl_array_clone(const struct gkyl_array* src)
{
  struct gkyl_array* arr = gkyl_malloc(sizeof(struct gkyl_array));
  arr->elemSz = src->elemSz;
  arr->size = src->size;
  arr->data = gkyl_calloc(arr->size, arr->elemSz);
  memcpy(arr->data, src->data, arr->size*arr->elemSz);
  arr->ref_count = (struct gkyl_ref_count) { array_free, 1 };
  return arr;
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
