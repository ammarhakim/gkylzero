#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include <gkyl_alloc.h>
#include <gkyl_array.h>
#include <gkyl_util.h>

// flags and corresponding bit-masks
enum array_flags { A_IS_CU_ARRAY };
static const uint32_t masks[] =
{ 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

// NV-GPU flags
#define SET_CU_ARRAY(flags) (flags) |= masks[A_IS_CU_ARRAY]
#define CLEAR_CU_ARRAY(flags) (flags) &= ~masks[A_IS_CU_ARRAY]
#define IS_CU_ARRAY(flags) (((flags) & masks[A_IS_CU_ARRAY]) != 0)

// size in bytes for various data-types
static const size_t array_elem_size[] = {
  [GKYL_INT] = sizeof(int),
  [GKYL_FLOAT] = sizeof(float),
  [GKYL_DOUBLE] = sizeof(double),
  [GKYL_USER] = 1,
};

static void
array_free(const struct gkyl_ref_count *ref)
{
  struct gkyl_array *arr = container_of(ref, struct gkyl_array, ref_count);
  gkyl_free(arr->data);
  gkyl_free(arr);
}

struct gkyl_array*
gkyl_array_new(enum gkyl_elem_type type, size_t ncomp, size_t size)
{
  struct gkyl_array* arr = gkyl_malloc(sizeof(struct gkyl_array));

  arr->type = type;
  arr->elemsz = array_elem_size[type];
  arr->ncomp = ncomp;
  arr->size = size;
  arr->flags = 0;
  arr->esznc = arr->elemsz*arr->ncomp;
  arr->ref_count = (struct gkyl_ref_count) { array_free, 1 };  
  arr->data = gkyl_calloc(arr->size, arr->esznc);
  
  return arr;
}

bool
gkyl_array_is_cu_dev(struct gkyl_array *const arr)
{
  return IS_CU_ARRAY(arr->flags);
}

struct gkyl_array*
gkyl_array_copy(struct gkyl_array* dest, const struct gkyl_array* src)
{
  assert(dest->esznc == src->esznc);
  
  long ncopy = src->size < dest->size ? src->size : dest->size;
  memcpy(dest->data, src->data, ncopy*src->esznc);
  return dest;
}

struct gkyl_array*
gkyl_array_clone(const struct gkyl_array* src)
{
  struct gkyl_array* arr = gkyl_malloc(sizeof(struct gkyl_array));
  
  arr->elemsz = src->elemsz;
  arr->ncomp = src->ncomp;
  arr->esznc = src->esznc;
  arr->size = src->size;
  arr->data = gkyl_calloc(arr->size, arr->esznc);
  memcpy(arr->data, src->data, arr->size*arr->esznc);
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
