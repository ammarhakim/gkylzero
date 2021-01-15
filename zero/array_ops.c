#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gkyl_array_ops.h>

// Compute number of elements of 'type' stored in array 'arr'
#define NELM(arr, type) ((arr->size*arr->elemsz)/sizeof(type))
// Compute number of components of 'type' stored in array 'arr'
#define NCOM(arr, type) ((arr->elemsz)/sizeof(type))

#define GKYL_ARRAY_CLEAR(type)                                  \
    struct gkyl_array*                                          \
    gkyl_array_clear_##type(struct gkyl_array* out, type val)   \
    {                                                           \
      type *out_d = out->data;                                  \
      for (size_t i=0; i<NELM(out,type); ++i)                   \
        out_d[i] = val;                                         \
      return out;                                               \
    }

GKYL_ARRAY_CLEAR(float)
GKYL_ARRAY_CLEAR(double)

#define GKYL_ARRAY_ACCUMULATE(type)                                     \
    struct gkyl_array*                                                  \
    gkyl_array_accumulate_##type(struct gkyl_array* out, type a,        \
      const struct gkyl_array* inp)                                     \
    {                                                                   \
      assert(out->size == inp->size && out->elemsz == inp->elemsz);     \
                                                                        \
      type *out_d = out->data;                                          \
      const type *inp_d = inp->data;                                    \
      for (size_t i=0; i<NELM(out,type); ++i)                           \
        out_d[i] += a*inp_d[i];                                         \
      return out;                                                       \
    }

GKYL_ARRAY_ACCUMULATE(float)
GKYL_ARRAY_ACCUMULATE(double)

#define GKYL_ARRAY_SET(type)                                            \
    struct gkyl_array*                                                  \
    gkyl_array_set_##type(struct gkyl_array* out, type a,               \
      const struct gkyl_array* inp)                                     \
    {                                                                   \
      assert(out->size == inp->size && out->elemsz == inp->elemsz);     \
                                                                        \
      type *out_d = out->data;                                          \
      const type *inp_d = inp->data;                                    \
      for (size_t i=0; i<NELM(out,type); ++i)                           \
        out_d[i] = a*inp_d[i];                                          \
      return out;                                                       \
    }

GKYL_ARRAY_SET(float)
GKYL_ARRAY_SET(double)

#define GKYL_ARRAY_ACCUMULATE_RANGE(type)                               \
    static void                                                         \
    array_acc1_##type(long n, long del, long outnc, long inpnc,         \
      type *restrict out, type a, type *restrict const inp)             \
    {                                                                   \
      for (long i=0; i<del; ++i)                                        \
        for (int c=0; c<n; ++c)                                         \
          out[i*outnc+c] += a*inp[i*inpnc+c];                           \
    }                                                                   \
                                                                        \
    struct gkyl_array*                                                  \
    gkyl_array_accumulate_range_##type(struct gkyl_array *out,          \
      type a, const struct gkyl_array *inp, const struct gkyl_range *range) \
    {                                                                   \
      assert(out->size == inp->size);                                   \
                                                                        \
      long outnc = NCOM(out, type), inpnc = NCOM(inp, type);            \
      long n = outnc<inpnc ? outnc : inpnc;                             \
                                                                        \
      struct gkyl_range_skip_iter skip;                                 \
      gkyl_range_skip_iter_init(&skip, range);                          \
                                                                        \
      struct gkyl_range_iter iter;                                      \
      gkyl_range_iter_init(&iter, &skip.range);                         \
                                                                        \
      long count = 0;                                                   \
      while (gkyl_range_iter_next(&iter)) {                             \
        long start = gkyl_range_idx(&skip.range, iter.idx);             \
        array_acc1_##type(n, skip.delta, outnc, inpnc,                  \
          gkyl_array_fetch(out, start), a, gkyl_array_fetch(inp, start)); \
      }                                                                 \
                                                                        \
      return out;                                                       \
    }

GKYL_ARRAY_ACCUMULATE_RANGE(float)
GKYL_ARRAY_ACCUMULATE_RANGE(double)

double
gkyl_array_reduce(struct gkyl_array *arr, enum gkyl_array_op op)
{
  double res, *arr_d = arr->data;

  switch (op) {
    case GKYL_MIN:
        res = DBL_MAX;
        for (size_t i=0; i<NELM(arr, double); ++i)
          res = fmin(res, arr_d[i]);
        break;
        
    case GKYL_MAX:
        res = DBL_MIN;
        for (size_t i=0; i<NELM(arr, double); ++i)
          res = fmax(res, arr_d[i]);
        break;
  }
  return res;
}

void
gkyl_array_copy_to_buffer(void *data, const struct gkyl_array *arr,
  const struct gkyl_range *range)
{
#define _F(loc) gkyl_array_fetch(arr, loc)  

  // construct skip iterator to allow copying (potentially) in chunks
  // rather than element by element
  struct gkyl_range_skip_iter skip;
  gkyl_range_skip_iter_init(&skip, range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skip.range);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&skip.range, iter.idx);
    memcpy(((char*) data) + arr->elemsz*skip.delta*count++, _F(start), arr->elemsz*skip.delta);
  }
  
#undef _F
}

void
gkyl_array_copy_from_buffer(struct gkyl_array *arr,
  const void *data, const struct gkyl_range *range)
{
#define _F(loc) gkyl_array_fetch(arr, loc)  

  // construct skip iterator to allow copying (potentially) in chunks
  // rather than element by element
  struct gkyl_range_skip_iter skip;
  gkyl_range_skip_iter_init(&skip, range);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &skip.range);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&skip.range, iter.idx);
    memcpy(_F(start), ((char*) data) + arr->elemsz*skip.delta*count++, arr->elemsz*skip.delta);
  }
  
#undef _F
}
