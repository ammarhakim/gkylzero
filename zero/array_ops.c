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

#define GKYL_ARRAY_REDUCE(type)                                         \
    type                                                                \
    gkyl_array_reduce_##type(const struct gkyl_array *arr, enum gkyl_array_op op, type *out) \
    {                                                                   \
      type res, *arr_d = arr->data;                                     \
                                                                        \
      switch (op) {                                                     \
        case GKYL_MIN:                                                  \
            res = DBL_MAX;                                              \
            for (size_t i=0; i<NELM(arr, type); ++i)                    \
              res = fmin(res, arr_d[i]);                                \
            break;                                                      \
                                                                        \
        case GKYL_MAX:                                                  \
            res = -DBL_MAX;                                             \
            for (size_t i=0; i<NELM(arr, type); ++i)                    \
              res = fmax(res, arr_d[i]);                                \
            break;                                                      \
      }                                                                 \
      if (out) *out = res;                                              \
      return res;                                                       \
    }                                                                  

GKYL_ARRAY_REDUCE(float)
GKYL_ARRAY_REDUCE(double)

// range based methods

#define GKYL_ARRAY_CLEAR_RANGE(type)                                    \
    static inline void                                                  \
    array_clear1_##type(long n, type *out, type val)                    \
    {                                                                   \
      for (int c=0; c<n; ++c)                                           \
        out[c] = val;                                                   \
    }                                                                   \
                                                                        \
    struct gkyl_array*                                                  \
    gkyl_array_clear_range_##type(struct gkyl_array *out,               \
      type val, const struct gkyl_range *range)                         \
    {                                                                   \
      long n = NCOM(out, type);                                         \
                                                                        \
      struct gkyl_range_iter iter;                                      \
      gkyl_range_iter_init(&iter, range);                               \
                                                                        \
      long count = 0;                                                   \
      while (gkyl_range_iter_next(&iter)) {                             \
        long start = gkyl_range_idx(range, iter.idx);                   \
        array_clear1_##type(n, gkyl_array_fetch(out, start), val);      \
      }                                                                 \
                                                                        \
      return out;                                                       \
    }

GKYL_ARRAY_CLEAR_RANGE(float)
GKYL_ARRAY_CLEAR_RANGE(double)

#define GKYL_ARRAY_ACCUMULATE_RANGE(type)                               \
    static inline void                                                  \
    array_acc1_##type(long n,                                           \
      type *restrict out, type a, const type *restrict inp)             \
    {                                                                   \
      for (int c=0; c<n; ++c)                                           \
        out[c] += a*inp[c];                                             \
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
      struct gkyl_range_iter iter;                                      \
      gkyl_range_iter_init(&iter, range);                               \
                                                                        \
      while (gkyl_range_iter_next(&iter)) {                             \
        long start = gkyl_range_idx(range, iter.idx);                   \
        array_acc1_##type(n,                                            \
          gkyl_array_fetch(out, start), a, gkyl_array_cfetch(inp, start)); \
      }                                                                 \
                                                                        \
      return out;                                                       \
    }

GKYL_ARRAY_ACCUMULATE_RANGE(float)
GKYL_ARRAY_ACCUMULATE_RANGE(double)

#define GKYL_ARRAY_SET_RANGE(type)                                      \
    static inline void                                                  \
    array_set1_##type(long n,                                           \
      type *restrict out, type a, const type *restrict inp)             \
    {                                                                   \
      for (int c=0; c<n; ++c)                                           \
        out[c] = a*inp[c];                                              \
    }                                                                   \
                                                                        \
    struct gkyl_array*                                                  \
    gkyl_array_set_range_##type(struct gkyl_array *out,                 \
      type a, const struct gkyl_array *inp, const struct gkyl_range *range) \
    {                                                                   \
      assert(out->size == inp->size);                                   \
                                                                        \
      long outnc = NCOM(out, type), inpnc = NCOM(inp, type);            \
      long n = outnc<inpnc ? outnc : inpnc;                             \
                                                                        \
      struct gkyl_range_iter iter;                                      \
      gkyl_range_iter_init(&iter, range);                               \
                                                                        \
      while (gkyl_range_iter_next(&iter)) {                             \
        long start = gkyl_range_idx(range, iter.idx);                   \
        array_set1_##type(n,                                            \
          gkyl_array_fetch(out, start), a, gkyl_array_cfetch(inp, start)); \
      }                                                                 \
                                                                        \
      return out;                                                       \
    }

GKYL_ARRAY_SET_RANGE(float)
GKYL_ARRAY_SET_RANGE(double)

#define GKYL_ARRAY_REDUCE_RANGE(type)                                   \
    void gkyl_array_reduce_range_##type(type *res,                      \
      const struct gkyl_array *arr, enum gkyl_array_op op, const struct gkyl_range *range) \
    {                                                                   \
      long n = NCOM(arr, type);                                         \
      struct gkyl_range_iter iter;                                      \
      gkyl_range_iter_init(&iter, range);                               \
                                                                        \
      switch (op) {                                                     \
        case GKYL_MIN:                                                  \
            for (long i=0; i<n; ++i) res[i] = DBL_MAX;                  \
                                                                        \
            while (gkyl_range_iter_next(&iter)) {                       \
              long start = gkyl_range_idx(range, iter.idx);             \
              const type *d = gkyl_array_cfetch(arr, start);            \
              for (long i=0; i<n; ++i)                                  \
                res[i] = fmin(res[i], d[i]);                            \
            }                                                           \
            break;                                                      \
        case GKYL_MAX:                                                  \
            for (long i=0; i<n; ++i) res[i] = -DBL_MAX;                 \
                                                                        \
            while (gkyl_range_iter_next(&iter)) {                       \
              long start = gkyl_range_idx(range, iter.idx);             \
              const type *d = gkyl_array_cfetch(arr, start);            \
              for (long i=0; i<n; ++i)                                  \
                res[i] = fmax(res[i], d[i]);                            \
            }                                                           \
            break;                                                      \
      }                                                                 \
    }                                                                  

GKYL_ARRAY_REDUCE_RANGE(double)
GKYL_ARRAY_REDUCE_RANGE(float)    

void
gkyl_array_copy_to_buffer(void *data, const struct gkyl_array *arr,
  const struct gkyl_range *range)
{
#define _F(loc) gkyl_array_cfetch(arr, loc)

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    memcpy(((char*) data) + arr->elemsz*count++, _F(start), arr->elemsz);
  }
  
#undef _F
}

void
gkyl_array_copy_from_buffer(struct gkyl_array *arr,
  const void *data, const struct gkyl_range *range)
{
#define _F(loc) gkyl_array_fetch(arr, loc)  

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, range);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(range, iter.idx);
    memcpy(_F(start), ((char*) data) + arr->elemsz*count++, arr->elemsz);
  }
  
#undef _F
}
