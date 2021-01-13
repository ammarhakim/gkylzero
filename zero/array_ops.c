#include <assert.h>
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
      for (unsigned i=0; i<NELM(out,type); ++i)                 \
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
      for (unsigned i=0; i<NELM(out,type); ++i)                         \
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
      for (unsigned i=0; i<NELM(out,type); ++i)                         \
        out_d[i] = a*inp_d[i];                                          \
      return out;                                                       \
    }

GKYL_ARRAY_SET(float)
GKYL_ARRAY_SET(double)

#define GKYL_NULL(type) static type _gkyl_null_##type(type x) { return x; }
GKYL_NULL(float)
GKYL_NULL(double)

#define GKYL_SQ(type) static type _gkyl_sq_##type(type x) { return x*x; }
GKYL_SQ(float)
GKYL_SQ(double)

#define GKYL_CUBE(type) static type _gkyl_cube_##type(type x) { return x*x*x; }
GKYL_CUBE(float)
GKYL_CUBE(double)

// List of double-precision unary operators
static struct { char *op; double (*f)(double x); } uniop_double[] = {
  { "cos", cos },
  { "cube", _gkyl_cube_double },
  { "sin", sin },
  { "square", _gkyl_sq_double },
  { "tan", tan },
};
// List of single-precision unary operators
static struct { char *op; float (*f)(float x); } uniop_float[] = {
  { "cos", cosf },
  { "cube", _gkyl_cube_float },
  { "sin", sinf },
  { "square", _gkyl_sq_float },
  { "tan", tanf },
};

#define FIND_UNIOP(type)                                                \
    static type (*find_uniop_##type(const char *op))(type)              \
    {                                                                   \
      size_t nv = sizeof(uniop_##type)/sizeof(uniop_##type[0]);         \
      for (unsigned i=0; i<nv; ++i)                                     \
        if (strcmp(op, uniop_double[i].op) == 0)                        \
          return uniop_##type[i].f;                                     \
      return _gkyl_null_##type;                                         \
    }

FIND_UNIOP(float)
FIND_UNIOP(double)

#define GKYL_ARRAY_UNIOP(type)                                          \
    struct gkyl_array*                                                  \
    gkyl_array_uniop_##type(const char *op, type a,                     \
      struct gkyl_array *out, type b, const struct gkyl_array *inp)     \
    {                                                                   \
      assert(out->size == inp->size && out->elemsz == inp->elemsz);     \
                                                                        \
      type (*f)(type) = find_uniop_##type(op);                          \
      assert(f != _gkyl_null_##type);                                   \
                                                                        \
      type *out_d = out->data;                                          \
      const type *inp_d = inp->data;                                    \
      for (unsigned i=0; i<NELM(out,type); ++i)                         \
        out_d[i] = a*out_d[i] + b*f(inp_d[i]);                          \
                                                                        \
      return out;                                                       \
    }

GKYL_ARRAY_UNIOP(float)
GKYL_ARRAY_UNIOP(double)

#define GKYL_ARRAY_ACCUMULATE_RANGE(type)                               \
    static void                                                         \
    array_acc1_##type(long n, long del, long outnc, long inpnc,         \
      type *restrict out, type a, type *restrict const inp)             \
    {                                                                   \
    for (long i=0; i<del; ++i)                                          \
      for (int c=0; c<n; ++c)                                           \
        out[i*outnc+c] += a*inp[i*inpnc+c];                             \
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
