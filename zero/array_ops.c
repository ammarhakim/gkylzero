#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gkyl_util.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>

struct gkyl_array*
gkyl_array_clear(struct gkyl_array* out, double val)
{
  assert(out->type == GKYL_DOUBLE);

  double *out_d = out->data;
  for (size_t i=0; i<NELM(out); ++i)
    out_d[i] = val;
  return out;
}

struct gkyl_array*
gkyl_array_accumulate(struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size && out->elemsz == inp->elemsz);

  double *out_d = out->data;
  const double *inp_d = inp->data;
  for (size_t i=0; i<NELM(out); ++i)
    out_d[i] += a*inp_d[i];
  return out;
}

struct gkyl_array*
gkyl_array_set(struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size && out->elemsz == inp->elemsz);

  double *out_d = out->data;
  const double *inp_d = inp->data;
  for (size_t i=0; i<NELM(out); ++i)
    out_d[i] = a*inp_d[i];
  return out;
}

struct gkyl_array*
gkyl_array_scale(struct gkyl_array* out, double a)
{
  return gkyl_array_set(out, a, out);
}

double
gkyl_array_reduce(const struct gkyl_array *arr, enum gkyl_array_op op, double *out)
{
  assert(arr->type == GKYL_DOUBLE);
  double res = 0, *arr_d = arr->data;

  switch (op) {
    case GKYL_MIN:
      res = DBL_MAX;
      for (size_t i=0; i<NELM(arr); ++i)
        res = fmin(res, arr_d[i]);
      break;

    case GKYL_MAX:
      res = -DBL_MAX;
      for (size_t i=0; i<NELM(arr); ++i)
        res = fmax(res, arr_d[i]);
      break;
  }
  if (out) *out = res;
  return res;
}

// range based methods
struct gkyl_array*
gkyl_array_clear_range(struct gkyl_array *out, double val, struct gkyl_range range)
{
  assert(out->type == GKYL_DOUBLE);
  long n = NCOM(out);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&range, iter.idx);
    array_clear1(n, gkyl_array_fetch(out, start), val);
  }

  return out;
}

struct gkyl_array*
gkyl_array_accumulate_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, struct gkyl_range range)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == inp->size);

  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&range, iter.idx);
    array_acc1(n,
      gkyl_array_fetch(out, start), a, gkyl_array_cfetch(inp, start));
  }

  return out;
}

struct gkyl_array*
gkyl_array_set_range(struct gkyl_array *out,
  double a, const struct gkyl_array *inp, struct gkyl_range range)
{
  assert(out->type == GKYL_DOUBLE && inp->type == GKYL_DOUBLE);
  assert(out->size == inp->size);

  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&range, iter.idx);
    array_set1(n,
      gkyl_array_fetch(out, start), a, gkyl_array_cfetch(inp, start));
  }

  return out;
}

struct gkyl_array* gkyl_array_scale_range(struct gkyl_array *out,
  double a, struct gkyl_range range)
{
  return gkyl_array_set_range(out, a, out, range);
}

void gkyl_array_reduce_range(double *res,
  const struct gkyl_array *arr, enum gkyl_array_op op, struct gkyl_range range)
{
  assert(arr->type == GKYL_DOUBLE);

  long n = NCOM(arr);
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  switch (op) {
    case GKYL_MIN:
      for (long i=0; i<n; ++i) res[i] = DBL_MAX;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(&range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        for (long i=0; i<n; ++i)
          res[i] = fmin(res[i], d[i]);
      }
      break;
    case GKYL_MAX:
      for (long i=0; i<n; ++i) res[i] = -DBL_MAX;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(&range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        for (long i=0; i<n; ++i)
          res[i] = fmax(res[i], d[i]);
      }
      break;
  }
}

struct gkyl_array*
gkyl_array_copy_range(struct gkyl_array *out,
  const struct gkyl_array *inp, struct gkyl_range range)
{
  assert(out->size == inp->size && out->elemsz == inp->elemsz);

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&range, iter.idx);
    memcpy(gkyl_array_fetch(out, start), gkyl_array_cfetch(inp, start), inp->esznc);
  }
  return out;
}

void
gkyl_array_copy_to_buffer(void *data, const struct gkyl_array *arr,
  struct gkyl_range range)
{
#define _F(loc) gkyl_array_cfetch(arr, loc)

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&range, iter.idx);
    memcpy(((char*) data) + arr->esznc*count++, _F(start), arr->esznc);
  }

#undef _F
}

void
gkyl_array_copy_from_buffer(struct gkyl_array *arr,
  const void *data, struct gkyl_range range)
{
#define _F(loc) gkyl_array_fetch(arr, loc)

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&range, iter.idx);
    memcpy(_F(start), ((char*) data) + arr->esznc*count++, arr->esznc);
  }

#undef _F
}
