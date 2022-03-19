#include "gkyl_range.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gkyl_util.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_array_reduce.h>

struct gkyl_array*
gkyl_array_clear(struct gkyl_array* out, double val)
{
  assert(out->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) {gkyl_array_clear_cu(out, val); return out; }
#endif

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

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out) && gkyl_array_is_cu_dev(inp)) { gkyl_array_accumulate_cu(out, a, inp); return out; }
#endif

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

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_set_cu(out, a, inp); return out; }
#endif

  double *out_d = out->data;
  const double *inp_d = inp->data;
  for (size_t i=0; i<NELM(out); ++i)
    out_d[i] = a*inp_d[i];
  return out;
}

struct gkyl_array*
gkyl_array_scale(struct gkyl_array* out, double a)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_scale_cu(out, a); return out; }
#endif

  return gkyl_array_set(out, a, out);
}

struct gkyl_array*
gkyl_array_scale_by_cell(struct gkyl_array* out, const struct gkyl_array* a)
{
  assert(out->type == GKYL_DOUBLE);
  assert(out->size == a->size && a->ncomp == 1);
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_scale_by_cell_cu(out, a); return out; }
#endif

  double *out_d = out->data;
  const double *a_d = a->data;
  for (size_t i=0; i<out->size; ++i)
    for (size_t c=0; c<out->ncomp; ++c)
      out_d[i*out->ncomp+c] = a_d[i]*out_d[i*out->ncomp+c];
  return out;
}

void 
gkyl_array_reduce(double *out, const struct gkyl_array *arr, enum gkyl_array_op op)
{
  assert(arr->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    switch (op) {
      case GKYL_MAX:
        gkyl_array_reduce_max_cu(out, arr);
        break;
      case GKYL_MIN:
        gkyl_array_reduce_min_cu(out, arr);
        break;
      case GKYL_SUM:
        gkyl_array_reduce_sum_cu(out, arr);
        break;
    }
    return;
  }
#endif

  long nc = NCOM(arr);
  double *arr_d = arr->data;

  switch (op) {
    case GKYL_MIN:
      for (long k=0; k<nc; ++k) out[k] = DBL_MAX;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        for (long k=0; k<nc; ++k)
          out[k] = fmin(out[k], d[k]);
      }
      break;

    case GKYL_MAX:
      for (long k=0; k<nc; ++k) out[k] = -DBL_MAX;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        for (long k=0; k<nc; ++k)
          out[k] = fmax(out[k], d[k]);
      }
      break;

    case GKYL_SUM:
      for (long k=0; k<nc; ++k) out[k] = 0;
      for (size_t i=0; i<arr->size; ++i) {
        const double *d = gkyl_array_cfetch(arr, i);
        for (long k=0; k<nc; ++k)
          out[k] += d[k];
      }
      break;
  }
}

// range based methods
struct gkyl_array*
gkyl_array_clear_range(struct gkyl_array *out, double val, struct gkyl_range range)
{
  assert(out->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_clear_range_cu(out, val, range); return out; }
#endif

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

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_accumulate_range_cu(out, a, inp, range); return out; }
#endif

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

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_set_range_cu(out, a, inp, range); return out; }
#endif

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

struct gkyl_array*
gkyl_array_scale_range(struct gkyl_array *out,
  double a, struct gkyl_range range)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_scale_range_cu(out, a, range); return out; }
#endif

  return gkyl_array_set_range(out, a, out, range);
}

void
gkyl_array_reduce_range(double *res,
  const struct gkyl_array *arr, enum gkyl_array_op op, struct gkyl_range range)
{
  assert(arr->type == GKYL_DOUBLE);

#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    switch (op) {
      case GKYL_MAX:
        gkyl_array_reduce_range_max_cu(res, arr, range);
        break;
      case GKYL_MIN:
        gkyl_array_reduce_range_min_cu(res, arr, range);
        break;
      case GKYL_SUM:
        gkyl_array_reduce_range_sum_cu(res, arr, range);
        break;
    }
    return;
  }
#endif

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
    case GKYL_SUM:
      for (long i=0; i<n; ++i) res[i] = 0;

      while (gkyl_range_iter_next(&iter)) {
        long start = gkyl_range_idx(&range, iter.idx);
        const double *d = gkyl_array_cfetch(arr, start);
        for (long i=0; i<n; ++i)
          res[i] += d[i];
      }
      break;
  }
}

struct gkyl_array*
gkyl_array_copy_range(struct gkyl_array *out,
  const struct gkyl_array *inp, struct gkyl_range range)
{
  assert(out->size == inp->size && out->elemsz == inp->elemsz);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_copy_range_cu(out, inp, range); return out; }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  while (gkyl_range_iter_next(&iter)) {
    long start = gkyl_range_idx(&range, iter.idx);
    memcpy(gkyl_array_fetch(out, start), gkyl_array_cfetch(inp, start), inp->esznc);
  }
  return out;
}

struct gkyl_array*
gkyl_array_copy_range_to_range(struct gkyl_array *out,
  const struct gkyl_array *inp, struct gkyl_range out_range, struct gkyl_range inp_range)
{
  assert(out->size == inp->size && out->elemsz == inp->elemsz);
  assert(out_range.volume == inp_range.volume);

#ifdef GKYL_HAVE_CUDA
  assert(gkyl_array_is_cu_dev(out)==gkyl_array_is_cu_dev(inp));
  if (gkyl_array_is_cu_dev(out)) { gkyl_array_copy_range_to_range_cu(out, inp, out_range, inp_range); return out; }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &out_range);

  while (gkyl_range_iter_next(&iter)) {
    long start_out = gkyl_range_idx(&out_range, iter.idx);
    long start_inp = gkyl_range_idx(&inp_range, iter.idx);
    memcpy(gkyl_array_fetch(out, start_out), gkyl_array_cfetch(inp, start_inp), inp->esznc);
  }
  return out;
}

void
gkyl_array_copy_to_buffer(void *data, const struct gkyl_array *arr,
  struct gkyl_range range)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) { gkyl_array_copy_to_buffer_cu(data, arr, range); return; }
#endif

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
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) { gkyl_array_copy_from_buffer_cu(arr, data, range); return; }
#endif

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

void
gkyl_array_copy_to_buffer_fn(void *data, const struct gkyl_array *arr,
  struct gkyl_range range, struct gkyl_array_copy_func *cf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) { gkyl_array_copy_to_buffer_fn_cu(data, arr, range, cf); return; }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  long count = 0;
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);

    const double *inp = gkyl_array_cfetch(arr, loc);
    double *out = flat_fetch(data, arr->esznc*count);
    cf->func(arr->ncomp, out, inp, cf->ctx);
    count += 1;
  }
}

void
gkyl_array_flip_copy_to_buffer_fn(void *data, const struct gkyl_array *arr,
  int dir, struct gkyl_range range, struct gkyl_array_copy_func *cf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    if (gkyl_array_is_cu_dev(arr)) { gkyl_array_flip_copy_to_buffer_fn_cu(data, arr, dir, range, cf); return; }
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &range);

  int fidx[GKYL_MAX_DIM]; // flipped index
  struct gkyl_range buff_range;
  gkyl_range_init(&buff_range, range.ndim, range.lower, range.upper);

  int uplo = range.upper[dir]+range.lower[dir];

  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&range, iter.idx);

    gkyl_copy_int_arr(range.ndim, iter.idx, fidx);
    fidx[dir] = uplo - iter.idx[dir];
    
    long count = gkyl_range_idx(&buff_range, fidx);

    const double *inp = gkyl_array_cfetch(arr, loc);
    double *out = flat_fetch(data, arr->esznc*count);
    cf->func(arr->ncomp, out, inp, cf->ctx);
  }  
}
