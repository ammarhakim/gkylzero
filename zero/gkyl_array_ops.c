#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gkyl_array_ops.h>

struct gkyl_array*
gkyl_array_clear_dbl(struct gkyl_array*out, double val)
{
  assert(out->elemSz==sizeof(double));

  double *out_d = out->data;
  for (unsigned i=0; i<out->size; ++i)
    out_d[i] = val;
  return out;
}
struct gkyl_array*
gkyl_array_clear_flt(struct gkyl_array*out, float val)
{
  assert(out->elemSz==sizeof(float));

  float *out_d = out->data;
  for (unsigned i=0; i<out->size; ++i)
    out_d[i] = val;
  return out;
}

struct gkyl_array*
gkyl_array_accumulate_dbl(struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  assert(out->elemSz==sizeof(double) && inp->elemSz==sizeof(double) &&
    out->size == inp->size);

  double *out_d = out->data;
  const double *inp_d = (double*) inp->data;
  for (unsigned i=0; i<out->size; ++i)
    out_d[i] += a*inp_d[i];
  return out;
}
struct gkyl_array*
gkyl_array_accumulate_flt(struct gkyl_array* out, float a,
  const struct gkyl_array* inp)
{
  assert(out->elemSz==sizeof(float) && inp->elemSz==sizeof(float) &&
    out->size == inp->size);

  float *out_d = out->data;
  const float *inp_d = inp->data;
  for (unsigned i=0; i<out->size; ++i)
    out_d[i] += a*inp_d[i];    
  return out;
}

struct gkyl_array*
gkyl_array_set_dbl(struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  assert(out->elemSz==sizeof(double) && inp->elemSz==sizeof(double) &&
    out->size == inp->size);

  double *out_d = out->data;
  const double *inp_d = inp->data;
  for (unsigned i=0; i<out->size; ++i)
    out_d[i] = a*inp_d[i];
  return out;
}
struct gkyl_array* gkyl_array_set_flt(struct gkyl_array* out,
  double a, const struct gkyl_array* inp)
{
  assert(out->elemSz==sizeof(float) && inp->elemSz==sizeof(float) &&
    out->size == inp->size);

  float *out_d = out->data;
  const float *inp_d = inp->data;
  for (unsigned i=0; i<out->size; ++i)
    out_d[i] = a*inp_d[i];
  return out;
}

static double
_gkyl_null_dbl(double x)
{
  return x;
}
static inline double
_gkyl_sq_dbl(double x)
{
  return x*x;
}
static inline double
_gkyl_cube_dbl(double x)
{
  return x*x*x;
}

// List of double-precision unary operators
static struct { char *op; double (*f)(double x); } uniop_funcs_dbl[] = {
  { "cos", cos },
  { "cube", _gkyl_cube_dbl },
  { "sin", sin },
  { "square", _gkyl_sq_dbl },
  { "tan", tan },
};

// Find function in uniop_dbl list, returning function pointer. [This
// very weird looking function signature reads: find_uniop_dbl is a
// function that takes a const char * and returns a pointer to a
// function with signature double (*)(double).]
double (*find_uniop_dbl(const char *op))(double)
{
  size_t nv = sizeof(uniop_funcs_dbl)/sizeof(uniop_funcs_dbl[0]);
  for (unsigned i=0; i<nv; ++i)
    if (strcmp(op, uniop_funcs_dbl[i].op) == 0)
      return uniop_funcs_dbl[i].f;
  return _gkyl_null_dbl;
}

struct gkyl_array*
gkyl_array_uniop_dbl(const char *op, double a,
  struct gkyl_array *out, double b, const struct gkyl_array *inp)
{
  assert(out->elemSz==sizeof(double) && inp->elemSz==sizeof(double) &&
    out->size == inp->size);
  
  double (*f)(double) = find_uniop_dbl(op);
  assert(f != _gkyl_null_dbl); // operator not found

  double *out_d = out->data;
  const double *inp_d = inp->data;
  for (unsigned i=0; i<out->size; ++i)
    out_d[i] = a*out_d[i] + b*f(inp_d[i]);

  return out;
}
