#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <gkyl_array_ops.h>

struct gkyl_array*
gkyl_array_clear_dbl(struct gkyl_array *out, double val)
{
  assert(out->elemSz==sizeof(double));

  double *out_d = out->data;
  for (unsigned i=0; i<out->size; ++i)
    out_d[i] = val;
  
  return out;
}
struct gkyl_array*
gkyl_array_clear_flt(struct gkyl_array *out, float val)
{
  assert(out->elemSz==sizeof(float));

  float *out_d = out->data;
  for (unsigned i=0; i<out->size; ++i)
    out_d[i] = val;
  
  return out;
}

struct gkyl_array*
gkyl_array_accumulate_dbl(struct gkyl_array * restrict out, double a,
  const struct gkyl_array * restrict inp)
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
gkyl_array_accumulate_flt(struct gkyl_array * restrict out, float a,
  const struct gkyl_array * restrict inp)
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
gkyl_array_set_dbl(struct gkyl_array * restrict out, double a,
  const struct gkyl_array * restrict inp)
{
  assert(out->elemSz==sizeof(double) && inp->elemSz==sizeof(double) &&
    out->size == inp->size);

  double *out_d = out->data;
  const double *inp_d = inp->data;
  for (unsigned i=0; i<out->size; ++i)
    out_d[i] = a*inp_d[i];
  
  return out;  
}
struct gkyl_array* gkyl_array_set_flt(struct gkyl_array *out,
  double a, const struct gkyl_array *inp)
{
  assert(out->elemSz==sizeof(float) && inp->elemSz==sizeof(float) &&
    out->size == inp->size);

  float *out_d = out->data;
  const float *inp_d = inp->data;
  for (unsigned i=0; i<out->size; ++i)
    out_d[i] = a*inp_d[i];
  
  return out;
}
