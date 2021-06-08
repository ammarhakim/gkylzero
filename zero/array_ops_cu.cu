#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_util.h>
}

__global__ void
gkyl_array_clear_cu_kernel(struct gkyl_array* out, double val)
{

  double *out_d = (double*) out->data;
  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x; 
      linc < NELM(out);
      linc += blockDim.x*gridDim.x)
  {
    out_d[linc] = val;
  }
}

__global__ void
gkyl_array_accumulate_cu_kernel(struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  double *out_d = (double*) out->data;
  const double *inp_d = (const double*) inp->data;
  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x; 
      linc < NELM(out);
      linc += blockDim.x*gridDim.x)
  {
    out_d[linc] += a*inp_d[linc];
  }
}

__global__ void
gkyl_array_set_cu_kernel(struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  double *out_d = (double*) out->data;
  const double *inp_d = (const double*) inp->data;
  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x; 
      linc < NELM(out);
      linc += blockDim.x*gridDim.x)
  {
    out_d[linc] = a*inp_d[linc];
  }
}

// Host-side wrappers for array operations
void
gkyl_array_clear_cu(struct gkyl_array* out, double val)
{
  gkyl_array_clear_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_device, val);
}

void
gkyl_array_accumulate_cu(struct gkyl_array* out, double a, const struct gkyl_array* inp)
{
  gkyl_array_accumulate_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_device, a, inp->on_device);
}

void
gkyl_array_set_cu(struct gkyl_array* out, double a, const struct gkyl_array* inp)
{
  gkyl_array_set_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_device, a, inp->on_device);
}

void
gkyl_array_scale_cu(struct gkyl_array* out, double a)
{
  gkyl_array_set_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_device, a, out->on_device);
}

// Range-based methods
// Range-based methods need to inverse index from linc to idx.
// Must use gkyl_sub_range_inv_idx so that linc=0 maps to idxc={0,0,...}
// since range can be a subrange.
// Then, convert back to a linear index on the super-range.
// This super range can include ghost cells and thus linear index will have
// have jumps over ghost cells.

__global__ void
gkyl_array_clear_range_cu_kernel(struct gkyl_array *out, double val, const struct gkyl_range range)
{
  long n = NCOM(out);
  int idx[GKYL_MAX_DIM];
  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x; 
      linc < range.volume;
      linc += blockDim.x*gridDim.x)
  {
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    array_clear1(n, (double*) gkyl_array_fetch(out, start), val);
  }
}

__global__ void
gkyl_array_accumulate_range_cu_kernel(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, const struct gkyl_range range)
{
  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;
  int idx[GKYL_MAX_DIM];
  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x; 
      linc < range.volume;
      linc += blockDim.x*gridDim.x)
  {
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    array_acc1(n,
      (double*) gkyl_array_fetch(out, start), a, (double*) gkyl_array_cfetch(inp, start));
  }
}

__global__ void
gkyl_array_set_range_cu_kernel(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, const struct gkyl_range range)
{
  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;
  int idx[GKYL_MAX_DIM];
  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x; 
      linc < range.volume;
      linc += blockDim.x*gridDim.x)
  {
    gkyl_sub_range_inv_idx(&range, linc, idx); 
    long start = gkyl_range_idx(&range, idx);
    array_set1(n,
      (double*) gkyl_array_fetch(out, start), a, (double*) gkyl_array_cfetch(inp, start));
  }
}

__global__ void 
gkyl_array_copy_range_cu_kernel(struct gkyl_array *out,
  const struct gkyl_array* inp, const struct gkyl_range range)
{
  int idx[GKYL_MAX_DIM];
  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x; 
      linc < range.volume;
      linc += blockDim.x*gridDim.x)
  {
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    memcpy((double*) gkyl_array_fetch(out, start), (double*) gkyl_array_cfetch(inp, start), inp->esznc);
  }
}

__global__ void 
gkyl_array_copy_to_buffer_cu_kernel(void *data, const struct gkyl_array *arr,
  const struct gkyl_range range)
{
  int idx[GKYL_MAX_DIM];
  long count = 0;
  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x; 
      linc < range.volume;
      linc += blockDim.x*gridDim.x)
  {
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    memcpy(((char*) data) + arr->esznc*count++, (const double*) gkyl_array_cfetch(arr, start), arr->esznc);
  }
}

__global__ void 
gkyl_array_copy_from_buffer_cu_kernel(struct gkyl_array *arr, const void *data,
  const struct gkyl_range range)
{
  int idx[GKYL_MAX_DIM];
  long count = 0;
  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x; 
      linc < range.volume;
      linc += blockDim.x*gridDim.x)
  {
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    memcpy((double*) gkyl_array_fetch(arr, start), ((char*) data) + arr->esznc*count++, arr->esznc);
  }
}

// Host-side wrappers for range-based array operations
void
gkyl_array_clear_range_cu(struct gkyl_array *out, double val, const struct gkyl_range range)
{
  gkyl_array_clear_range_cu_kernel<<<range.nblocks, range.nthreads>>>(out->on_device, val, range);
}

void
gkyl_array_accumulate_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, const struct gkyl_range range)
{
  gkyl_array_accumulate_range_cu_kernel<<<range.nblocks, range.nthreads>>>(out->on_device, a, inp->on_device, range);
}

void
gkyl_array_set_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, const struct gkyl_range range)
{
  gkyl_array_set_range_cu_kernel<<<range.nblocks, range.nthreads>>>(out->on_device, a, inp->on_device, range);
}

void
gkyl_array_scale_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_range range)
{
  gkyl_array_set_range_cu_kernel<<<range.nblocks, range.nthreads>>>(out->on_device, a, out->on_device, range);
}

void
gkyl_array_copy_range_cu(struct gkyl_array *out,
  const struct gkyl_array* inp, const struct gkyl_range range)
{
  gkyl_array_copy_range_cu_kernel<<<range.nblocks, range.nthreads>>>(out->on_device, inp->on_device, range);
}

void 
gkyl_array_copy_to_buffer_cu(void *data, 
  const struct gkyl_array *arr, const struct gkyl_range range)
{
  gkyl_array_copy_to_buffer_cu_kernel<<<range.nblocks, range.nthreads>>>(data, arr->on_device, range);
}

void 
gkyl_array_copy_from_buffer_cu(struct gkyl_array *arr,
  const void *data, const struct gkyl_range range)
{
  gkyl_array_copy_from_buffer_cu_kernel<<<range.nblocks, range.nthreads>>>(arr->on_device, data, range);
}
