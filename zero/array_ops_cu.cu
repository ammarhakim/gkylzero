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
gkyl_array_clear_range_cu_kernel(struct gkyl_array *out, double val, struct gkyl_range range)
{
  long n = NCOM(out);
  int idx[GKYL_MAX_DIM];
  for(unsigned long linc = blockIdx.x;
      linc < range.volume;
      linc += gridDim.x)
  {
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    double* out_d = (double*) gkyl_array_fetch(out, start);
    if(threadIdx.x < n)
      out_d[threadIdx.x] = val;
  }
}

__global__ void
gkyl_array_accumulate_range_cu_kernel(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, struct gkyl_range range)
{
  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;
  int idx[GKYL_MAX_DIM];
  for(unsigned long linc = blockIdx.x;
      linc < range.volume;
      linc += gridDim.x)
  {
    gkyl_sub_range_inv_idx(&range, linc, idx);
    long start = gkyl_range_idx(&range, idx);
    double* out_d = (double*) gkyl_array_fetch(out, start);
    double* inp_d = (double*) gkyl_array_cfetch(inp, start);
    if(threadIdx.x < n)
      out_d[threadIdx.x] += a*inp_d[threadIdx.x];
  }
}

__global__ void
gkyl_array_set_range_cu_kernel(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, struct gkyl_range range)
{
  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;
  int idx[GKYL_MAX_DIM];
  for(unsigned long linc = blockIdx.x;
      linc < range.volume;
      linc += gridDim.x)
  {
    gkyl_sub_range_inv_idx(&range, linc, idx); 
    long start = gkyl_range_idx(&range, idx);
    double* out_d = (double*) gkyl_array_fetch(out, start);
    double* inp_d = (double*) gkyl_array_cfetch(inp, start);
    if(threadIdx.x < n)
      out_d[threadIdx.x] = a*inp_d[threadIdx.x];
  }
}

__global__ void 
gkyl_array_copy_range_cu_kernel(struct gkyl_array *out, const struct gkyl_array* inp,
  struct gkyl_range out_range, struct gkyl_range inp_range)
{
  int idx_out[GKYL_MAX_DIM], idx_inp[GKYL_MAX_DIM];
  for(unsigned long linc = blockIdx.x;
      linc < out_range.volume;
      linc += gridDim.x) {

    gkyl_sub_range_inv_idx(&out_range, linc, idx_out);
    gkyl_sub_range_inv_idx(&inp_range, linc, idx_inp);
    long start_out = gkyl_range_idx(&out_range, idx_out);
    long start_inp = gkyl_range_idx(&inp_range, idx_inp);
    double *out_data = (double*) gkyl_array_fetch(out, start_out);
    const double *inp_data = (const double*) gkyl_array_cfetch(inp, start_inp);
    if(threadIdx.x < NCOM(out))
      out_data[threadIdx.x] = inp_data[threadIdx.x];
  }
}

__global__ void 
gkyl_array_copy_to_buffer_cu_kernel(void *data, const struct gkyl_array *arr,
  struct gkyl_range range, struct gkyl_range data_range)
{
  int idx[GKYL_MAX_DIM], idx_d[1];
  double *d_data = (double*) data;
  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < range.volume;
      linc += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&range, linc, idx);
    // Since data range is just a 1D index, can just use range_inv_idx
    gkyl_range_inv_idx(&data_range, linc, idx_d);
    long start = gkyl_range_idx(&range, idx);
    long start_d = gkyl_range_idx(&data_range, idx_d);
    const double *arr_data = (const double*) gkyl_array_cfetch(arr, start);
    for (unsigned i = 0; i < NCOM(arr); ++i)
      d_data[i+start_d*NCOM(arr)] = arr_data[i];
  }
}

__global__ void 
gkyl_array_copy_from_buffer_cu_kernel(struct gkyl_array *arr, const void *data,
  struct gkyl_range range, struct gkyl_range data_range)
{
  int idx[GKYL_MAX_DIM], idx_d[1];
  const double *d_data = (const double*) data;
  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < range.volume;
      linc += blockDim.x*gridDim.x) {

    gkyl_sub_range_inv_idx(&range, linc, idx);
    // Since data range is just a 1D index, can just use range_inv_idx
    gkyl_range_inv_idx(&data_range, linc, idx_d);
    long start = gkyl_range_idx(&range, idx);
    long start_d = gkyl_range_idx(&data_range, idx_d);
    double *arr_data = (double*) gkyl_array_fetch(arr, start);
    for (unsigned i = 0; i < NCOM(arr); ++i)
      arr_data[i] = d_data[i+start_d*NCOM(arr)];
  }
}

// Host-side wrappers for range-based array operations
void
gkyl_array_clear_range_cu(struct gkyl_array *out, double val, struct gkyl_range range)
{
  int nthreads = out->ncomp;
  int nblocks = range.volume;
  gkyl_array_clear_range_cu_kernel<<<nblocks, nthreads>>>(out->on_device, val, range);
}

void
gkyl_array_accumulate_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, struct gkyl_range range)
{
  int nthreads = out->ncomp;
  int nblocks = range.volume;
  gkyl_array_accumulate_range_cu_kernel<<<nblocks, nthreads>>>(out->on_device, a, inp->on_device, range);
}

void
gkyl_array_set_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, struct gkyl_range range)
{
  int nthreads = out->ncomp;
  int nblocks = range.volume;
  gkyl_array_set_range_cu_kernel<<<nblocks, nthreads>>>(out->on_device, a, inp->on_device, range);
}

void
gkyl_array_scale_range_cu(struct gkyl_array *out,
  double a, struct gkyl_range range)
{
  int nthreads = out->ncomp;
  int nblocks = range.volume;
  gkyl_array_set_range_cu_kernel<<<nblocks, nthreads>>>(out->on_device, a, out->on_device, range);
}

void
gkyl_array_copy_range_cu(struct gkyl_array *out,
  const struct gkyl_array *inp, struct gkyl_range range)
{
  int nthreads = inp->ncomp;
  int nblocks = range.volume;
  gkyl_array_copy_range_cu_kernel<<<nblocks, nthreads>>>(out->on_device, inp->on_device, range, range);
}

void
gkyl_array_copy_range_to_range_cu(struct gkyl_array *out,
  const struct gkyl_array *inp, struct gkyl_range out_range, struct gkyl_range inp_range)
{
  gkyl_array_copy_range_cu_kernel<<<out_range.nblocks, out_range.nthreads>>>(out->on_device, inp->on_device, out_range, inp_range);
}

void 
gkyl_array_copy_to_buffer_cu(void *data, 
  const struct gkyl_array *arr, struct gkyl_range range)
{
  gkyl_range data_range;
  int shape[1] = {range.volume};
  gkyl_range_init_from_shape(&data_range, 1, shape);
  gkyl_array_copy_to_buffer_cu_kernel<<<range.nblocks, range.nthreads>>>(data, arr->on_device, range, data_range);
}

void 
gkyl_array_copy_from_buffer_cu(struct gkyl_array *arr,
  const void *data, struct gkyl_range range)
{
  gkyl_range data_range;
  int shape[1] = {range.volume};
  gkyl_range_init_from_shape(&data_range, 1, shape);
  gkyl_array_copy_from_buffer_cu_kernel<<<range.nblocks, range.nthreads>>>(arr->on_device, data, range, data_range);
}
