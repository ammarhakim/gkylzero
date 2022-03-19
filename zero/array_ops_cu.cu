/* -*- c++ -*- */

#include <math.h>
#include <time.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_util.h>
}

// start ID for use in various loops
#define START_ID (threadIdx.x + blockIdx.x*blockDim.x)

// NOTE: This is duplicated in dg_bin_ops_cu. Should be cleaned up 01/05/22
static void
gkyl_get_array_range_kernel_launch_dims(dim3* dimGrid, dim3* dimBlock, gkyl_range range, int ncomp)
{
  int volume = range.volume;
  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  int ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;
  dimBlock->x = min(ncomp*ac1, GKYL_DEFAULT_NUM_THREADS);
  dimGrid->x = gkyl_int_div_up(ncomp*ac1, dimBlock->x);

  dimBlock->y = gkyl_int_div_up(GKYL_DEFAULT_NUM_THREADS, ncomp*ac1);
  dimGrid->y = gkyl_int_div_up(volume, ac1*dimBlock->y);
}

__global__ void
gkyl_array_clear_cu_kernel(struct gkyl_array* out, double val)
{
  double *out_d = (double*) out->data;
  for (unsigned long linc = START_ID; linc < NELM(out); linc += blockDim.x*gridDim.x)
    out_d[linc] = val;
}

__global__ void
gkyl_array_accumulate_cu_kernel(struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  double *out_d = (double*) out->data;
  const double *inp_d = (const double*) inp->data;
  for (unsigned long linc = START_ID; linc < NELM(out); linc += blockDim.x*gridDim.x)
    out_d[linc] += a*inp_d[linc];
}

__global__ void
gkyl_array_set_cu_kernel(struct gkyl_array* out, double a,
  const struct gkyl_array* inp)
{
  double *out_d = (double*) out->data;
  const double *inp_d = (const double*) inp->data;
  for (unsigned long linc = START_ID; linc < NELM(out); linc += blockDim.x*gridDim.x)
    out_d[linc] = a*inp_d[linc];
}

__global__ void
gkyl_array_scale_by_cell_cu_kernel(struct gkyl_array* out, const struct gkyl_array* a)
{
  double *out_d = (double*) out->data;
  const double *a_d = (double*) a->data;
  for (unsigned long linc = START_ID; linc < NELM(out); linc += blockDim.x*gridDim.x)
    out_d[linc] = a_d[linc/out->ncomp]*out_d[linc];
} 

// Host-side wrappers for array operations
void
gkyl_array_clear_cu(struct gkyl_array* out, double val)
{
  gkyl_array_clear_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_dev, val);
}

void
gkyl_array_accumulate_cu(struct gkyl_array* out, double a, const struct gkyl_array* inp)
{
  gkyl_array_accumulate_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_dev, a, inp->on_dev);
}

void
gkyl_array_set_cu(struct gkyl_array* out, double a, const struct gkyl_array* inp)
{
  gkyl_array_set_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_dev, a, inp->on_dev);
}

void
gkyl_array_scale_cu(struct gkyl_array* out, double a)
{
  gkyl_array_set_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_dev, a, out->on_dev);
}

void
gkyl_array_scale_by_cell_cu(struct gkyl_array* out, const struct gkyl_array* a)
{
  gkyl_array_scale_by_cell_cu_kernel<<<out->nblocks, out->nthreads>>>(out->on_dev, a->on_dev);
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
  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc1 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
  // linc2 = idx2 + ac2*idx3 + ...
  for (unsigned long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
      linc2 < range.volume/ac1;
      linc2 += gridDim.y*blockDim.y)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc2.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&range, ac1*linc2, idx);
    long start = gkyl_range_idx(&range, idx);
    
    double* out_d = (double*) gkyl_array_fetch(out, start);
    // do operation on contiguous data block
    if (linc1 < n*ac1)
      out_d[linc1] = val;
  }
}

__global__ void
gkyl_array_accumulate_range_cu_kernel(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, struct gkyl_range range)
{
  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;
  int idx[GKYL_MAX_DIM];

  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc1 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
  long c = linc1 % n;
  long idx1 = linc1 / n;
  // get corresponding linc1 index for inp and out 
  // (one of these will not be contiguous if outnc!=inpnc)
  long linc1_in = c + inpnc*idx1; 
  long linc1_out = c + outnc*idx1; 
  // linc2 = idx2 + ac2*idx3 + ...
  for (unsigned long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
      linc2 < range.volume/ac1;
      linc2 += gridDim.y*blockDim.y)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc2.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&range, ac1*linc2, idx);
    long start = gkyl_range_idx(&range, idx);
    
    double* out_d = (double*) gkyl_array_fetch(out, start);
    const double* inp_d = (const double*) gkyl_array_cfetch(inp, start);
    // do operation on contiguous data block
    if (linc1 < n*ac1)
      out_d[linc1_out] += a*inp_d[linc1_in];
  }
}

__global__ void
gkyl_array_set_range_cu_kernel(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, struct gkyl_range range)
{
  long outnc = NCOM(out), inpnc = NCOM(inp);
  long n = outnc<inpnc ? outnc : inpnc;
  int idx[GKYL_MAX_DIM];

  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc1 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
  long c = linc1 % n;
  long idx1 = linc1 / n;
  // get corresponding linc1 index for inp and out 
  // (one of these will not be contiguous if outnc!=inpnc)
  long linc1_in = c + inpnc*idx1; 
  long linc1_out = c + outnc*idx1; 
  // linc2 = idx2 + ac2*idx3 + ...
  for (unsigned long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
      linc2 < range.volume/ac1;
      linc2 += gridDim.y*blockDim.y)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc2.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&range, ac1*linc2, idx);
    long start = gkyl_range_idx(&range, idx);
    
    double* out_d = (double*) gkyl_array_fetch(out, start);
    const double* inp_d = (const double*) gkyl_array_cfetch(inp, start);
    // do operation on contiguous data block
    if (linc1 < n*ac1)
      out_d[linc1_out] = a*inp_d[linc1_in];
  }
}

__global__ void 
gkyl_array_copy_range_cu_kernel(struct gkyl_array *out, const struct gkyl_array* inp,
  struct gkyl_range out_range, struct gkyl_range inp_range)
{
  int idx_out[GKYL_MAX_DIM], idx_inp[GKYL_MAX_DIM];
  long n = NCOM(out); // assume ncomp_in == ncomp_out
  int ndim = inp_range.ndim;
  // ac1 = size of last dimension of inp_range (fastest moving dimension)
  long ac1 = inp_range.iac[ndim-1] > 0 ? inp_range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc1 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
  // linc2 = idx2 + ac2*idx3 + ...
  for (unsigned long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
      linc2 < inp_range.volume/ac1;
      linc2 += gridDim.y*blockDim.y)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc2.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    // NOTE: the above necessarily applies only to inp_range
    gkyl_sub_range_inv_idx(&out_range, ac1*linc2, idx_out);
    gkyl_sub_range_inv_idx(&inp_range, ac1*linc2, idx_inp);
    long start_out = gkyl_range_idx(&out_range, idx_out);
    long start_inp = gkyl_range_idx(&inp_range, idx_inp);
    
    double* out_d = (double*) gkyl_array_fetch(out, start_out);
    const double* inp_d = (const double*) gkyl_array_cfetch(inp, start_inp);
    // do operation on contiguous data block
    if (linc1 < n*ac1)
      out_d[linc1] = inp_d[linc1];
  }
}

__global__ void 
gkyl_array_copy_to_buffer_cu_kernel(void *data, const struct gkyl_array *arr,
  struct gkyl_range range)
{
  double *d_data = (double*) data;
  int idx[GKYL_MAX_DIM];
  long n = NCOM(arr); // assume ncomp_in == ncomp_out
  int ndim = range.ndim;
  // ac1 = size of last dimension of range (fastest moving dimension)
  long ac1 = range.iac[ndim-1] > 0 ? range.iac[ndim-1] : 1;

  // 2D thread grid
  // linc1 = c + n*idx1 (contiguous data, including component index c, with idx1 = 0,.., ac1-1)
  long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
  // linc2 = idx2 + ac2*idx3 + ...
  for (unsigned long linc2 = threadIdx.y + blockIdx.y*blockDim.y;
      linc2 < range.volume/ac1;
      linc2 += gridDim.y*blockDim.y)
  {
    // full linear cell index (not including components) is 
    // idx1 + ac1*idx2 + ac1*ac2*idx3 + ... = idx1 + ac1*linc2.
    // we want to find the start linear index of each contiguous data block, 
    // which corresponds to idx1 = 0. 
    // so linear index of start of contiguous block is ac1*linc2.
    gkyl_sub_range_inv_idx(&range, ac1*linc2, idx);
    long start = gkyl_range_idx(&range, idx);
    
    const double* arr_d = (const double*) gkyl_array_cfetch(arr, start);
    // read from contiguous data block
    if (linc1 < n*ac1)
      d_data[linc1 + n*ac1*linc2] = arr_d[linc1];
  }
}

__global__ void 
gkyl_array_copy_from_buffer_cu_kernel(struct gkyl_array *arr, const void *data,
  struct gkyl_range range)
{
  int idx[GKYL_MAX_DIM];
  const double *d_data = (const double*) data;
  long n = NCOM(arr);
  
  // since input data is just a linear array, just stream through data linearly
  // linc = c + n*idx1 + n*ac1*idx2 + ...
  for (unsigned long linc = START_ID; linc < range.volume*n; linc += blockDim.x*gridDim.x) {
    int c = linc % n;
    long linc2 = linc / n; // = idx1 + ac1*idx2 + ...
    gkyl_sub_range_inv_idx(&range, linc2, idx);
    long start = gkyl_range_idx(&range, idx);
    
    double *arr_data = (double*) gkyl_array_fetch(arr, start);
    arr_data[c] = d_data[linc];
  }
}

__global__ void 
gkyl_array_copy_to_buffer_fn_cu_kernel(void *data, const struct gkyl_array *arr,
  struct gkyl_range range, struct gkyl_array_copy_func *cf)
{
  long count = 0;
  int idx[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume; linc1 += blockDim.x*gridDim.x) {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    // convert back to a linear index on the super-range (with ghost cells)
    // linc will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&range, idx);

    const double *inp = (const double*) gkyl_array_cfetch(arr, linc);
    double *out = (double*) flat_fetch(data, arr->esznc*count);
    cf->func(arr->ncomp, out, inp, cf->ctx);

    count += 1;
  }
}

__global__ void 
gkyl_array_flip_copy_to_buffer_fn_cu_kernel(void *data, const struct gkyl_array *arr,
  int dir, struct gkyl_range range, struct gkyl_range buff_range,
  struct gkyl_array_copy_func *cf)
{
  int idx[GKYL_MAX_DIM];
  int fidx[GKYL_MAX_DIM]; // flipped index

  int uplo = range.upper[dir]+range.lower[dir];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x*blockDim.x;
      linc1 < range.volume; linc1 += blockDim.x*gridDim.x) {
    // inverse index from linc1 to idxc
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idxc={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&range, linc1, idx);

    gkyl_copy_int_arr(range.ndim, idx, fidx);
    fidx[dir] = uplo - idx[dir];

    // convert back to a linear index on the super-range (with ghost cells)
    // linc and flipped linc (flinc) will have jumps in it to jump over ghost cells
    long linc = gkyl_range_idx(&range, idx);
    long flinc = gkyl_range_idx(&buff_range, fidx);

    const double *inp = (const double*) gkyl_array_cfetch(arr, linc);
    double *out = (double*) flat_fetch(data, arr->esznc*flinc);
    cf->func(arr->ncomp, out, inp, cf->ctx);
  }
}

// Host-side wrappers for range-based array operations
void
gkyl_array_clear_range_cu(struct gkyl_array *out, double val, struct gkyl_range range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, range, out->ncomp);

  gkyl_array_clear_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev, val, range);
}

void
gkyl_array_accumulate_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, struct gkyl_range range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, range, min(out->ncomp, inp->ncomp));

  gkyl_array_accumulate_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev, a, inp->on_dev, range);
}

void
gkyl_array_set_range_cu(struct gkyl_array *out,
  double a, const struct gkyl_array* inp, struct gkyl_range range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, range, min(out->ncomp, inp->ncomp));

  gkyl_array_set_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev,
    a, inp->on_dev, range);
}

void
gkyl_array_scale_range_cu(struct gkyl_array *out,
  double a, struct gkyl_range range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, range, out->ncomp);

  gkyl_array_set_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev,
    a, out->on_dev, range);
}

void
gkyl_array_copy_range_cu(struct gkyl_array *out,
  const struct gkyl_array *inp, struct gkyl_range range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, range, out->ncomp);

  gkyl_array_copy_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev,
    inp->on_dev, range, range);
}

void
gkyl_array_copy_range_to_range_cu(struct gkyl_array *out,
  const struct gkyl_array *inp, struct gkyl_range out_range, struct gkyl_range inp_range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, inp_range, out->ncomp);

  gkyl_array_copy_range_cu_kernel<<<dimGrid, dimBlock>>>(out->on_dev,
    inp->on_dev, out_range, inp_range);
}

void 
gkyl_array_copy_to_buffer_cu(void *data, 
  const struct gkyl_array *arr, struct gkyl_range range)
{
  dim3 dimGrid, dimBlock;
  gkyl_get_array_range_kernel_launch_dims(&dimGrid, &dimBlock, range, arr->ncomp);

  gkyl_array_copy_to_buffer_cu_kernel<<<dimGrid, dimBlock>>>(data,
    arr->on_dev, range);
}

void 
gkyl_array_copy_from_buffer_cu(struct gkyl_array *arr,
  const void *data, struct gkyl_range range)
{
  int nelem = range.volume*arr->ncomp;
  int nthreads = GKYL_DEFAULT_NUM_THREADS;
  int nblocks = gkyl_int_div_up(nelem, nthreads);
  gkyl_array_copy_from_buffer_cu_kernel<<<nblocks, nthreads>>>(arr->on_dev,
    data, range);
}

void
gkyl_array_copy_to_buffer_fn_cu(void *data, const struct gkyl_array *arr,
  struct gkyl_range range, struct gkyl_array_copy_func *cf)
{
  int nblocks = range.nblocks;
  int nthreads = range.nthreads;

  gkyl_array_copy_to_buffer_fn_cu_kernel<<<nblocks, nthreads>>>(data, arr, range, cf);
}

void
gkyl_array_flip_copy_to_buffer_fn_cu(void *data, const struct gkyl_array *arr,
  int dir, struct gkyl_range range, struct gkyl_array_copy_func *cf)
{
  int nblocks = range.nblocks;
  int nthreads = range.nthreads;

  struct gkyl_range buff_range;
  gkyl_range_init(&buff_range, range.ndim, range.lower, range.upper);
  
  gkyl_array_flip_copy_to_buffer_fn_cu_kernel<<<nblocks, nthreads>>>(data, arr, dir, range,
    buff_range, cf);
}
