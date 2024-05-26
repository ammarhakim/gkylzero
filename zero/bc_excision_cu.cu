/* -*- c++ -*- */

extern "C" {
#include <gkyl_bc_excision.h>
#include <gkyl_bc_excision_priv.h>
}

__global__ static void
gkyl_bc_excision_advance_cu_ker(int tan_dir, int tan_cellsD2, int num_basis,
  const struct gkyl_range ghost_r, const struct gkyl_range buff_r, 
  const struct gkyl_array *ghost_buffer, struct gkyl_array *distf)
{
  int idx[GKYL_MAX_DIM];
  int sidx[GKYL_MAX_DIM]; // Index shifted along the direction tangential to the boundary.

  for (unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < ghost_r.volume; linc += blockDim.x*gridDim.x) {

    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange
    gkyl_sub_range_inv_idx(&ghost_r, linc, idx);

    gkyl_copy_int_arr(ghost_r.ndim, idx, sidx);
    // Shift the index in the direction tangential to the boundary so that we place
    // f from a ghost cell in the ghost cell on the opposite side of the excision.
    int tan_idx = idx[tan_dir];
    sidx[tan_dir] = tan_idx > tan_cellsD2 ? tan_idx-tan_cellsD2 : tan_idx+tan_cellsD2;

    long buff_loc = gkyl_range_idx(&buff_r, sidx);
    long out_loc = gkyl_range_idx(&ghost_r, idx);

    const double *buff = (const double*) gkyl_array_cfetch(ghost_buffer, buff_loc);
    double *out = (double*) gkyl_array_fetch(distf, out_loc);

    for (int i=0; i<num_basis; i++)
      out[i] = buff[i];
  }
}

void
gkyl_bc_excision_advance_cu(const struct gkyl_bc_excision *up, const struct gkyl_array *ghost_buffer,
  struct gkyl_array *distf)
{
  if (up->ghost_r.volume > 0) {
    int nblocks = up->ghost_r.nblocks, nthreads = up->ghost_r.nthreads;

    gkyl_bc_excision_advance_cu_ker<<<nblocks, nthreads>>>(up->tan_dir, up->tan_dir_num_cellsD2,
      up->num_basis, up->ghost_r, up->buff_r, ghost_buffer->on_dev, distf->on_dev);
  }
}
