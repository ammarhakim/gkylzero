/* -*- c++ -*- */

extern "C" {
#include <gkyl_skin_surf_from_ghost.h>
#include <gkyl_skin_surf_from_ghost_priv.h>
}

// CUDA kernel to set device pointers to kernel that computes the reflected f.
__global__ static void
gkyl_skin_surf_from_ghost_set_cu_ker_ptrs(const struct gkyl_basis *basis,
  enum gkyl_edge_loc edge, struct gkyl_skin_surf_from_ghost_kernels *kers)
{
  int dim = basis->ndim;
  enum gkyl_basis_type b_type = basis->b_type;
  int poly_order = basis->poly_order;

  switch (b_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kers->ghost_to_skin = ser_skin_surf_from_ghost_list[edge].list[dim-2].kernels[poly_order-1];
      break;
    default:
      assert(false);
  }
};

void
gkyl_skin_surf_from_ghost_choose_kernel_cu(const struct gkyl_basis *basis,
  enum gkyl_edge_loc edge, struct gkyl_skin_surf_from_ghost_kernels *kers)
{
  gkyl_skin_surf_from_ghost_set_cu_ker_ptrs<<<1,1>>>(basis, edge, kers);
}

__global__ static void
gkyl_skin_surf_from_ghost_advance_cu_ker(int cdim, int dir, const struct gkyl_range skin_r, const struct gkyl_range ghost_r,
  const struct gkyl_basis *basis, struct gkyl_array *phi, struct gkyl_skin_surf_from_ghost_kernels *kers)
{
  int fidx[GKYL_MAX_DIM]; // Flipped index.
  int pidx[GKYL_MAX_DIM];

  int pdim = skin_r.ndim;
  int vpar_dir = cdim;
  int uplo = skin_r.upper[vpar_dir]+skin_r.lower[vpar_dir];

  for(unsigned long linc = threadIdx.x + blockIdx.x*blockDim.x;
      linc < skin_r.volume; linc += blockDim.x*gridDim.x) {

    // inverse index from linc1 to idx
    // must use gkyl_sub_range_inv_idx so that linc1=0 maps to idx={1,1,...}
    // since update_range is a subrange

    gkyl_sub_range_inv_idx(&skin_r, linc, pidx);

    gkyl_copy_int_arr(pdim, pidx, fidx);
    fidx[vpar_dir] = uplo - pidx[vpar_dir];
    // Turn this skin fidx into a ghost fidx.
    fidx[dir] = ghost_r.lower[dir];

    long skin_loc = gkyl_range_idx(&skin_r, pidx);
    long ghost_loc = gkyl_range_idx(&ghost_r, fidx);

    const double *inp = (const double*) gkyl_array_cfetch(phi, ghost_loc);
    double *out = (double*) gkyl_array_fetch(phi, skin_loc);

    // Write something like out = inp
    kers->ghost_to_skin(inp,out);
  }
}

void
gkyl_skin_surf_from_ghost_advance_cu(const struct gkyl_skin_surf_from_ghost *up, struct gkyl_array *phi)
{
  if (up->skin_r->volume > 0) {
    int nblocks = up->skin_r->nblocks, nthreads = up->skin_r->nthreads;

    gkyl_skin_surf_from_ghost_advance_cu_ker<<<nblocks, nthreads>>>(up->cdim, up->dir, *up->skin_r, *up->ghost_r, 
      up->basis, phi->on_dev, up->kernels_cu);
  }
}
