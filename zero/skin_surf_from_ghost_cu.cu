/* -*- c++ -*- */

extern "C" {
#include <gkyl_skin_surf_from_ghost.h>
#include <gkyl_skin_surf_from_ghost_priv.h>
}

// CUDA kernel to set device pointers to the kernel that transfers ghost cell values to skin cells.
__global__ static void
gkyl_skin_surf_from_ghost_set_cu_ker_ptrs(const struct gkyl_basis *basis,
  enum gkyl_edge_loc edge, struct gkyl_skin_surf_from_ghost_kernels *kers)
{
  // Get the dimension and basis type information from the provided basis object.
  int dim = basis->ndim; // Spatial dimension.
  enum gkyl_basis_type b_type = basis->b_type; // Basis function type.
  int poly_order = basis->poly_order; // Polynomial order of the basis.

  // Select the appropriate kernel based on the basis type and edge location.
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      // Set the ghost_to_skin kernel for the chosen basis and dimensionality.
      kers->ghost_to_skin = ser_skin_surf_from_ghost_list[edge].list[dim-1].kernels[poly_order-1];
      break;
    default:
      // If an unsupported basis type is encountered, assert failure.
      assert(false);
  }
};

// Function to launch a CUDA kernel that selects the appropriate kernel on the GPU.
void
gkyl_skin_surf_from_ghost_choose_kernel_cu(const struct gkyl_basis *basis,
  enum gkyl_edge_loc edge, struct gkyl_skin_surf_from_ghost_kernels *kers)
{
  // Launch the kernel with a single thread to set the kernel pointers.
  gkyl_skin_surf_from_ghost_set_cu_ker_ptrs<<<1,1>>>(basis, edge, kers);
}

// CUDA kernel to copy ghost cell values to the adjacent skin (boundary) cells on the GPU.
__global__ static void
gkyl_skin_surf_from_ghost_advance_cu_ker(int cdim, int dir, const struct gkyl_range skin_r, const struct gkyl_range ghost_r,
  const struct gkyl_basis *basis, struct gkyl_array *field, struct gkyl_skin_surf_from_ghost_kernels *kers)
{
  int fidx[GKYL_MAX_DIM]; // Flipped index to map ghost to skin.
  int pidx[GKYL_MAX_DIM]; // Original index in the skin range.

  int pdim = skin_r.ndim; // Problem dimensionality.
  int vpar_dir = cdim; // Parallel direction index.
  int uplo = skin_r.upper[vpar_dir] + skin_r.lower[vpar_dir]; // Boundary offset for reflection.

  // Loop over all points in the skin range using CUDA threads.
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x;
       linc < skin_r.volume; linc += blockDim.x * gridDim.x) {

    // Convert the linear index to a multi-dimensional index in the skin range.
    gkyl_sub_range_inv_idx(&skin_r, linc, pidx);

    // Copy the current index and flip the direction to the ghost side.
    gkyl_copy_int_arr(pdim, pidx, fidx);
    fidx[vpar_dir] = uplo - pidx[vpar_dir]; // Reflect across the boundary.
    fidx[dir] = ghost_r.lower[dir]; // Set the correct ghost cell index.

    // Compute the linear indices for both skin and ghost locations.
    long skin_loc = gkyl_range_idx(&skin_r, pidx);
    long ghost_loc = gkyl_range_idx(&ghost_r, fidx);

    // Fetch the values from the ghost region and copy them to the skin region.
    const double *inp = (const double*) gkyl_array_cfetch(field, ghost_loc);
    double *out = (double*) gkyl_array_fetch(field, skin_loc);

    // Apply the ghost_to_skin kernel to transfer the ghost values to the skin cells.
    kers->ghost_to_skin(inp, out);
  }
}

// Function to launch the CUDA kernel that performs the ghost-to-skin value transfer on the GPU.
void
gkyl_skin_surf_from_ghost_advance_cu(const struct gkyl_skin_surf_from_ghost *up, struct gkyl_array *field)
{
  // Only proceed if the skin range has a non-zero volume (i.e., there are skin cells to update).
  if (up->skin_r->volume > 0) {
    int nblocks = up->skin_r->nblocks, nthreads = up->skin_r->nthreads; // CUDA grid configuration.

    // Launch the CUDA kernel to advance the ghost-to-skin update.
    gkyl_skin_surf_from_ghost_advance_cu_ker<<<nblocks, nthreads>>>(up->cdim, up->dir, *up->skin_r, *up->ghost_r, 
      up->basis, field->on_dev, up->kernels_cu);
  }
}