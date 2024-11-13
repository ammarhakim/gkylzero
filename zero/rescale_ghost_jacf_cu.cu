/* -*- c++ -*- */

extern "C" {
#include <gkyl_rescale_ghost_jacf.h>
#include <gkyl_rescale_ghost_jacf_priv.h>
#include <gkyl_array_ops.h>
#include <float.h>
}

// CUDA kernel to set device pointers to kernels.
__global__ static void
gkyl_pghost_cdm_set_cu_ker_ptrs(struct gkyl_rescale_ghost_jacf_kernels *kernels,
  const struct gkyl_basis *cbasis, struct gkyl_basis pbasis)
{
  int cdim = cbasis->ndim;
  enum gkyl_basis_type cbasis_type = cbasis->b_type;
  int cpoly_order = cbasis->poly_order;
  switch (cbasis_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->conf_inv_op = choose_ser_inv_kern(cdim, cpoly_order);
      kernels->conf_mul_op = choose_ser_mul_kern(cdim, cpoly_order);
      break;
    default:
      assert(false);
      break;
  }

  int pdim = pbasis.ndim;
  enum gkyl_basis_type pbasis_type = pbasis.b_type;
  int ppoly_order = pbasis.poly_order;
  switch (pbasis_type) {
    case GKYL_BASIS_MODAL_GKHYBRID:
    case GKYL_BASIS_MODAL_SERENDIPITY:
      kernels->conf_phase_mul_op = choose_mul_conf_phase_kern(pbasis_type, cdim, pdim-cdim, ppoly_order);
      break;
    default:
      assert(false);
      break;
  }
};

void
pghost_cdm_choose_kernel_cu(struct gkyl_rescale_ghost_jacf_kernels *kernels,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis)
{
  gkyl_pghost_cdm_set_cu_ker_ptrs<<<1,1>>>(kernels, cbasis, *pbasis);
}

// CUDA kernel to copy ghost cell values to the adjacent skin (boundary) cells on the GPU.
__global__ static void
gkyl_rescale_ghost_jacf_advance_cu_ker(struct gkyl_rescale_ghost_jacf_kernels *kers,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_basis *conf_basis,
  const struct gkyl_range conf_skin_r, const struct gkyl_range conf_ghost_r,
  const struct gkyl_range phase_ghost_r, const struct gkyl_array *jac, struct gkyl_array *jf)
{
  int sidx[GKYL_MAX_DIM]; // skin idx
  int gidx[GKYL_MAX_DIM]; // ghost idx

  int cdim = conf_skin_r.ndim;

  // Loop over all points in the skin range using CUDA threads.
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x;
       linc < phase_ghost_r.volume; linc += blockDim.x * gridDim.x) {

    // Convert the linear index to a multi-dimensional index in the skin range.
    gkyl_sub_range_inv_idx(&phase_ghost_r, linc, gidx);

    // Get skin cell corresponding to this ghost cell.
    gkyl_copy_int_arr(cdim, gidx, sidx);
    sidx[dir] = edge == GKYL_LOWER_EDGE? gidx[dir]+1 : gidx[dir]-1;

    // Compute the linear indices for both skin and ghost locations.
    long clinidx_skin = gkyl_range_idx(&conf_skin_r, sidx);
    long clinidx_ghost = gkyl_range_idx(&conf_ghost_r, gidx);

    const double *jacskin_c = (const double *) gkyl_array_cfetch(jac, clinidx_skin);
    const double *jacghost_c = (const double *) gkyl_array_cfetch(jac, clinidx_ghost);

    // Calculate the reciprocal of jac in the ghost cell.
    double jacghost_inv_c[20]; // MF 2024/10/30: hardcoded to 3x p=2 ser for now.
    kers->conf_inv_op(jacghost_c, jacghost_inv_c);

    // Flip the jactor in the skin cell in the direction of the boundary. 
    double jacskin_flipped_c[20]; // MF 2024/10/30: hardcoded to 3x p=2 ser for now.
    conf_basis->flip_odd_sign(dir, jacskin_c, jacskin_flipped_c);

    // Multiply the ghost cell jf by the reciprocal of j_ghost and by
    // the flipped j_skin.
    long plinidx_ghost = gkyl_range_idx(&phase_ghost_r, gidx);
    double *jf_c = (double *) gkyl_array_fetch(jf, plinidx_ghost);

    kers->conf_phase_mul_op(jacghost_inv_c, jf_c, jf_c);
    kers->conf_phase_mul_op(jacskin_flipped_c, jf_c, jf_c);
  }
}

// Function to launch the CUDA kernel that performs the ghost-to-skin value transfer on the GPU.
void
gkyl_rescale_ghost_jacf_advance_cu(const struct gkyl_rescale_ghost_jacf *up,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_range *conf_skin_r, const struct gkyl_range *conf_ghost_r,
  const struct gkyl_range *phase_ghost_r, const struct gkyl_array *jac, struct gkyl_array *jf)
{
  // Only proceed if the skin range has a non-zero volume (i.e., there are skin cells to update).
  int nblocks = phase_ghost_r->nblocks, nthreads = phase_ghost_r->nthreads; // CUDA grid configuration.

  // Launch the CUDA kernel to advance the ghost-to-skin update.
  gkyl_rescale_ghost_jacf_advance_cu_ker<<<nblocks, nthreads>>>(up->kernels, dir, edge,
    up->conf_basis, *conf_skin_r, *conf_ghost_r, *phase_ghost_r, jac->on_dev, jf->on_dev);
}
