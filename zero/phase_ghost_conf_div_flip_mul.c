#include <gkyl_phase_ghost_conf_div_flip_mul.h>
#include <gkyl_phase_ghost_conf_div_flip_mul_priv.h>
#include <gkyl_alloc.h>
#include <assert.h>

struct gkyl_phase_ghost_conf_div_flip_mul*
gkyl_phase_ghost_conf_div_flip_mul_new(const struct gkyl_basis *conf_basis,
  const struct gkyl_basis *phase_basis, bool use_gpu)
{

  // Allocate space for new updater.
  struct gkyl_phase_ghost_conf_div_flip_mul *up = gkyl_malloc(sizeof(*up));

  up->use_gpu = use_gpu;
  up->conf_basis = conf_basis;

  int poly_order = conf_basis->poly_order;
  assert(poly_order == 1); // MF 2024/10/30: Because of the inv_op below.

  // Choose the kernel that does the skin surf from ghost copy
  if (!use_gpu)
    up->kernels = gkyl_malloc(sizeof(struct gkyl_phase_ghost_conf_div_flip_mul_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels = gkyl_cu_malloc(sizeof(struct gkyl_phase_ghost_conf_div_flip_mul_kernels));
#endif

  phase_ghost_conf_div_flip_mul_choose_kernel(up->kernels, conf_basis, phase_basis, use_gpu);

  return up;
}

void
gkyl_phase_ghost_conf_div_flip_mul_advance(const struct gkyl_phase_ghost_conf_div_flip_mul *up,
  int dir, enum gkyl_edge_loc edge, const struct gkyl_range *conf_skin_r, const struct gkyl_range *conf_ghost_r,
  const struct gkyl_range *phase_ghost_r, const struct gkyl_array *jac, struct gkyl_array *jf)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_phase_ghost_conf_div_flip_mul_advance_cu(up, dir, edge,
      conf_skin_r, conf_ghost_r, phase_ghost_r, jac, jf);
    return;
  }
#endif

  int sidx[GKYL_MAX_DIM]; // Skin index.
  int cdim = conf_ghost_r->ndim; 

  struct gkyl_range vel_rng;
  struct gkyl_range_iter conf_iter, vel_iter;

  int rem_dir[GKYL_MAX_DIM] = {0};
  for (int d=0; d<cdim; ++d) rem_dir[d] = 1;

  gkyl_range_iter_init(&conf_iter, conf_ghost_r);
  while (gkyl_range_iter_next(&conf_iter)) {

    // Get skin cell corresponding to this ghost cell.
    gkyl_copy_int_arr(cdim, conf_iter.idx, sidx);
    sidx[dir] = edge == GKYL_LOWER_EDGE? conf_iter.idx[dir]+1 : conf_iter.idx[dir]-1; 

    long clinidx_skin  = gkyl_range_idx(conf_skin_r, sidx);
    long clinidx_ghost = gkyl_range_idx(conf_ghost_r, conf_iter.idx);

    const double *jacskin_c = gkyl_array_cfetch(jac, clinidx_skin);
    const double *jacghost_c = gkyl_array_cfetch(jac, clinidx_ghost);

    // Calculate the reciprocal of jac in the ghost cell.
    double jacghost_inv_c[up->conf_basis->num_basis];
    up->kernels->conf_inv_op(jacghost_c, jacghost_inv_c);

    // Flip the jactor in the skin cell in the direction of the boundary. 
    double jacskin_flipped_c[up->conf_basis->num_basis];
    up->conf_basis->flip_odd_sign(dir, jacskin_c, jacskin_flipped_c);

    // Multiply the ghost cell jf by the reciprocal of j_ghost and by
    // the flipped j_skin.
    gkyl_range_deflate(&vel_rng, phase_ghost_r, rem_dir, conf_iter.idx);
    gkyl_range_iter_no_split_init(&vel_iter, &vel_rng);
    while (gkyl_range_iter_next(&vel_iter)) {
      long plinidx_ghost = gkyl_range_idx(&vel_rng, vel_iter.idx);
      double *jf_c = gkyl_array_fetch(jf, plinidx_ghost);

      up->kernels->conf_phase_mul_op(jacghost_inv_c, jf_c, jf_c);
      up->kernels->conf_phase_mul_op(jacskin_flipped_c, jf_c, jf_c);
    }
  }
}

void
gkyl_phase_ghost_conf_div_flip_mul_release(struct gkyl_phase_ghost_conf_div_flip_mul *up)
{
  // Release memory associated with this updater.
  if (!up->use_gpu)
    gkyl_free(up->kernels);
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu)
    gkyl_cu_free(up->kernels);
#endif
  gkyl_free(up);
}
