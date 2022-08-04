#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_bc_sheath_gyrokinetic_priv.h>
#include <gkyl_alloc.h>
#include <assert.h>

struct gkyl_bc_sheath_gyrokinetic*
gkyl_bc_sheath_gyrokinetic_new(int dir, enum gkyl_edge_loc edge, const struct gkyl_range *local_range_ext,
  const int *num_ghosts, const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid, int cdim,
  double q2Dm, bool use_gpu)
{

  // Allocate space for new updater.
  struct gkyl_bc_sheath_gyrokinetic *up = gkyl_malloc(sizeof(struct gkyl_bc_sheath_gyrokinetic));

  up->dir = dir;
  up->cdim = cdim;
  up->edge = edge;
  up->use_gpu = use_gpu;
  up->q2Dm = q2Dm;
  up->basis = basis;
  up->grid = grid;

  // Create the skin/ghost ranges.
  gkyl_skin_ghost_ranges(&up->skin_r, &up->ghost_r, dir, edge,
                         local_range_ext, num_ghosts);

  // Need the configuration space range to index into phi.
  int rlo[3] = {0}, rup[3] = {0};
  for (int d=0; d<cdim; d++) {
    rlo[d] = local_range_ext->lower[d];
    rup[d] = local_range_ext->upper[d];
  }
  gkyl_sub_range_init(&up->conf_r, local_range_ext, rlo, rup);

  // Choose the kernel that does the reflection/no reflection/partial
  // reflection.
  up->ker_reflectedf = bc_gksheath_choose_reflectedf_kernel(basis->ndim, basis->b_type, basis->poly_order, edge);
  assert(up->ker_reflectedf);

  return up;
}

/* Modeled after gkyl_array_flip_copy_to_buffer_fn */
void
bc_sheath_gyrokinetic_advance(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, struct gkyl_array *distf)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(arr)) {
    bc_sheath_gyrokinetic_advance_cu(up, phi, phi_wall, distf);
    return;
  }
#endif

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->skin_r);

  int fidx[GKYL_MAX_DIM]; // Flipped index.
  double xc[GKYL_MAX_DIM];

  int vpar_dir = up->cdim;
  double dvpar = up->grid->dx[vpar_dir];
  int uplo = up->skin_r.upper[vpar_dir]+up->skin_r.lower[vpar_dir];

  while (gkyl_range_iter_next(&iter)) {

    gkyl_copy_int_arr(up->skin_r.ndim, iter.idx, fidx);
    fidx[vpar_dir] = uplo - iter.idx[vpar_dir];
    // Turn this skin fidx into a ghost fidx.
    fidx[up->dir] = up->ghost_r.lower[up->dir];

    gkyl_rect_grid_cell_center(up->grid, iter.idx, xc);
    double vpar_c = xc[vpar_dir];
    double vparAbsSq_lo = vpar_c > 0.? pow(vpar_c-0.5*dvpar,2) : pow(vpar_c+0.5*dvpar,2);
    double vparAbsSq_up = vpar_c > 0.? pow(vpar_c+0.5*dvpar,2) : pow(vpar_c-0.5*dvpar,2);

    long skin_loc = gkyl_range_idx(&up->skin_r, iter.idx);
    long ghost_loc = gkyl_range_idx(&up->ghost_r, fidx);

    const double *inp = gkyl_array_cfetch(distf, skin_loc);
    double *out = gkyl_array_fetch(distf, ghost_loc);

    long conf_loc = gkyl_range_idx(&up->conf_r, iter.idx);
    const double *phi_p = gkyl_array_cfetch(phi, conf_loc);
    const double *phi_wall_p = gkyl_array_cfetch(phi_wall, conf_loc);

    // Calculate reflected distribution function fhat.
    // note: reflected distribution can be
    // 1) fhat=0 (no reflection, i.e. absorb),
    // 2) fhat=f (full reflection)
    // 3) fhat=c*f (partial reflection)
    double fhat[up->basis->num_basis];
    up->ker_reflectedf(vpar_c, dvpar, vparAbsSq_lo, vparAbsSq_up, up->q2Dm, phi_p, phi_wall_p, inp, fhat);

    // Reflect fhat into skin cells.
    bc_gksheath_reflect(up->dir, up->basis, up->cdim, out, fhat);
  }
}

void gkyl_bc_sheath_gyrokinetic_release(struct gkyl_bc_sheath_gyrokinetic *up)
{
  // Release memory associated with array_copy_func.
//  if (up->use_gpu) {
//  }
  // Release updater memory.
  gkyl_free(up);
}
