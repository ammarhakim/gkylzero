#include <gkyl_bc_sheath_gyrokinetic.h>
#include <gkyl_bc_sheath_gyrokinetic_priv.h>
#include <gkyl_alloc.h>
#include <assert.h>

struct gkyl_bc_sheath_gyrokinetic*
gkyl_bc_sheath_gyrokinetic_new(int dir, enum gkyl_edge_loc edge, const struct gkyl_basis *basis,
  const struct gkyl_range *skin_r, const struct gkyl_range *ghost_r,
  const struct gkyl_rect_grid *grid, int cdim, double q2Dm, bool use_gpu)
{

  // Allocate space for new updater.
  struct gkyl_bc_sheath_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  up->dir = dir;
  up->cdim = cdim;
  up->edge = edge;
  up->use_gpu = use_gpu;
  up->q2Dm = q2Dm;
  up->basis = basis;
  up->grid = grid;
  up->skin_r = skin_r;
  up->ghost_r = ghost_r;

  // Choose the kernel that does the reflection/no reflection/partial
  // reflection.
  up->kernels = gkyl_malloc(sizeof(struct gkyl_bc_sheath_gyrokinetic_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_bc_sheath_gyrokinetic_kernels));
    gkyl_bc_gksheath_choose_reflectedf_kernel_cu(basis, edge, up->kernels_cu);
  } else {
    up->kernels->reflectedf = bc_gksheath_choose_reflectedf_kernel(basis, edge);
    assert(up->kernels->reflectedf);
    up->kernels_cu = up->kernels;
  }
#else
  up->kernels->reflectedf = bc_gksheath_choose_reflectedf_kernel(basis, edge);
  assert(up->kernels->reflectedf);
  up->kernels_cu = up->kernels;
#endif

  return up;
}

/* Modeled after gkyl_array_flip_copy_to_buffer_fn */
void
gkyl_bc_sheath_gyrokinetic_advance(const struct gkyl_bc_sheath_gyrokinetic *up, const struct gkyl_array *phi,
  const struct gkyl_array *phi_wall, struct gkyl_array *distf, const struct gkyl_range *conf_r)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_bc_sheath_gyrokinetic_advance_cu(up, phi, phi_wall, distf, conf_r);
    return;
  }
#endif

  int fidx[GKYL_MAX_DIM]; // Flipped index.
  int cidx[3];
  double xc[GKYL_MAX_DIM];

  int vpar_dir = up->cdim;
  double dvpar = up->grid->dx[vpar_dir];
  double dvparD2 = dvpar*0.5;
  int uplo = up->skin_r->upper[vpar_dir]+up->skin_r->lower[vpar_dir];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, up->skin_r);
  while (gkyl_range_iter_next(&iter)) {

    gkyl_copy_int_arr(up->skin_r->ndim, iter.idx, fidx);
    fidx[vpar_dir] = uplo - iter.idx[vpar_dir];
    // Turn this skin fidx into a ghost fidx.
    fidx[up->dir] = up->ghost_r->lower[up->dir];

    gkyl_rect_grid_cell_center(up->grid, iter.idx, xc);
    double vpar_c = xc[vpar_dir];
    double vparAbsSq_lo = vpar_c > 0.? pow(vpar_c-dvparD2,2) : pow(vpar_c+dvparD2,2);
    double vparAbsSq_up = vpar_c > 0.? pow(vpar_c+dvparD2,2) : pow(vpar_c-dvparD2,2);

    long skin_loc = gkyl_range_idx(up->skin_r, iter.idx);
    long ghost_loc = gkyl_range_idx(up->ghost_r, fidx);

    const double *inp = (const double*) gkyl_array_cfetch(distf, skin_loc);
    double *out = (double*) gkyl_array_fetch(distf, ghost_loc);

    for (int d=0; d<up->cdim; d++) cidx[d] = iter.idx[d]; 
    long conf_loc = gkyl_range_idx(conf_r, cidx);
    const double *phi_p = (const double*) gkyl_array_cfetch(phi, conf_loc);
    const double *phi_wall_p = (const double*) gkyl_array_cfetch(phi_wall, conf_loc);

    // Calculate reflected distribution function fhat.
    // note: reflected distribution can be
    // 1) fhat=0 (no reflection, i.e. absorb),
    // 2) fhat=f (full reflection)
    // 3) fhat=c*f (partial reflection)
    double fhat[up->basis->num_basis];
    up->kernels->reflectedf(vpar_c, dvpar, vparAbsSq_lo, vparAbsSq_up, up->q2Dm, phi_p, phi_wall_p, inp, fhat);

    // Reflect fhat into skin cells.
    bc_gksheath_reflect(up->dir, up->basis, up->cdim, out, fhat);
  }
}

void gkyl_bc_sheath_gyrokinetic_release(struct gkyl_bc_sheath_gyrokinetic *up)
{
  // Release memory associated with this updater.
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->kernels_cu);
  }
#endif
  gkyl_free(up->kernels);
  gkyl_free(up);
}
