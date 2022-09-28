#include <gkyl_sheath_rarefaction_pot.h>
#include <gkyl_sheath_rarefaction_pot_priv.h>
#include <gkyl_alloc.h>
#include <assert.h>

struct gkyl_sheath_rarefaction_pot*
gkyl_sheath_rarefaction_pot_new(enum gkyl_edge_loc edge, const struct gkyl_range *local_range_ext,
  const int *num_ghosts, const struct gkyl_basis *basis, const struct gkyl_rect_grid *grid,
  double elem_charge, double mass_e, double mass_i, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_sheath_rarefaction_pot *up = gkyl_malloc(sizeof(struct gkyl_sheath_rarefaction_pot));

  up->edge = edge;
  up->use_gpu = use_gpu;
  up->basis = basis;
  up->grid = grid;
  up->elem_charge = elem_charge;
  up->mass_e = mass_e;
  up->mass_i = mass_i;

  int ndim = local_range_ext->ndim;
  assert((ndim==1) || (ndim==3)); // Only meant for 1x and 3x quasineutral sims.

  int dir = ndim-1;

  // Create the skin/ghost ranges.
  struct gkyl_range ghost_r;
  gkyl_skin_ghost_ranges(&up->skin_r, &ghost_r, dir, edge, local_range_ext, num_ghosts);

  // Choose the kernel modifies the potential at the boundary.
  up->kernels = gkyl_malloc(sizeof(struct gkyl_sheath_rarefaction_pot_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_sheath_rarefaction_pot_kernels));
    gkyl_sheath_rarepot_choose_phimod_kernel_cu(basis, edge, up->kernels_cu);
  } else {
    up->kernels->phimod = sheath_rarepot_choose_phimod_kernel(basis, edge);
    assert(up->kernels->phimod);
    up->kernels_cu = up->kernels;
  }
#else
  up->kernels->phimod = sheath_rarepot_choose_phimod_kernel(basis, edge);
  assert(up->kernels->phimod);
  up->kernels_cu = up->kernels;
#endif

  return up;
}

void
gkyl_sheath_rarefaction_pot_advance(const struct gkyl_sheath_rarefaction_pot *up,
  const struct gkyl_array *moms_e, const struct gkyl_array *m2par_e,
  const struct gkyl_array *moms_i, const struct gkyl_array *m2par_i,
  const struct gkyl_array *phi_wall, struct gkyl_array *phi)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_sheath_rarefaction_pot_advance_cu(up, phi, phi_wall, distf);
    return;
  }
#endif

    printf("\n");
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, &up->skin_r);
  while (gkyl_range_iter_next(&iter)) {
    long loc = gkyl_range_idx(&up->skin_r, iter.idx);

    const double *momse_p = (const double*) gkyl_array_cfetch(moms_e, loc);
    const double *m2pare_p = (const double*) gkyl_array_cfetch(m2par_e, loc);
    const double *momsi_p = (const double*) gkyl_array_cfetch(moms_i, loc);
    const double *m2pari_p = (const double*) gkyl_array_cfetch(m2par_i, loc);
    const double *phiwall_p = (const double*) gkyl_array_cfetch(phi_wall, loc);
    double *phi_p = (double*) gkyl_array_fetch(phi, loc);

    // Modify the potential at the boundary to account for rarefaction wave.
    up->kernels->phimod(up->elem_charge, up->mass_e, momse_p, m2pare_p, up->mass_i, momsi_p, m2pari_p, phiwall_p, phi_p);
  }
}

void gkyl_sheath_rarefaction_pot_release(struct gkyl_sheath_rarefaction_pot *up)
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
