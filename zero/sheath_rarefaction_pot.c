#include <gkyl_sheath_rarefaction_pot.h>
#include <gkyl_sheath_rarefaction_pot_priv.h>
#include <gkyl_alloc.h>
#include <assert.h>

struct gkyl_sheath_rarefaction_pot*
gkyl_sheath_rarefaction_pot_new(enum gkyl_edge_loc edge, const struct gkyl_basis *basis,
  double elem_charge, double mass_e, double mass_i, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_sheath_rarefaction_pot *up = gkyl_malloc(sizeof(struct gkyl_sheath_rarefaction_pot));

  up->elem_charge = elem_charge;
  up->mass_e = mass_e;
  up->mass_i = mass_i;
  up->use_gpu = use_gpu;

  // Choose the kernel modifies the potential at the boundary.
  if (!use_gpu)
    up->kernels = gkyl_malloc(sizeof(struct gkyl_sheath_rarefaction_pot_kernels));

#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels = gkyl_cu_malloc(sizeof(struct gkyl_sheath_rarefaction_pot_kernels));
#endif

  sheath_rarepot_choose_phimod_kernel(basis, edge, up->kernels, use_gpu);

  return up;
}

void
gkyl_sheath_rarefaction_pot_advance(const struct gkyl_sheath_rarefaction_pot *up,
  const struct gkyl_range *skin_range, const struct gkyl_range *surf_range,
  const struct gkyl_array *moms_e, const struct gkyl_array *moms_i,
  const struct gkyl_array *phi_wall, struct gkyl_array *phi)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_sheath_rarefaction_pot_advance_cu(up, skin_range, surf_range, moms_e, moms_i, phi_wall, phi);
    return;
  }
#endif

  int idx_surf[GKYL_MAX_CDIM];

  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, skin_range);
  while (gkyl_range_iter_next(&iter)) {
    idx_surf[0] = 1;
    for (int d=0; d<skin_range->ndim-1; d++)
      idx_surf[d] = iter.idx[d];

    long linidx_vol = gkyl_range_idx(skin_range, iter.idx);
    long linidx_surf = gkyl_range_idx(surf_range, idx_surf);

    const double *phiwall_p = (const double*) gkyl_array_cfetch(phi_wall, linidx_surf);

    const double *momse_p = (const double*) gkyl_array_cfetch(moms_e, linidx_vol);
    const double *momsi_p = (const double*) gkyl_array_cfetch(moms_i, linidx_vol);
    double *phi_p = (double*) gkyl_array_fetch(phi, linidx_vol);

    // Modify the potential at the boundary to account for rarefaction wave.
    up->kernels->phimod(up->elem_charge, up->mass_e, momse_p, up->mass_i, momsi_p, phiwall_p, phi_p);
  }
}

void gkyl_sheath_rarefaction_pot_release(struct gkyl_sheath_rarefaction_pot *up)
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
