#include <gkyl_positivity_shift_gyrokinetic.h>
#include <gkyl_positivity_shift_gyrokinetic_priv.h>
#include <gkyl_alloc.h>
#include <gkyl_array_ops.h>
#include <float.h>

struct gkyl_positivity_shift_gyrokinetic*
gkyl_positivity_shift_gyrokinetic_new(struct gkyl_basis pbasis, struct gkyl_rect_grid grid,
  double mass, const struct gk_geometry *gk_geom, const struct gkyl_velocity_map *vel_map, bool use_gpu)
{
  // Allocate space for new updater.
  struct gkyl_positivity_shift_gyrokinetic *up = gkyl_malloc(sizeof(*up));

  up->ffloor_fac = 0.0; // ffloor will be set to max(f)*ffloor_fac.
  up->grid = grid;
  up->num_basis = pbasis.num_basis;
  up->mass = mass;
  up->gk_geom = gkyl_gk_geometry_acquire(gk_geom);
  up->vel_map = gkyl_velocity_map_acquire(vel_map);
  up->use_gpu = use_gpu;
  up->cellav_fac = 1./pow(sqrt(2.),pbasis.ndim);

  if (!use_gpu) {
    up->kernels = gkyl_malloc(sizeof(struct gkyl_positivity_shift_gyrokinetic_kernels));

    up->ffloor = gkyl_malloc(sizeof(double[1]));
    up->ffloor[0] = 0.0;  // Gets updated after 1st call to _advance.
  }
#ifdef GKYL_HAVE_CUDA
  if (use_gpu) {
    up->kernels = gkyl_cu_malloc(sizeof(struct gkyl_positivity_shift_gyrokinetic_kernels));

    up->ffloor = gkyl_cu_malloc(sizeof(double[1]));
    double ffloor_zero[] = {0.};  // Gets updated after 1st call to _advance.
    gkyl_cu_memcpy(up->ffloor, ffloor_zero, sizeof(double[1]), GKYL_CU_MEMCPY_H2D);
  }
#endif

  // Choose kernels that shift f and compute int moms of Deltaf.
  pos_shift_gk_choose_shift_kernel(up->kernels, pbasis, use_gpu);

  return up;
}

void
gkyl_positivity_shift_gyrokinetic_advance(gkyl_positivity_shift_gyrokinetic* up,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  struct gkyl_array *GKYL_RESTRICT distf, struct gkyl_array *GKYL_RESTRICT mom)
{
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_positivity_shift_gyrokinetic_advance_cu(up, phase_rng, conf_rng, distf, mom);
    return;
  }
#endif

  gkyl_array_clear_range(mom, 0.0, conf_rng);

  int pidx[GKYL_MAX_DIM];
  double distf_max = -DBL_MAX;
  int pdim = up->grid.ndim;
  int cdim = up->gk_geom->grid.ndim;

  struct gkyl_range_iter phase_iter;
  gkyl_range_iter_init(&phase_iter, phase_rng);
  while (gkyl_range_iter_next(&phase_iter)) {
    long plinidx = gkyl_range_idx(phase_rng, phase_iter.idx);

    double *distf_d = gkyl_array_fetch(distf, plinidx);

    double Deltaf[up->num_basis];
    bool shiftedf = up->kernels->shift(up->ffloor[0], distf_d, Deltaf);

    distf_max = fmax(distf_max, distf_d[0]); 

    if (shiftedf) {
      int idx_vel[2];
      for (int d=cdim; d<pdim; d++) idx_vel[d-cdim] = phase_iter.idx[d];

      long clinidx = gkyl_range_idx(conf_rng, phase_iter.idx);
      long vlinidx = gkyl_range_idx(&up->vel_map->local_vel, idx_vel);

      up->kernels->int_mom(up->grid.dx,
        gkyl_array_cfetch(up->vel_map->vmap, vlinidx), up->mass,
        gkyl_array_cfetch(up->gk_geom->bmag, clinidx),
        Deltaf, gkyl_array_fetch(mom, clinidx));
    }
  }

  up->ffloor[0] = up->ffloor_fac * distf_max * up->cellav_fac; 
}

void
gkyl_positivity_shift_gyrokinetic_release(gkyl_positivity_shift_gyrokinetic* up)
{
  // Release memory associated with this updater.
  gkyl_gk_geometry_release(up->gk_geom);
  gkyl_velocity_map_release(up->vel_map);
  if (!up->use_gpu) {
    gkyl_free(up->ffloor);
    gkyl_free(up->kernels);
  }
#ifdef GKYL_HAVE_CUDA
  if (up->use_gpu) {
    gkyl_cu_free(up->ffloor);
    gkyl_cu_free(up->kernels);
  }
#endif
  gkyl_free(up);
}
