#include <gkyl_ambi_bolt_potential.h>
#include <gkyl_ambi_bolt_potential_priv.h>
#include <gkyl_alloc.h>

gkyl_ambi_bolt_potential*
gkyl_ambi_bolt_potential_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *basis,
  const struct gkyl_range *local_range_ext, const int *num_ghosts,
  double mass_e, double charge_e, double temp_e, bool use_gpu)
{
  struct gkyl_ambi_bolt_potential *up = gkyl_malloc(sizeof(struct gkyl_ambi_bolt_potential));

  up->ndim = basis->ndim;
  up->num_basis = basis->num_basis;
  up->use_gpu = use_gpu;

  int sheath_dir = up->ndim;

  up->dz = grid->dx[sheath_dir];
  up->mass_e = mass_e;
  up->charge_e = charge_e;
  up->temp_e = temp_e;

  // Create ranges over the skin cells along the field line.
  struct gkyl_range ghost_r;
  gkyl_skin_ghost_ranges(&up->loskin_r, &ghost_r, sheath_dir, GKYL_LOWER_EDGE, local_range_ext, num_ghosts);
  gkyl_skin_ghost_ranges(&up->upskin_r, &ghost_r, sheath_dir, GKYL_UPPER_EDGE, local_range_ext, num_ghosts);

  up->kernels = gkyl_malloc(sizeof(struct gkyl_ambi_bolt_potential_kernels));
#ifdef GKYL_HAVE_CUDA
  if (use_gpu)
    up->kernels_cu = gkyl_cu_malloc(sizeof(struct gkyl_ambi_bolt_potential_kernels));
  else
    up->kernels_cu = up->kernels;
#else
  up->kernels_cu = up->kernels;
#endif

  // Select sheath_calc and phi_calc kernels.
  ambi_bolt_potential_choose_kernels(basis, up->kernels);

  return up;
}

void
gkyl_ambi_bolt_potential_sheath_calc(struct gkyl_ambi_bolt_potential *up, enum gkyl_edge_loc edge,
  struct gkyl_array *jacob_geo_inv, struct gkyl_array *gammai, struct gkyl_array *m0i,
  struct gkyl_array *sheath_vals)
{
  unsigned int keridx = (edge == GKYL_LOWER_EDGE) ? 0 : 1;

  int sonly_idx[GKYL_MAX_CDIM]; // Index in skin-only array/grid.

  struct gkyl_range *skin_r = (edge == GKYL_LOWER_EDGE) ? &up->loskin_r : &up->upskin_r;
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, skin_r);
  while (gkyl_range_iter_next(&iter)) {

    long skin_loc = gkyl_range_idx(skin_r, iter.idx);

    gkyl_copy_int_arr(skin_r->ndim, iter.idx, sonly_idx);
    sonly_idx[up->ndim] = 0;

    long sonly_loc = gkyl_range_idx(skin_r, sonly_idx);

    const double *jacinv_p = (const double*) gkyl_array_cfetch(jacob_geo_inv, skin_loc);
    const double *gammai_p = (const double*) gkyl_array_cfetch(gammai, skin_loc);
    const double *m0i_p = (const double*) gkyl_array_cfetch(m0i, skin_loc);

    double *out_p = (double*) gkyl_array_cfetch(sheath_vals, sonly_loc);

    up->kernels->sheath_calc[keridx](up->dz, up->charge_e, up->mass_e, up->temp_e, jacinv_p, gammai_p, m0i_p, out_p);
  
  }
}

void
gkyl_ambi_bolt_potential_release(gkyl_ambi_bolt_potential *up)
{
  gkyl_free(up);
}
