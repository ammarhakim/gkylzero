#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_pkpm_fluid_flux_surf.h>
#include <gkyl_dg_pkpm_fluid_flux_surf_priv.h>
#include <gkyl_util.h>

gkyl_dg_pkpm_fluid_flux_surf *
gkyl_dg_pkpm_fluid_flux_surf_inew(const struct gkyl_dg_pkpm_fluid_flux_surf_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_pkpm_fluid_flux_surf_cu_dev_inew(inp);
  }
#endif
  struct gkyl_dg_pkpm_fluid_flux_surf *up = gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim;
  up->conf_grid = *inp->conf_grid;
  assert(inp->phase_grid);
  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->use_gpu = inp->use_gpu;
  // p=1: the surface basis of a face is the configuration basis with one direction removed.
  up->num_surf_basis = inp->conf_basis->num_basis / 2;
  up->num_face_nodes = 1;
  for (int d = 1; d < cdim; ++d) {
    up->num_face_nodes *= 5; // the transverse nodes of the kinetic face flux
  }
  assert(inp->pos_map);
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;

  for (int d = 0; d < cdim; ++d) {
    up->fluid_flux_accum[d] = choose_pkpm_fluid_flux_accum_kern(
      inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order, d
    );
    up->fluid_flux[d] = choose_pkpm_fluid_flux_kern(
      inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order, d
    );
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host

  return up;
}

void
gkyl_dg_pkpm_fluid_flux_surf_advance(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *conf_range_ext, const struct gkyl_range *phase_range,
  const struct gkyl_range *phase_range_ext, const struct gkyl_array *bvar, const struct gkyl_array *u,
  const struct gkyl_array *xi, const struct gkyl_array *theta, const struct gkyl_array *fluid,
  const struct gkyl_array *f, struct gkyl_array *cflrate, struct gkyl_array *fluid_flux_surf
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(fluid_flux_surf)) {
    return gkyl_dg_pkpm_fluid_flux_surf_advance_cu(
      up, conf_range, conf_range_ext, phase_range, phase_range_ext, bvar, u, xi, theta, fluid, f, cflrate, fluid_flux_surf
    );
  }
#endif
  int cdim = up->cdim;
  int nflux = 5 * up->num_surf_basis; // per direction: [F_rho, F_x, F_y, F_z, F_E]
  int idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);

  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(conf_range_ext, iter.idx);
    double *flux = gkyl_array_fetch(fluid_flux_surf, cidx);
    double *cflrate_d = gkyl_array_fetch(cflrate, cidx);

    for (int dir = 0; dir < cdim; ++dir) {
      // Lower face: the left neighbor (interior or ghost) is the l side.
      gkyl_copy_int_arr(cdim, iter.idx, idx_l);
      idx_l[dir] = idx_l[dir] - 1;
      cflrate_d[0] += pkpm_fluid_flux_face(
        up, dir, idx_l, iter.idx, conf_range_ext, phase_range, phase_range_ext, bvar, u, xi, theta, fluid, f,
        &flux[dir * nflux]
      );

      // Upper domain boundary: the face owned by the upper ghost cell, with the
      // current cell as the l side and the ghost as the r side.
      if (iter.idx[dir] == conf_range->upper[dir]) {
        gkyl_copy_int_arr(cdim, iter.idx, idx_r);
        idx_r[dir] = idx_r[dir] + 1;
        double *flux_r = gkyl_array_fetch(fluid_flux_surf, gkyl_range_idx(conf_range_ext, idx_r));
        cflrate_d[0] += pkpm_fluid_flux_face(
          up, dir, iter.idx, idx_r, conf_range_ext, phase_range, phase_range_ext, bvar, u, xi, theta, fluid, f,
          &flux_r[dir * nflux]
        );
      }
    }
  }
}

int
gkyl_dg_pkpm_fluid_flux_surf_num_surf_basis(const struct gkyl_dg_pkpm_fluid_flux_surf *up)
{
  return up->num_surf_basis;
}

void
gkyl_dg_pkpm_fluid_flux_surf_release(struct gkyl_dg_pkpm_fluid_flux_surf *up)
{
  gkyl_vlasov_position_map_release(up->pos_map);
  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_pkpm_fluid_flux_surf *
gkyl_dg_pkpm_fluid_flux_surf_cu_dev_inew(const struct gkyl_dg_pkpm_fluid_flux_surf_inp *inp)
{
  assert(false);
  return 0;
}

void
gkyl_dg_pkpm_fluid_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *conf_range_ext, const struct gkyl_range *phase_range,
  const struct gkyl_range *phase_range_ext, const struct gkyl_array *bvar, const struct gkyl_array *u,
  const struct gkyl_array *xi, const struct gkyl_array *theta, const struct gkyl_array *fluid,
  const struct gkyl_array *f, struct gkyl_array *cflrate, struct gkyl_array *fluid_flux_surf
)
{
  assert(false);
}

#endif
