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
  up->cdim = cdim;
  up->use_gpu = inp->use_gpu;
  // p=1: the surface basis of a face is the configuration basis with one direction removed.
  up->num_surf_basis = inp->conf_basis->num_basis / 2;
  up->wv_eqn = gkyl_wv_eqn_acquire(inp->wv_eqn);
  up->geom = gkyl_wave_geom_acquire(inp->geom);

  for (int d = 0; d < cdim; ++d) {
    up->fluid_flux[d] = choose_pkpm_fluid_flux_kern(
      inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order, d
    );
    up->fluid_flux_absorb_accum[d] = choose_pkpm_fluid_flux_absorb_accum_kern(
      inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order, d
    );
    up->fluid_flux_absorb[d] = choose_pkpm_fluid_flux_absorb_kern(
      inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order, d
    );
  }
  up->has_phase_grid = inp->phase_grid != 0;
  if (up->has_phase_grid) {
    up->phase_grid = *inp->phase_grid;
  }
  up->num_face_nodes = 1;
  for (int d = 1; d < cdim; ++d) {
    up->num_face_nodes *= 3;
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host

  return up;
}

void
gkyl_dg_pkpm_fluid_flux_surf_advance(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *conf_range_ext, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *fluid,
  struct gkyl_array *cflrate, struct gkyl_array *fluid_flux_surf
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(fluid_flux_surf)) {
    return gkyl_dg_pkpm_fluid_flux_surf_advance_cu(
      up, conf_range, conf_range_ext, vlasov_pkpm_moms, bvar, u, fluid, cflrate, fluid_flux_surf
    );
  }
#endif
  int cdim = up->cdim;
  int nflux = 5 * up->num_surf_basis; // per direction: [F_x, F_y, F_z, F_rho, F_E]
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM], xc_r[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(cdim, iter.idx, idx);
    long cidx = gkyl_range_idx(conf_range, idx);
    gkyl_rect_grid_cell_center(&up->conf_grid, idx, xc);
    const struct gkyl_wave_cell_geom *geom_c = gkyl_wave_geom_get(up->geom, idx);

    const double *moms_c = gkyl_array_cfetch(vlasov_pkpm_moms, cidx);
    const double *b_c = gkyl_array_cfetch(bvar, cidx);
    const double *u_c = gkyl_array_cfetch(u, cidx);
    const double *fluid_c = gkyl_array_cfetch(fluid, cidx);
    double *cflrate_d = gkyl_array_fetch(cflrate, cidx);
    double *flux = gkyl_array_fetch(fluid_flux_surf, cidx);

    for (int dir = 0; dir < cdim; ++dir) {
      // Lower face: the left neighbor (interior or ghost) is the l side.
      gkyl_copy_int_arr(cdim, iter.idx, idx_l);
      idx_l[dir] = idx_l[dir] - 1;
      long cidx_l = gkyl_range_idx(conf_range, idx_l);
      cflrate_d[0] += up->fluid_flux[dir](
        up->wv_eqn, geom_c, xc, up->conf_grid.dx, gkyl_array_cfetch(vlasov_pkpm_moms, cidx_l), moms_c,
        gkyl_array_cfetch(bvar, cidx_l), b_c, gkyl_array_cfetch(u, cidx_l), u_c,
        gkyl_array_cfetch(fluid, cidx_l), fluid_c, &flux[dir * nflux]
      );

      // Upper domain boundary: the face owned by the upper ghost cell, with the
      // current cell as the l side and the ghost as the r side.
      if (idx[dir] == conf_range->upper[dir]) {
        gkyl_copy_int_arr(cdim, iter.idx, idx_r);
        idx_r[dir] = idx_r[dir] + 1;
        long cidx_r = gkyl_range_idx(conf_range_ext, idx_r);
        gkyl_rect_grid_cell_center(&up->conf_grid, idx_r, xc_r);
        const struct gkyl_wave_cell_geom *geom_r = gkyl_wave_geom_get(up->geom, idx_r);
        double *flux_r = gkyl_array_fetch(fluid_flux_surf, cidx_r);
        cflrate_d[0] += up->fluid_flux[dir](
          up->wv_eqn, geom_r, xc_r, up->conf_grid.dx, moms_c, gkyl_array_cfetch(vlasov_pkpm_moms, cidx_r),
          b_c, gkyl_array_cfetch(bvar, cidx_r), u_c, gkyl_array_cfetch(u, cidx_r), fluid_c,
          gkyl_array_cfetch(fluid, cidx_r), &flux_r[dir * nflux]
        );
      }
    }
  }
}

void
gkyl_dg_pkpm_fluid_flux_surf_absorb(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, int dir, enum gkyl_edge_loc edge,
  const struct gkyl_range *conf_skin, const struct gkyl_range *conf_range_ext,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *fluid, const struct gkyl_array *f, struct gkyl_array *fluid_flux_surf
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(fluid_flux_surf)) {
    return gkyl_dg_pkpm_fluid_flux_surf_absorb_cu(
      up, dir, edge, conf_skin, conf_range_ext, phase_range, phase_range_ext, bvar, u, xi, fluid, f,
      fluid_flux_surf
    );
  }
#endif
  assert(up->has_phase_grid);
  int cdim = up->cdim;
  int nflux = 5 * up->num_surf_basis;
  int edge_sign = edge == GKYL_LOWER_EDGE ? -1 : 1;
  int idx_f[GKYL_MAX_DIM];
  double I[4 * 9];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_skin);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(conf_range_ext, iter.idx);
    const double *b_c = gkyl_array_cfetch(bvar, cidx);
    const double *u_c = gkyl_array_cfetch(u, cidx);
    const double *xi_c = gkyl_array_cfetch(xi, cidx);
    pkpm_fluid_flux_absorb_moments(
      up, dir, edge_sign, iter.idx, phase_range, phase_range_ext, b_c, u_c, xi_c, f, I
    );
    // The face is the lower face of the skin cell (lower boundary) or of the ghost cell
    // beyond it (upper boundary).
    gkyl_copy_int_arr(cdim, iter.idx, idx_f);
    if (edge == GKYL_UPPER_EDGE) {
      idx_f[dir] += 1;
    }
    double *flux = gkyl_array_fetch(fluid_flux_surf, gkyl_range_idx(conf_range_ext, idx_f));
    up->fluid_flux_absorb[dir](
      edge_sign, b_c, u_c, xi_c, gkyl_array_cfetch(fluid, cidx), I, &flux[dir * nflux]
    );
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
  gkyl_wv_eqn_release(up->wv_eqn);
  gkyl_wave_geom_release(up->geom);
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
  const struct gkyl_range *conf_range_ext, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *fluid,
  struct gkyl_array *cflrate, struct gkyl_array *fluid_flux_surf
)
{
  assert(false);
}

void
gkyl_dg_pkpm_fluid_flux_surf_absorb_cu(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, int dir, enum gkyl_edge_loc edge,
  const struct gkyl_range *conf_skin, const struct gkyl_range *conf_range_ext,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *fluid, const struct gkyl_array *f, struct gkyl_array *fluid_flux_surf
)
{
  assert(false);
}

#endif
