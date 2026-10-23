#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_pkpm_conf_flux_surf.h>
#include <gkyl_dg_pkpm_conf_flux_surf_priv.h>
#include <gkyl_util.h>

gkyl_dg_pkpm_conf_flux_surf *
gkyl_dg_pkpm_conf_flux_surf_inew(const struct gkyl_dg_pkpm_conf_flux_surf_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_pkpm_conf_flux_surf_cu_dev_inew(inp);
  }
#endif
  struct gkyl_dg_pkpm_conf_flux_surf *up = gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim;
  assert(pdim == cdim + 1);
  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->pdim = pdim;
  up->use_gpu = inp->use_gpu;
  // p=1 configuration space: the surface basis of a configuration face is the
  // phase basis with one p=1 configuration direction removed.
  up->num_surf_basis = inp->phase_basis->num_basis / 2;
  up->num_face_coef = PKPM_FC_NUM;
  for (int d = 1; d < cdim; ++d) {
    up->num_face_coef *= 5; // the transverse nodes of the face
  }
  assert(inp->pos_map);
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;

  for (int d = 0; d < cdim; ++d) {
    up->face_coef[d] = choose_pkpm_conf_face_coef_kern(
      inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order, d
    );
    up->conf_flux[d] = choose_pkpm_conf_flux_kern(
      inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order, d
    );
  }

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host

  return up;
}

void
gkyl_dg_pkpm_conf_flux_surf_face_coef(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *theta, struct gkyl_array *face_coef
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(face_coef)) {
    return gkyl_dg_pkpm_conf_flux_surf_face_coef_cu(up, conf_range, bvar, u, xi, theta, face_coef);
  }
#endif
  int cdim = up->cdim, nfc = up->num_face_coef;
  int idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, conf_range);
  while (gkyl_range_iter_next(&iter)) {
    long cidx = gkyl_range_idx(conf_range, iter.idx);
    double *fc = gkyl_array_fetch(face_coef, cidx);
    for (int dir = 0; dir < cdim; ++dir) {
      // Lower face: the left neighbor (interior or ghost) is the l side.
      gkyl_copy_int_arr(cdim, iter.idx, idx_l);
      idx_l[dir] = idx_l[dir] - 1;
      long cidx_l = gkyl_range_idx(conf_range, idx_l);
      up->face_coef[dir](
        gkyl_array_cfetch(bvar, cidx_l), gkyl_array_cfetch(bvar, cidx), gkyl_array_cfetch(u, cidx_l),
        gkyl_array_cfetch(u, cidx), gkyl_array_cfetch(xi, cidx_l), gkyl_array_cfetch(xi, cidx),
        gkyl_array_cfetch(theta, cidx_l), gkyl_array_cfetch(theta, cidx), &fc[dir * nfc]
      );
      // Upper domain boundary: the face owned by the upper ghost cell.
      if (iter.idx[dir] == conf_range->upper[dir]) {
        gkyl_copy_int_arr(cdim, iter.idx, idx_r);
        idx_r[dir] = idx_r[dir] + 1;
        long cidx_r = gkyl_range_idx(conf_range, idx_r);
        double *fc_r = gkyl_array_fetch(face_coef, cidx_r);
        up->face_coef[dir](
          gkyl_array_cfetch(bvar, cidx), gkyl_array_cfetch(bvar, cidx_r), gkyl_array_cfetch(u, cidx),
          gkyl_array_cfetch(u, cidx_r), gkyl_array_cfetch(xi, cidx), gkyl_array_cfetch(xi, cidx_r),
          gkyl_array_cfetch(theta, cidx), gkyl_array_cfetch(theta, cidx_r), &fc_r[dir * nfc]
        );
      }
    }
  }
}

void
gkyl_dg_pkpm_conf_flux_surf_advance(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *face_coef, const struct gkyl_array *fin, struct gkyl_array *cflrate,
  struct gkyl_array *conf_flux_surf
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(conf_flux_surf)) {
    return gkyl_dg_pkpm_conf_flux_surf_advance_cu(
      up, conf_range, phase_range, phase_range_ext, face_coef, fin, cflrate, conf_flux_surf
    );
  }
#endif
  int pdim = up->pdim;
  int cdim = up->cdim;
  int nflux = 8 * up->num_surf_basis; // per direction: the eight blocks
  int nfc = up->num_face_coef;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM], xc_r[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idx);
    long cidx = gkyl_range_idx(conf_range, idx);
    long pidx = gkyl_range_idx(phase_range, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    const double *fc_c = gkyl_array_cfetch(face_coef, cidx);
    const double *jac_c = gkyl_array_cfetch(up->jacob_pos, cidx);
    const double *f_c = gkyl_array_cfetch(fin, pidx);
    double *cflrate_d = gkyl_array_fetch(cflrate, pidx);
    double *flux = gkyl_array_fetch(conf_flux_surf, pidx);

    for (int dir = 0; dir < cdim; ++dir) {
      // Lower face: the left neighbor (interior or ghost) is the l side.
      gkyl_copy_int_arr(pdim, iter.idx, idx_l);
      idx_l[dir] = idx_l[dir] - 1;
      long cidx_l = gkyl_range_idx(conf_range, idx_l);
      long pidx_l = gkyl_range_idx(phase_range, idx_l);
      cflrate_d[0] += up->conf_flux[dir](
        xc, up->phase_grid.dx, gkyl_array_cfetch(up->jacob_pos, cidx_l), jac_c, &fc_c[dir * nfc],
        gkyl_array_cfetch(fin, pidx_l), f_c, &flux[dir * nflux]
      );

      // Upper domain boundary: the face owned by the upper ghost cell, with the
      // current cell as the l side and the ghost as the r side (both traces are
      // genuine: the ghost holds b, u and f after boundary conditions).
      if (idx[dir] == phase_range->upper[dir]) {
        gkyl_copy_int_arr(pdim, iter.idx, idx_r);
        idx_r[dir] = idx_r[dir] + 1;
        long cidx_r = gkyl_range_idx(conf_range, idx_r);
        long pidx_r = gkyl_range_idx(phase_range_ext, idx_r);
        gkyl_rect_grid_cell_center(&up->phase_grid, idx_r, xc_r);
        double *flux_r = gkyl_array_fetch(conf_flux_surf, pidx_r);
        const double *fc_r = gkyl_array_cfetch(face_coef, cidx_r);
        // The ghost face bounds the current cell: its speed enters this cell's cflrate.
        cflrate_d[0] += up->conf_flux[dir](
          xc_r, up->phase_grid.dx, jac_c, gkyl_array_cfetch(up->jacob_pos, cidx_r), &fc_r[dir * nfc],
          f_c, gkyl_array_cfetch(fin, pidx_r), &flux_r[dir * nflux]
        );
      }
    }
  }
}

int
gkyl_dg_pkpm_conf_flux_surf_face_coef_ncomp(const struct gkyl_dg_pkpm_conf_flux_surf *up)
{
  return up->cdim * up->num_face_coef;
}

int
gkyl_dg_pkpm_conf_flux_surf_num_surf_basis(const struct gkyl_dg_pkpm_conf_flux_surf *up)
{
  return up->num_surf_basis;
}

void
gkyl_dg_pkpm_conf_flux_surf_release(struct gkyl_dg_pkpm_conf_flux_surf *up)
{
  gkyl_vlasov_position_map_release(up->pos_map);
  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_pkpm_conf_flux_surf *
gkyl_dg_pkpm_conf_flux_surf_cu_dev_inew(const struct gkyl_dg_pkpm_conf_flux_surf_inp *inp)
{
  assert(false);
  return 0;
}

void
gkyl_dg_pkpm_conf_flux_surf_face_coef_cu(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *theta, struct gkyl_array *face_coef
)
{
  assert(false);
}

void
gkyl_dg_pkpm_conf_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *face_coef, const struct gkyl_array *fin, struct gkyl_array *cflrate,
  struct gkyl_array *conf_flux_surf
)
{
  assert(false);
}

#endif
