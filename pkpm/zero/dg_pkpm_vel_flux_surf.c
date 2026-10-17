#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_pkpm_vel_flux_surf.h>
#include <gkyl_dg_pkpm_vel_flux_surf_priv.h>
#include <gkyl_util.h>

gkyl_dg_pkpm_vel_flux_surf *
gkyl_dg_pkpm_vel_flux_surf_inew(const struct gkyl_dg_pkpm_vel_flux_surf_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_pkpm_vel_flux_surf_cu_dev_inew(inp);
  }
#endif
  struct gkyl_dg_pkpm_vel_flux_surf *up = gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim;
  assert(pdim == cdim + 1);
  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->pdim = pdim;
  up->use_gpu = inp->use_gpu;
  // The surface basis of the v_par face is the p=1 configuration basis.
  up->num_surf_basis = inp->conf_basis->num_basis;
  up->vel_flux =
    choose_pkpm_vel_flux_kern(inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order);

  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host

  return up;
}

void
gkyl_dg_pkpm_vel_flux_surf_advance(
  struct gkyl_dg_pkpm_vel_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *bvar, const struct gkyl_array *grad_b,
  const struct gkyl_array *u, const struct gkyl_array *grad_u, const struct gkyl_array *theta, const struct gkyl_array *a_P,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *vel_flux_surf
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(vel_flux_surf)) {
    return gkyl_dg_pkpm_vel_flux_surf_advance_cu(
      up, conf_range, phase_range, bvar, grad_b, u, grad_u, theta, a_P, fin, cflrate, vel_flux_surf
    );
  }
#endif
  int pdim = up->pdim;
  int cdim = up->cdim;
  int nflux = 2 * up->num_surf_basis;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];
  struct gkyl_range_iter iter;
  gkyl_range_iter_init(&iter, phase_range);

  while (gkyl_range_iter_next(&iter)) {
    gkyl_copy_int_arr(pdim, iter.idx, idx);
    long cidx = gkyl_range_idx(conf_range, idx);
    long pidx = gkyl_range_idx(phase_range, idx);
    double *flux = gkyl_array_fetch(vel_flux_surf, pidx);

    // The lower face of the lowest cell is the velocity domain edge: zero total flux.
    if (idx[cdim] == phase_range->lower[cdim]) {
      for (int k = 0; k < nflux; ++k) {
        flux[k] = 0.0;
      }
      continue;
    }

    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);
    gkyl_copy_int_arr(pdim, iter.idx, idx_l);
    idx_l[cdim] = idx_l[cdim] - 1;
    long pidx_l = gkyl_range_idx(phase_range, idx_l);
    double *cflrate_d = gkyl_array_fetch(cflrate, pidx);
    cflrate_d[0] += up->vel_flux(
      xc, up->phase_grid.dx, gkyl_array_cfetch(bvar, cidx), gkyl_array_cfetch(grad_b, cidx),
      gkyl_array_cfetch(u, cidx), gkyl_array_cfetch(grad_u, cidx), gkyl_array_cfetch(theta, cidx), gkyl_array_cfetch(a_P, cidx),
      gkyl_array_cfetch(fin, pidx_l), gkyl_array_cfetch(fin, pidx), flux
    );
  }
}

int
gkyl_dg_pkpm_vel_flux_surf_num_surf_basis(const struct gkyl_dg_pkpm_vel_flux_surf *up)
{
  return up->num_surf_basis;
}

void
gkyl_dg_pkpm_vel_flux_surf_release(struct gkyl_dg_pkpm_vel_flux_surf *up)
{
  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_pkpm_vel_flux_surf *
gkyl_dg_pkpm_vel_flux_surf_cu_dev_inew(const struct gkyl_dg_pkpm_vel_flux_surf_inp *inp)
{
  assert(false);
  return 0;
}

void
gkyl_dg_pkpm_vel_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_vel_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *bvar, const struct gkyl_array *grad_b,
  const struct gkyl_array *u, const struct gkyl_array *grad_u, const struct gkyl_array *theta, const struct gkyl_array *a_P,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *vel_flux_surf
)
{
  assert(false);
}

#endif
