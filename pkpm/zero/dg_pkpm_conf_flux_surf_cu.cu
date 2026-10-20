/* -*- c++ -*- */

#include <math.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_pkpm_conf_flux_surf.h>
#include <gkyl_dg_pkpm_conf_flux_surf_priv.h>
#include <gkyl_util.h>
}

// One thread per interior phase-space cell: the cell's lower faces in every
// configuration direction, plus the upper ghost-owned face at the domain edge
// (written by exactly one interior cell per direction, so no atomics are needed).
__global__ void
gkyl_dg_pkpm_conf_flux_surf_advance_cu_kernel(
  struct gkyl_dg_pkpm_conf_flux_surf *up, struct gkyl_range conf_range,
  struct gkyl_range phase_range, struct gkyl_range phase_range_ext,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *theta, const struct gkyl_array *fin, struct gkyl_array *cflrate,
  struct gkyl_array *conf_flux_surf
)
{
  int pdim = up->pdim;
  int cdim = up->cdim;
  int nflux = 8 * up->num_surf_basis;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM], xc_r[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < phase_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    long pidx = gkyl_range_idx(&phase_range, idx);
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);

    const double *b_c = (const double *)gkyl_array_cfetch(bvar, cidx);
    const double *u_c = (const double *)gkyl_array_cfetch(u, cidx);
    const double *xi_c = (const double *)gkyl_array_cfetch(xi, cidx);
    const double *theta_c = (const double *)gkyl_array_cfetch(theta, cidx);
    const double *jac_c = (const double *)gkyl_array_cfetch(up->jacob_pos, cidx);
    const double *f_c = (const double *)gkyl_array_cfetch(fin, pidx);
    double *cflrate_d = (double *)gkyl_array_fetch(cflrate, pidx);
    double *flux = (double *)gkyl_array_fetch(conf_flux_surf, pidx);

    for (int dir = 0; dir < cdim; ++dir) {
      gkyl_copy_int_arr(pdim, idx, idx_l);
      idx_l[dir] = idx_l[dir] - 1;
      long cidx_l = gkyl_range_idx(&conf_range, idx_l);
      long pidx_l = gkyl_range_idx(&phase_range, idx_l);
      cflrate_d[0] += up->conf_flux[dir](
        xc, up->phase_grid.dx, (const double *)gkyl_array_cfetch(up->jacob_pos, cidx_l), jac_c,
        (const double *)gkyl_array_cfetch(bvar, cidx_l), b_c,
        (const double *)gkyl_array_cfetch(u, cidx_l), u_c,
        (const double *)gkyl_array_cfetch(xi, cidx_l), xi_c,
        (const double *)gkyl_array_cfetch(theta, cidx_l), theta_c,
        (const double *)gkyl_array_cfetch(fin, pidx_l), f_c, &flux[dir * nflux]
      );

      if (idx[dir] == phase_range.upper[dir]) {
        gkyl_copy_int_arr(pdim, idx, idx_r);
        idx_r[dir] = idx_r[dir] + 1;
        long cidx_r = gkyl_range_idx(&conf_range, idx_r);
        long pidx_r = gkyl_range_idx(&phase_range_ext, idx_r);
        gkyl_rect_grid_cell_center(&up->phase_grid, idx_r, xc_r);
        double *flux_r = (double *)gkyl_array_fetch(conf_flux_surf, pidx_r);
        cflrate_d[0] += up->conf_flux[dir](
          xc_r, up->phase_grid.dx, jac_c, (const double *)gkyl_array_cfetch(up->jacob_pos, cidx_r),
          b_c, (const double *)gkyl_array_cfetch(bvar, cidx_r), u_c,
          (const double *)gkyl_array_cfetch(u, cidx_r), xi_c,
          (const double *)gkyl_array_cfetch(xi, cidx_r), theta_c,
          (const double *)gkyl_array_cfetch(theta, cidx_r), f_c,
          (const double *)gkyl_array_cfetch(fin, pidx_r), &flux_r[dir * nflux]
        );
      }
    }
  }
}

void
gkyl_dg_pkpm_conf_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_conf_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_range *phase_range_ext,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *theta, const struct gkyl_array *fin, struct gkyl_array *cflrate,
  struct gkyl_array *conf_flux_surf
)
{
  gkyl_dg_pkpm_conf_flux_surf_advance_cu_kernel<<<phase_range->nblocks, phase_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, *phase_range_ext, bvar->on_dev, u->on_dev,
    xi->on_dev, theta->on_dev, fin->on_dev, cflrate->on_dev, conf_flux_surf->on_dev
  );
}

__global__ static void
dg_pkpm_conf_flux_surf_set_cu_dev_ptrs(
  struct gkyl_dg_pkpm_conf_flux_surf *up, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  for (int d = 0; d < cdim; ++d) {
    up->conf_flux[d] = choose_pkpm_conf_flux_kern(b_type, cdim, poly_order, d);
  }
}

struct gkyl_dg_pkpm_conf_flux_surf *
gkyl_dg_pkpm_conf_flux_surf_cu_dev_inew(const struct gkyl_dg_pkpm_conf_flux_surf_inp *inp)
{
  struct gkyl_dg_pkpm_conf_flux_surf *up =
    (struct gkyl_dg_pkpm_conf_flux_surf *)gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim;
  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->pdim = pdim;
  up->use_gpu = inp->use_gpu;
  up->num_surf_basis = inp->phase_basis->num_basis / 2;
  // The device copy needs the device pointer of the Jacobian array; the host object keeps the
  // acquired map for its own lifetime management.
  assert(inp->pos_map);
  up->pos_map = 0;
  up->jacob_pos = inp->pos_map->jacob_pos->on_dev;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_pkpm_conf_flux_surf *up_cu =
    (struct gkyl_dg_pkpm_conf_flux_surf *)gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(*up), GKYL_CU_MEMCPY_H2D);

  dg_pkpm_conf_flux_surf_set_cu_dev_ptrs<<<1, 1>>>(
    up_cu, inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order
  );

  up->on_dev = up_cu;
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;
  return up;
}
