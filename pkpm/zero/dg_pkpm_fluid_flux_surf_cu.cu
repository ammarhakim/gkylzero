/* -*- c++ -*- */

#include <math.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_pkpm_fluid_flux_surf.h>
#include <gkyl_dg_pkpm_fluid_flux_surf_priv.h>
#include <gkyl_util.h>
}

// One thread per interior configuration cell: its lower faces in every direction,
// plus the upper ghost-owned face at the domain edge.
__global__ void
gkyl_dg_pkpm_fluid_flux_surf_advance_cu_kernel(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, struct gkyl_range conf_range,
  struct gkyl_range conf_range_ext, const struct gkyl_array *vlasov_pkpm_moms, const struct gkyl_array *bvar,
  const struct gkyl_array *u, const struct gkyl_array *fluid, struct gkyl_array *cflrate,
  struct gkyl_array *fluid_flux_surf
)
{
  int cdim = up->cdim;
  int nflux = 5 * up->num_surf_basis;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM], xc_r[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    gkyl_rect_grid_cell_center(&up->conf_grid, idx, xc);
    const struct gkyl_wave_cell_geom *geom_c = gkyl_wave_geom_get(up->geom, idx);

    const double *moms_c = (const double *)gkyl_array_cfetch(vlasov_pkpm_moms, cidx);
    const double *b_c = (const double *)gkyl_array_cfetch(bvar, cidx);
    const double *u_c = (const double *)gkyl_array_cfetch(u, cidx);
    const double *fluid_c = (const double *)gkyl_array_cfetch(fluid, cidx);
    double *cflrate_d = (double *)gkyl_array_fetch(cflrate, cidx);
    double *flux = (double *)gkyl_array_fetch(fluid_flux_surf, cidx);

    for (int dir = 0; dir < cdim; ++dir) {
      gkyl_copy_int_arr(cdim, idx, idx_l);
      idx_l[dir] = idx_l[dir] - 1;
      long cidx_l = gkyl_range_idx(&conf_range, idx_l);
      cflrate_d[0] += up->fluid_flux[dir](
        up->wv_eqn, geom_c, xc, up->conf_grid.dx,
        (const double *)gkyl_array_cfetch(up->jacob_pos, cidx_l), (const double *)gkyl_array_cfetch(up->jacob_pos, cidx),
        (const double *)gkyl_array_cfetch(vlasov_pkpm_moms, cidx_l), moms_c,
        (const double *)gkyl_array_cfetch(bvar, cidx_l), b_c,
        (const double *)gkyl_array_cfetch(u, cidx_l), u_c,
        (const double *)gkyl_array_cfetch(fluid, cidx_l), fluid_c, &flux[dir * nflux]
      );

      if (idx[dir] == conf_range.upper[dir]) {
        gkyl_copy_int_arr(cdim, idx, idx_r);
        idx_r[dir] = idx_r[dir] + 1;
        long cidx_r = gkyl_range_idx(&conf_range_ext, idx_r);
        gkyl_rect_grid_cell_center(&up->conf_grid, idx_r, xc_r);
        const struct gkyl_wave_cell_geom *geom_r = gkyl_wave_geom_get(up->geom, idx_r);
        double *flux_r = (double *)gkyl_array_fetch(fluid_flux_surf, cidx_r);
        cflrate_d[0] += up->fluid_flux[dir](
          up->wv_eqn, geom_r, xc_r, up->conf_grid.dx,
          (const double *)gkyl_array_cfetch(up->jacob_pos, cidx), (const double *)gkyl_array_cfetch(up->jacob_pos, cidx_r),
          moms_c,
          (const double *)gkyl_array_cfetch(vlasov_pkpm_moms, cidx_r), b_c,
          (const double *)gkyl_array_cfetch(bvar, cidx_r), u_c,
          (const double *)gkyl_array_cfetch(u, cidx_r), fluid_c,
          (const double *)gkyl_array_cfetch(fluid, cidx_r), &flux_r[dir * nflux]
        );
      }
    }
  }
}

void
gkyl_dg_pkpm_fluid_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *conf_range_ext, const struct gkyl_array *vlasov_pkpm_moms,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *fluid,
  struct gkyl_array *cflrate, struct gkyl_array *fluid_flux_surf
)
{
  gkyl_dg_pkpm_fluid_flux_surf_advance_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, *conf_range_ext, vlasov_pkpm_moms->on_dev, bvar->on_dev, u->on_dev,
    fluid->on_dev, cflrate->on_dev, fluid_flux_surf->on_dev
  );
}

__global__ static void
dg_pkpm_fluid_flux_surf_set_cu_dev_ptrs(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  for (int d = 0; d < cdim; ++d) {
    up->fluid_flux[d] = choose_pkpm_fluid_flux_kern(b_type, cdim, poly_order, d);
    up->fluid_flux_absorb_accum[d] = choose_pkpm_fluid_flux_absorb_accum_kern(b_type, cdim, poly_order, d);
    up->fluid_flux_absorb[d] = choose_pkpm_fluid_flux_absorb_kern(b_type, cdim, poly_order, d);
  }
}

struct gkyl_dg_pkpm_fluid_flux_surf *
gkyl_dg_pkpm_fluid_flux_surf_cu_dev_inew(const struct gkyl_dg_pkpm_fluid_flux_surf_inp *inp)
{
  struct gkyl_dg_pkpm_fluid_flux_surf *up =
    (struct gkyl_dg_pkpm_fluid_flux_surf *)gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim;
  up->conf_grid = *inp->conf_grid;
  up->cdim = cdim;
  up->use_gpu = inp->use_gpu;
  up->num_surf_basis = inp->conf_basis->num_basis / 2;
  up->has_phase_grid = inp->phase_grid != 0;
  if (up->has_phase_grid) {
    up->phase_grid = *inp->phase_grid;
  }
  up->num_face_nodes = 1;
  for (int d = 1; d < cdim; ++d) {
    up->num_face_nodes *= 3;
  }

  // The device copy needs the device pointers of the wave equation and geometry;
  // the host object keeps the host pointers for its own lifetime management.
  struct gkyl_wv_eqn *eqn = gkyl_wv_eqn_acquire(inp->wv_eqn);
  struct gkyl_wave_geom *geom = gkyl_wave_geom_acquire(inp->geom);
  up->wv_eqn = eqn->on_dev;
  up->geom = geom->on_dev;
  assert(inp->pos_map);
  up->pos_map = 0;
  up->jacob_pos = inp->pos_map->jacob_pos->on_dev;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_pkpm_fluid_flux_surf *up_cu =
    (struct gkyl_dg_pkpm_fluid_flux_surf *)gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(*up), GKYL_CU_MEMCPY_H2D);

  dg_pkpm_fluid_flux_surf_set_cu_dev_ptrs<<<1, 1>>>(
    up_cu, inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order
  );

  up->on_dev = up_cu;
  up->wv_eqn = eqn;
  up->geom = geom;
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;
  return up;
}

__global__ static void
gkyl_dg_pkpm_fluid_flux_surf_absorb_cu_kernel(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, int dir, int edge_sign, struct gkyl_range conf_skin,
  struct gkyl_range conf_range_ext, struct gkyl_range phase_range, struct gkyl_range phase_range_ext,
  const struct gkyl_array *bvar, const struct gkyl_array *u, const struct gkyl_array *xi,
  const struct gkyl_array *fluid, const struct gkyl_array *f, struct gkyl_array *fluid_flux_surf
)
{
  int cdim = up->cdim;
  int nflux = 5 * up->num_surf_basis;
  int idx[GKYL_MAX_DIM], idx_f[GKYL_MAX_DIM];
  double I[4 * 9];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_skin.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_skin, linc1, idx);
    long cidx = gkyl_range_idx(&conf_range_ext, idx);
    const double *b_c = (const double *)gkyl_array_cfetch(bvar, cidx);
    const double *u_c = (const double *)gkyl_array_cfetch(u, cidx);
    const double *xi_c = (const double *)gkyl_array_cfetch(xi, cidx);
    pkpm_fluid_flux_absorb_moments(up, dir, edge_sign, idx, &phase_range, &phase_range_ext, b_c, u_c, xi_c, f, I);
    gkyl_copy_int_arr(cdim, idx, idx_f);
    if (edge_sign > 0) {
      idx_f[dir] += 1;
    }
    double *flux = (double *)gkyl_array_fetch(fluid_flux_surf, gkyl_range_idx(&conf_range_ext, idx_f));
    up->fluid_flux_absorb[dir](
      edge_sign, (const double *)gkyl_array_cfetch(up->jacob_pos, cidx), b_c, u_c, xi_c,
      (const double *)gkyl_array_cfetch(fluid, cidx), I, &flux[dir * nflux]
    );
  }
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
  assert(up->has_phase_grid);
  int edge_sign = edge == GKYL_LOWER_EDGE ? -1 : 1;
  gkyl_dg_pkpm_fluid_flux_surf_absorb_cu_kernel<<<conf_skin->nblocks, conf_skin->nthreads>>>(
    up->on_dev, dir, edge_sign, *conf_skin, *conf_range_ext, *phase_range, *phase_range_ext, bvar->on_dev,
    u->on_dev, xi->on_dev, fluid->on_dev, f->on_dev, fluid_flux_surf->on_dev
  );
}
