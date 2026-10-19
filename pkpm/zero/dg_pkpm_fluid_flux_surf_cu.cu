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

// One thread per interior configuration cell: its lower faces in every direction, plus
// the upper ghost-owned face at the domain edge. The velocity cells are summed inside.
__global__ void
gkyl_dg_pkpm_fluid_flux_surf_advance_cu_kernel(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, struct gkyl_range conf_range, struct gkyl_range conf_range_ext,
  struct gkyl_range phase_range, struct gkyl_range phase_range_ext, const struct gkyl_array *bvar,
  const struct gkyl_array *u, const struct gkyl_array *xi, const struct gkyl_array *theta,
  const struct gkyl_array *fluid, const struct gkyl_array *f, struct gkyl_array *cflrate,
  struct gkyl_array *fluid_flux_surf
)
{
  int cdim = up->cdim;
  int nflux = 5 * up->num_surf_basis;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM], idx_r[GKYL_MAX_DIM];
  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < conf_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc1, idx);
    long cidx = gkyl_range_idx(&conf_range_ext, idx);
    double *flux = (double *)gkyl_array_fetch(fluid_flux_surf, cidx);
    double *cflrate_d = (double *)gkyl_array_fetch(cflrate, cidx);
    for (int dir = 0; dir < cdim; ++dir) {
      gkyl_copy_int_arr(cdim, idx, idx_l);
      idx_l[dir] = idx_l[dir] - 1;
      cflrate_d[0] += pkpm_fluid_flux_face(
        up, dir, idx_l, idx, &conf_range_ext, &phase_range, &phase_range_ext, bvar, u, xi, theta, fluid, f,
        &flux[dir * nflux]
      );
      if (idx[dir] == conf_range.upper[dir]) {
        gkyl_copy_int_arr(cdim, idx, idx_r);
        idx_r[dir] = idx_r[dir] + 1;
        double *flux_r = (double *)gkyl_array_fetch(fluid_flux_surf, gkyl_range_idx(&conf_range_ext, idx_r));
        cflrate_d[0] += pkpm_fluid_flux_face(
          up, dir, idx, idx_r, &conf_range_ext, &phase_range, &phase_range_ext, bvar, u, xi, theta, fluid, f,
          &flux_r[dir * nflux]
        );
      }
    }
  }
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
  gkyl_dg_pkpm_fluid_flux_surf_advance_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, *conf_range_ext, *phase_range, *phase_range_ext, bvar->on_dev, u->on_dev,
    xi->on_dev, theta->on_dev, fluid->on_dev, f->on_dev, cflrate->on_dev, fluid_flux_surf->on_dev
  );
}

__global__ static void
dg_pkpm_fluid_flux_surf_set_cu_dev_ptrs(
  struct gkyl_dg_pkpm_fluid_flux_surf *up, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  for (int d = 0; d < cdim; ++d) {
    up->fluid_flux_accum[d] = choose_pkpm_fluid_flux_accum_kern(b_type, cdim, poly_order, d);
    up->fluid_flux[d] = choose_pkpm_fluid_flux_kern(b_type, cdim, poly_order, d);
  }
}

struct gkyl_dg_pkpm_fluid_flux_surf *
gkyl_dg_pkpm_fluid_flux_surf_cu_dev_inew(const struct gkyl_dg_pkpm_fluid_flux_surf_inp *inp)
{
  struct gkyl_dg_pkpm_fluid_flux_surf *up =
    (struct gkyl_dg_pkpm_fluid_flux_surf *)gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim;
  up->conf_grid = *inp->conf_grid;
  assert(inp->phase_grid);
  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->use_gpu = inp->use_gpu;
  up->num_surf_basis = inp->conf_basis->num_basis / 2;
  up->num_face_nodes = 1;
  for (int d = 1; d < cdim; ++d) {
    up->num_face_nodes *= 3;
  }
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
  // the host object keeps the host pointers for its own lifetime management
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;

  return up;
}
