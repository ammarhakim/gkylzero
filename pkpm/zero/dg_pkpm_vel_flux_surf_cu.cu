/* -*- c++ -*- */

#include <math.h>

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_array_ops.h>
#include <gkyl_array_ops_priv.h>
#include <gkyl_dg_pkpm_vel_flux_surf.h>
#include <gkyl_dg_pkpm_vel_flux_surf_priv.h>
#include <gkyl_util.h>
}

// One thread per phase-space cell: its lower v_par face (zero at the lower
// velocity domain edge).
__global__ void
gkyl_dg_pkpm_vel_flux_surf_advance_cu_kernel(
  struct gkyl_dg_pkpm_vel_flux_surf *up, struct gkyl_range conf_range,
  struct gkyl_range phase_range, const struct gkyl_array *bvar, const struct gkyl_array *grad_b,
  const struct gkyl_array *u, const struct gkyl_array *grad_u, const struct gkyl_array *xi, const struct gkyl_array *grad_xi,
  const struct gkyl_array *theta, const struct gkyl_array *alpha_0, const struct gkyl_array *c_h,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *vel_flux_surf
)
{
  int pdim = up->pdim;
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM], idx_l[GKYL_MAX_DIM];
  double xc[GKYL_MAX_DIM];

  for (unsigned long linc1 = threadIdx.x + blockIdx.x * blockDim.x; linc1 < phase_range.volume;
       linc1 += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&phase_range, linc1, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    long pidx = gkyl_range_idx(&phase_range, idx);
    if (idx[cdim] == phase_range.lower[cdim]) {
      continue;
    }
    gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);
    gkyl_copy_int_arr(pdim, idx, idx_l);
    idx_l[cdim] = idx_l[cdim] - 1;
    long pidx_l = gkyl_range_idx(&phase_range, idx_l);
    double *cflrate_d = (double *)gkyl_array_fetch(cflrate, pidx);
    cflrate_d[0] += up->vel_flux(
      xc, up->phase_grid.dx, (const double *)gkyl_array_cfetch(bvar, cidx),
      (const double *)gkyl_array_cfetch(grad_b, cidx),
      (const double *)gkyl_array_cfetch(u, cidx), (const double *)gkyl_array_cfetch(grad_u, cidx),
      (const double *)gkyl_array_cfetch(xi, cidx), (const double *)gkyl_array_cfetch(grad_xi, cidx),
      (const double *)gkyl_array_cfetch(theta, cidx), (const double *)gkyl_array_cfetch(alpha_0, cidx),
      (const double *)gkyl_array_cfetch(c_h, cidx), (const double *)gkyl_array_cfetch(fin, pidx_l),
      (const double *)gkyl_array_cfetch(fin, pidx), (double *)gkyl_array_fetch(vel_flux_surf, pidx)
    );
  }
}



void
gkyl_dg_pkpm_vel_flux_surf_advance_cu(
  struct gkyl_dg_pkpm_vel_flux_surf *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *bvar, const struct gkyl_array *grad_b,
  const struct gkyl_array *u, const struct gkyl_array *grad_u, const struct gkyl_array *xi, const struct gkyl_array *grad_xi,
  const struct gkyl_array *theta, const struct gkyl_array *alpha_0, const struct gkyl_array *c_h,
  const struct gkyl_array *fin, struct gkyl_array *cflrate, struct gkyl_array *vel_flux_surf
)
{
  gkyl_dg_pkpm_vel_flux_surf_advance_cu_kernel<<<phase_range->nblocks, phase_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, bvar->on_dev, grad_b->on_dev, u->on_dev, grad_u->on_dev,
    xi->on_dev, grad_xi->on_dev, theta->on_dev,
    alpha_0->on_dev, c_h->on_dev, fin->on_dev, cflrate->on_dev, vel_flux_surf->on_dev
  );
}



__global__ static void
dg_pkpm_vel_flux_surf_set_cu_dev_ptrs(
  struct gkyl_dg_pkpm_vel_flux_surf *up, enum gkyl_basis_type b_type, int cdim, int poly_order
)
{
  up->vel_flux = choose_pkpm_vel_flux_kern(b_type, cdim, poly_order);
}

struct gkyl_dg_pkpm_vel_flux_surf *
gkyl_dg_pkpm_vel_flux_surf_cu_dev_inew(const struct gkyl_dg_pkpm_vel_flux_surf_inp *inp)
{
  struct gkyl_dg_pkpm_vel_flux_surf *up =
    (struct gkyl_dg_pkpm_vel_flux_surf *)gkyl_malloc(sizeof(*up));

  int cdim = inp->conf_basis->ndim, pdim = inp->phase_basis->ndim;
  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->pdim = pdim;
  up->use_gpu = inp->use_gpu;
  up->num_surf_basis = inp->conf_basis->num_basis;

  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);

  struct gkyl_dg_pkpm_vel_flux_surf *up_cu =
    (struct gkyl_dg_pkpm_vel_flux_surf *)gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(*up), GKYL_CU_MEMCPY_H2D);

  dg_pkpm_vel_flux_surf_set_cu_dev_ptrs<<<1, 1>>>(
    up_cu, inp->conf_basis->b_type, cdim, inp->conf_basis->poly_order
  );

  up->on_dev = up_cu;
  return up;
}
