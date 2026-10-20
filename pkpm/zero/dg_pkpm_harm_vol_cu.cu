/* -*- c++ -*- */

extern "C" {
#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_dg_pkpm_harm_vol_priv.h>
}

// One thread per configuration cell: its coefficients (local memory), then its velocity cells.
__global__ static void
gkyl_dg_pkpm_harm_vol_advance_cu_kernel(
  const struct gkyl_dg_pkpm_harm_vol *up, struct gkyl_range conf_range, struct gkyl_range phase_range,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *xi, const struct gkyl_array *grad_xi,
  const struct gkyl_array *grad_ln_n, const struct gkyl_array *theta, const struct gkyl_array *dtb,
  const struct gkyl_array *em, const struct gkyl_array *fin, struct gkyl_array *cflrate,
  struct gkyl_array *rhs
)
{
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM];
  struct pkpm_hm_coeffs hc;
  for (unsigned long linc = threadIdx.x + blockIdx.x * blockDim.x; linc < conf_range.volume;
       linc += gridDim.x * blockDim.x) {
    gkyl_sub_range_inv_idx(&conf_range, linc, idx);
    long cidx = gkyl_range_idx(&conf_range, idx);
    const double *bvar_d = (const double *)gkyl_array_cfetch(bvar, cidx);
    const double *grad_b_d = (const double *)gkyl_array_cfetch(grad_b, cidx);
    const double *u_d = (const double *)gkyl_array_cfetch(u, cidx);
    const double *grad_u_d = (const double *)gkyl_array_cfetch(grad_u, cidx);
    const double *xi_d = (const double *)gkyl_array_cfetch(xi, cidx);
    const double *grad_xi_d = (const double *)gkyl_array_cfetch(grad_xi, cidx);
    const double *gln_d = (const double *)gkyl_array_cfetch(grad_ln_n, cidx);
    const double *theta_d = (const double *)gkyl_array_cfetch(theta, cidx);
    const double *dtb_d = (const double *)gkyl_array_cfetch(dtb, cidx);
    const double *em_d = (const double *)gkyl_array_cfetch(em, cidx);
    const double *jac_d = (const double *)gkyl_array_cfetch(up->jacob_pos, cidx);
    pkpm_hm_build(up, bvar_d, grad_b_d, u_d, grad_u_d, xi_d, grad_xi_d, gln_d, theta_d, dtb_d, em_d, &hc);
    for (int v = phase_range.lower[cdim]; v <= phase_range.upper[cdim]; ++v) {
      idx[cdim] = v;
      long pidx = gkyl_range_idx(&phase_range, idx);
      double xc[GKYL_MAX_DIM];
      gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);
      double dw = up->phase_grid.dx[cdim];
      pkpm_hm_apply(up, &hc, xc[cdim], dw, (const double *)gkyl_array_cfetch(fin, pidx),
        (double *)gkyl_array_fetch(rhs, pidx));
      pkpm_hm_apply_cf(up, &hc, jac_d, (const double *)gkyl_array_cfetch(fin, pidx),
        (double *)gkyl_array_fetch(rhs, pidx));
      double *cfl = (double *)gkyl_array_fetch(cflrate, pidx);
      cfl[0] += pkpm_hm_rate(up, bvar_d, grad_b_d, u_d, grad_u_d, xi_d, grad_xi_d, gln_d, theta_d,
        dtb_d, em_d, xc[cdim] - 0.5 * dw, xc[cdim] + 0.5 * dw);
    }
  }
}

void
gkyl_dg_pkpm_harm_vol_advance_cu(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *bvar, const struct gkyl_array *grad_b,
  const struct gkyl_array *u, const struct gkyl_array *grad_u, const struct gkyl_array *xi,
  const struct gkyl_array *grad_xi, const struct gkyl_array *grad_ln_n, const struct gkyl_array *theta,
  const struct gkyl_array *dtB_over_magB, const struct gkyl_array *em, const struct gkyl_array *fin,
  struct gkyl_array *cflrate, struct gkyl_array *rhs
)
{
  gkyl_dg_pkpm_harm_vol_advance_cu_kernel<<<conf_range->nblocks, conf_range->nthreads>>>(
    up->on_dev, *conf_range, *phase_range, bvar->on_dev, grad_b->on_dev, u->on_dev, grad_u->on_dev,
    xi->on_dev, grad_xi->on_dev, grad_ln_n->on_dev, theta->on_dev, dtB_over_magB->on_dev, em->on_dev,
    fin->on_dev, cflrate->on_dev, rhs->on_dev
  );
}

struct gkyl_dg_pkpm_harm_vol *
gkyl_dg_pkpm_harm_vol_cu_dev_inew(const struct gkyl_dg_pkpm_harm_vol_inp *inp)
{
  struct gkyl_dg_pkpm_harm_vol *up = (struct gkyl_dg_pkpm_harm_vol *)gkyl_malloc(sizeof(*up));
  gkyl_dg_pkpm_harm_vol_setup(up, inp);
  up->flags = 0;
  GKYL_SET_CU_ALLOC(up->flags);
  // The device copy needs the device pointer of the Jacobian array; the host object keeps the
  // acquired map for its own lifetime management.
  up->jacob_pos = inp->pos_map->jacob_pos->on_dev;
  struct gkyl_dg_pkpm_harm_vol *up_cu =
    (struct gkyl_dg_pkpm_harm_vol *)gkyl_cu_malloc(sizeof(*up_cu));
  gkyl_cu_memcpy(up_cu, up, sizeof(*up), GKYL_CU_MEMCPY_H2D);
  up->on_dev = up_cu;
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->jacob_pos = inp->pos_map->jacob_pos;
  return up;
}
