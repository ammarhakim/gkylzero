#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_dg_pkpm_harm_vol_priv.h>

// Basis maps, degree <= 2 modes and the nonzero configuration triples of the tables.
void
gkyl_dg_pkpm_harm_vol_setup(struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_dg_pkpm_harm_vol_inp *inp)
{
  int cdim = inp->conf_basis->ndim;
  assert(inp->conf_basis->poly_order == 1 && inp->conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY);
  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->nc = inp->conf_basis->num_basis;
  up->np = inp->phase_basis->num_basis;
  up->qbym = inp->qbym;
  up->use_gpu = inp->use_gpu;
  assert(inp->pos_map);
  up->pos_map = 0;
  up->jacob_pos = inp->pos_map->jacob_pos;
  const int *cdeg = cdim == 1 ? pkpm_harm_1x1v_p1_conf_deg :
    cdim == 2 ? pkpm_harm_2x1v_p1_conf_deg : pkpm_harm_3x1v_p1_conf_deg;
  const int *pconf = cdim == 1 ? pkpm_harm_1x1v_p1_phase_conf :
    cdim == 2 ? pkpm_harm_2x1v_p1_phase_conf : pkpm_harm_3x1v_p1_phase_conf;
  const int *pwdeg = cdim == 1 ? pkpm_harm_1x1v_p1_phase_wdeg :
    cdim == 2 ? pkpm_harm_2x1v_p1_phase_wdeg : pkpm_harm_3x1v_p1_phase_wdeg;
  for (int k = 0; k < 3 * up->nc; ++k) {
    up->conf_deg[k] = cdeg[k];
  }
  for (int p = 0; p < up->np; ++p) {
    up->phase_of[pconf[p]][pwdeg[p]] = p;
  }
  // degree <= 2 modes, k = k0 + 3 k1 + 9 k2
  up->n2 = 1;
  up->n5 = 1;
  for (int d = 0; d < cdim; ++d) {
    up->n2 *= 3;
    up->n5 *= 5;
  }
  int kd[PKPM_HM_N2][3];
  for (int k = 0; k < up->n2; ++k) {
    kd[k][0] = k % 3;
    kd[k][1] = cdim > 1 ? (k / 3) % 3 : 0;
    kd[k][2] = cdim > 2 ? k / 9 : 0;
    up->p2_idx[k] = kd[k][0] + 5 * kd[k][1] + 25 * kd[k][2];
  }
  // nonzero int phi_m Q_k phi_l = prod_d t1(m_d, k_d, l_d)
  up->ntrip = 0;
  for (int m = 0; m < up->nc; ++m) {
    for (int k = 0; k < up->n2; ++k) {
      for (int l = 0; l < up->nc; ++l) {
        double v = 1.0;
        for (int d = 0; d < cdim; ++d) {
          v *= pkpm_harm_t1[(up->conf_deg[3 * m + d] * 3 + kd[k][d]) * 2 + up->conf_deg[3 * l + d]];
        }
        if (v != 0.0) {
          assert(up->ntrip < PKPM_HM_NTRIP);
          up->trip_m[up->ntrip] = m;
          up->trip_k[up->ntrip] = k;
          up->trip_l[up->ntrip] = l;
          up->trip_v[up->ntrip] = v;
          up->ntrip += 1;
        }
      }
    }
  }
  // nonzero int d_d phi_m Q_k phi_l = sqrt 3 prod_e t1(m'_e, k_e, l_e), m' = m with m'_d = 0
  // (d_d Phat_1 = sqrt 3 Phat_0 on the reference cell)
  for (int d = 0; d < 3; ++d) {
    up->ntripd[d] = 0;
  }
  for (int d = 0; d < cdim; ++d) {
    for (int m = 0; m < up->nc; ++m) {
      if (up->conf_deg[3 * m + d] != 1) {
        continue;
      }
      for (int k = 0; k < up->n2; ++k) {
        for (int l = 0; l < up->nc; ++l) {
          double v = 1.7320508075688772;
          for (int e = 0; e < cdim; ++e) {
            int me = e == d ? 0 : up->conf_deg[3 * m + e];
            v *= pkpm_harm_t1[(me * 3 + kd[k][e]) * 2 + up->conf_deg[3 * l + e]];
          }
          if (v != 0.0) {
            assert(up->ntripd[d] < PKPM_HM_NTRIPD);
            up->tripd_m[d][up->ntripd[d]] = m;
            up->tripd_k[d][up->ntripd[d]] = k;
            up->tripd_l[d][up->ntripd[d]] = l;
            up->tripd_v[d][up->ntripd[d]] = v;
            up->ntripd[d] += 1;
          }
        }
      }
    }
  }
}

struct gkyl_dg_pkpm_harm_vol *
gkyl_dg_pkpm_harm_vol_inew(const struct gkyl_dg_pkpm_harm_vol_inp *inp)
{
#ifdef GKYL_HAVE_CUDA
  if (inp->use_gpu) {
    return gkyl_dg_pkpm_harm_vol_cu_dev_inew(inp);
  }
#endif
  struct gkyl_dg_pkpm_harm_vol *up = gkyl_malloc(sizeof(*up));
  gkyl_dg_pkpm_harm_vol_setup(up, inp);
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  return up;
}

void
gkyl_dg_pkpm_harm_vol_advance(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *bvar, const struct gkyl_array *grad_b,
  const struct gkyl_array *u, const struct gkyl_array *grad_u, const struct gkyl_array *xi,
  const struct gkyl_array *grad_xi, const struct gkyl_array *grad_ln_n, const struct gkyl_array *theta,
  const struct gkyl_array *dtB_over_magB, const struct gkyl_array *em, const struct gkyl_array *fin,
  struct gkyl_array *cflrate, struct gkyl_array *rhs
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(rhs)) {
    gkyl_dg_pkpm_harm_vol_advance_cu(
      up, conf_range, phase_range, bvar, grad_b, u, grad_u, xi, grad_xi, grad_ln_n, theta,
      dtB_over_magB, em, fin, cflrate, rhs
    );
    return;
  }
#endif
  int cdim = up->cdim;
  struct pkpm_hm_coeffs *hc = gkyl_malloc(sizeof(*hc));
  int idx[GKYL_MAX_DIM];
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, conf_range);
  while (gkyl_range_iter_next(&citer)) {
    long cidx = gkyl_range_idx(conf_range, citer.idx);
    const double *bvar_d = gkyl_array_cfetch(bvar, cidx), *grad_b_d = gkyl_array_cfetch(grad_b, cidx);
    const double *u_d = gkyl_array_cfetch(u, cidx), *grad_u_d = gkyl_array_cfetch(grad_u, cidx);
    const double *xi_d = gkyl_array_cfetch(xi, cidx), *grad_xi_d = gkyl_array_cfetch(grad_xi, cidx);
    const double *gln_d = gkyl_array_cfetch(grad_ln_n, cidx), *theta_d = gkyl_array_cfetch(theta, cidx);
    const double *dtb_d = gkyl_array_cfetch(dtB_over_magB, cidx), *em_d = gkyl_array_cfetch(em, cidx);
    const double *jac_d = gkyl_array_cfetch(up->jacob_pos, cidx);
    pkpm_hm_build(up, bvar_d, grad_b_d, u_d, grad_u_d, xi_d, grad_xi_d, gln_d, theta_d, dtb_d, em_d, hc);
    gkyl_copy_int_arr(cdim, citer.idx, idx);
    for (int v = phase_range->lower[cdim]; v <= phase_range->upper[cdim]; ++v) {
      idx[cdim] = v;
      long pidx = gkyl_range_idx(phase_range, idx);
      double xc[GKYL_MAX_DIM];
      gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);
      double dw = up->phase_grid.dx[cdim];
      pkpm_hm_apply(up, hc, xc[cdim], dw, gkyl_array_cfetch(fin, pidx), gkyl_array_fetch(rhs, pidx));
      pkpm_hm_apply_cf(up, hc, jac_d, gkyl_array_cfetch(fin, pidx), gkyl_array_fetch(rhs, pidx));
      double *cfl = gkyl_array_fetch(cflrate, pidx);
      cfl[0] += pkpm_hm_rate(up, bvar_d, grad_b_d, u_d, grad_u_d, xi_d, grad_xi_d, gln_d, theta_d,
        dtb_d, em_d, xc[cdim] - 0.5 * dw, xc[cdim] + 0.5 * dw);
    }
  }
  gkyl_free(hc);
}

void
gkyl_dg_pkpm_harm_vol_release(struct gkyl_dg_pkpm_harm_vol *up)
{
  gkyl_vlasov_position_map_release(up->pos_map);
  if (GKYL_IS_CU_ALLOC(up->flags)) {
    gkyl_cu_free(up->on_dev);
  }
  gkyl_free(up);
}

#ifndef GKYL_HAVE_CUDA

struct gkyl_dg_pkpm_harm_vol *
gkyl_dg_pkpm_harm_vol_cu_dev_inew(const struct gkyl_dg_pkpm_harm_vol_inp *inp)
{
  assert(false);
  return 0;
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
  assert(false);
}

#endif
