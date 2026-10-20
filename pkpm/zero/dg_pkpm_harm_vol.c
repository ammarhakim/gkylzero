#include <assert.h>

#include <gkyl_alloc.h>
#include <gkyl_alloc_flags_priv.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_dg_pkpm_harm_vol_priv.h>

// Add the slot term gain * (slot) to the stored-operator entry [o, i, kind, j] (kind 0 source,
// 1 w flux; gkyl_pkpm_harm_ops.h).
static void
pkpm_hm_ent_term(struct gkyl_dg_pkpm_harm_vol *up, int o, int i, int kind, int j, int slot, double gain)
{
  int e = 0;
  while (e < PKPM_HO_NENT && !(pkpm_ho_ent[e][0] == o && pkpm_ho_ent[e][1] == i &&
    pkpm_ho_ent[e][2] == kind && pkpm_ho_ent[e][3] == j)) {
    ++e;
  }
  assert(e < PKPM_HO_NENT);
  int t = up->ent_nterm[e];
  assert(t < PKPM_HM_ENTTERM);
  up->ent_slot[e][t] = slot;
  up->ent_gain[e][t] = gain;
  up->ent_nterm[e] = t + 1;
}

void
gkyl_dg_pkpm_harm_vol_setup(struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_dg_pkpm_harm_vol_inp *inp)
{
  int cdim = inp->conf_basis->ndim;
  assert(inp->conf_basis->poly_order == 1 && inp->conf_basis->b_type == GKYL_BASIS_MODAL_SERENDIPITY);
  up->phase_grid = *inp->phase_grid;
  up->cdim = cdim;
  up->nc = inp->conf_basis->num_basis;
  up->qbym = inp->qbym;
  up->use_gpu = inp->use_gpu;
  assert(inp->pos_map);
  up->pos_map = 0;
  up->jacob_pos = inp->pos_map->jacob_pos;
  const int *cdeg = cdim == 1 ? pkpm_harm_1x1v_p1_conf_deg :
    cdim == 2 ? pkpm_harm_2x1v_p1_conf_deg : pkpm_harm_3x1v_p1_conf_deg;
  for (int k = 0; k < 3 * up->nc; ++k) {
    up->conf_deg[k] = cdeg[k];
  }
  for (int q = 0; q < up->nc; ++q) {
    up->conf_of[cdeg[3 * q] + 2 * cdeg[3 * q + 1] + 4 * cdeg[3 * q + 2]] = q;
  }
  up->nslot = PKPM_HS_NSLOT(cdim);
  up->nop = PKPM_HO_NOP(cdim);
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
    // Phat_0(0) = 1/sqrt 2, Phat_1(0) = 0, Phat_2(0) = -sqrt(5/2)/2
    const double c1[3] = {0.7071067811865475, 0.0, -0.7905694150420949};
    up->cen[k] = 1.0;
    for (int d = 0; d < cdim; ++d) {
      up->cen[k] *= c1[kd[k][d]];
    }
  }
  // the slot terms of the stored-operator entries (the w-flux couplings of F and G with themselves
  // are the configuration-pass kernel's)
  for (int e = 0; e < PKPM_HO_NENT; ++e) {
    up->ent_nterm[e] = 0;
  }
  for (int b = 2; b < 8; ++b) {
    pkpm_hm_ent_term(up, b, b, 1, 1, PKPM_HS_KW, 1.0);
    pkpm_hm_ent_term(up, b, b, 1, 2, PKPM_HS_BG, 1.0);
  }
  for (int a = 0; a < 3; ++a) {
    for (int j = 0; j < 2; ++j) {
      pkpm_hm_ent_term(up, 0, 2 + a, 1, j, PKPM_HS_TF(a, j), 1.0);
      pkpm_hm_ent_term(up, 1, 5 + a, 1, j, PKPM_HS_TG(a, j), 1.0);
      pkpm_hm_ent_term(up, 2 + a, 1, 1, j, PKPM_HS_TM(a, j), 1.0);
      pkpm_hm_ent_term(up, 5 + a, 1, 1, j, PKPM_HS_TNG(a, j), 1.0);
      pkpm_hm_ent_term(up, 5 + a, 0, 1, j, PKPM_HS_TNF(a, j), 1.0);
    }
    for (int c = 0; c < 3; ++c) {
      pkpm_hm_ent_term(up, 2 + a, 5 + c, 1, 0, PKPM_HS_RMN(a, c), 1.0);
      pkpm_hm_ent_term(up, 5 + a, 5 + c, 1, 0, PKPM_HS_RNN(a, c), 1.0);
      pkpm_hm_ent_term(up, 5 + a, 2 + c, 1, 0, PKPM_HS_RNM(a, c), 1.0);
    }
    pkpm_hm_ent_term(up, 0, 2 + a, 0, 0, PKPM_HS_SFM(a), 1.0);
    pkpm_hm_ent_term(up, 1, 5 + a, 0, 0, PKPM_HS_SGN(a), 1.0);
    pkpm_hm_ent_term(up, 2 + a, 1, 0, 0, PKPM_HS_SMG(a), 1.0);
    pkpm_hm_ent_term(up, 5 + a, 1, 0, 0, PKPM_HS_SNG0(a), 1.0);
    pkpm_hm_ent_term(up, 5 + a, 0, 0, 0, PKPM_HS_SNF(a), 1.0);
    for (int j = 1; j < 3; ++j) {
      pkpm_hm_ent_term(up, 1, 2 + a, 0, j, PKPM_HS_SGM(a, j), 1.0);
      pkpm_hm_ent_term(up, 2 + a, 0, 0, j, PKPM_HS_SMF(a, j), 1.0);
      pkpm_hm_ent_term(up, 5 + a, 1, 0, j, PKPM_HS_SNG1(a, j), 1.0);
    }
    for (int c = 0; c < 3; ++c) {
      for (int j = 0; j < 2; ++j) {
        pkpm_hm_ent_term(up, 2 + a, 2 + c, 0, j, PKPM_HS_SMM(a, c, j), 1.0);
        pkpm_hm_ent_term(up, 5 + a, 5 + c, 0, j, PKPM_HS_SNN(a, c, j), 1.0);
      }
    }
  }
  for (int e = 0; e < PKPM_HO_NENT; ++e) {
    assert(up->ent_nterm[e] > 0);
  }
  // the source pairs of the rate estimate
  up->nsrc = 0;
  for (int e = 0; e < PKPM_HO_NENT; ++e) {
    if (pkpm_ho_ent[e][2] != 0) {
      continue;
    }
    int p = 0;
    while (p < up->nsrc && !(up->src_o[p] == pkpm_ho_ent[e][0] && up->src_i[p] == pkpm_ho_ent[e][1])) {
      ++p;
    }
    if (p == up->nsrc) {
      assert(up->nsrc < PKPM_HM_MAXSRC);
      up->src_o[p] = pkpm_ho_ent[e][0];
      up->src_i[p] = pkpm_ho_ent[e][1];
      up->src_ent[p][0] = up->src_ent[p][1] = up->src_ent[p][2] = -1;
      up->nsrc += 1;
    }
    up->src_ent[p][pkpm_ho_ent[e][3]] = e;
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
  up->mat_kern = ser_pkpm_harm_mat_kernels[up->cdim - 1];
  up->vol_kern = ser_pkpm_harm_vol_kernels[up->cdim - 1];
  up->pos_map = gkyl_vlasov_position_map_acquire(inp->pos_map);
  up->flags = 0;
  GKYL_CLEAR_CU_ALLOC(up->flags);
  up->on_dev = up; // self-reference on host
  return up;
}

int
gkyl_dg_pkpm_harm_vol_num_coef(const struct gkyl_dg_pkpm_harm_vol *up)
{
  return up->nslot * up->n2;
}

void
gkyl_dg_pkpm_harm_vol_coef(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *xi, const struct gkyl_array *grad_xi,
  const struct gkyl_array *grad_ln_n, const struct gkyl_array *theta,
  const struct gkyl_array *dtB_over_magB, const struct gkyl_array *em, struct gkyl_array *coef
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(coef)) {
    gkyl_dg_pkpm_harm_vol_coef_cu(
      up, conf_range, bvar, grad_b, u, grad_u, xi, grad_xi, grad_ln_n, theta, dtB_over_magB, em, coef
    );
    return;
  }
#endif
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, conf_range);
  while (gkyl_range_iter_next(&citer)) {
    long cidx = gkyl_range_idx(conf_range, citer.idx);
    pkpm_hm_build(up, gkyl_array_cfetch(bvar, cidx), gkyl_array_cfetch(grad_b, cidx),
      gkyl_array_cfetch(u, cidx), gkyl_array_cfetch(grad_u, cidx), gkyl_array_cfetch(xi, cidx),
      gkyl_array_cfetch(grad_xi, cidx), gkyl_array_cfetch(grad_ln_n, cidx),
      gkyl_array_cfetch(theta, cidx), gkyl_array_cfetch(dtB_over_magB, cidx),
      gkyl_array_cfetch(em, cidx), gkyl_array_fetch(coef, cidx));
  }
}

void
gkyl_dg_pkpm_harm_vol_advance(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *fin,
  struct gkyl_array *cflrate, struct gkyl_array *rhs
)
{
#ifdef GKYL_HAVE_CUDA
  if (gkyl_array_is_cu_dev(rhs)) {
    gkyl_dg_pkpm_harm_vol_advance_cu(up, conf_range, phase_range, coef, fin, cflrate, rhs);
    return;
  }
#endif
  int cdim = up->cdim;
  int idx[GKYL_MAX_DIM];
  double ho[PKPM_HO_NOP_MAX * PKPM_HM_N2], mat[PKPM_HM_MAT], poly[PKPM_HM_MAXSRC][3];
  double dw = up->phase_grid.dx[cdim];
  struct gkyl_range_iter citer;
  gkyl_range_iter_init(&citer, conf_range);
  while (gkyl_range_iter_next(&citer)) {
    long cidx = gkyl_range_idx(conf_range, citer.idx);
    // the stored operator of the configuration cell, its matrices and its source rate
    pkpm_hm_pack(up, gkyl_array_cfetch(coef, cidx), ho);
    up->mat_kern(up->phase_grid.dx, gkyl_array_cfetch(up->jacob_pos, cidx), ho, mat);
    pkpm_hm_rate_poly(up, ho, poly);
    gkyl_copy_int_arr(cdim, citer.idx, idx);
    for (int v = phase_range->lower[cdim]; v <= phase_range->upper[cdim]; ++v) {
      idx[cdim] = v;
      long pidx = gkyl_range_idx(phase_range, idx);
      double xc[GKYL_MAX_DIM];
      gkyl_rect_grid_cell_center(&up->phase_grid, idx, xc);
      up->vol_kern(xc, up->phase_grid.dx, mat, gkyl_array_cfetch(fin, pidx), gkyl_array_fetch(rhs, pidx));
      double *cfl = gkyl_array_fetch(cflrate, pidx);
      cfl[0] += pkpm_hm_rate(up, poly, xc[cdim] - 0.5 * dw, xc[cdim] + 0.5 * dw);
    }
  }
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
gkyl_dg_pkpm_harm_vol_coef_cu(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_array *bvar, const struct gkyl_array *grad_b, const struct gkyl_array *u,
  const struct gkyl_array *grad_u, const struct gkyl_array *xi, const struct gkyl_array *grad_xi,
  const struct gkyl_array *grad_ln_n, const struct gkyl_array *theta,
  const struct gkyl_array *dtB_over_magB, const struct gkyl_array *em, struct gkyl_array *coef
)
{
  assert(false);
}

void
gkyl_dg_pkpm_harm_vol_advance_cu(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_range *conf_range,
  const struct gkyl_range *phase_range, const struct gkyl_array *coef, const struct gkyl_array *fin,
  struct gkyl_array *cflrate, struct gkyl_array *rhs
)
{
  assert(false);
}

#endif
