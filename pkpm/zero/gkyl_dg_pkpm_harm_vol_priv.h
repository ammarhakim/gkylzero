// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_pkpm_harm_modal_tables.h>
#include <gkyl_pkpm_harm_ops.h>
#include <gkyl_pkpm_harm_slots.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>

// Tensor orthonormal-Legendre expansions in configuration space: mode (n0, n1, n2) at
// n0 + 5 n1 + 25 n2, degrees up to 4 per direction (n_d = 0 beyond cdim).
#define PKPM_HM_N 125
// Projections onto tensor degree <= 2 per direction (3^cdim modes, k0 + 3 k1 + 9 k2).
#define PKPM_HM_N2 27
// Slot terms of a stored-operator entry (the stored slots packed into the operator,
// gkyl_pkpm_harm_ops.h) and (output, input) block pairs of the source-rate estimate.
#define PKPM_HM_ENTTERM 4
#define PKPM_HM_MAXSRC 64
// The nc x nc matrices of the stored operator of a configuration cell (pkpm_harm_mat_*).
#define PKPM_HM_MAT (PKPM_HO_NOP_MAX * 64)

// Generated kernels: the matrices of the stored operator of a configuration cell, and the volume
// terms of a phase cell from them.
typedef void (*pkpm_harm_mat_t)(
  const double *dxv, const double *jacob_pos, const double *ho, double *GKYL_RESTRICT mat
);
typedef void (*pkpm_harm_vol_t)(
  const double *w, const double *dxv, const double *mat, const double *f, double *GKYL_RESTRICT out
);

struct gkyl_dg_pkpm_harm_vol {
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  int cdim; // Configuration-space dimensions.
  int nc; // Configuration basis size.
  int n2; // Number of degree <= 2 tensor modes (3^cdim).
  int n5; // Number of expansion slots in use (5^cdim; the higher ones belong to absent directions).
  int nslot; // Number of stored coefficient slots (gkyl_pkpm_harm_slots.h).
  int nop; // Number of entries of the stored operator (gkyl_pkpm_harm_ops.h).
  double qbym; // Charge to mass ratio (gyration).
  int conf_deg[8 * 3]; // 1D Legendre degrees of each configuration basis function.
  int conf_of[8]; // Configuration basis function of each degree pattern s0 + 2 s1 + 4 s2.
  int p2_idx[PKPM_HM_N2]; // Expansion index of each degree <= 2 mode.
  double cen[PKPM_HM_N2]; // Value of each degree <= 2 mode at the cell center.
  // the slot terms of each operator entry: entry e = sum_t ent_gain[e][t] (slot ent_slot[e][t])
  int ent_nterm[PKPM_HO_NENT], ent_slot[PKPM_HO_NENT][PKPM_HM_ENTTERM];
  double ent_gain[PKPM_HO_NENT][PKPM_HM_ENTTERM];
  // the (output, input) source pairs of the rate estimate: the source entry of each power of w
  // (-1 if absent); the gyration is not rated (it is the exponential part of the time integrator)
  int nsrc, src_o[PKPM_HM_MAXSRC], src_i[PKPM_HM_MAXSRC], src_ent[PKPM_HM_MAXSRC][3];
  pkpm_harm_mat_t mat_kern; // matrices of the stored operator of a configuration cell
  pkpm_harm_vol_t vol_kern; // volume terms of a phase cell
  const struct gkyl_vlasov_position_map *pos_map; // Position map (acquired; host pointer, 0 on the device copy).
  const struct gkyl_array *jacob_pos; // Per-direction Jacobians (borrowed from pos_map; device pointer on device).
  bool use_gpu; // Whether the updater runs on the GPU.
  uint32_t flags;
  struct gkyl_dg_pkpm_harm_vol *on_dev; // pointer to itself or device data
};

// Fill the updater (sizes, basis maps, degree <= 2 modes, the slot terms of the operator entries
// and the source pairs of the rate) from the input and the generated tables; host function shared
// by the host and device constructors (the kernel pointers are set by each).
void gkyl_dg_pkpm_harm_vol_setup(struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_dg_pkpm_harm_vol_inp *inp);

// Generated kernels (p=1 serendipity; indexed by cdim-1).
GKYL_CU_D static const pkpm_harm_mat_t ser_pkpm_harm_mat_kernels[] = {
  pkpm_harm_mat_1x_ser_p1,
  pkpm_harm_mat_2x_ser_p1,
  pkpm_harm_mat_3x_ser_p1,
};
GKYL_CU_D static const pkpm_harm_vol_t ser_pkpm_harm_vol_kernels[] = {
  pkpm_harm_vol_1x1v_ser_p1,
  pkpm_harm_vol_2x1v_ser_p1,
  pkpm_harm_vol_3x1v_ser_p1,
};

GKYL_CU_DH static inline void
pkpm_hm_zero(const struct gkyl_dg_pkpm_harm_vol *up, double *e)
{
  for (int n = 0; n < up->n5; ++n) {
    e[n] = 0.0;
  }
}

// out += a x (whole arrays).
GKYL_CU_DH static inline void
pkpm_hm_axpy(const struct gkyl_dg_pkpm_harm_vol *up, double *out, double a, const double *x)
{
  for (int n = 0; n < up->n5; ++n) {
    out[n] += a * x[n];
  }
}

// Expansion of a p=1 configuration field given by its coefficients in the gkyl basis.
GKYL_CU_DH static inline void
pkpm_hm_p1(const struct gkyl_dg_pkpm_harm_vol *up, const double *f, double *out)
{
  pkpm_hm_zero(up, out);
  for (int q = 0; q < up->nc; ++q) {
    const int *s = &up->conf_deg[3 * q];
    out[s[0] + 5 * s[1] + 25 * s[2]] = f[q];
  }
}

// One direction of a p=1 Legendre factor: out = Phat_s(x_d) in on the box of limits olim
// (s = 0: Phat_0 = 1/sqrt 2, olim_d = min(lim_d, wout); s = 1: Phat_1 Phat_n = up[n] Phat_{n+1} +
// dn[n] Phat_{n-1}, olim_d = min(lim_d + 1, wout)); in is read on its box lim.
GKYL_CU_DH static inline void
pkpm_hm_dir(const int *lim, int d, int s, int wout, const double *in, double *out, int *olim)
{
  const int stride[3] = {1, 5, 25};
  for (int e = 0; e < 3; ++e) {
    olim[e] = lim[e];
  }
  int ld = lim[d] + s;
  olim[d] = ld < wout ? ld : wout;
  for (int n2 = 0; n2 <= olim[2]; ++n2) {
    for (int n1 = 0; n1 <= olim[1]; ++n1) {
      for (int n0 = 0; n0 <= olim[0]; ++n0) {
        int n[3] = {n0, n1, n2}, idx = n0 + 5 * n1 + 25 * n2, m = n[d];
        double v;
        if (s == 0) {
          v = 0.7071067811865475 * in[idx];
        } else {
          v = 0.0;
          if (m >= 1) {
            v += pkpm_harm_xp_up[m - 1] * in[idx - stride[d]];
          }
          if (m + 1 <= lim[d]) {
            v += pkpm_harm_xp_dn[m + 1] * in[idx + stride[d]];
          }
        }
        out[idx] = v;
      }
    }
  }
}

// out = c0 x0 + c1 x1 on the union of the boxes l0, l1 (each read on its own box, zero outside);
// ol receives the union.
GKYL_CU_DH static inline void
pkpm_hm_comb(
  double c0, const double *x0, const int *l0, double c1, const double *x1, const int *l1, double *out,
  int *ol
)
{
  for (int e = 0; e < 3; ++e) {
    ol[e] = l0[e] > l1[e] ? l0[e] : l1[e];
  }
  for (int n2 = 0; n2 <= ol[2]; ++n2) {
    for (int n1 = 0; n1 <= ol[1]; ++n1) {
      for (int n0 = 0; n0 <= ol[0]; ++n0) {
        int idx = n0 + 5 * n1 + 25 * n2;
        double v = 0.0;
        if (n0 <= l0[0] && n1 <= l0[1] && n2 <= l0[2]) {
          v += c0 * x0[idx];
        }
        if (n0 <= l1[0] && n1 <= l1[1] && n2 <= l1[2]) {
          v += c1 * x1[idx];
        }
        out[idx] = v;
      }
    }
  }
}

// out += coef P_wout[a f]: a an expansion exact through degree wa per direction, f a p=1
// configuration field (gkyl basis coefficients). f = sum_s f_s prod_d Phat_{s_d}(x_d) is applied
// direction by direction with its coefficients contracted as early as possible (sum
// factorization): two 1D operations in the first direction, then one per remaining coefficient
// pattern. The result is exact through degree min(wout, wa - 1 + [a complete]) as the caller plans.
GKYL_CU_DH static inline void
pkpm_hm_mul_acc(
  const struct gkyl_dg_pkpm_harm_vol *up, const double *a, int wa, const double *f, double coef,
  double *out, int wout
)
{
  int cdim = up->cdim, ns = 1 << cdim;
  double fs[8];
  bool any = false;
  for (int s = 0; s < ns; ++s) {
    fs[s] = f[up->conf_of[s]];
    any = any || fs[s] != 0.0;
  }
  if (!any) {
    return;
  }
  int lim[3] = {0, 0, 0};
  for (int d = 0; d < cdim; ++d) {
    lim[d] = wa;
  }
  // first direction: a/sqrt 2 and Phat_1(x_0) a
  double A0[PKPM_HM_N], A1[PKPM_HM_N];
  int l0[3], l1[3];
  pkpm_hm_dir(lim, 0, 0, wout, a, A0, l0);
  pkpm_hm_dir(lim, 0, 1, wout, a, A1, l1);
  double R[PKPM_HM_N];
  int lr[3];
  if (cdim == 1) {
    pkpm_hm_comb(fs[0], A0, l0, fs[1], A1, l1, R, lr);
  } else {
    // contract the first direction, apply the second
    double C[2][PKPM_HM_N], D0[PKPM_HM_N], D1[PKPM_HM_N], E[2][PKPM_HM_N];
    int lc[2][3], ld0[3], ld1[3], le[2][3];
    int nz = cdim > 2 ? 2 : 1;
    for (int s2 = 0; s2 < nz; ++s2) {
      for (int s1 = 0; s1 < 2; ++s1) {
        pkpm_hm_comb(fs[2 * s1 + 4 * s2], A0, l0, fs[1 + 2 * s1 + 4 * s2], A1, l1, C[s1], lc[s1]);
      }
      pkpm_hm_dir(lc[0], 1, 0, wout, C[0], D0, ld0);
      pkpm_hm_dir(lc[1], 1, 1, wout, C[1], D1, ld1);
      pkpm_hm_comb(1.0, D0, ld0, 1.0, D1, ld1, E[s2], le[s2]);
    }
    if (cdim == 2) {
      for (int e = 0; e < 3; ++e) {
        lr[e] = le[0][e];
      }
      for (int n = 0; n < up->n5; ++n) {
        R[n] = E[0][n];
      }
    } else {
      // third direction
      pkpm_hm_dir(le[0], 2, 0, wout, E[0], D0, ld0);
      pkpm_hm_dir(le[1], 2, 1, wout, E[1], D1, ld1);
      pkpm_hm_comb(1.0, D0, ld0, 1.0, D1, ld1, R, lr);
    }
  }
  for (int n2 = 0; n2 <= lr[2]; ++n2) {
    for (int n1 = 0; n1 <= lr[1]; ++n1) {
      for (int n0 = 0; n0 <= lr[0]; ++n0) {
        int idx = n0 + 5 * n1 + 25 * n2;
        out[idx] += coef * R[idx];
      }
    }
  }
}

// out = P_wout[a f].
GKYL_CU_DH static inline void
pkpm_hm_mul(
  const struct gkyl_dg_pkpm_harm_vol *up, const double *a, int wa, const double *f, double *out,
  int wout
)
{
  pkpm_hm_zero(up, out);
  pkpm_hm_mul_acc(up, a, wa, f, 1.0, out, wout);
}

// Slot s of the cell's coefficient array += coef (degree <= 2 projection of e).
GKYL_CU_DH static inline void
pkpm_hm_put(const struct gkyl_dg_pkpm_harm_vol *up, const double *e, double coef, double *cell, int s)
{
  double *dst = &cell[s * up->n2];
  for (int k = 0; k < up->n2; ++k) {
    dst[k] += coef * e[up->p2_idx[k]];
  }
}

// v_a - b_a (b.v) for expansions v (complete through degree wv) into windows wout: out[a] (3
// expansions), with bv scratch.
GKYL_CU_DH static inline void
pkpm_hm_perp(
  const struct gkyl_dg_pkpm_harm_vol *up, const double *const *b, double (*v)[PKPM_HM_N], int wv,
  double (*out)[PKPM_HM_N], int wout
)
{
  double bv[PKPM_HM_N];
  pkpm_hm_zero(up, bv);
  for (int a = 0; a < 3; ++a) {
    pkpm_hm_mul_acc(up, v[a], wv, b[a], 1.0, bv, wv + 1);
  }
  for (int a = 0; a < 3; ++a) {
    for (int n = 0; n < up->n5; ++n) {
      out[a][n] = v[a][n];
    }
    pkpm_hm_mul_acc(up, bv, wv + 1, b[a], -1.0, out[a], wout);
  }
}

/**
 * Coefficients of the first-harmonic terms of one configuration cell into its slots
 * (gkyl_pkpm_harm_slots.h; definitions of gkyl_pkpm_harm_coeffs.h, plan 17.11; all products exact,
 * each coefficient kept as its projection onto degree <= 2 per direction):
 *   w flux: alpha_k = -K w - (b.grad xi) w^2; F: t.M0 + beta G; G: t.N + beta H;
 *           M0: G t + R_B1 N; N: H t + R_B1 O
 *   source: F: -ell_xi.M0; G: f_perp.M0 - (ell_xi + 2g).N;
 *           M0: f_perp F - w b.ell_xi M0 + L_perp M0 - G I_perp(ell_xi + g) - b(C.M0 + beta G);
 *           N: 2 G f_perp - w b.ell_xi N + L_perp N - R_E1 N - H I_perp(ell_xi + 3g) - b(C.N + beta H)
 *   (the gyration of L_perp as the field B), configuration-flux couplings,
 * with H = 4 theta G - 2 theta^2 F, O = 6 theta N - 6 theta^2 M0. cell holds nslot n2 doubles.
 */
GKYL_CU_DH static inline void
pkpm_hm_build(
  const struct gkyl_dg_pkpm_harm_vol *up, const double *bvar, const double *grad_b,
  const double *u, const double *grad_u, const double *xi, const double *grad_xi,
  const double *grad_ln_n, const double *theta, const double *dtb, const double *em,
  double *cell
)
{
  int cdim = up->cdim, nc = up->nc;
  const double *b[3], *uu[3], *dd[3];
  for (int a = 0; a < 3; ++a) {
    b[a] = &bvar[a * nc];
    uu[a] = &u[a * nc];
    dd[a] = &dtb[a * nc];
  }
  double zero_p1[8] = {0.0};
  // U[a][c] = d_c u_a, Db[a][c] = d_c b_a, g[c] = d_c xi, l[c] = d_c ln n (zero beyond cdim)
  const double *U[3][3], *Db[3][3], *g[3], *l[3];
  for (int c = 0; c < 3; ++c) {
    g[c] = c < cdim ? &grad_xi[(c * 3) * nc] : zero_p1;
    l[c] = c < cdim ? &grad_ln_n[c * nc] : zero_p1;
    for (int a = 0; a < 3; ++a) {
      U[a][c] = c < cdim ? &grad_u[(c * 3 + a) * nc] : zero_p1;
      Db[a][c] = c < cdim ? &grad_b[(c * 3 + a) * nc] : zero_p1;
    }
  }
  for (int k = 0; k < up->nslot * up->n2; ++k) {
    cell[k] = 0.0;
  }

  double e1[PKPM_HM_N], e2[PKPM_HM_N], e3[PKPM_HM_N];

  // --- C = C0 + w C1: C0 = I_perp d + (u.grad) b, C1 = xi (b.grad) b (degree 3, complete)
  double C0[3][PKPM_HM_N], C1[3][PKPM_HM_N];
  {
    double dv[3][PKPM_HM_N];
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_p1(up, dd[a], dv[a]);
    }
    pkpm_hm_perp(up, b, dv, 1, C0, 3);
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_zero(up, e1);
      for (int c = 0; c < cdim; ++c) {
        pkpm_hm_p1(up, Db[a][c], e2);
        pkpm_hm_mul_acc(up, e2, 1, uu[c], 1.0, C0[a], 2);
        pkpm_hm_mul_acc(up, e2, 1, b[c], 1.0, e1, 2);
      }
      pkpm_hm_mul(up, e1, 2, xi, C1[a], 3);
    }
  }
  // --- U b, U^T b (degree 2), K = b.U b (3), b.g (2), beta = xi div b (2)
  double Ub[3][PKPM_HM_N], UTb[3][PKPM_HM_N], K[PKPM_HM_N], bg[PKPM_HM_N], beta[PKPM_HM_N];
  for (int a = 0; a < 3; ++a) {
    pkpm_hm_zero(up, Ub[a]);
    pkpm_hm_zero(up, UTb[a]);
  }
  for (int a = 0; a < 3; ++a) {
    for (int c = 0; c < cdim; ++c) {
      pkpm_hm_p1(up, U[a][c], e1);
      pkpm_hm_mul_acc(up, e1, 1, b[c], 1.0, Ub[a], 2);
      pkpm_hm_mul_acc(up, e1, 1, b[a], 1.0, UTb[c], 2);
    }
  }
  pkpm_hm_zero(up, K);
  for (int a = 0; a < 3; ++a) {
    pkpm_hm_mul_acc(up, Ub[a], 2, b[a], 1.0, K, 3);
  }
  pkpm_hm_zero(up, bg);
  pkpm_hm_zero(up, e1);
  for (int c = 0; c < cdim; ++c) {
    pkpm_hm_p1(up, g[c], e2);
    pkpm_hm_mul_acc(up, e2, 1, b[c], 1.0, bg, 2);
    pkpm_hm_p1(up, Db[c][c], e2);
    pkpm_hm_axpy(up, e1, 1.0, e2);
  }
  pkpm_hm_mul(up, e1, 1, xi, beta, 2);

  // --- t = t0 + w t1 = I_perp(C - U^T b - w g): window 4 (exact through 4)
  double t0[3][PKPM_HM_N], t1[3][PKPM_HM_N];
  {
    double X[3][PKPM_HM_N];
    for (int a = 0; a < 3; ++a) {
      for (int n = 0; n < up->n5; ++n) {
        X[a][n] = C0[a][n] - UTb[a][n];
      }
    }
    pkpm_hm_perp(up, b, X, 3, t0, 4);
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_p1(up, g[a], e1);
      for (int n = 0; n < up->n5; ++n) {
        X[a][n] = C1[a][n] - e1[n];
      }
    }
    pkpm_hm_perp(up, b, X, 3, t1, 4);
  }
  // --- R_B1 = (1/2)(beta I + I_perp (B_1 + B_1^T)), B_1 = xi grad b: degree 4, complete
  double PS[3][3][PKPM_HM_N];
  {
    double S[3][PKPM_HM_N], col[3][PKPM_HM_N];
    for (int c = 0; c < 3; ++c) {
      for (int a = 0; a < 3; ++a) {
        pkpm_hm_p1(up, Db[a][c], e1);
        pkpm_hm_p1(up, Db[c][a], e2);
        pkpm_hm_axpy(up, e1, 1.0, e2);
        pkpm_hm_mul(up, e1, 1, xi, S[a], 2);
      }
      pkpm_hm_perp(up, b, S, 2, col, 4);
      for (int a = 0; a < 3; ++a) {
        for (int n = 0; n < up->n5; ++n) {
          PS[a][c][n] = col[a][n];
        }
      }
    }
  }

  // ===== w fluxes =====
  // alpha_k = -w K - w^2 b.g on every block
  pkpm_hm_put(up, K, -1.0, cell, PKPM_HS_KW);
  pkpm_hm_put(up, bg, -1.0, cell, PKPM_HS_BG);
  // two-field couplings: F <- beta G, G <- beta H = 4 theta beta G - 2 theta^2 beta F
  pkpm_hm_put(up, beta, 1.0, cell, PKPM_HS_BF);
  pkpm_hm_mul(up, beta, 2, theta, e1, 3); // theta beta (complete)
  pkpm_hm_put(up, e1, 4.0, cell, PKPM_HS_BGG);
  pkpm_hm_mul(up, e1, 3, theta, e2, 2);
  pkpm_hm_put(up, e2, -2.0, cell, PKPM_HS_BGF);
  for (int a = 0; a < 3; ++a) {
    double *tj[2] = {t0[a], t1[a]};
    for (int j = 0; j < 2; ++j) {
      // F: t.M0; G: t.N; M0: G t
      pkpm_hm_put(up, tj[j], 1.0, cell, PKPM_HS_TF(a, j));
      pkpm_hm_put(up, tj[j], 1.0, cell, PKPM_HS_TG(a, j));
      pkpm_hm_put(up, tj[j], 1.0, cell, PKPM_HS_TM(a, j));
      // N: H t = (4 theta G - 2 theta^2 F) t
      pkpm_hm_mul(up, tj[j], 4, theta, e1, 3);
      pkpm_hm_put(up, e1, 4.0, cell, PKPM_HS_TNG(a, j));
      pkpm_hm_mul(up, e1, 3, theta, e2, 2);
      pkpm_hm_put(up, e2, -2.0, cell, PKPM_HS_TNF(a, j));
    }
    // M0: R_B1 N; N: R_B1 O = 6 theta R N - 6 theta^2 R M0
    for (int c = 0; c < 3; ++c) {
      for (int n = 0; n < up->n5; ++n) {
        e3[n] = 0.5 * PS[a][c][n] + (a == c ? 0.5 * beta[n] : 0.0);
      }
      pkpm_hm_put(up, e3, 1.0, cell, PKPM_HS_RMN(a, c));
      pkpm_hm_mul(up, e3, 4, theta, e1, 3);
      pkpm_hm_put(up, e1, 6.0, cell, PKPM_HS_RNN(a, c));
      pkpm_hm_mul(up, e1, 3, theta, e2, 2);
      pkpm_hm_put(up, e2, -6.0, cell, PKPM_HS_RNM(a, c));
    }
  }

  // ===== sources =====
  // f_perp = w fp1 + w^2 fp2 = -w I_perp(U b + C): window 2
  double fp1[3][PKPM_HM_N], fp2[3][PKPM_HM_N];
  {
    double Z[3][PKPM_HM_N];
    for (int a = 0; a < 3; ++a) {
      for (int n = 0; n < up->n5; ++n) {
        Z[a][n] = -(Ub[a][n] + C0[a][n]);
      }
    }
    pkpm_hm_perp(up, b, Z, 3, fp1, 2);
    for (int a = 0; a < 3; ++a) {
      for (int n = 0; n < up->n5; ++n) {
        Z[a][n] = -C1[a][n];
      }
    }
    pkpm_hm_perp(up, b, Z, 3, fp2, 2);
  }
  // ell_xi = xi grad ln n (2), sigma_0 = -w b.ell_xi (3)
  double xil[3][PKPM_HM_N], xibl[PKPM_HM_N];
  pkpm_hm_zero(up, xibl);
  for (int a = 0; a < 3; ++a) {
    pkpm_hm_p1(up, l[a], e1);
    pkpm_hm_mul(up, e1, 1, xi, xil[a], 2);
    pkpm_hm_mul_acc(up, xil[a], 2, b[a], 1.0, xibl, 3);
  }
  // I_perp(ell_xi + g) and I_perp(ell_xi + 3 g): degree 4, complete
  double plg[3][PKPM_HM_N], pl3g[3][PKPM_HM_N];
  {
    double V[3][PKPM_HM_N], V3[3][PKPM_HM_N];
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_p1(up, g[a], e1);
      for (int n = 0; n < up->n5; ++n) {
        V[a][n] = xil[a][n] + e1[n];
        V3[a][n] = xil[a][n] + 3.0 * e1[n];
      }
    }
    pkpm_hm_perp(up, b, V, 2, plg, 4);
    pkpm_hm_perp(up, b, V3, 2, pl3g, 4);
  }
  for (int a = 0; a < 3; ++a) {
    // F: -ell_xi.M0; G: f_perp.M0 - (ell_xi + 2 g).N
    pkpm_hm_put(up, xil[a], -1.0, cell, PKPM_HS_SFM(a));
    pkpm_hm_put(up, fp1[a], 1.0, cell, PKPM_HS_SGM(a, 1));
    pkpm_hm_put(up, fp2[a], 1.0, cell, PKPM_HS_SGM(a, 2));
    pkpm_hm_put(up, xil[a], -1.0, cell, PKPM_HS_SGN(a));
    pkpm_hm_p1(up, g[a], e1);
    pkpm_hm_put(up, e1, -2.0, cell, PKPM_HS_SGN(a));
    // M0: f_perp F; N: 2 G f_perp
    pkpm_hm_put(up, fp1[a], 1.0, cell, PKPM_HS_SMF(a, 1));
    pkpm_hm_put(up, fp2[a], 1.0, cell, PKPM_HS_SMF(a, 2));
    pkpm_hm_put(up, fp1[a], 2.0, cell, PKPM_HS_SNG1(a, 1));
    pkpm_hm_put(up, fp2[a], 2.0, cell, PKPM_HS_SNG1(a, 2));
    // M0: -G I_perp(ell_xi + g) - b beta G
    pkpm_hm_mul(up, beta, 2, b[a], e3, 3); // b_a beta (complete)
    pkpm_hm_put(up, plg[a], -1.0, cell, PKPM_HS_SMG(a));
    pkpm_hm_put(up, e3, -1.0, cell, PKPM_HS_SMG(a));
    // N: -H I_perp(ell_xi + 3 g) - b beta H, H = 4 theta G - 2 theta^2 F
    for (int n = 0; n < up->n5; ++n) {
      e2[n] = pl3g[a][n] + e3[n]; // I_perp(ell_xi + 3 g) + b beta: degree 4, complete
    }
    pkpm_hm_mul(up, e2, 4, theta, e1, 3);
    pkpm_hm_put(up, e1, -4.0, cell, PKPM_HS_SNG0(a));
    pkpm_hm_mul(up, e1, 3, theta, e2, 2);
    pkpm_hm_put(up, e2, 2.0, cell, PKPM_HS_SNF(a));
    // diagonal of M0 and N: sigma_0 = -w b.ell_xi, -w b.g from L_perp; N also -R_E1's
    // (1/2)(eta + 2 w b.g) = (1/2)(div u - K) + (1/2) w beta + 2 w b.g
    pkpm_hm_put(up, xibl, -1.0, cell, PKPM_HS_SMM(a, a, 1));
    pkpm_hm_put(up, bg, -1.0, cell, PKPM_HS_SMM(a, a, 1));
    pkpm_hm_put(up, xibl, -1.0, cell, PKPM_HS_SNN(a, a, 1));
    pkpm_hm_put(up, bg, -3.0, cell, PKPM_HS_SNN(a, a, 1));
    pkpm_hm_put(up, beta, -0.5, cell, PKPM_HS_SNN(a, a, 1));
    pkpm_hm_put(up, K, 0.5, cell, PKPM_HS_SNN(a, a, 0));
    for (int c = 0; c < cdim; ++c) {
      pkpm_hm_p1(up, U[c][c], e1);
      pkpm_hm_put(up, e1, -0.5, cell, PKPM_HS_SNN(a, a, 0));
    }
  }
  // Matrix couplings M0 <- M0 and N <- N: L_perp = -I_perp U - w I_perp B_1 (+ the gyration),
  // -b C. (both), and for N -(1/2) I_perp(U + U^T) - (1/2) w I_perp(B_1 + B_1^T)
  for (int c = 0; c < 3; ++c) {
    double col[3][PKPM_HM_N], pcol[3][PKPM_HM_N];
    // -I_perp U (column c): degree 3
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_p1(up, U[a][c], col[a]);
    }
    pkpm_hm_perp(up, b, col, 1, pcol, 2);
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_put(up, pcol[a], -1.0, cell, PKPM_HS_SMM(a, c, 0));
      pkpm_hm_put(up, pcol[a], -1.0, cell, PKPM_HS_SNN(a, c, 0));
    }
    // -(1/2) I_perp(U + U^T) for N
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_p1(up, U[a][c], col[a]);
      pkpm_hm_p1(up, U[c][a], e1);
      pkpm_hm_axpy(up, col[a], 1.0, e1);
    }
    pkpm_hm_perp(up, b, col, 1, pcol, 2);
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_put(up, pcol[a], -0.5, cell, PKPM_HS_SNN(a, c, 0));
    }
    // -w I_perp B_1 (column c): B_1 = xi grad b, degree 4
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_p1(up, Db[a][c], e1);
      pkpm_hm_mul(up, e1, 1, xi, col[a], 2);
    }
    pkpm_hm_perp(up, b, col, 2, pcol, 2);
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_put(up, pcol[a], -1.0, cell, PKPM_HS_SMM(a, c, 1));
      pkpm_hm_put(up, pcol[a], -1.0, cell, PKPM_HS_SNN(a, c, 1));
      // -(1/2) w I_perp(B_1 + B_1^T) for N
      pkpm_hm_put(up, PS[a][c], -0.5, cell, PKPM_HS_SNN(a, c, 1));
    }
    // -b_a C_c (powers 0 and 1)
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_mul(up, C0[c], 3, b[a], e1, 2);
      pkpm_hm_put(up, e1, -1.0, cell, PKPM_HS_SMM(a, c, 0));
      pkpm_hm_put(up, e1, -1.0, cell, PKPM_HS_SNN(a, c, 0));
      pkpm_hm_mul(up, C1[c], 3, b[a], e1, 2);
      pkpm_hm_put(up, e1, -1.0, cell, PKPM_HS_SMM(a, c, 1));
      pkpm_hm_put(up, e1, -1.0, cell, PKPM_HS_SNN(a, c, 1));
    }
  }
  // gyration (q/m) v x B: the field B (the pair terms carry q/m and the cross product)
  for (int e = 0; e < 3; ++e) {
    pkpm_hm_p1(up, &em[(3 + e) * nc], e1);
    pkpm_hm_put(up, e1, 1.0, cell, PKPM_HS_GYR(e));
  }
  // configuration-space flux couplings in direction d: xi e_d (F <- M0, G <- N), xi (e_d - b b_d),
  // 4 xi theta (e_d - b b_d) and -2 xi theta^2 (e_d - b b_d) (M0 <- G, N <- G, N <- F)
  for (int d = 0; d < cdim; ++d) {
    pkpm_hm_p1(up, xi, e1);
    pkpm_hm_put(up, e1, 1.0, cell, PKPM_HS_CF(d, 0));
    pkpm_hm_put(up, e1, 1.0, cell, PKPM_HS_CF(d, 1));
    pkpm_hm_p1(up, b[d], e2);
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_mul(up, e2, 1, b[a], e3, 2); // b_a b_d (complete)
      for (int n = 0; n < up->n5; ++n) {
        e3[n] = -e3[n];
      }
      if (a == d) {
        double one = 1.0; // the constant 1: (sqrt 2)^cdim Phat_0 ... Phat_0
        for (int c = 0; c < cdim; ++c) {
          one *= 1.4142135623730951;
        }
        e3[0] += one;
      }
      pkpm_hm_mul(up, e3, 2, xi, e1, 3); // xi (e_d - b b_d)_a (complete)
      pkpm_hm_put(up, e1, 1.0, cell, PKPM_HS_CF(d, 2 + a));
      double ex[PKPM_HM_N];
      pkpm_hm_mul(up, e1, 3, theta, ex, 3); // exact through 3
      pkpm_hm_put(up, ex, 4.0, cell, PKPM_HS_CF(d, 5 + a));
      pkpm_hm_mul(up, ex, 3, theta, e1, 2);
      pkpm_hm_put(up, e1, -2.0, cell, PKPM_HS_CF(d, 8 + a));
    }
  }
}

// The stored operator of a configuration cell (gkyl_pkpm_harm_ops.h) from its stored slots: every
// entry the sum of its slot terms, the gyration entries qbym B_e, the configuration-flux couplings.
GKYL_CU_DH static inline void
pkpm_hm_pack(const struct gkyl_dg_pkpm_harm_vol *up, const double *cell, double *ho)
{
  int n2 = up->n2;
  for (int e = 0; e < PKPM_HO_NENT; ++e) {
    double *out = &ho[e * n2];
    for (int k = 0; k < n2; ++k) {
      out[k] = 0.0;
    }
    for (int t = 0; t < up->ent_nterm[e]; ++t) {
      const double *c = &cell[up->ent_slot[e][t] * n2];
      double g = up->ent_gain[e][t];
      for (int k = 0; k < n2; ++k) {
        out[k] += g * c[k];
      }
    }
  }
  for (int e = 0; e < 3; ++e) {
    for (int k = 0; k < n2; ++k) {
      ho[PKPM_HO_GYR(e) * n2 + k] = up->qbym * cell[PKPM_HS_GYR(e) * n2 + k];
    }
  }
  for (int d = 0; d < up->cdim; ++d) {
    for (int q = 0; q < 11; ++q) {
      for (int k = 0; k < n2; ++k) {
        ho[PKPM_HO_CF(d, q) * n2 + k] = cell[PKPM_HS_CF(d, q) * n2 + k];
      }
    }
  }
}

// The source polynomials of the rate estimate at the cell center: poly[p][j] = the source entry of
// pair p with power j at the center (0 where absent).
GKYL_CU_DH static inline void
pkpm_hm_rate_poly(const struct gkyl_dg_pkpm_harm_vol *up, const double *ho, double (*poly)[3])
{
  int n2 = up->n2;
  for (int p = 0; p < up->nsrc; ++p) {
    for (int j = 0; j < 3; ++j) {
      double v = 0.0;
      if (up->src_ent[p][j] >= 0) {
        const double *c = &ho[up->src_ent[p][j] * n2];
        for (int k = 0; k < n2; ++k) {
          v += up->cen[k] * c[k];
        }
      }
      poly[p][j] = v;
    }
  }
}

// Rate estimate of the known first-harmonic sources of a phase cell (time-step bound): the largest
// absolute row sum (Gershgorin) of the source matrix at the configuration cell center, at the two
// w edges of the cell.
GKYL_CU_DH static inline double
pkpm_hm_rate(const struct gkyl_dg_pkpm_harm_vol *up, double (*poly)[3], double wlo, double wup)
{
  double rate = 0.0;
  for (int e = 0; e < 2; ++e) {
    double w = e == 0 ? wlo : wup, row[8] = {0.0};
    for (int p = 0; p < up->nsrc; ++p) {
      row[up->src_o[p]] += fabs(poly[p][0] + w * (poly[p][1] + w * poly[p][2]));
    }
    for (int o = 0; o < 8; ++o) {
      rate = fmax(rate, row[o]);
    }
  }
  return rate;
}
