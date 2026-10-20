// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_pkpm_harm_coeffs.h>
#include <gkyl_pkpm_harm_modal_tables.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>

// Tensor orthonormal-Legendre expansions in configuration space: mode (n0, n1, n2) at
// n0 + 5 n1 + 25 n2, degrees up to 4 per direction (n_d = 0 beyond cdim).
#define PKPM_HM_N 125
// Projections onto tensor degree <= 2 per direction (3^cdim modes, k0 + 3 k1 + 9 k2).
#define PKPM_HM_N2 27
// Nonzero (test, coefficient, distribution) configuration triples of int phi_m Q_k phi_l.
#define PKPM_HM_NTRIP 125
// Nonzero triples of int d_d phi_m Q_k phi_l per direction (phi_m linear in x_d).
#define PKPM_HM_NTRIPD 50

struct gkyl_dg_pkpm_harm_vol {
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  int cdim; // Configuration-space dimensions.
  int nc, np; // Configuration and phase basis sizes.
  int n2; // Number of degree <= 2 tensor modes (3^cdim).
  int n5; // Number of expansion slots in use (5^cdim; the higher ones belong to absent directions).
  double qbym; // Charge to mass ratio (gyration).
  int conf_deg[8 * 3]; // 1D Legendre degrees of each configuration basis function.
  int phase_of[8][3]; // Phase basis index of (configuration basis function, Legendre degree in w).
  int p2_idx[PKPM_HM_N2]; // Expansion index of each degree <= 2 mode.
  int ntrip; // Number of nonzero configuration triples.
  int trip_m[PKPM_HM_NTRIP], trip_k[PKPM_HM_NTRIP], trip_l[PKPM_HM_NTRIP];
  double trip_v[PKPM_HM_NTRIP]; // int phi_m Q_k phi_l over the reference cell.
  int ntripd[3]; // Number of nonzero derivative triples per direction.
  int tripd_m[3][PKPM_HM_NTRIPD], tripd_k[3][PKPM_HM_NTRIPD], tripd_l[3][PKPM_HM_NTRIPD];
  double tripd_v[3][PKPM_HM_NTRIPD]; // int d_d phi_m Q_k phi_l over the reference cell.
  const struct gkyl_vlasov_position_map *pos_map; // Position map (acquired; host pointer, 0 on the device copy).
  const struct gkyl_array *jacob_pos; // Per-direction Jacobians (borrowed from pos_map; device pointer on device).
  bool use_gpu; // Whether the updater runs on the GPU.
  uint32_t flags;
  struct gkyl_dg_pkpm_harm_vol *on_dev; // pointer to itself or device data
};

// The (output block, input block) pairs with nonzero coefficients.
#define PKPM_HM_NPAIR 57
GKYL_CU_TABLE __attribute__((unused)) static const int pkpm_hm_pairs[PKPM_HM_NPAIR][2] = {
  {0, 2}, {0, 3}, {0, 4}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}, {1, 7},
  {2, 0}, {3, 0}, {4, 0}, {2, 1}, {3, 1}, {4, 1},
  {2, 2}, {2, 3}, {2, 4}, {3, 2}, {3, 3}, {3, 4}, {4, 2}, {4, 3}, {4, 4},
  {2, 5}, {2, 6}, {2, 7}, {3, 5}, {3, 6}, {3, 7}, {4, 5}, {4, 6}, {4, 7},
  {5, 0}, {6, 0}, {7, 0}, {5, 1}, {6, 1}, {7, 1},
  {5, 2}, {5, 3}, {5, 4}, {6, 2}, {6, 3}, {6, 4}, {7, 2}, {7, 3}, {7, 4},
  {5, 5}, {5, 6}, {5, 7}, {6, 5}, {6, 6}, {6, 7}, {7, 5}, {7, 6}, {7, 7},
};

// Coupling pairs (output block, input block) of the configuration-space fluxes (plan 17.14):
// F <- xi M0_d, G <- xi N_d, M0_a <- xi G (e_d - b b_d)_a, N_a <- xi H (e_d - b b_d)_a with
// H = 4 theta G - 2 theta^2 F.
#define PKPM_HM_NCF 15
GKYL_CU_TABLE __attribute__((unused)) static const int pkpm_hm_cf_pairs[PKPM_HM_NCF][2] = {
  {0, 2}, {0, 3}, {0, 4}, {1, 5}, {1, 6}, {1, 7}, {2, 1}, {3, 1}, {4, 1},
  {5, 1}, {6, 1}, {7, 1}, {5, 0}, {6, 0}, {7, 0},
};

// Coefficients of the first-harmonic terms of one configuration cell: the known w flux
// Phi_o = sum_i phi[o][i] U_i and source S_o = sum_i src[o][i] U_i, each coefficient a polynomial
// sum_j w^j c_j(x) whose c_j is kept as its exact projection onto degree <= 2 per direction (all
// that the integrals against (test function) x (distribution), degree 2 per direction, see), and
// the configuration-space flux couplings cf[d][p] of the pairs pkpm_hm_cf_pairs in direction d
// (w independent, projected likewise).
struct pkpm_hm_coeffs {
  double phi[8][8][3][PKPM_HM_N2];
  double src[8][8][3][PKPM_HM_N2];
  double cf[3][PKPM_HM_NCF][PKPM_HM_N2];
  int jmax[PKPM_HM_NPAIR]; // Highest w power with a nonzero coefficient per pair (-1: none).
};

// Fill the updater (sizes, basis maps, degree <= 2 modes, nonzero configuration triples) from
// the input and the generated tables; host function shared by the host and device constructors.
void gkyl_dg_pkpm_harm_vol_setup(struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_dg_pkpm_harm_vol_inp *inp);


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

// out += coef P_wout[a f]: a an expansion exact through degree wa per direction, f a p=1
// configuration field (gkyl basis coefficients). Each basis function of f is a product of 1D
// Legendre factors, applied direction by direction (Phat_0 = 1/sqrt 2; Phat_1 by the recurrence);
// the result is exact through degree min(wout, wa - 1 + [a complete]) as the caller plans.
GKYL_CU_DH static inline void
pkpm_hm_mul_acc(
  const struct gkyl_dg_pkpm_harm_vol *up, const double *a, int wa, const double *f, double coef,
  double *out, int wout
)
{
  int cdim = up->cdim;
  const int stride[3] = {1, 5, 25};
  double t1[PKPM_HM_N], t2[PKPM_HM_N];
  for (int q = 0; q < up->nc; ++q) {
    if (f[q] == 0.0) {
      continue;
    }
    int lim[3] = {0, 0, 0};
    for (int d = 0; d < cdim; ++d) {
      lim[d] = wa;
    }
    for (int n = 0; n < up->n5; ++n) {
      t1[n] = a[n];
    }
    for (int d = 0; d < cdim; ++d) {
      int s = up->conf_deg[3 * q + d];
      int nlim[3] = {lim[0], lim[1], lim[2]};
      nlim[d] = lim[d] + s < wout ? lim[d] + s : wout;
      for (int n2 = 0; n2 <= nlim[2]; ++n2) {
        for (int n1 = 0; n1 <= nlim[1]; ++n1) {
          for (int n0 = 0; n0 <= nlim[0]; ++n0) {
            int n[3] = {n0, n1, n2}, idx = n0 + 5 * n1 + 25 * n2, m = n[d];
            double v = 0.0;
            if (s == 0) {
              v = m <= lim[d] ? 0.7071067811865475 * t1[idx] : 0.0;
            } else {
              if (m >= 1 && m - 1 <= lim[d]) {
                v += pkpm_harm_xp_up[m - 1] * t1[idx - stride[d]];
              }
              if (m + 1 <= lim[d]) {
                v += pkpm_harm_xp_dn[m + 1] * t1[idx + stride[d]];
              }
            }
            t2[idx] = v;
          }
        }
      }
      for (int n = 0; n < up->n5; ++n) {
        t1[n] = 0.0;
      }
      for (int n2 = 0; n2 <= nlim[2]; ++n2) {
        for (int n1 = 0; n1 <= nlim[1]; ++n1) {
          for (int n0 = 0; n0 <= nlim[0]; ++n0) {
            int idx = n0 + 5 * n1 + 25 * n2;
            t1[idx] = t2[idx];
          }
        }
      }
      lim[d] = nlim[d];
    }
    // directions not touched by f keep the window wa: truncate them to wout
    for (int n2 = 0; n2 <= lim[2] && n2 <= (cdim > 2 ? wout : 0); ++n2) {
      for (int n1 = 0; n1 <= lim[1] && n1 <= (cdim > 1 ? wout : 0); ++n1) {
        for (int n0 = 0; n0 <= lim[0] && n0 <= wout; ++n0) {
          int idx = n0 + 5 * n1 + 25 * n2;
          out[idx] += coef * f[q] * t1[idx];
        }
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

// dst[k] += coef (degree <= 2 projection of e)[k].
GKYL_CU_DH static inline void
pkpm_hm_p2_acc(const struct gkyl_dg_pkpm_harm_vol *up, const double *e, double coef, double *dst)
{
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
 * Coefficients of the first-harmonic terms of one configuration cell (definitions of
 * gkyl_pkpm_harm_coeffs.h, plan 17.11; all products exact, each coefficient kept as its projection
 * onto degree <= 2 per direction):
 *   flux:   F: t.M0; G: t.N; M0: alpha_k M0 + G t + R_B1 N; N: alpha_k N + H t + R_B1 O
 *   source: F: -ell_xi.M0; G: f_perp.M0 - (ell_xi + 2g).N;
 *           M0: f_perp F - w b.ell_xi M0 + L_perp M0 - G I_perp(ell_xi + g) - b(C.M0 + beta G);
 *           N: 2 G f_perp - w b.ell_xi N + L_perp N - R_E1 N - H I_perp(ell_xi + 3g) - b(C.N + beta H)
 * with H = 4 theta G - 2 theta^2 F, O = 6 theta N - 6 theta^2 M0.
 */
GKYL_CU_DH static inline void
pkpm_hm_build(
  const struct gkyl_dg_pkpm_harm_vol *up, const double *bvar, const double *grad_b,
  const double *u, const double *grad_u, const double *xi, const double *grad_xi,
  const double *grad_ln_n, const double *theta, const double *dtb, const double *em,
  struct pkpm_hm_coeffs *hc
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

  for (int o = 0; o < 8; ++o) {
    for (int i = 0; i < 8; ++i) {
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < PKPM_HM_N2; ++k) {
          hc->phi[o][i][j][k] = 0.0;
          hc->src[o][i][j][k] = 0.0;
        }
      }
    }
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
  for (int a = 0; a < 3; ++a) {
    // F: t.M0; G: t.N
    pkpm_hm_p2_acc(up, t0[a], 1.0, hc->phi[0][2 + a][0]);
    pkpm_hm_p2_acc(up, t1[a], 1.0, hc->phi[0][2 + a][1]);
    pkpm_hm_p2_acc(up, t0[a], 1.0, hc->phi[1][5 + a][0]);
    pkpm_hm_p2_acc(up, t1[a], 1.0, hc->phi[1][5 + a][1]);
    // alpha_k = -w K - w^2 b.g on the M0 and N diagonals
    pkpm_hm_p2_acc(up, K, -1.0, hc->phi[2 + a][2 + a][1]);
    pkpm_hm_p2_acc(up, bg, -1.0, hc->phi[2 + a][2 + a][2]);
    pkpm_hm_p2_acc(up, K, -1.0, hc->phi[5 + a][5 + a][1]);
    pkpm_hm_p2_acc(up, bg, -1.0, hc->phi[5 + a][5 + a][2]);
    // M0: G t
    pkpm_hm_p2_acc(up, t0[a], 1.0, hc->phi[2 + a][1][0]);
    pkpm_hm_p2_acc(up, t1[a], 1.0, hc->phi[2 + a][1][1]);
    // N: H t = (4 theta G - 2 theta^2 F) t
    pkpm_hm_mul(up, t0[a], 4, theta, e1, 3);
    pkpm_hm_p2_acc(up, e1, 4.0, hc->phi[5 + a][1][0]);
    pkpm_hm_mul(up, e1, 3, theta, e2, 2);
    pkpm_hm_p2_acc(up, e2, -2.0, hc->phi[5 + a][0][0]);
    pkpm_hm_mul(up, t1[a], 4, theta, e1, 3);
    pkpm_hm_p2_acc(up, e1, 4.0, hc->phi[5 + a][1][1]);
    pkpm_hm_mul(up, e1, 3, theta, e2, 2);
    pkpm_hm_p2_acc(up, e2, -2.0, hc->phi[5 + a][0][1]);
    // M0: R_B1 N; N: R_B1 O = 6 theta R N - 6 theta^2 R M0
    for (int c = 0; c < 3; ++c) {
      for (int n = 0; n < up->n5; ++n) {
        e3[n] = 0.5 * PS[a][c][n] + (a == c ? 0.5 * beta[n] : 0.0);
      }
      pkpm_hm_p2_acc(up, e3, 1.0, hc->phi[2 + a][5 + c][0]);
      pkpm_hm_mul(up, e3, 4, theta, e1, 3);
      pkpm_hm_p2_acc(up, e1, 6.0, hc->phi[5 + a][5 + c][0]);
      pkpm_hm_mul(up, e1, 3, theta, e2, 2);
      pkpm_hm_p2_acc(up, e2, -6.0, hc->phi[5 + a][2 + c][0]);
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
    pkpm_hm_p2_acc(up, xil[a], -1.0, hc->src[0][2 + a][0]);
    pkpm_hm_p2_acc(up, fp1[a], 1.0, hc->src[1][2 + a][1]);
    pkpm_hm_p2_acc(up, fp2[a], 1.0, hc->src[1][2 + a][2]);
    pkpm_hm_p2_acc(up, xil[a], -1.0, hc->src[1][5 + a][0]);
    pkpm_hm_p1(up, g[a], e1);
    pkpm_hm_p2_acc(up, e1, -2.0, hc->src[1][5 + a][0]);
    // M0: f_perp F; N: 2 G f_perp
    pkpm_hm_p2_acc(up, fp1[a], 1.0, hc->src[2 + a][0][1]);
    pkpm_hm_p2_acc(up, fp2[a], 1.0, hc->src[2 + a][0][2]);
    pkpm_hm_p2_acc(up, fp1[a], 2.0, hc->src[5 + a][1][1]);
    pkpm_hm_p2_acc(up, fp2[a], 2.0, hc->src[5 + a][1][2]);
    // M0: -G I_perp(ell_xi + g) - b beta G
    pkpm_hm_mul(up, beta, 2, b[a], e3, 3); // b_a beta (complete)
    pkpm_hm_p2_acc(up, plg[a], -1.0, hc->src[2 + a][1][0]);
    pkpm_hm_p2_acc(up, e3, -1.0, hc->src[2 + a][1][0]);
    // N: -H I_perp(ell_xi + 3 g) - b beta H, H = 4 theta G - 2 theta^2 F
    for (int n = 0; n < up->n5; ++n) {
      e2[n] = pl3g[a][n] + e3[n]; // I_perp(ell_xi + 3 g) + b beta: degree 4, complete
    }
    pkpm_hm_mul(up, e2, 4, theta, e1, 3);
    pkpm_hm_p2_acc(up, e1, -4.0, hc->src[5 + a][1][0]);
    pkpm_hm_mul(up, e1, 3, theta, e2, 2);
    pkpm_hm_p2_acc(up, e2, 2.0, hc->src[5 + a][0][0]);
    // diagonal of M0 and N: sigma_0 = -w b.ell_xi, -w b.g from L_perp; N also -R_E1's
    // (1/2)(eta + 2 w b.g) = (1/2)(div u - K) + (1/2) w beta + 2 w b.g
    pkpm_hm_p2_acc(up, xibl, -1.0, hc->src[2 + a][2 + a][1]);
    pkpm_hm_p2_acc(up, bg, -1.0, hc->src[2 + a][2 + a][1]);
    pkpm_hm_p2_acc(up, xibl, -1.0, hc->src[5 + a][5 + a][1]);
    pkpm_hm_p2_acc(up, bg, -3.0, hc->src[5 + a][5 + a][1]);
    pkpm_hm_p2_acc(up, beta, -0.5, hc->src[5 + a][5 + a][1]);
    pkpm_hm_p2_acc(up, K, 0.5, hc->src[5 + a][5 + a][0]);
    for (int c = 0; c < cdim; ++c) {
      pkpm_hm_p1(up, U[c][c], e1);
      pkpm_hm_p2_acc(up, e1, -0.5, hc->src[5 + a][5 + a][0]);
    }
  }
  // Matrix couplings M0 <- M0 and N <- N: L_perp = -I_perp U - w I_perp B_1 + (q/m)(. x B),
  // -b C. (both), and for N -(1/2) I_perp(U + U^T) - (1/2) w I_perp(B_1 + B_1^T)
  for (int c = 0; c < 3; ++c) {
    double col[3][PKPM_HM_N], pcol[3][PKPM_HM_N];
    // -I_perp U (column c): degree 3
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_p1(up, U[a][c], col[a]);
    }
    pkpm_hm_perp(up, b, col, 1, pcol, 2);
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_p2_acc(up, pcol[a], -1.0, hc->src[2 + a][2 + c][0]);
      pkpm_hm_p2_acc(up, pcol[a], -1.0, hc->src[5 + a][5 + c][0]);
    }
    // -(1/2) I_perp(U + U^T) for N
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_p1(up, U[a][c], col[a]);
      pkpm_hm_p1(up, U[c][a], e1);
      pkpm_hm_axpy(up, col[a], 1.0, e1);
    }
    pkpm_hm_perp(up, b, col, 1, pcol, 2);
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_p2_acc(up, pcol[a], -0.5, hc->src[5 + a][5 + c][0]);
    }
    // -w I_perp B_1 (column c): B_1 = xi grad b, degree 4
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_p1(up, Db[a][c], e1);
      pkpm_hm_mul(up, e1, 1, xi, col[a], 2);
    }
    pkpm_hm_perp(up, b, col, 2, pcol, 2);
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_p2_acc(up, pcol[a], -1.0, hc->src[2 + a][2 + c][1]);
      pkpm_hm_p2_acc(up, pcol[a], -1.0, hc->src[5 + a][5 + c][1]);
      // -(1/2) w I_perp(B_1 + B_1^T) for N
      pkpm_hm_p2_acc(up, PS[a][c], -0.5, hc->src[5 + a][5 + c][1]);
    }
    // -b_a C_c (powers 0 and 1)
    for (int a = 0; a < 3; ++a) {
      pkpm_hm_mul(up, C0[c], 3, b[a], e1, 2);
      pkpm_hm_p2_acc(up, e1, -1.0, hc->src[2 + a][2 + c][0]);
      pkpm_hm_p2_acc(up, e1, -1.0, hc->src[5 + a][5 + c][0]);
      pkpm_hm_mul(up, C1[c], 3, b[a], e1, 2);
      pkpm_hm_p2_acc(up, e1, -1.0, hc->src[2 + a][2 + c][1]);
      pkpm_hm_p2_acc(up, e1, -1.0, hc->src[5 + a][5 + c][1]);
    }
  }
  // configuration-space flux couplings in direction d: xi e_d (F <- M0, G <- N), xi (e_d - b b_d),
  // 4 xi theta (e_d - b b_d) and -2 xi theta^2 (e_d - b b_d) (M0 <- G, N <- G, N <- F)
  for (int d = 0; d < 3; ++d) {
    for (int p = 0; p < PKPM_HM_NCF; ++p) {
      for (int k = 0; k < PKPM_HM_N2; ++k) {
        hc->cf[d][p][k] = 0.0;
      }
    }
  }
  for (int d = 0; d < cdim; ++d) {
    pkpm_hm_p1(up, xi, e1);
    pkpm_hm_p2_acc(up, e1, 1.0, hc->cf[d][d]);
    pkpm_hm_p2_acc(up, e1, 1.0, hc->cf[d][3 + d]);
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
      pkpm_hm_p2_acc(up, e1, 1.0, hc->cf[d][6 + a]);
      double ex[PKPM_HM_N];
      pkpm_hm_mul(up, e1, 3, theta, ex, 3); // exact through 3
      pkpm_hm_p2_acc(up, ex, 4.0, hc->cf[d][9 + a]);
      pkpm_hm_mul(up, ex, 3, theta, e1, 2);
      pkpm_hm_p2_acc(up, e1, -2.0, hc->cf[d][12 + a]);
    }
  }
  // gyration (q/m) v x B: (v x B)_a = eps_{a c e} v_c B_e
  for (int a = 0; a < 3; ++a) {
    int c1 = (a + 1) % 3, c2 = (a + 2) % 3;
    // (v x B)_a = v_{a+1} B_{a+2} - v_{a+2} B_{a+1}
    pkpm_hm_p1(up, &em[(3 + c2) * nc], e1);
    pkpm_hm_p2_acc(up, e1, up->qbym, hc->src[2 + a][2 + c1][0]);
    pkpm_hm_p2_acc(up, e1, up->qbym, hc->src[5 + a][5 + c1][0]);
    pkpm_hm_p1(up, &em[(3 + c1) * nc], e1);
    pkpm_hm_p2_acc(up, e1, -up->qbym, hc->src[2 + a][2 + c2][0]);
    pkpm_hm_p2_acc(up, e1, -up->qbym, hc->src[5 + a][5 + c2][0]);
  }
  // highest w power in use per pair (the apply skips the others)
  for (int p = 0; p < PKPM_HM_NPAIR; ++p) {
    int o = pkpm_hm_pairs[p][0], i = pkpm_hm_pairs[p][1];
    hc->jmax[p] = -1;
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < up->n2; ++k) {
        if (hc->src[o][i][j][k] != 0.0 || hc->phi[o][i][j][k] != 0.0) {
          hc->jmax[p] = j;
        }
      }
    }
  }
}

/**
 * Volume integrals of one phase cell (center w_c, size dw in w) from the coefficients of its
 * configuration cell: out[o] += int psi_m S_o + dv1 int d_eta psi_m Phi_o, exact. The w part of
 * each coefficient sum_j w^j c_j with w = w_c + (dw/2) eta becomes sum_r (binomial weights) eta^r
 * against the 1D tables mw (int Phat_a eta^r Phat_b) and nw (int Phat_a' eta^r Phat_b); the
 * configuration part is the sum over the nonzero triples int phi_m Q_k phi_l.
 */
GKYL_CU_DH static inline void
pkpm_hm_apply(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct pkpm_hm_coeffs *hc, double wc, double dw,
  const double *f, double *out
)
{
  int nc = up->nc, np = up->np, n2 = up->n2;
  double h = 0.5 * dw, dv1 = 2.0 / dw;
  const double bw[3][3] = {{1.0, 0.0, 0.0}, {wc, h, 0.0}, {wc * wc, 2.0 * wc * h, h * h}};
  // distribution blocks contracted with the w tables: UM[i][r][mw][l], UN likewise
  double UM[8][3][3][8], UN[8][3][3][8];
  for (int i = 0; i < 8; ++i) {
    for (int r = 0; r < 3; ++r) {
      for (int mw = 0; mw < 3; ++mw) {
        for (int lc = 0; lc < nc; ++lc) {
          double sm = 0.0, sn = 0.0;
          for (int lw = 0; lw < 3; ++lw) {
            double ul = f[i * np + up->phase_of[lc][lw]];
            sm += pkpm_harm_mw[(r * 3 + mw) * 3 + lw] * ul;
            sn += pkpm_harm_nw[(r * 3 + mw) * 3 + lw] * ul;
          }
          UM[i][r][mw][lc] = sm;
          UN[i][r][mw][lc] = dv1 * sn;
        }
      }
    }
  }
  for (int p = 0; p < PKPM_HM_NPAIR; ++p) {
    int o = pkpm_hm_pairs[p][0], i = pkpm_hm_pairs[p][1], jm = hc->jmax[p];
    if (jm < 0) {
      continue;
    }
    double eS[3][PKPM_HM_N2], eF[3][PKPM_HM_N2];
    for (int r = 0; r <= jm; ++r) {
      for (int k = 0; k < n2; ++k) {
        double ss = 0.0, sf = 0.0;
        for (int j = r; j <= jm; ++j) {
          ss += bw[j][r] * hc->src[o][i][j][k];
          sf += bw[j][r] * hc->phi[o][i][j][k];
        }
        eS[r][k] = ss;
        eF[r][k] = sf;
      }
    }
    for (int t = 0; t < up->ntrip; ++t) {
      int mc = up->trip_m[t], k = up->trip_k[t], lc = up->trip_l[t];
      for (int mw = 0; mw < 3; ++mw) {
        double v = 0.0;
        for (int r = 0; r <= jm; ++r) {
          v += eS[r][k] * UM[i][r][mw][lc] + eF[r][k] * UN[i][r][mw][lc];
        }
        out[o * np + up->phase_of[mc][mw]] += up->trip_v[t] * v;
      }
    }
  }
}

/**
 * Volume integrals of the configuration-space flux couplings of one phase cell: out[o] +=
 * sum_d 2/(dx_d J_d) int d_d psi_m cf[d] U_i for the pairs (o, i) of pkpm_hm_cf_pairs; the w part
 * is the identity (the couplings are w independent), the configuration part the derivative
 * triples int d_d phi_m Q_k phi_l. jac holds the per-direction Jacobians J_d of the cell.
 */
GKYL_CU_DH static inline void
pkpm_hm_apply_cf(
  const struct gkyl_dg_pkpm_harm_vol *up, const struct pkpm_hm_coeffs *hc, const double *jac,
  const double *f, double *out
)
{
  int np = up->np;
  for (int d = 0; d < up->cdim; ++d) {
    double dfac = 2.0 / (up->phase_grid.dx[d] * jac[2 * d]);
    for (int p = 0; p < PKPM_HM_NCF; ++p) {
      int o = pkpm_hm_cf_pairs[p][0], i = pkpm_hm_cf_pairs[p][1];
      const double *c = hc->cf[d][p];
      for (int t = 0; t < up->ntripd[d]; ++t) {
        double cv = c[up->tripd_k[d][t]];
        if (cv == 0.0) {
          continue;
        }
        cv *= dfac * up->tripd_v[d][t];
        int mc = up->tripd_m[d][t], lc = up->tripd_l[d][t];
        for (int mw = 0; mw < 3; ++mw) {
          out[o * np + up->phase_of[mc][mw]] += cv * f[i * np + up->phase_of[lc][mw]];
        }
      }
    }
  }
}

// Rate estimate of the known first-harmonic sources of a phase cell (time-step bound): the
// estimator of gkyl_pkpm_harm_coeffs.h at the configuration cell center and the two w edges.
GKYL_CU_DH static inline double
pkpm_hm_rate(
  const struct gkyl_dg_pkpm_harm_vol *up, const double *bvar, const double *grad_b,
  const double *u, const double *grad_u, const double *xi, const double *grad_xi,
  const double *grad_ln_n, const double *theta, const double *dtb, const double *em, double wlo,
  double wup
)
{
  int cdim = up->cdim, nc = up->nc;
  double c0 = 1.0;
  for (int d = 0; d < cdim; ++d) {
    c0 *= 0.7071067811865475;
  }
  double b[3], uq[3], dq[3], B[3], gu[9], gb[9], gx[3], gl[3];
  for (int a = 0; a < 3; ++a) {
    b[a] = c0 * bvar[a * nc];
    uq[a] = c0 * u[a * nc];
    dq[a] = c0 * dtb[a * nc];
    B[a] = c0 * em[(3 + a) * nc];
  }
  for (int d = 0; d < cdim; ++d) {
    gx[d] = c0 * grad_xi[(d * 3) * nc];
    gl[d] = c0 * grad_ln_n[d * nc];
    for (int a = 0; a < 3; ++a) {
      gu[d * 3 + a] = c0 * grad_u[(d * 3 + a) * nc];
      gb[d * 3 + a] = c0 * grad_b[(d * 3 + a) * nc];
    }
  }
  struct pkpm_harm_coeffs hc;
  pkpm_harm_coeffs_set(&hc, cdim, b, uq, c0 * xi[0], c0 * theta[0], dq, gu, gb, gx, gl, up->qbym, B);
  return fmax(pkpm_harm_source_rate(&hc, wlo), pkpm_harm_source_rate(&hc, wup));
}
