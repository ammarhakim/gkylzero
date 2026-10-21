// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_dg_pkpm_harm_exp.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_util.h>

// Coefficients of one configuration node: (q/m) B_n, then the 3x3 matrices (row major) of the
// functions of z = h L_n applied to a vector: e^z, phi_1(z), the three final weights
// w_1 = phi_1 - 3 phi_2 + 4 phi_3, w_2 = -phi_2 + 4 phi_3, w_3 = 4 phi_2 - 8 phi_3, e^(z/2), phi_1(z/2).
// On the plane perpendicular to b_n, J v = v x b_n acts as i, so a function c(z) is
// Re c I_perp + Im c J; along b_n it is c(0) (the SSP-RK3 weights). These are the coefficients of the
// gyration alone; a phase node with a gyrophase-mixing rate nu > 0 (z = h(-nu + i |Omega|) on the
// perpendicular plane) builds the matrices it needs on the fly.
enum {
  PKPM_HX_OM = 0, // (q/m) B_n (3)
  PKPM_HX_E = 3,
  PKPM_HX_P1 = 12,
  PKPM_HX_W1 = 21,
  PKPM_HX_W2 = 30,
  PKPM_HX_W3 = 39,
  PKPM_HX_EH = 48,
  PKPM_HX_P1H = 57,
  PKPM_HX_NC = 66 // per configuration node
};

struct gkyl_dg_pkpm_harm_exp {
  int cdim; // Configuration-space dimensions.
  int num_conf_nodes; // 2^cdim.
  int num_phase_nodes; // 3 * 2^cdim.
  double qbym; // Charge to mass ratio.
  double mix; // Gyrophase-mixing rate in units of |Omega| (0: off).
  double half_dw; // Half the velocity cell width.
  bool use_gpu; // Whether the updater runs on the GPU.
  uint32_t flags;
  struct gkyl_dg_pkpm_harm_exp *on_dev; // pointer to itself or device data
};

// phi_k(z), k = 0..3, of a complex z: Taylor series for |z| < 1 (truncated where the next term is
// below round-off), the recurrence phi_(k+1) = (phi_k - 1/k!)/z otherwise.
GKYL_CU_DH static inline void
pkpm_hx_phis(double zr, double zi, double *pr, double *pi)
{
  const double ifac[24] = {
    1.0, 1.0, 1.0 / 2.0, 1.0 / 6.0, 1.0 / 24.0, 1.0 / 120.0, 1.0 / 720.0, 1.0 / 5040.0,
    1.0 / 40320.0, 1.0 / 362880.0, 1.0 / 3628800.0, 1.0 / 39916800.0, 1.0 / 479001600.0,
    1.0 / 6227020800.0, 1.0 / 87178291200.0, 1.0 / 1307674368000.0, 1.0 / 20922789888000.0,
    1.0 / 355687428096000.0, 1.0 / 6402373705728000.0, 1.0 / 121645100408832000.0,
    1.0 / 2432902008176640000.0, 1.0 / 51090942171709440000.0,
    1.0 / 1124000727777607680000.0, 1.0 / 25852016738884976640000.0
  };
  double az2 = zr * zr + zi * zi;
  if (az2 < 1.0) {
    // highest power kept: |z|^(J+1)/(J+1)! below 1e-17
    int J = az2 < 1e-6 ? 4 : (az2 < 1e-4 ? 6 : (az2 < 1e-2 ? 9 : (az2 < 0.25 ? 14 : 20)));
    for (int k = 0; k < 4; ++k) {
      // Horner on sum_(j=0)^J z^j/(j+k)!
      double ar = ifac[J + k], ai = 0.0;
      for (int j = J - 1; j >= 0; --j) {
        double tr = ar * zr - ai * zi + ifac[j + k];
        ai = ar * zi + ai * zr;
        ar = tr;
      }
      pr[k] = ar;
      pi[k] = ai;
    }
  }
  else {
    double ez = exp(zr);
    pr[0] = ez * cos(zi);
    pi[0] = ez * sin(zi);
    double iz = 1.0 / az2; // 1/z = conj(z)/|z|^2
    for (int k = 0; k < 3; ++k) {
      double ar = pr[k] - ifac[k], ai = pi[k];
      pr[k + 1] = (ar * zr + ai * zi) * iz;
      pi[k + 1] = (ai * zr - ar * zi) * iz;
    }
  }
}

// The matrix of c_par along oh plus (cr, ci) = c(z) on the perpendicular plane.
GKYL_CU_DH static inline void
pkpm_hx_node_mat(const double *oh, double cp, double cr, double ci, double *M)
{
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) {
      M[3 * a + b] = (cp - cr) * oh[a] * oh[b] + (a == b ? cr : 0.0);
    }
  }
  // + ci (v x oh)
  M[1] += ci * oh[2];
  M[2] -= ci * oh[1];
  M[3] -= ci * oh[2];
  M[5] += ci * oh[0];
  M[6] += ci * oh[1];
  M[7] -= ci * oh[0];
}

// |(q/m) B_n| and its direction (zero where B_n = 0).
GKYL_CU_DH static inline double
pkpm_hx_axis(const double *om, double *oh)
{
  double mag = sqrt(om[0] * om[0] + om[1] * om[1] + om[2] * om[2]);
  for (int a = 0; a < 3; ++a) {
    oh[a] = mag > 0.0 ? om[a] / mag : 0.0;
  }
  return mag;
}

// Coefficients of one configuration node from (q/m) B_n and the step h (the gyration alone).
GKYL_CU_DH static inline void
pkpm_hx_node_coef(const double *om, double h, double *c)
{
  double oh[3];
  double mag = pkpm_hx_axis(om, oh);
  for (int a = 0; a < 3; ++a) {
    c[PKPM_HX_OM + a] = om[a];
  }
  double pr[4], pi[4];
  pkpm_hx_phis(0.0, h * mag, pr, pi);
  pkpm_hx_node_mat(oh, 1.0, pr[0], pi[0], &c[PKPM_HX_E]);
  pkpm_hx_node_mat(oh, 1.0, pr[1], pi[1], &c[PKPM_HX_P1]);
  pkpm_hx_node_mat(oh, 1.0 / 6.0, pr[1] - 3.0 * pr[2] + 4.0 * pr[3], pi[1] - 3.0 * pi[2] + 4.0 * pi[3], &c[PKPM_HX_W1]);
  pkpm_hx_node_mat(oh, 1.0 / 6.0, -pr[2] + 4.0 * pr[3], -pi[2] + 4.0 * pi[3], &c[PKPM_HX_W2]);
  pkpm_hx_node_mat(oh, 2.0 / 3.0, 4.0 * pr[2] - 8.0 * pr[3], 4.0 * pi[2] - 8.0 * pi[3], &c[PKPM_HX_W3]);
  pkpm_hx_phis(0.0, 0.5 * h * mag, pr, pi);
  pkpm_hx_node_mat(oh, 1.0, pr[0], pi[0], &c[PKPM_HX_EH]);
  pkpm_hx_node_mat(oh, 1.0, pr[1], pi[1], &c[PKPM_HX_P1H]);
}

// The matrices of the slots k0 and k1 (k1 < 0: none) of a node damped at the rate nu on the
// perpendicular plane: z = h(-nu + i |Omega|), or z/2 for the slots e^(z/2) and phi_1(z/2) (k0 and k1
// share the argument). The component along b_n is not damped.
GKYL_CU_DH static inline void
pkpm_hx_damped_mats(const double *c, double nu, double h, int k0, int k1, double *M0, double *M1)
{
  double oh[3];
  double mag = pkpm_hx_axis(&c[PKPM_HX_OM], oh);
  double s = (k0 == PKPM_HX_EH || k0 == PKPM_HX_P1H) ? 0.5 * h : h;
  double pr[4], pi[4];
  pkpm_hx_phis(-s * nu, s * mag, pr, pi);
  for (int i = 0; i < 2; ++i) {
    int k = i == 0 ? k0 : k1;
    double *M = i == 0 ? M0 : M1;
    if (k < 0) {
      continue;
    }
    switch (k) {
      case PKPM_HX_E:
      case PKPM_HX_EH:
        pkpm_hx_node_mat(oh, 1.0, pr[0], pi[0], M);
        break;
      case PKPM_HX_P1:
      case PKPM_HX_P1H:
        pkpm_hx_node_mat(oh, 1.0, pr[1], pi[1], M);
        break;
      case PKPM_HX_W1:
        pkpm_hx_node_mat(oh, 1.0 / 6.0, pr[1] - 3.0 * pr[2] + 4.0 * pr[3], pi[1] - 3.0 * pi[2] + 4.0 * pi[3], M);
        break;
      case PKPM_HX_W2:
        pkpm_hx_node_mat(oh, 1.0 / 6.0, -pr[2] + 4.0 * pr[3], -pi[2] + 4.0 * pi[3], M);
        break;
      default: // PKPM_HX_W3
        pkpm_hx_node_mat(oh, 2.0 / 3.0, 4.0 * pr[2] - 8.0 * pr[3], 4.0 * pi[2] - 8.0 * pi[3], M);
        break;
    }
  }
}

// out += s M v.
GKYL_CU_DH static inline void
pkpm_hx_mv(const double *M, double s, const double *v, double *out)
{
  for (int a = 0; a < 3; ++a) {
    out[a] += s * (M[3 * a] * v[0] + M[3 * a + 1] * v[1] + M[3 * a + 2] * v[2]);
  }
}

// Gyrophase-mixing rate of one phase node: mix |Omega| (1 - 1/r)_+, with r the larger of the
// realizability ratios r_M = sqrt2 |M0_perp|/sqrt(F G) and r_N = |N_perp| sqrt(F)/G^(3/2), which are
// <= 1 when the one-harmonic reconstruction is non-negative (r = infinity where F or G <= 0). Capped at
// ln(r)/h, the rate that takes the unforced state to the bound in one step.
GKYL_CU_DH static inline double
pkpm_hx_rate(double mix, double mag, double h, double F, double G, const double *m0, const double *n)
{
  double rate = mix * mag;
  if (rate <= 0.0) {
    return 0.0;
  }
  if (F <= 0.0 || G <= 0.0) {
    return rate;
  }
  double am = sqrt(m0[0] * m0[0] + m0[1] * m0[1] + m0[2] * m0[2]);
  double an = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  double r = fmax(sqrt(2.0) * am / sqrt(F * G), an * sqrt(F) / (G * sqrt(G)));
  if (r <= 1.0) {
    return 0.0;
  }
  return fmin(rate * (1.0 - 1.0 / r), log(r) / h);
}

GKYL_CU_DH static inline void
pkpm_hx_to_nodes(int cdim, const double *f, double *fn)
{
  switch (cdim) {
    case 1: pkpm_harm_exp_to_nodes_1x1v_ser_p1(f, fn); break;
    case 2: pkpm_harm_exp_to_nodes_2x1v_ser_p1(f, fn); break;
    case 3: pkpm_harm_exp_to_nodes_3x1v_ser_p1(f, fn); break;
  }
}

GKYL_CU_DH static inline void
pkpm_hx_to_modal(int cdim, const double *fn, double *f)
{
  switch (cdim) {
    case 1: pkpm_harm_exp_to_modal_1x1v_ser_p1(fn, f); break;
    case 2: pkpm_harm_exp_to_modal_2x1v_ser_p1(fn, f); break;
    case 3: pkpm_harm_exp_to_modal_3x1v_ser_p1(fn, f); break;
  }
}

GKYL_CU_DH static inline void
pkpm_hx_fg_to_nodes(int cdim, const double *f, double *fgn)
{
  switch (cdim) {
    case 1: pkpm_harm_exp_fg_to_nodes_1x1v_ser_p1(f, fgn); break;
    case 2: pkpm_harm_exp_fg_to_nodes_2x1v_ser_p1(f, fgn); break;
    case 3: pkpm_harm_exp_fg_to_nodes_3x1v_ser_p1(f, fgn); break;
  }
}

GKYL_CU_DH static inline void
pkpm_hx_wint(int cdim, const double *f, double *wi)
{
  switch (cdim) {
    case 1: pkpm_harm_exp_wint_1x1v_ser_p1(f, wi); break;
    case 2: pkpm_harm_exp_wint_2x1v_ser_p1(f, wi); break;
    case 3: pkpm_harm_exp_wint_3x1v_ser_p1(f, wi); break;
  }
}

// Coefficients of one configuration cell.
GKYL_CU_DH static inline void
pkpm_hx_coef_cell(const struct gkyl_dg_pkpm_harm_exp *up, double h, const double *em, double *coef)
{
  double Bn[8 * 3];
  switch (up->cdim) {
    case 1: pkpm_harm_exp_conf_nodes_1x_ser_p1(em, Bn); break;
    case 2: pkpm_harm_exp_conf_nodes_2x_ser_p1(em, Bn); break;
    case 3: pkpm_harm_exp_conf_nodes_3x_ser_p1(em, Bn); break;
  }
  for (int n = 0; n < up->num_conf_nodes; ++n) {
    double om[3] = {up->qbym * Bn[3 * n], up->qbym * Bn[3 * n + 1], up->qbym * Bn[3 * n + 2]};
    pkpm_hx_node_coef(om, h, &coef[n * PKPM_HX_NC]);
  }
}

// Stage accumulation of one phase cell: N = rhs - U x (q/m) B_n at the nodes, then S and A. At stage
// 1 (U the state at the start of the step) it first sets the mixing rate nu of every node, and act to
// 1 when a node of the cell is damped (act may be shared by the cells of a column: every writer
// writes 1).
GKYL_CU_DH static inline void
pkpm_hx_accum_cell(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const double *coef, const double *fin,
  const double *rhs, double *nu, double *act, double *S, double *A
)
{
  double Un[24 * 6], Rn[24 * 6];
  pkpm_hx_to_nodes(up->cdim, fin, Un);
  pkpm_hx_to_nodes(up->cdim, rhs, Rn);
  if (stage == 1 && up->mix > 0.0) {
    double fgn[24 * 2];
    pkpm_hx_fg_to_nodes(up->cdim, fin, fgn);
    for (int p = 0; p < up->num_phase_nodes; ++p) {
      double oh[3];
      double mag = pkpm_hx_axis(&coef[(p / 3) * PKPM_HX_NC + PKPM_HX_OM], oh);
      nu[p] = pkpm_hx_rate(up->mix, mag, h, fgn[2 * p], fgn[2 * p + 1], &Un[6 * p], &Un[6 * p + 3]);
      if (nu[p] > 0.0) {
        act[0] = 1.0;
      }
    }
  }
  int wslot = stage == 1 ? PKPM_HX_W1 : (stage == 2 ? PKPM_HX_W2 : PKPM_HX_W3);
  for (int p = 0; p < up->num_phase_nodes; ++p) {
    const double *c = &coef[(p / 3) * PKPM_HX_NC];
    const double *om = &c[PKPM_HX_OM];
    double Md[9];
    const double *Mw = &c[wslot];
    if (up->mix > 0.0 && nu[p] > 0.0) {
      pkpm_hx_damped_mats(c, nu[p], h, wslot, -1, Md, 0);
      Mw = Md;
    }
    for (int v = 0; v < 2; ++v) {
      const double *u = &Un[6 * p + 3 * v];
      double *r = &Rn[6 * p + 3 * v];
      r[0] -= u[1] * om[2] - u[2] * om[1];
      r[1] -= u[2] * om[0] - u[0] * om[2];
      r[2] -= u[0] * om[1] - u[1] * om[0];
      double *s = &S[6 * p + 3 * v], *a = &A[6 * p + 3 * v];
      if (stage == 1) {
        for (int i = 0; i < 3; ++i) {
          s[i] = r[i];
          a[i] = 0.0;
        }
      }
      else if (stage == 2) {
        for (int i = 0; i < 3; ++i) {
          s[i] += r[i];
        }
      }
      pkpm_hx_mv(Mw, h, r, a);
    }
  }
}

// Stage output of one phase cell in the Shu-Osher form of the driver (see gkyl_dg_pkpm_harm_exp_state).
GKYL_CU_DH static inline void
pkpm_hx_state_cell(
  const struct gkyl_dg_pkpm_harm_exp *up, int stage, double h, const double *coef, const double *nu,
  const double *f0, const double *S, const double *A, double *out
)
{
  double Mn[24 * 6], Xn[24 * 6];
  pkpm_hx_to_nodes(up->cdim, f0, Mn);
  int k0 = stage == 2 ? PKPM_HX_EH : PKPM_HX_E;
  int k1 = stage == 1 ? PKPM_HX_P1 : (stage == 2 ? PKPM_HX_P1H : -1);
  for (int p = 0; p < up->num_phase_nodes; ++p) {
    const double *c = &coef[(p / 3) * PKPM_HX_NC];
    double Md0[9], Md1[9];
    const double *M0 = &c[k0], *M1 = k1 >= 0 ? &c[k1] : 0;
    if (up->mix > 0.0 && nu[p] > 0.0) {
      pkpm_hx_damped_mats(c, nu[p], h, k0, k1, Md0, Md1);
      M0 = Md0;
      M1 = Md1;
    }
    for (int v = 0; v < 2; ++v) {
      const double *m = &Mn[6 * p + 3 * v];
      double *x = &Xn[6 * p + 3 * v];
      x[0] = x[1] = x[2] = 0.0;
      if (stage == 1) {
        pkpm_hx_mv(M0, 1.0, m, x);
        pkpm_hx_mv(M1, h, &S[6 * p + 3 * v], x);
      }
      else if (stage == 2) {
        // X_2 = 4 U3 - 3 M
        pkpm_hx_mv(M0, 4.0, m, x);
        pkpm_hx_mv(M1, h, &S[6 * p + 3 * v], x);
        for (int i = 0; i < 3; ++i) {
          x[i] -= 3.0 * m[i];
        }
      }
      else {
        // X_3 = (3 M^(n+1) - M)/2
        pkpm_hx_mv(M0, 1.5, m, x);
        for (int i = 0; i < 3; ++i) {
          x[i] += 1.5 * A[6 * p + 3 * v + i] - 0.5 * m[i];
        }
      }
    }
  }
  pkpm_hx_to_modal(up->cdim, Xn, out);
}

// The integrals over w of F and M0_perp at the configuration nodes of one column (the phase cells
// f[0..num_cells-1] of a configuration cell, stride apart): col[node*4 + (0: F, 1-3: M0_perp)].
GKYL_CU_DH static inline void
pkpm_hx_col_add(const struct gkyl_dg_pkpm_harm_exp *up, const double *f, double *col)
{
  double wi[8 * 4];
  pkpm_hx_wint(up->cdim, f, wi);
  for (int i = 0; i < 4 * up->num_conf_nodes; ++i) {
    col[i] += up->half_dw * wi[i];
  }
}

// Momentum restoration of one phase cell: with m = int M0_perp dw / int F dw at the configuration
// nodes, M0_perp -= F m and N_perp -= 2 G m (the first-harmonic change of a perpendicular frame shift).
GKYL_CU_DH static inline void
pkpm_hx_restore_cell(const struct gkyl_dg_pkpm_harm_exp *up, const double *col, double *out)
{
  double Xn[24 * 6], fgn[24 * 2];
  pkpm_hx_to_nodes(up->cdim, out, Xn);
  pkpm_hx_fg_to_nodes(up->cdim, out, fgn);
  for (int p = 0; p < up->num_phase_nodes; ++p) {
    const double *cn = &col[(p / 3) * 4];
    if (cn[0] <= 0.0) {
      continue;
    }
    for (int a = 0; a < 3; ++a) {
      double m = cn[1 + a] / cn[0];
      Xn[6 * p + a] -= fgn[2 * p] * m;
      Xn[6 * p + 3 + a] -= 2.0 * fgn[2 * p + 1] * m;
    }
  }
  pkpm_hx_to_modal(up->cdim, Xn, out);
}
