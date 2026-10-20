#pragma once

#include <gkyl_util.h>

// Layout of the stored first-harmonic coefficients of a configuration cell (plan 17.15): every
// coefficient of the known terms of the eight blocks [F, G, M0_perp (x, y, z), N_perp (x, y, z)]
// is a polynomial sum_j w^j c_j(x), and each c_j is stored as its exact projection onto the tensor
// Legendre modes of degree <= 2 per configuration direction (3^cdim modes, k = k0 + 3 k1 + 9 k2).
// That projection is all that any Galerkin integral of the solver sees: in the volume and on the
// w faces the coefficient multiplies (test function) x (distribution trace), of degree <= 2 per
// configuration direction. A slot is one c_j; slot s occupies [s n2, (s + 1) n2) of the cell's
// coefficient array, n2 = 3^cdim.
//
// w flux (alpha_k U + T(U), the coefficients of d_w terms):
//   alpha_k = -K w - (b.grad xi) w^2 on every block (two shared slots);
//   T: F <- t.M0, G <- beta H + t.N (the two-field beta terms separately), M0 <- G t + R_B1 N,
//      N <- H t + R_B1 O with H = 4 theta G - 2 theta^2 F, O = 6 theta N - 6 theta^2 M0.
// Sources (o, i, j): F <- -ell_xi.M0; G <- f_perp.M0 - (ell_xi + 2 g).N; M0 <- f_perp F
//   + sigma_0 M0 + L_perp M0 - G I_perp(ell_xi + g) - b(C.M0 + beta G); N <- 2 G f_perp
//   + sigma_0 N + L_perp N - R_E1 N - H I_perp(ell_xi + 3 g) - b(C.N + beta H); the gyration
//   (q/m) v x B is kept as the field B itself (three slots).
// Configuration-flux couplings in direction d (w independent): F <- xi M0_d, G <- xi N_d,
//   M0 <- xi G (e_d - b b_d), N <- xi H (e_d - b b_d).

// w flux: shared diagonal alpha_k
#define PKPM_HS_KW 0 // -K (w^1)
#define PKPM_HS_BG 1 // -b.grad xi (w^2)
// w flux: first-harmonic couplings T (a, c = 0, 1, 2; j the power of w)
#define PKPM_HS_TF(a, j) (2 + 2 * (a) + (j)) // F <- M0_a: t_a (j = 0, 1)
#define PKPM_HS_TG(a, j) (8 + 2 * (a) + (j)) // G <- N_a: t_a
#define PKPM_HS_TM(a, j) (14 + 2 * (a) + (j)) // M0_a <- G: t_a
#define PKPM_HS_TNG(a, j) (20 + 2 * (a) + (j)) // N_a <- G: 4 theta t_a
#define PKPM_HS_TNF(a, j) (26 + 2 * (a) + (j)) // N_a <- F: -2 theta^2 t_a
#define PKPM_HS_RMN(a, c) (32 + 3 * (a) + (c)) // M0_a <- N_c: (R_B1)_ac (w^0)
#define PKPM_HS_RNN(a, c) (41 + 3 * (a) + (c)) // N_a <- N_c: 6 theta (R_B1)_ac
#define PKPM_HS_RNM(a, c) (50 + 3 * (a) + (c)) // N_a <- M0_c: -6 theta^2 (R_B1)_ac
// w flux: the two-field couplings of F and G (w^0)
#define PKPM_HS_BF 59 // F <- G: beta = xi div b
#define PKPM_HS_BGG 60 // G <- G: 4 theta beta
#define PKPM_HS_BGF 61 // G <- F: -2 theta^2 beta
// sources
#define PKPM_HS_SFM(a) (62 + (a)) // F <- M0_a (w^0): -ell_xi
#define PKPM_HS_SGM(a, j) (65 + 2 * (a) + (j) - 1) // G <- M0_a (w^1, w^2): f_perp
#define PKPM_HS_SGN(a) (71 + (a)) // G <- N_a (w^0): -(ell_xi + 2 g)
#define PKPM_HS_SMF(a, j) (74 + 2 * (a) + (j) - 1) // M0_a <- F (w^1, w^2): f_perp
#define PKPM_HS_SNG1(a, j) (80 + 2 * (a) + (j) - 1) // N_a <- G (w^1, w^2): 2 f_perp
#define PKPM_HS_SMG(a) (86 + (a)) // M0_a <- G (w^0): -I_perp(ell_xi + g) - b beta
#define PKPM_HS_SNG0(a) (89 + (a)) // N_a <- G (w^0): -4 theta (I_perp(ell_xi + 3 g) + b beta)
#define PKPM_HS_SNF(a) (92 + (a)) // N_a <- F (w^0): 2 theta^2 (I_perp(ell_xi + 3 g) + b beta)
#define PKPM_HS_SMM(a, c, j) (95 + 6 * (a) + 2 * (c) + (j)) // M0_a <- M0_c (w^0, w^1), no gyration
#define PKPM_HS_SNN(a, c, j) (113 + 6 * (a) + 2 * (c) + (j)) // N_a <- N_c (w^0, w^1), no gyration
#define PKPM_HS_GYR(e) (131 + (e)) // B_e of the gyration (q/m) v x B (w^0)
// configuration-flux couplings in direction d: q = 0 F <- M0_d, 1 G <- N_d, 2 + a M0_a <- G,
// 5 + a N_a <- G, 8 + a N_a <- F
#define PKPM_HS_CF(d, q) (134 + 11 * (d) + (q))
// number of slots and of degree <= 2 modes of a configuration dimension
#define PKPM_HS_NSLOT(cdim) (134 + 11 * (cdim))
#define PKPM_HS_NSLOT_MAX 167
#define PKPM_HS_W_END 62 // the w-flux slots are [0, 62)
#define PKPM_HS_SRC_END 134 // the source slots are [62, 134)

// Node values of the w flux at a configuration point (the w-face and closure kernels): the w-flux
// slots [0, PKPM_HS_W_END) evaluated there, then theta, b (three components), alpha_0 and c.
#define PKPM_HS_NODE_THETA 62
#define PKPM_HS_NODE_B 63
#define PKPM_HS_NODE_A0 66
#define PKPM_HS_NODE_CH 67
#define PKPM_HS_NODE_NV 68

