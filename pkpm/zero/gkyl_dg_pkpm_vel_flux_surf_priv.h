// Private header: not for direct use
#pragma once

#include <assert.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_pkpm_harm_slots.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>

// Largest number of w-face nodes (5^cdim).
#define PKPM_VF_MAXNODE 125

typedef double (*pkpm_vel_flux_t)(
  const double *w, const double *dxv, const double *nv, const double *f_l, const double *f_r,
  double *GKYL_RESTRICT flux
);

struct gkyl_dg_pkpm_vel_flux_surf {
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  int cdim; // Configuration-space dimensions.
  int pdim; // Phase-space dimensions.
  int num_surf_basis; // Surface modal basis size per component.
  int nc; // Configuration basis size.
  int n2; // Degree <= 2 tensor modes of the stored coefficients (3^cdim).
  int nnode; // Face nodes (5^cdim).
  double e1[5][3]; // The 1D Legendre modes of degree <= 2 at the 1D face nodes.
  double ev[PKPM_VF_MAXNODE][8]; // The p=1 configuration basis at the face nodes.
  bool use_gpu; // bool to determine if on GPU.
  pkpm_vel_flux_t vel_flux; // Flux kernel on the v_par face.

  uint32_t flags;
  struct gkyl_dg_pkpm_vel_flux_surf *on_dev; // pointer to itself or device data
};

// Sizes and the node tables of the updater (host function shared by the host and device
// constructors).
void gkyl_dg_pkpm_vel_flux_surf_setup(
  struct gkyl_dg_pkpm_vel_flux_surf *up, const struct gkyl_dg_pkpm_vel_flux_surf_inp *inp
);

// v_par flux kernels (p=1 serendipity only; indexed by cdim-1)
GKYL_CU_D static const pkpm_vel_flux_t ser_pkpm_vel_flux_kernels[] = {
  pkpm_vel_flux_vpar_1x1v_ser_p1,
  pkpm_vel_flux_vpar_2x1v_ser_p1,
  pkpm_vel_flux_vpar_3x1v_ser_p1,
};

GKYL_CU_D
static pkpm_vel_flux_t
choose_pkpm_vel_flux_kern(enum gkyl_basis_type b_type, int cdim, int poly_order)
{
  assert(poly_order == 1);
  switch (b_type) {
    case GKYL_BASIS_MODAL_SERENDIPITY:
      return ser_pkpm_vel_flux_kernels[cdim - 1];
      break;
    default:
      assert(false);
      break;
  }
  return 0;
}

// Node values of the w flux of one configuration cell at the face nodes (PKPM_HS_NODE_NV per node):
// the w-flux slots of the stored coefficients, theta, b, alpha_0 and c. The face nodes are the tensor
// product of the 1D five-point nodes, node i = i_{cdim-1} + 5 i_{cdim-2} + 25 i_{cdim-3}, and the
// slots are evaluated there direction by direction (sum factorization; mode k = k_0 + 3 k_1 + 9 k_2).
GKYL_CU_DH static inline void
pkpm_vf_node_values(
  const struct gkyl_dg_pkpm_vel_flux_surf *up, const double *coef, const double *theta,
  const double *bvar, const double *alpha_0, const double *c_h, double *nv
)
{
  int cdim = up->cdim, nc = up->nc, n2 = up->n2;
  const int nn = PKPM_HS_NODE_NV;
  for (int s = 0; s < PKPM_HS_W_END; ++s) {
    const double *c = &coef[s * n2];
    if (cdim == 1) {
      for (int i0 = 0; i0 < 5; ++i0) {
        nv[i0 * nn + s] = up->e1[i0][0] * c[0] + up->e1[i0][1] * c[1] + up->e1[i0][2] * c[2];
      }
    }
    else if (cdim == 2) {
      double a[5][3];
      for (int i0 = 0; i0 < 5; ++i0) {
        for (int k1 = 0; k1 < 3; ++k1) {
          a[i0][k1] = up->e1[i0][0] * c[3 * k1] + up->e1[i0][1] * c[3 * k1 + 1] + up->e1[i0][2] * c[3 * k1 + 2];
        }
      }
      for (int i0 = 0; i0 < 5; ++i0) {
        for (int i1 = 0; i1 < 5; ++i1) {
          nv[(5 * i0 + i1) * nn + s] = up->e1[i1][0] * a[i0][0] + up->e1[i1][1] * a[i0][1] + up->e1[i1][2] * a[i0][2];
        }
      }
    }
    else {
      double a[5][3][3], b[5][5][3];
      for (int i0 = 0; i0 < 5; ++i0) {
        for (int k = 0; k < 9; ++k) {
          a[i0][k % 3][k / 3] = up->e1[i0][0] * c[3 * k] + up->e1[i0][1] * c[3 * k + 1] + up->e1[i0][2] * c[3 * k + 2];
        }
      }
      for (int i0 = 0; i0 < 5; ++i0) {
        for (int i1 = 0; i1 < 5; ++i1) {
          for (int k2 = 0; k2 < 3; ++k2) {
            b[i0][i1][k2] = up->e1[i1][0] * a[i0][0][k2] + up->e1[i1][1] * a[i0][1][k2] + up->e1[i1][2] * a[i0][2][k2];
          }
        }
      }
      for (int i0 = 0; i0 < 5; ++i0) {
        for (int i1 = 0; i1 < 5; ++i1) {
          for (int i2 = 0; i2 < 5; ++i2) {
            nv[(25 * i0 + 5 * i1 + i2) * nn + s] =
              up->e1[i2][0] * b[i0][i1][0] + up->e1[i2][1] * b[i0][i1][1] + up->e1[i2][2] * b[i0][i1][2];
          }
        }
      }
    }
  }
  for (int i = 0; i < up->nnode; ++i) {
    double *v = &nv[i * nn];
    const double *ev = up->ev[i];
    double th = 0.0, b0 = 0.0, b1 = 0.0, b2 = 0.0, a0 = 0.0, ch = 0.0;
    for (int a = 0; a < nc; ++a) {
      th += ev[a] * theta[a];
      b0 += ev[a] * bvar[a];
      b1 += ev[a] * bvar[nc + a];
      b2 += ev[a] * bvar[2 * nc + a];
      a0 += ev[a] * alpha_0[a];
      ch += ev[a] * c_h[a];
    }
    v[PKPM_HS_NODE_THETA] = th;
    v[PKPM_HS_NODE_B] = b0;
    v[PKPM_HS_NODE_B + 1] = b1;
    v[PKPM_HS_NODE_B + 2] = b2;
    v[PKPM_HS_NODE_A0] = a0;
    v[PKPM_HS_NODE_CH] = ch;
  }
}
