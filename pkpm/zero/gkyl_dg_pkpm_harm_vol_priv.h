// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_dg_pkpm_harm_vol.h>
#include <gkyl_pkpm_harm_ops.h>
#include <gkyl_pkpm_kinetic_kernels.h>
#include <gkyl_range.h>
#include <gkyl_rect_grid.h>
#include <gkyl_util.h>
#include <gkyl_vlasov_position_map.h>

// Tensor Legendre modes of degree <= 2 per configuration direction (3^cdim, k0 + 3 k1 + 9 k2).
#define PKPM_HM_N2 27
// (Output, input) block pairs of the source-rate estimate.
#define PKPM_HM_MAXSRC 64
// The nc x nc matrices of the stored operator of a configuration cell (pkpm_harm_mat_*).
#define PKPM_HM_MAT (PKPM_HO_NOP_MAX * 64)

// Generated kernels: the stored operator of a configuration cell from the p=1 fields, its
// matrices, and the volume terms of a phase cell from them.
typedef void (*pkpm_harm_build_t)(
  double qbym, const double *bvar, const double *grad_b, const double *u, const double *grad_u,
  const double *xi, const double *grad_xi, const double *grad_ln_n, const double *theta,
  const double *dtb, const double *em, double *GKYL_RESTRICT ho
);
typedef void (*pkpm_harm_mat_t)(
  const double *dxv, const double *jacob_pos, const double *ho, double *GKYL_RESTRICT mat
);
typedef void (*pkpm_harm_vol_t)(
  const double *w, const double *dxv, const double *mat, const double *f, double *GKYL_RESTRICT out
);

struct gkyl_dg_pkpm_harm_vol {
  struct gkyl_rect_grid phase_grid; // Phase-space grid.
  int cdim; // Configuration-space dimensions.
  int n2; // Number of degree <= 2 tensor modes (3^cdim).
  int nop; // Number of entries of the stored operator (gkyl_pkpm_harm_ops.h).
  double qbym; // Charge to mass ratio (gyration).
  double cen[PKPM_HM_N2]; // Value of each degree <= 2 mode at the cell center.
  // the (output, input) source pairs of the rate estimate: the source entry of each power of w
  // (-1 if absent); the gyration is not rated (it is the exponential part of the time integrator)
  int nsrc, src_o[PKPM_HM_MAXSRC], src_i[PKPM_HM_MAXSRC], src_ent[PKPM_HM_MAXSRC][3];
  pkpm_harm_build_t build_kern; // stored operator of a configuration cell
  pkpm_harm_mat_t mat_kern; // matrices of the stored operator of a configuration cell
  pkpm_harm_vol_t vol_kern; // volume terms of a phase cell
  const struct gkyl_vlasov_position_map *pos_map; // Position map (acquired; host pointer, 0 on the device copy).
  const struct gkyl_array *jacob_pos; // Per-direction Jacobians (borrowed from pos_map; device pointer on device).
  bool use_gpu; // Whether the updater runs on the GPU.
  uint32_t flags;
  struct gkyl_dg_pkpm_harm_vol *on_dev; // pointer to itself or device data
};

// Fill the updater (sizes, the degree <= 2 modes at the cell center, the source pairs of the
// rate) from the input and the generated tables; host function shared by the host and device
// constructors (the kernel pointers are set by each).
void gkyl_dg_pkpm_harm_vol_setup(struct gkyl_dg_pkpm_harm_vol *up, const struct gkyl_dg_pkpm_harm_vol_inp *inp);

// Generated kernels (p=1 serendipity; indexed by cdim-1).
GKYL_CU_D static const pkpm_harm_build_t ser_pkpm_harm_build_kernels[] = {
  pkpm_harm_build_1x_ser_p1,
  pkpm_harm_build_2x_ser_p1,
  pkpm_harm_build_3x_ser_p1,
};
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
