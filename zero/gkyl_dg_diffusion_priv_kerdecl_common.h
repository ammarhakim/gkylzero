#pragma once

#include <gkyl_dg_eqn.h>

// Common objects used in private dg_diffusion header files.

// Types for various kernels
typedef double (*diffusion_surf_t)(const double *w, const double *dx,
  const double *coeff, const double *ql, const double *qc, const double *qr,
  double* GKYL_RESTRICT out);

typedef double (*diffusion_boundary_surf_t)(const double *w, const double *dx,
  const double *coeff, int edge, const double *fSkin, const double *fEdge,
  double* GKYL_RESTRICT out);

struct dg_diffusion {
  struct gkyl_dg_eqn eqn;
  diffusion_surf_t surf[GKYL_MAX_CDIM];
  diffusion_boundary_surf_t boundary_surf[GKYL_MAX_CDIM];
  struct gkyl_range conf_range;
  struct gkyl_dg_diffusion_auxfields auxfields;
  bool const_coeff;
  bool diff_in_dir[GKYL_MAX_CDIM];
  int num_equations, num_basis;
};

#define _cfD(idx) diffusion->const_coeff? (const double *) gkyl_array_cfetch(diffusion->auxfields.D, 0) : (const double *) gkyl_array_cfetch(diffusion->auxfields.D, gkyl_range_idx(&diffusion->conf_range, idx))

// for use in kernel tables
typedef struct { vol_termf_t kernels[7]; } gkyl_dg_diffusion_vol_kern_list_diffdir;
typedef struct { gkyl_dg_diffusion_vol_kern_list_diffdir list[2]; } gkyl_dg_diffusion_vol_kern_list_polyOrder;
typedef struct { gkyl_dg_diffusion_vol_kern_list_polyOrder list[3]; } gkyl_dg_diffusion_vol_kern_list;

// for use in kernel tables
typedef struct { diffusion_surf_t kernels[2]; } gkyl_dg_diffusion_surf_kernels_polyOrder;
typedef struct { gkyl_dg_diffusion_surf_kernels_polyOrder list[6]; } gkyl_dg_diffusion_surf_kern_list;

typedef struct { diffusion_boundary_surf_t kernels[2]; } gkyl_dg_diffusion_boundary_surf_kernels_polyOrder;
typedef struct { gkyl_dg_diffusion_boundary_surf_kernels_polyOrder list[6]; } gkyl_dg_diffusion_boundary_surf_kern_list;

// ............... Homogeneous (constant) diffusion coefficient ............... //

// Serendipity volume kernels
// Need to be separated like this for GPU build

// 1x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_1x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_1x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_1x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_1x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 1x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_1x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_1x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_1x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_1x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 1x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_1x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_1x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 2x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 2x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 2x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 3x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 3x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 3x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}

// Volume kernel list.
GKYL_CU_D
static const gkyl_dg_diffusion_vol_kern_list ser_vol_kernels_constcoeff[] = {
  // 1x
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_dg_diffusion_order2_vol_1x_ser_p1_constcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_order2_vol_1x_ser_p2_constcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
      // 4th order diffusion.
      {.list={
          {ker_dg_diffusion_order4_vol_1x_ser_p1_constcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_order4_vol_1x_ser_p2_constcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
      // 6th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_order6_vol_1x_ser_p2_constcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
    },
  },
  // 2x
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsx,ker_dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsy,ker_dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsxy,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsx,ker_dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsy,ker_dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsxy,NULL,NULL,NULL,NULL},},
      },
      // 4th order diffusion.
      {.list={
          {ker_dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsx,ker_dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsy,ker_dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsxy,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsx,ker_dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsy,ker_dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsxy,NULL,NULL,NULL,NULL},},
      },
      // 6th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsx,ker_dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsy,ker_dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsxy,NULL,NULL,NULL,NULL},},
      },
    },
  },
  // 3x
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsx,ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsy,ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxy,ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsz,ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxz,ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsyz,ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxyz},
          {ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsx,ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsy,ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxy,ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsz,ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxz,ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsyz,ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxyz},},
      },
      // 4th order diffusion.
      {.list={
          {ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsx,ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsy,ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxy,ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsz,ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxz,ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsyz,ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxyz},
          {ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsx,ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsy,ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxy,ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsz,ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxz,ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsyz,ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxyz},},
      },
      // 6th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsx,ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsy,ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxy,ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsz,ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxz,ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsyz,ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxyz},},
      },
    },
  },
};

// ............... Inhomogeneous (spatially varying) diffusion coefficient ............... //

// Serendipity volume kernels
// Need to be separated like this for GPU build

// 1x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_1x_ser_p1_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_1x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_1x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_1x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 1x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_1x_ser_p1_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_1x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_1x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_1x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 1x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_1x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_1x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 2x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 2x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 2x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 3x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 3x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
// 3x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxy(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  return dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxyz(xc, dx, _cfD(idx), qIn, qRhsOut);
}

// Volume kernel list.
GKYL_CU_D
static const gkyl_dg_diffusion_vol_kern_list ser_vol_kernels_varcoeff[] = {
  // 1x
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_dg_diffusion_order2_vol_1x_ser_p1_varcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_order2_vol_1x_ser_p2_varcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
      // 4th order diffusion.
      {.list={
          {ker_dg_diffusion_order4_vol_1x_ser_p1_varcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_order4_vol_1x_ser_p2_varcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
      // 6th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_order6_vol_1x_ser_p2_varcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
    },
  },
  // 2x
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsx,ker_dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsy,ker_dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsxy,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsx,ker_dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsy,ker_dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsxy,NULL,NULL,NULL,NULL},},
      },
      // 4th order diffusion.
      {.list={
          {ker_dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsx,ker_dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsy,ker_dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsxy,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsx,ker_dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsy,ker_dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsxy,NULL,NULL,NULL,NULL},},
      },
      // 6th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsx,ker_dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsy,ker_dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsxy,NULL,NULL,NULL,NULL},},
      },
    },
  },
  // 3x
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsx,ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsy,ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxy,ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsz,ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxz,ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsyz,ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxyz},
          {ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsx,ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsy,ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxy,ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsz,ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxz,ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsyz,ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxyz},},
      },
      // 4th order diffusion.
      {.list={
          {ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsx,ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsy,ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxy,ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsz,ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxz,ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsyz,ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxyz},
          {ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsx,ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsy,ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxy,ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsz,ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxz,ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsyz,ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxyz},},
      },
      // 6th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsx,ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsy,ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxy,ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsz,ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxz,ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsyz,ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxyz},},
      },
    },
  },
};

