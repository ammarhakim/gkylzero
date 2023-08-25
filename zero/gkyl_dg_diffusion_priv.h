#pragma once

#include <gkyl_dg_diffusion_kernels.h>
#include <gkyl_ref_count.h>
#include <gkyl_dg_eqn.h>

// private header for use in diffusion DG equation object creation
// functions

// Types for various kernels
typedef double (*diffusion_surf_t)(const double *w, const double *dx,
  const double *coeff, const double *ql, const double *qc, const double *qr,
  double* GKYL_RESTRICT out);

typedef double (*diffusion_boundary_surf_t)(const double *w, const double *dx,
  const double *coeff, int edge, const double *fSkin, const double *fEdge,
  double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vol_termf_t kernels[7]; } gkyl_dg_diffusion_vol_kern_list_diffdir;
typedef struct { gkyl_dg_diffusion_vol_kern_list_diffdir list[2]; } gkyl_dg_diffusion_vol_kern_list_polyOrder;
typedef struct { gkyl_dg_diffusion_vol_kern_list_polyOrder list[3]; } gkyl_dg_diffusion_vol_kern_list;

typedef struct { diffusion_surf_t kernels[2]; } gkyl_dg_diffusion_surf_kernels_polyOrder;
typedef struct { gkyl_dg_diffusion_surf_kernels_polyOrder list[3]; } gkyl_dg_diffusion_surf_kern_list;

typedef struct { diffusion_boundary_surf_t kernels[2]; } gkyl_dg_diffusion_boundary_surf_kernels_polyOrder;
typedef struct { gkyl_dg_diffusion_boundary_surf_kernels_polyOrder list[3]; } gkyl_dg_diffusion_boundary_surf_kern_list;

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

#define _cfD(loc) diffusion->const_coeff? gkyl_array_cfetch(diffusion->auxfields.D, 0) : gkyl_array_cfetch(diffusion->auxfields.D, loc)

#define CKVOL(lst,cdim,diff_order,poly_order,diffdir_linidx) lst[cdim-1].list[diff_order/2-1].list[poly_order-1].kernels[diffdir_linidx]

#define CKSURF(lst,diff_order,cdim,poly_order) lst[diff_order/2-1].list[cdim-1].kernels[poly_order-1]

// ............... Homogeneous (constant) diffusion coefficient ............... //

// Serendipity volume kernels
// Need to be separated like this for GPU build

// 1x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_1x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_1x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_1x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_1x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 1x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_1x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_1x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_1x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_1x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 1x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_1x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_1x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 2x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_2x_ser_p1_constcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_2x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 2x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_2x_ser_p1_constcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_2x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 2x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_2x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 3x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p1_constcoeff_diffdirsxyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p2_constcoeff_diffdirsxyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 3x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p1_constcoeff_diffdirsxyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p2_constcoeff_diffdirsxyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 3x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_3x_ser_p2_constcoeff_diffdirsxyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
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

// Surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_surfx_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_surfx_1x_ser_p1_constcoeff, dg_diffusion_order2_surfx_1x_ser_p2_constcoeff },
      { dg_diffusion_order2_surfx_2x_ser_p1_constcoeff, dg_diffusion_order2_surfx_2x_ser_p2_constcoeff },
      { dg_diffusion_order2_surfx_3x_ser_p1_constcoeff, dg_diffusion_order2_surfx_3x_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_surfx_1x_ser_p1_constcoeff, dg_diffusion_order4_surfx_1x_ser_p2_constcoeff },
      { dg_diffusion_order4_surfx_2x_ser_p1_constcoeff, dg_diffusion_order4_surfx_2x_ser_p2_constcoeff },
      { dg_diffusion_order4_surfx_3x_ser_p1_constcoeff, dg_diffusion_order4_surfx_3x_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_surfx_1x_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_surfx_2x_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_surfx_3x_ser_p2_constcoeff },},
  },
};
// Surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_surfy_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { dg_diffusion_order2_surfy_2x_ser_p1_constcoeff, dg_diffusion_order2_surfy_2x_ser_p2_constcoeff },
      { dg_diffusion_order2_surfy_3x_ser_p1_constcoeff, dg_diffusion_order2_surfy_3x_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { dg_diffusion_order4_surfy_2x_ser_p1_constcoeff, dg_diffusion_order4_surfy_2x_ser_p2_constcoeff },
      { dg_diffusion_order4_surfy_3x_ser_p1_constcoeff, dg_diffusion_order4_surfy_3x_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, dg_diffusion_order6_surfy_2x_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_surfy_3x_ser_p2_constcoeff },},
  },
};
// Surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_surfz_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_surfz_3x_ser_p1_constcoeff, dg_diffusion_order2_surfz_3x_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_surfz_3x_ser_p1_constcoeff, dg_diffusion_order4_surfz_3x_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_surfz_3x_ser_p2_constcoeff },},
  },
};

// Boundary surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_boundary_surfx_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_boundary_surfx_1x_ser_p1_constcoeff, dg_diffusion_order2_boundary_surfx_1x_ser_p2_constcoeff },
      { dg_diffusion_order2_boundary_surfx_2x_ser_p1_constcoeff, dg_diffusion_order2_boundary_surfx_2x_ser_p2_constcoeff },
      { dg_diffusion_order2_boundary_surfx_3x_ser_p1_constcoeff, dg_diffusion_order2_boundary_surfx_3x_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_boundary_surfx_1x_ser_p1_constcoeff, dg_diffusion_order4_boundary_surfx_1x_ser_p2_constcoeff },
      { dg_diffusion_order4_boundary_surfx_2x_ser_p1_constcoeff, dg_diffusion_order4_boundary_surfx_2x_ser_p2_constcoeff },
      { dg_diffusion_order4_boundary_surfx_3x_ser_p1_constcoeff, dg_diffusion_order4_boundary_surfx_3x_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_boundary_surfx_1x_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_boundary_surfx_2x_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_boundary_surfx_3x_ser_p2_constcoeff },},
  },
};
// Boundary surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_boundary_surfy_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { dg_diffusion_order2_boundary_surfy_2x_ser_p1_constcoeff, dg_diffusion_order2_boundary_surfy_2x_ser_p2_constcoeff },
      { dg_diffusion_order2_boundary_surfy_3x_ser_p1_constcoeff, dg_diffusion_order2_boundary_surfy_3x_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { dg_diffusion_order4_boundary_surfy_2x_ser_p1_constcoeff, dg_diffusion_order4_boundary_surfy_2x_ser_p2_constcoeff },
      { dg_diffusion_order4_boundary_surfy_3x_ser_p1_constcoeff, dg_diffusion_order4_boundary_surfy_3x_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, dg_diffusion_order6_boundary_surfy_2x_ser_p2_constcoeff },
      { NULL, dg_diffusion_order6_boundary_surfy_3x_ser_p2_constcoeff },},
  },
};
// Boundary surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_boundary_surfz_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_boundary_surfz_3x_ser_p1_constcoeff, dg_diffusion_order2_boundary_surfz_3x_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_boundary_surfz_3x_ser_p1_constcoeff, dg_diffusion_order4_boundary_surfz_3x_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_boundary_surfz_3x_ser_p2_constcoeff },},
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
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_1x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_1x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_1x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 1x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_1x_ser_p1_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_1x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_1x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_1x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 1x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_1x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_1x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 2x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_2x_ser_p1_varcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_2x_ser_p2_varcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 2x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_2x_ser_p1_varcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_2x_ser_p2_varcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 2x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_2x_ser_p2_varcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 3x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p1_varcoeff_diffdirsxyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order2_vol_3x_ser_p2_varcoeff_diffdirsxyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 3x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p1_varcoeff_diffdirsxyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order4_vol_3x_ser_p2_varcoeff_diffdirsxyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
// 3x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxy(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  long cidx = gkyl_range_idx(&diffusion->conf_range, idx);
  return dg_diffusion_order6_vol_3x_ser_p2_varcoeff_diffdirsxyz(xc, dx, _cfD(cidx), qIn, qRhsOut);
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

// Surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_surfx_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_surfx_1x_ser_p1_varcoeff, dg_diffusion_order2_surfx_1x_ser_p2_varcoeff },
      { dg_diffusion_order2_surfx_2x_ser_p1_varcoeff, dg_diffusion_order2_surfx_2x_ser_p2_varcoeff },
      { dg_diffusion_order2_surfx_3x_ser_p1_varcoeff, dg_diffusion_order2_surfx_3x_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_surfx_1x_ser_p1_varcoeff, dg_diffusion_order4_surfx_1x_ser_p2_varcoeff },
      { dg_diffusion_order4_surfx_2x_ser_p1_varcoeff, dg_diffusion_order4_surfx_2x_ser_p2_varcoeff },
      { dg_diffusion_order4_surfx_3x_ser_p1_varcoeff, dg_diffusion_order4_surfx_3x_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_surfx_1x_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_surfx_2x_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_surfx_3x_ser_p2_varcoeff },},
  },
};
// Surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_surfy_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { dg_diffusion_order2_surfy_2x_ser_p1_varcoeff, dg_diffusion_order2_surfy_2x_ser_p2_varcoeff },
      { dg_diffusion_order2_surfy_3x_ser_p1_varcoeff, dg_diffusion_order2_surfy_3x_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { dg_diffusion_order4_surfy_2x_ser_p1_varcoeff, dg_diffusion_order4_surfy_2x_ser_p2_varcoeff },
      { dg_diffusion_order4_surfy_3x_ser_p1_varcoeff, dg_diffusion_order4_surfy_3x_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, dg_diffusion_order6_surfy_2x_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_surfy_3x_ser_p2_varcoeff },},
  },
};
// Surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_surf_kern_list ser_surfz_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_surfz_3x_ser_p1_varcoeff, dg_diffusion_order2_surfz_3x_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_surfz_3x_ser_p1_varcoeff, dg_diffusion_order4_surfz_3x_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_surfz_3x_ser_p2_varcoeff },},
  },
};

// Boundary surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_boundary_surfx_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_order2_boundary_surfx_1x_ser_p1_varcoeff, dg_diffusion_order2_boundary_surfx_1x_ser_p2_varcoeff },
      { dg_diffusion_order2_boundary_surfx_2x_ser_p1_varcoeff, dg_diffusion_order2_boundary_surfx_2x_ser_p2_varcoeff },
      { dg_diffusion_order2_boundary_surfx_3x_ser_p1_varcoeff, dg_diffusion_order2_boundary_surfx_3x_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_order4_boundary_surfx_1x_ser_p1_varcoeff, dg_diffusion_order4_boundary_surfx_1x_ser_p2_varcoeff },
      { dg_diffusion_order4_boundary_surfx_2x_ser_p1_varcoeff, dg_diffusion_order4_boundary_surfx_2x_ser_p2_varcoeff },
      { dg_diffusion_order4_boundary_surfx_3x_ser_p1_varcoeff, dg_diffusion_order4_boundary_surfx_3x_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_order6_boundary_surfx_1x_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_boundary_surfx_2x_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_boundary_surfx_3x_ser_p2_varcoeff },},
  },
};
// Boundary surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_boundary_surfy_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { dg_diffusion_order2_boundary_surfy_2x_ser_p1_varcoeff, dg_diffusion_order2_boundary_surfy_2x_ser_p2_varcoeff },
      { dg_diffusion_order2_boundary_surfy_3x_ser_p1_varcoeff, dg_diffusion_order2_boundary_surfy_3x_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { dg_diffusion_order4_boundary_surfy_2x_ser_p1_varcoeff, dg_diffusion_order4_boundary_surfy_2x_ser_p2_varcoeff },
      { dg_diffusion_order4_boundary_surfy_3x_ser_p1_varcoeff, dg_diffusion_order4_boundary_surfy_3x_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, dg_diffusion_order6_boundary_surfy_2x_ser_p2_varcoeff },
      { NULL, dg_diffusion_order6_boundary_surfy_3x_ser_p2_varcoeff },},
  },
};
// Boundary surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_boundary_surf_kern_list ser_boundary_surfz_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order2_boundary_surfz_3x_ser_p1_varcoeff, dg_diffusion_order2_boundary_surfz_3x_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_order4_boundary_surfz_3x_ser_p1_varcoeff, dg_diffusion_order4_boundary_surfz_3x_ser_p2_varcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_order6_boundary_surfz_3x_ser_p2_varcoeff },},
  },
};

/**
 * Free diffusion equation object
 *
 * @param ref Reference counter for constant diffusion equation
 */
void gkyl_diffusion_free(const struct gkyl_ref_count* ref);

GKYL_CU_D static double surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR,
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  if (diffusion->diff_in_dir[dir]) {
    long cidx = gkyl_range_idx(&diffusion->conf_range, idxC);
    for (size_t c=0; c<diffusion->num_equations; c++) {
      int off = c*diffusion->num_basis;
      diffusion->surf[dir](xcC, dxC, _cfD(cidx), qInL+off, qInC+off, qInR+off, qRhsOut+off);
    }
  }
  return 0.;  // CFL frequency computed in volume term.
}

GKYL_CU_D static double boundary_surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcEdge, const double* xcSkin, const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{ 
  struct dg_diffusion* diffusion = container_of(eqn, struct dg_diffusion, eqn);
  
  if (diffusion->diff_in_dir[dir]) {
    long cidx = gkyl_range_idx(&diffusion->conf_range, idxSkin);
    for (size_t c=0; c<diffusion->num_equations; c++) {
      int off = c*diffusion->num_basis;
      diffusion->boundary_surf[dir](xcSkin, dxSkin, _cfD(cidx), edge, qInSkin+off, qInEdge+off, qRhsOut+off);
    }
  }
  return 0.;  // CFL frequency computed in volume term.
}

#undef _cfD
