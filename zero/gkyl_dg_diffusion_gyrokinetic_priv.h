#pragma once

#include <gkyl_dg_diffusion_gyrokinetic.h>
#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>
#include <gkyl_ref_count.h>

// private header for use in diffusion DG equation object creation
// functions

static inline int diffdirs_linidx(const bool *isdirdiff, int cdim) {
  // Compute the linear index into the array of volume kernels (one
  // kernel for each combination of diffusive directions).
  bool diff_in_dir[GKYL_MAX_CDIM];
  if (isdirdiff)
    for (int d=0; d<cdim; d++) diff_in_dir[d] = isdirdiff[d];
  else
    for (int d=0; d<cdim; d++) diff_in_dir[d] = true;

  // Linear index into list of volume kernels.
  int dirs_bin_key[] = {1,2,4,8,16,32}; // Binary: 000001, 000010, 000100, 001000, 010000, 100000.
  int dirs_linidx = 0; // Binary 000000.
  for (int d=0; d<cdim; d++) {
     if (diff_in_dir[d]) dirs_linidx = dirs_linidx | dirs_bin_key[d];
  }
  dirs_linidx -= 1;
  return dirs_linidx;
}

// Types for various kernels
typedef double (*diffusion_surf_t)(const double *w, const double *dx,
  const double *coeff, const double *jacobgeo_inv,
  const double *ql, const double *qc, const double *qr, double* GKYL_RESTRICT out);

typedef double (*diffusion_boundary_surf_t)(const double *w, const double *dx,
  const double *coeff, const double *jacobgeo_inv, int edge,
  const double *qSkin, const double *qEdge, double* GKYL_RESTRICT out);

struct dg_diffusion_gyrokinetic {
  struct gkyl_dg_eqn eqn;
  diffusion_surf_t surf[GKYL_MAX_CDIM];
  diffusion_boundary_surf_t boundary_surf[GKYL_MAX_CDIM];
  struct gkyl_range diff_range;
  struct gkyl_dg_diffusion_gyrokinetic_auxfields auxfields;
  bool const_coeff;
  bool diff_in_dir[GKYL_MAX_CDIM];
  double skip_cell_thresh;
  int num_basis;
};

#define _cfD(idx) diffusion->const_coeff? (const double *) gkyl_array_cfetch(diffusion->auxfields.D, 0) : (const double *) gkyl_array_cfetch(diffusion->auxfields.D, gkyl_range_idx(&diffusion->diff_range, idx))

#define _cfJacInv(idx) (const double *) gkyl_array_cfetch(diffusion->auxfields.jacobgeo_inv, gkyl_range_idx(&diffusion->diff_range, idx))

// for use in kernel tables
typedef struct { vol_termf_t kernels[7]; } gkyl_dg_diffusion_gyrokinetic_vol_kern_list_diffdir;
typedef struct { gkyl_dg_diffusion_gyrokinetic_vol_kern_list_diffdir list[2]; } gkyl_dg_diffusion_gyrokinetic_vol_kern_list_polyOrder;
typedef struct { gkyl_dg_diffusion_gyrokinetic_vol_kern_list_polyOrder list[3]; } gkyl_dg_diffusion_gyrokinetic_vol_kern_list;

// for use in kernel tables
typedef struct { diffusion_surf_t kernels[2]; } gkyl_dg_diffusion_gyrokinetic_surf_kernels_polyOrder;
typedef struct { gkyl_dg_diffusion_gyrokinetic_surf_kernels_polyOrder list[6]; } gkyl_dg_diffusion_gyrokinetic_surf_kern_list;

typedef struct { diffusion_boundary_surf_t kernels[2]; } gkyl_dg_diffusion_gyrokinetic_boundary_surf_kernels_polyOrder;
typedef struct { gkyl_dg_diffusion_gyrokinetic_boundary_surf_kernels_polyOrder list[6]; } gkyl_dg_diffusion_gyrokinetic_boundary_surf_kern_list;

// ............... Homogeneous (constant) diffusion coefficient ............... //

// Serendipity volume kernels
// Need to be separated like this for GPU build

// 1x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_1x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_1x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_1x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_1x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
// 1x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_1x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_1x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_1x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_1x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
// 1x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order6_vol_1x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order6_vol_1x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
// 2x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_constcoeff_diffdirsy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_constcoeff_diffdirsxy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
// 2x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_2x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_2x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_2x_ser_p1_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_2x_ser_p1_constcoeff_diffdirsy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_2x_ser_p1_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_2x_ser_p1_constcoeff_diffdirsxy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_2x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_2x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_2x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_2x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_2x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_2x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
// 2x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order6_vol_2x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order6_vol_2x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order6_vol_2x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order6_vol_2x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order6_vol_2x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order6_vol_2x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
// 3x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsxy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsxz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsyz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsxyz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsxz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsyz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsxyz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
// 3x 4th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsxy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsxz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsyz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsxyz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsxz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsyz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsxyz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
// 3x 6th order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsxy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsxz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsyz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsxyz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}

// Volume kernel list.
GKYL_CU_D
static const gkyl_dg_diffusion_gyrokinetic_vol_kern_list ser_vol_kernels_constcoeff[] = {
  // 1x
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_dg_diffusion_gyrokinetic_order2_vol_1x_ser_p1_constcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_gyrokinetic_order2_vol_1x_ser_p2_constcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
      // 4th order diffusion.
      {.list={
          {ker_dg_diffusion_gyrokinetic_order4_vol_1x_ser_p1_constcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_gyrokinetic_order4_vol_1x_ser_p2_constcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
      // 6th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_gyrokinetic_order6_vol_1x_ser_p2_constcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
    },
  },
  // 2x
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_constcoeff_diffdirsx,ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_constcoeff_diffdirsy,ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_constcoeff_diffdirsxy,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_constcoeff_diffdirsx,ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_constcoeff_diffdirsy,ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_constcoeff_diffdirsxy,NULL,NULL,NULL,NULL},},
      },
      // 4th order diffusion.
      {.list={
          {ker_dg_diffusion_gyrokinetic_order4_vol_2x_ser_p1_constcoeff_diffdirsx,ker_dg_diffusion_gyrokinetic_order4_vol_2x_ser_p1_constcoeff_diffdirsy,ker_dg_diffusion_gyrokinetic_order4_vol_2x_ser_p1_constcoeff_diffdirsxy,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_gyrokinetic_order4_vol_2x_ser_p2_constcoeff_diffdirsx,ker_dg_diffusion_gyrokinetic_order4_vol_2x_ser_p2_constcoeff_diffdirsy,ker_dg_diffusion_gyrokinetic_order4_vol_2x_ser_p2_constcoeff_diffdirsxy,NULL,NULL,NULL,NULL},},
      },
      // 6th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_gyrokinetic_order6_vol_2x_ser_p2_constcoeff_diffdirsx,ker_dg_diffusion_gyrokinetic_order6_vol_2x_ser_p2_constcoeff_diffdirsy,ker_dg_diffusion_gyrokinetic_order6_vol_2x_ser_p2_constcoeff_diffdirsxy,NULL,NULL,NULL,NULL},},
      },
    },
  },
  // 3x
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsx,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsy,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsxy,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsz,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsxz,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsyz,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_constcoeff_diffdirsxyz},
          {ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsx,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsy,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsxy,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsz,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsxz,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsyz,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_constcoeff_diffdirsxyz},},
      },
      // 4th order diffusion.
      {.list={
          {ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsx,ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsy,ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsxy,ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsz,ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsxz,ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsyz,ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p1_constcoeff_diffdirsxyz},
          {ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsx,ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsy,ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsxy,ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsz,ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsxz,ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsyz,ker_dg_diffusion_gyrokinetic_order4_vol_3x_ser_p2_constcoeff_diffdirsxyz},},
      },
      // 6th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsx,ker_dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsy,ker_dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsxy,ker_dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsz,ker_dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsxz,ker_dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsyz,ker_dg_diffusion_gyrokinetic_order6_vol_3x_ser_p2_constcoeff_diffdirsxyz},},
      },
    },
  },
};

// Surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_gyrokinetic_surf_kern_list ser_gyrokinetic_surfx_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_gyrokinetic_order2_surfx_1x1v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order2_surfx_1x1v_ser_p2_constcoeff },
      { dg_diffusion_gyrokinetic_order2_surfx_1x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order2_surfx_1x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_surfx_2x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order2_surfx_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_surfx_3x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order2_surfx_3x2v_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_gyrokinetic_order4_surfx_1x1v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order4_surfx_1x1v_ser_p2_constcoeff },
      { dg_diffusion_gyrokinetic_order4_surfx_1x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order4_surfx_1x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order4_surfx_2x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order4_surfx_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order4_surfx_3x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order4_surfx_3x2v_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_gyrokinetic_order6_surfx_1x1v_ser_p2_constcoeff },
      { NULL, dg_diffusion_gyrokinetic_order6_surfx_1x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_gyrokinetic_order6_surfx_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_gyrokinetic_order6_surfx_3x2v_ser_p2_constcoeff },},
  },
};
// Surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_gyrokinetic_surf_kern_list ser_gyrokinetic_surfy_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_surfy_2x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order2_surfy_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_surfy_3x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order2_surfy_3x2v_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order4_surfy_2x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order4_surfy_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order4_surfy_3x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order4_surfy_3x2v_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_gyrokinetic_order6_surfy_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_gyrokinetic_order6_surfy_3x2v_ser_p2_constcoeff },},
  },
};
// Surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_gyrokinetic_surf_kern_list ser_gyrokinetic_surfz_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_surfz_3x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order2_surfz_3x2v_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order4_surfz_3x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order4_surfz_3x2v_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_gyrokinetic_order6_surfz_3x2v_ser_p2_constcoeff },},
  },
};

// Boundary surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_gyrokinetic_boundary_surf_kern_list ser_gyrokinetic_boundary_surfx_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_gyrokinetic_order2_boundary_surfx_1x1v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order2_boundary_surfx_1x1v_ser_p2_constcoeff },
      { dg_diffusion_gyrokinetic_order2_boundary_surfx_1x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order2_boundary_surfx_1x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_boundary_surfx_2x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order2_boundary_surfx_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_boundary_surfx_3x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order2_boundary_surfx_3x2v_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { dg_diffusion_gyrokinetic_order4_boundary_surfx_1x1v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order4_boundary_surfx_1x1v_ser_p2_constcoeff },
      { dg_diffusion_gyrokinetic_order4_boundary_surfx_1x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order4_boundary_surfx_1x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order4_boundary_surfx_2x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order4_boundary_surfx_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order4_boundary_surfx_3x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order4_boundary_surfx_3x2v_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, dg_diffusion_gyrokinetic_order6_boundary_surfx_1x1v_ser_p2_constcoeff },
      { NULL, dg_diffusion_gyrokinetic_order6_boundary_surfx_1x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_gyrokinetic_order6_boundary_surfx_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_gyrokinetic_order6_boundary_surfx_3x2v_ser_p2_constcoeff },},
  },
};
// Boundary surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_gyrokinetic_boundary_surf_kern_list ser_gyrokinetic_boundary_surfy_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_boundary_surfy_2x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order2_boundary_surfy_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_boundary_surfy_3x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order2_boundary_surfy_3x2v_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order4_boundary_surfy_2x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order4_boundary_surfy_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order4_boundary_surfy_3x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order4_boundary_surfy_3x2v_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_gyrokinetic_order6_boundary_surfy_2x2v_ser_p2_constcoeff },
      { NULL, NULL },
      { NULL, dg_diffusion_gyrokinetic_order6_boundary_surfy_3x2v_ser_p2_constcoeff },},
  },
};
// Boundary surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_gyrokinetic_boundary_surf_kern_list ser_gyrokinetic_boundary_surfz_kernels_constcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_boundary_surfz_3x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order2_boundary_surfz_3x2v_ser_p2_constcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order4_boundary_surfz_3x2v_ser_p1_constcoeff, dg_diffusion_gyrokinetic_order4_boundary_surfz_3x2v_ser_p2_constcoeff },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, dg_diffusion_gyrokinetic_order6_boundary_surfz_3x2v_ser_p2_constcoeff },},
  },
};

// ............... Inhomogeneous (spatially varying) diffusion coefficient ............... //

// Serendipity volume kernels
// Need to be separated like this for GPU build

// 1x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_1x_ser_p1_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_1x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_1x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_1x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
// 2x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_varcoeff_diffdirsy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_varcoeff_diffdirsxy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_varcoeff_diffdirsy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_varcoeff_diffdirsxy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
// 3x 2nd order diffusion.
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsxy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsxz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsyz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsxyz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsx(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsx(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsxy(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsxy(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsxz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsxz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsyz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}
GKYL_CU_DH static double ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsxyz(const struct gkyl_dg_eqn *eqn,
  const double* xc, const double* dx, const int* idx, const double* qIn, double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  return dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsxyz(xc, dx, _cfD(idx), _cfJacInv(idx), qIn, qRhsOut);
}

// Volume kernel list.
GKYL_CU_D
static const gkyl_dg_diffusion_gyrokinetic_vol_kern_list ser_vol_kernels_varcoeff[] = {
  // 1x
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_dg_diffusion_gyrokinetic_order2_vol_1x_ser_p1_varcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_gyrokinetic_order2_vol_1x_ser_p2_varcoeff_diffdirsx,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
      // 4th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
      // 6th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
    },
  },
  // 2x
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_varcoeff_diffdirsx,ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_varcoeff_diffdirsy,ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p1_varcoeff_diffdirsxy,NULL,NULL,NULL,NULL},
          {ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_varcoeff_diffdirsx,ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_varcoeff_diffdirsy,ker_dg_diffusion_gyrokinetic_order2_vol_2x_ser_p2_varcoeff_diffdirsxy,NULL,NULL,NULL,NULL},},
      },
      // 4th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
      // 6th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
    },
  },
  // 3x
  {.list={
      // 2nd order diffusion.
      {.list={
          {ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsx,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsy,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsxy,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsz,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsxz,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsyz,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p1_varcoeff_diffdirsxyz},
          {ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsx,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsy,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsxy,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsz,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsxz,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsyz,ker_dg_diffusion_gyrokinetic_order2_vol_3x_ser_p2_varcoeff_diffdirsxyz},},
      },
      // 4th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
      // 6th order diffusion.
      {.list={
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},
          {NULL,NULL,NULL,NULL,NULL,NULL,NULL},},
      },
    },
  },
};

// Surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_gyrokinetic_surf_kern_list ser_gyrokinetic_surfx_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_gyrokinetic_order2_surfx_1x1v_ser_p1_varcoeff, dg_diffusion_gyrokinetic_order2_surfx_1x1v_ser_p2_varcoeff },
      { dg_diffusion_gyrokinetic_order2_surfx_1x2v_ser_p1_varcoeff, dg_diffusion_gyrokinetic_order2_surfx_1x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_surfx_2x2v_ser_p1_varcoeff, dg_diffusion_gyrokinetic_order2_surfx_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_surfx_3x2v_ser_p1_varcoeff, dg_diffusion_gyrokinetic_order2_surfx_3x2v_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },},
  },
};
// Surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_gyrokinetic_surf_kern_list ser_gyrokinetic_surfy_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_surfy_2x2v_ser_p1_varcoeff, dg_diffusion_gyrokinetic_order2_surfy_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_surfy_3x2v_ser_p1_varcoeff, dg_diffusion_gyrokinetic_order2_surfy_3x2v_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },},
  },
};
// Surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_gyrokinetic_surf_kern_list ser_gyrokinetic_surfz_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_surfz_3x2v_ser_p1_varcoeff, dg_diffusion_gyrokinetic_order2_surfz_3x2v_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },},
  },
};

// Boundary surface kernel list: x-direction
GKYL_CU_D static const gkyl_dg_diffusion_gyrokinetic_boundary_surf_kern_list ser_gyrokinetic_boundary_surfx_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { dg_diffusion_gyrokinetic_order2_boundary_surfx_1x1v_ser_p1_varcoeff, dg_diffusion_gyrokinetic_order2_boundary_surfx_1x1v_ser_p2_varcoeff },
      { dg_diffusion_gyrokinetic_order2_boundary_surfx_1x2v_ser_p1_varcoeff, dg_diffusion_gyrokinetic_order2_boundary_surfx_1x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_boundary_surfx_2x2v_ser_p1_varcoeff, dg_diffusion_gyrokinetic_order2_boundary_surfx_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_boundary_surfx_3x2v_ser_p1_varcoeff, dg_diffusion_gyrokinetic_order2_boundary_surfx_3x2v_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },},
  },
};
// Boundary surface kernel list: y-direction
GKYL_CU_D static const gkyl_dg_diffusion_gyrokinetic_boundary_surf_kern_list ser_gyrokinetic_boundary_surfy_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_boundary_surfy_2x2v_ser_p1_varcoeff, dg_diffusion_gyrokinetic_order2_boundary_surfy_2x2v_ser_p2_varcoeff },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_boundary_surfy_3x2v_ser_p1_varcoeff, dg_diffusion_gyrokinetic_order2_boundary_surfy_3x2v_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },},
  },
};
// Boundary surface kernel list: z-direction
GKYL_CU_D static const gkyl_dg_diffusion_gyrokinetic_boundary_surf_kern_list ser_gyrokinetic_boundary_surfz_kernels_varcoeff[] = {
  // 2nd order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { dg_diffusion_gyrokinetic_order2_boundary_surfz_3x2v_ser_p1_varcoeff, dg_diffusion_gyrokinetic_order2_boundary_surfz_3x2v_ser_p2_varcoeff },},
  },
  // 4th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },},
  },
  // 6th order diffusion.
  {.list= {
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },
      { NULL, NULL },},
  },
};

#define SURFKERIDX(cdim,vdim) (cdim-1+vdim-1)*2-(vdim-1)

// Macro for choosing volume and surface kernels.
#define CKVOL(lst,cdim,diff_order,poly_order,diffdir_linidx) lst[cdim-1].list[diff_order/2-1].list[poly_order-1].kernels[diffdir_linidx]
#define CKSURF(lst,diff_order,cdim,vdim,poly_order) lst[diff_order/2-1].list[SURFKERIDX(cdim,vdim)].kernels[poly_order-1]

GKYL_CU_D static double surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcL, const double* xcC, const double* xcR, 
  const double* dxL, const double* dxC, const double* dxR,
  const int* idxL, const int* idxC, const int* idxR,
  const double* qInL, const double* qInC, const double* qInR,
  double* GKYL_RESTRICT qRhsOut)
{
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);

  if (fabs(qInL[0]) < diffusion->skip_cell_thresh && fabs(qInC[0]) < diffusion->skip_cell_thresh && fabs(qInR[0]) < diffusion->skip_cell_thresh) {
    return 0.;
  }
  
  if (diffusion->diff_in_dir[dir])
    diffusion->surf[dir](xcC, dxC, _cfD(idxC), _cfJacInv(idxC), qInL, qInC, qInR, qRhsOut);

  return 0.;  // CFL frequency computed in volume term.
}

GKYL_CU_D static double boundary_surf(const struct gkyl_dg_eqn* eqn, int dir,
  const double* xcEdge, const double* xcSkin, const double* dxEdge, const double* dxSkin,
  const int* idxEdge, const int* idxSkin, const int edge,
  const double* qInEdge, const double* qInSkin, double* GKYL_RESTRICT qRhsOut)
{ 
  struct dg_diffusion_gyrokinetic* diffusion = container_of(eqn, struct dg_diffusion_gyrokinetic, eqn);
  
  if (fabs(qInEdge[0]) < diffusion->skip_cell_thresh && fabs(qInSkin[0]) < diffusion->skip_cell_thresh) {
    return 0.;
  }
  if (diffusion->diff_in_dir[dir])
    diffusion->boundary_surf[dir](xcSkin, dxSkin, _cfD(idxSkin), _cfJacInv(idxSkin), edge, qInSkin, qInEdge, qRhsOut);

  return 0.;  // CFL frequency computed in volume term.
}

#undef _cfD
#undef _cfJacInv

/**
 * Free diffusion equation object
 *
 * @param ref Reference counter for constant diffusion equation
 */
void gkyl_dg_diffusion_gyrokinetic_free(const struct gkyl_ref_count* ref);

#ifdef GKYL_HAVE_CUDA
/**
 * Create a new gyrokinetic diffusion equation object on the device.
 *
 * @param basis Basis functions of the equation system.
 * @param cbasis Configuration space basis.
 * @param is_diff_constant If diffusion coefficient spatially constant.
 * @param diff_in_dir Whether to apply diffusion in each direction.
 * @param diff_order Diffusion order.
 * @param diff_range Range object to index the diffusion coefficient.
 * @param skip_cell_threshold Threshold which to skip cells
 * @return Pointer to diffusion equation object
 */
struct gkyl_dg_eqn*
gkyl_dg_diffusion_gyrokinetic_cu_dev_new(const struct gkyl_basis *basis, const struct gkyl_basis *cbasis,
  bool is_diff_const, const bool *diff_in_dir, int diff_order, const struct gkyl_range *diff_range, double skip_cell_threshold);
#endif
