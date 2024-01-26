// Private header: not for direct use
#pragma once

#include <math.h>

#include <gkyl_array.h>
#include <gkyl_basis.h>
#include <gkyl_mom_cross_bgk_gyrokinetic_kernels.h>
#include <gkyl_mat.h>
#include <gkyl_range.h>
#include <gkyl_util.h>
#include <assert.h>

// Function pointer type for cross moments calculation
typedef void (*gyrokinetic_mom_cross_bgk_t)(const double beta, 
  const double m_self, const double *moms_self, 
  const double m_other, const double *moms_other, 
  const double *nu_sr, const double *nu_rs, double *moms_cross);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below.
GKYL_CU_D
static struct { int vdim[3]; } cv_index[] = {
  {-1, -1, -1}, // 0x makes no sense
  {-1,  0,  1}, // 1x kernel indices
  {-1, -1,  2}, // 2x kernel indices
  {-1, -1,  3}, // 3x kernel indices
};

// For use in kernel tables
typedef struct { gyrokinetic_mom_cross_bgk_t kernels[3]; } gkyl_mom_cross_bgk_gyrokinetic_kern_list;

// Cross moments kernel list
GKYL_CU_D
static const gkyl_mom_cross_bgk_gyrokinetic_kern_list ser_mom_cross_gyrokinetic_list[] = {
  // 1x kernels 
  { NULL, gyrokinetic_mom_cross_bgk_1x1v_ser_p1, gyrokinetic_mom_cross_bgk_1x1v_ser_p2}, //0
  { NULL, gyrokinetic_mom_cross_bgk_1x2v_ser_p1, gyrokinetic_mom_cross_bgk_1x2v_ser_p2}, //1
  // 2x kernels
  { NULL, gyrokinetic_mom_cross_bgk_2x2v_ser_p1, NULL}, // no gyrokinetic_mom_cross_bgk_2x2v_ser_p2 due to the lack of gkyl_basis_ser_2x_p2_inv.h //2
  // 3x kernels
  { NULL, gyrokinetic_mom_cross_bgk_3x2v_ser_p1, NULL}, // no gyrokinetic_mom_cross_bgk_2x2v_ser_p2 due to the lack of gkyl_basis_ser_3x_p2_inv.h //3
};

GKYL_CU_D
static gyrokinetic_mom_cross_bgk_t
choose_mom_cross_bgk_gyrokinetic_kern(int cdim, int vdim, int poly_order){
  return ser_mom_cross_gyrokinetic_list[cv_index[cdim].vdim[vdim]].kernels[poly_order];
}

// Primary struct in this updater
struct gkyl_mom_cross_bgk_gyrokinetic {
  bool use_gpu;
  gyrokinetic_mom_cross_bgk_t mom_cross_calc; // a pointer to cross moments kernel

  struct gkyl_mom_cross_bgk_gyrokinetic *on_dev; 
};
