#pragma once

// Private header, not for direct use in user code

#include <gkyl_vlasov_lbo_mom_kernels.h>

typedef void (*lbo_momf_t)(const int *idx, const int *atLower, const double *vBoundary,
  const double *dxv, const double *fIn, double* GKYL_RESTRICT out);

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { lbo_momf_t kernels[3]; } gkyl_lbo_mom_kern_list;

//
// Serendipity basis kernels
//

// f boundary integral moment correction kernel lists
GKYL_CU_D
static const gkyl_lbo_mom_kern_list ser_boundary_integral_f_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_integral_1x1v_ser_f_p1, vlasov_boundary_integral_1x1v_ser_f_p2 }, // 0
  { NULL, vlasov_boundary_integral_1x2v_ser_f_p1, vlasov_boundary_integral_1x2v_ser_f_p2 }, // 1
  { NULL, vlasov_boundary_integral_1x3v_ser_f_p1, vlasov_boundary_integral_1x3v_ser_f_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_integral_2x2v_ser_f_p1, vlasov_boundary_integral_2x2v_ser_f_p2 }, // 3
  { NULL, vlasov_boundary_integral_2x3v_ser_f_p1, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_integral_3x3v_ser_f_p1, NULL }, // 5
};

// v*f boundary integral moment correction kernel lists
GKYL_CU_D
static const gkyl_lbo_mom_kern_list ser_boundary_integral_vf_kernels[] = {
  // 1x kernels
  { NULL, vlasov_boundary_integral_1x1v_ser_vf_p1, vlasov_boundary_integral_1x1v_ser_vf_p2 }, // 0
  { NULL, vlasov_boundary_integral_1x2v_ser_vf_p1, vlasov_boundary_integral_1x2v_ser_vf_p2 }, // 1
  { NULL, vlasov_boundary_integral_1x3v_ser_vf_p1, vlasov_boundary_integral_1x3v_ser_vf_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_boundary_integral_2x2v_ser_vf_p1, vlasov_boundary_integral_2x2v_ser_vf_p2 }, // 3
  { NULL, vlasov_boundary_integral_2x3v_ser_vf_p1, NULL }, // 4
  // 3x kernels
  { NULL, vlasov_boundary_integral_3x3v_ser_vf_p1, NULL }, // 5
};

struct lbo_mom_type {
  struct gkyl_mom_type momt;
  lbo_momf_t kernel; // moment calculation kernel
  const double *vBoundary;
  const int *atLower;
};

GKYL_CU_D
static void
kernel(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out)
{
  struct lbo_mom_type *mom_bcorr = container_of(momt, struct lbo_mom_type, momt);

  return mom_bcorr->kernel(idx, mom_bcorr->atLower, mom_bcorr->vBoundary, dx, f, out);
}
