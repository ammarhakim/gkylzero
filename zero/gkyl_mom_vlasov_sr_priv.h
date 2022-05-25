#pragma once

// Private header, not for direct use in user code

#include <gkyl_array.h>
#include <gkyl_mom_type.h>
#include <gkyl_ref_count.h>
#include <gkyl_mom_vlasov_kernels.h>
#include <gkyl_range.h>
#include <gkyl_util.h>

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

typedef void (*vlasov_sr_momf_t)(const double *xc, const double *dx,
  const int *idx, const double *p_over_gamma, const double *fIn, double* GKYL_RESTRICT out);

// for use in kernel tables
typedef struct { vlasov_sr_momf_t kernels[3]; } gkyl_mom_kern_list;

//
// Serendipity basis kernels
//

// M0 kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ser_m0_kernels[] = {
  // 1x kernels
  { NULL, vlasov_sr_M0_1x1v_ser_p1, vlasov_sr_M0_1x1v_ser_p2 }, // 0
  { NULL, vlasov_sr_M0_1x2v_ser_p1, vlasov_sr_M0_1x2v_ser_p2 }, // 1
  { NULL, vlasov_sr_M0_1x3v_ser_p1, vlasov_sr_M0_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_sr_M0_2x2v_ser_p1, vlasov_sr_M0_2x2v_ser_p2 }, // 3
  { NULL, vlasov_sr_M0_2x3v_ser_p1, vlasov_sr_M0_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_sr_M0_3x3v_ser_p1, NULL                  }, // 5
};

// M1i kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ser_m1i_kernels[] = {
  // 1x kernels
  { NULL, vlasov_sr_M1i_1x1v_ser_p1, vlasov_sr_M1i_1x1v_ser_p2 }, // 0
  { NULL, vlasov_sr_M1i_1x2v_ser_p1, vlasov_sr_M1i_1x2v_ser_p2 }, // 1
  { NULL, vlasov_sr_M1i_1x3v_ser_p1, vlasov_sr_M1i_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_sr_M1i_2x2v_ser_p1, vlasov_sr_M1i_2x2v_ser_p2 }, // 3
  { NULL, vlasov_sr_M1i_2x3v_ser_p1, vlasov_sr_M1i_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_sr_M1i_3x3v_ser_p1, NULL                   }, // 5
};

// Integrated moments kernel list
GKYL_CU_D
static const gkyl_mom_kern_list ser_int_mom_kernels[] = {
  // 1x kernels
  { NULL, vlasov_sr_int_mom_1x1v_ser_p1, vlasov_sr_int_mom_1x1v_ser_p2 }, // 0
  { NULL, vlasov_sr_int_mom_1x2v_ser_p1, vlasov_sr_int_mom_1x2v_ser_p2 }, // 1
  { NULL, vlasov_sr_int_mom_1x3v_ser_p1, vlasov_sr_int_mom_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, vlasov_sr_int_mom_2x2v_ser_p1, vlasov_sr_int_mom_2x2v_ser_p2 }, // 3
  { NULL, vlasov_sr_int_mom_2x3v_ser_p1, vlasov_sr_int_mom_2x3v_ser_p2 }, // 4
  // 3x kernels
  { NULL, vlasov_sr_int_mom_3x3v_ser_p1, NULL                     }, // 5
};

//
// Tensor-product basis kernels
//

// // M0 kernel list
// GKYL_CU_D
// static const gkyl_mom_kern_list ten_m0_kernels[] = {
//   // 1x kernels
//   { NULL, vlasov_sr_M0_1x1v_ser_p1, vlasov_sr_M0_1x1v_tensor_p2 }, // 0
//   { NULL, vlasov_sr_M0_1x2v_ser_p1, vlasov_sr_M0_1x2v_tensor_p2 }, // 1
//   { NULL, vlasov_sr_M0_1x3v_ser_p1, vlasov_sr_M0_1x3v_tensor_p2 }, // 2
//   // 2x kernels
//   { NULL, vlasov_sr_M0_2x2v_ser_p1, vlasov_sr_M0_2x2v_tensor_p2 }, // 3
//   { NULL, vlasov_sr_M0_2x3v_ser_p1, NULL                  }, // 4
//   // 3x kernels
//   { NULL, vlasov_sr_M0_3x3v_ser_p1, NULL                  }, // 5
// };

// // M1i kernel list
// GKYL_CU_D
// static const gkyl_mom_kern_list ten_m1i_kernels[] = {
//   // 1x kernels
//   { NULL, vlasov_sr_M1i_1x1v_ser_p1, vlasov_sr_M1i_1x1v_tensor_p2 }, // 0
//   { NULL, vlasov_sr_M1i_1x2v_ser_p1, vlasov_sr_M1i_1x2v_tensor_p2 }, // 1
//   { NULL, vlasov_sr_M1i_1x3v_ser_p1, vlasov_sr_M1i_1x3v_tensor_p2 }, // 2
//   // 2x kernels
//   { NULL, vlasov_sr_M1i_2x2v_ser_p1, vlasov_sr_M1i_2x2v_tensor_p2 }, // 3
//   { NULL, vlasov_sr_M1i_2x3v_ser_p1, NULL                   }, // 4
//   // 3x kernels
//   { NULL, vlasov_sr_M1i_3x3v_ser_p1, NULL                   }, // 5
// };

struct mom_type_vlasov_sr {
  struct gkyl_mom_type momt;
  vlasov_sr_momf_t kernel; // moment calculation kernel
  struct gkyl_range vel_range; // velocity space range
  struct gkyl_mom_vlasov_sr_auxfields auxfields; // Auxiliary fields.
};

/**
 * Free moment object.
 *
 * @param ref Reference counter for moment to free
 */
void gkyl_mom_vm_sr_free(const struct gkyl_ref_count *ref);

GKYL_CU_D
static void
kernel(const struct gkyl_mom_type *momt, const double *xc, const double *dx,
  const int *idx, const double *f, double* out, void *param)
{
  struct mom_type_vlasov_sr *mom_vm_sr = container_of(momt, struct mom_type_vlasov_sr, momt);

  int cdim = mom_vm_sr->momt.cdim;
  int pdim = mom_vm_sr->momt.pdim;
  int idx_vel[GKYL_MAX_DIM];
  for (int i=0; i<pdim-cdim; ++i)
    idx_vel[i] = idx[cdim+i];

  long vidx = gkyl_range_idx(&mom_vm_sr->vel_range, idx_vel);

  return mom_vm_sr->kernel(xc, dx, idx, 
    (const double*) gkyl_array_cfetch(mom_vm_sr->auxfields.p_over_gamma, vidx),
    f, out);
}
