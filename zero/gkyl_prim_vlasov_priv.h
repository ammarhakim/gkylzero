// Private header: not for direct use
#pragma once

#include <gkyl_prim_vlasov_kernels.h>
#include <gkyl_mat.h>
#include <gkyl_util.h>

// The cv_index[cd].vdim[vd] is used to index the various list of
// kernels below
static struct { int vdim[4]; } cv_index[] = {
  {-1, -1, -1, -1}, // 0x makes no sense
  {-1,  0,  1,  2}, // 1x kernel indices
  {-1, -1,  3,  4}, // 2x kernel indices
  {-1, -1, -1,  5}, // 3x kernel indices  
};

// for use in kernel tables
typedef struct { self_prim_t kernels[3]; } gkyl_prim_vlasov_kern_list;

//
// Serendipity basis kernels
//

// self primitive moment kernel list
GKYL_CU_D
static const gkyl_prim_vlasov_kern_list ser_self_prim_kernels[] = {
  // 1x kernels
  { NULL, NULL, vlasov_self_prim_moments_1x1v_ser_p2 }, // 0
  { NULL, NULL, vlasov_self_prim_moments_1x2v_ser_p2 }, // 1
  { NULL, NULL, vlasov_self_prim_moments_1x3v_ser_p2 }, // 2
  // 2x kernels
  { NULL, NULL, vlasov_self_prim_moments_2x2v_ser_p2 }, // 3
  { NULL, NULL, NULL }, // 4
  // 3x kernels
  { NULL, NULL, NULL }, // 5
};
