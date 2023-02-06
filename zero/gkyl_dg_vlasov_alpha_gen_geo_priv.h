#pragma once

// Private header, not for direct use in user code

#include <gkyl_dg_vlasov_alpha_gen_geo_kernels.h>

typedef void (*dg_vlasov_alpha_gen_geof_t)(const double *w, const double *dxv, const double *tvComp, const double *gij, double* GKYL_RESTRICT alpha_geo);

// for use in kernel tables
typedef struct { dg_vlasov_alpha_gen_geof_t kernels[3]; } gkyl_vlasov_alpha_gen_geo_kern_list;

//
// Serendipity basis kernels
// 

// Alpha gen geo kernel list
GKYL_CU_D
static const gkyl_vlasov_alpha_gen_geo_kern_list ser_vlasov_alpha_gen_geo_kernels[] = {
  { NULL, vlasov_alpha_gen_geo_3x3v_ser_p1, NULL } // 0 only 3x3v kernels exist
};
