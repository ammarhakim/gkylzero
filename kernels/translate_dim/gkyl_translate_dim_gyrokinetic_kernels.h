#pragma once

#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH void translate_dim_gyrokinetic_2x2v_ser_p1_from_1x2v_p1(const double *flow, double *fout);

GKYL_CU_DH void translate_dim_gyrokinetic_3x2v_ser_p1_from_1x2v_p1(const double *flow, double *fout);
GKYL_CU_DH void translate_dim_gyrokinetic_3x2v_ser_p1_from_2x2v_p1(const double *flow, double *fout);


EXTERN_C_END
