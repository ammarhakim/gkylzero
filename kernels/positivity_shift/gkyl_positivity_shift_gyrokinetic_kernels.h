#pragma once

#include <gkyl_util.h>

EXTERN_C_BEG

GKYL_CU_DH bool positivity_shift_gyrokinetic_1x1v_ser_p1(double ffloor, double *distf, double *Deltaf);

GKYL_CU_DH bool positivity_shift_gyrokinetic_1x2v_ser_p1(double ffloor, double *distf, double *Deltaf);

GKYL_CU_DH bool positivity_shift_gyrokinetic_2x2v_ser_p1(double ffloor, double *distf, double *Deltaf);

GKYL_CU_DH bool positivity_shift_gyrokinetic_3x2v_ser_p1(double ffloor, double *distf, double *Deltaf);


EXTERN_C_END
