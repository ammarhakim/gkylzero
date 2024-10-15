#pragma once

#include <gkyl_util.h>
#include <math.h> 

EXTERN_C_BEG

void dg_interpolate_gyrokinetic_1x1v_ser_p1(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
void dg_interpolate_gyrokinetic_1x2v_ser_p1(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
void dg_interpolate_gyrokinetic_2x2v_ser_p1(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);
void dg_interpolate_gyrokinetic_3x2v_ser_p1(const double *wDo, const double *wTar, const double *dxDo, const double *dxTar, const double *fldDo, double *fldTar);



EXTERN_C_END
