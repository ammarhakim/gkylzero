#pragma once

#include <gkyl_util.h>
#include <math.h>

EXTERN_C_BEG

GKYL_CU_DH void gkyl_array_average_1x_ser_p1_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_y_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_x_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_yz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_xz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_xy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_z_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_y_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_x_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);

EXTERN_C_END