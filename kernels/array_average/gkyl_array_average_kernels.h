#pragma once

#include <gkyl_util.h>
EXTERN_C_BEG

GKYL_CU_DH void gkyl_array_average_1x_ser_p1_avgx_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgx_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgxy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgx_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgxy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgxz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgyz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgxyz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_1x_ser_p2_avgx_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p2_avgx_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p2_avgy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_2x_ser_p2_avgxy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p2_avgx_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p2_avgy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p2_avgz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p2_avgxy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p2_avgxz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p2_avgyz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);
GKYL_CU_DH void gkyl_array_average_3x_ser_p2_avgxyz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out);

EXTERN_C_END