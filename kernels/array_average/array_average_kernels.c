#include <gkyl_array_average_kernels.h>

GKYL_CU_DH void gkyl_array_average_1x_ser_p1_avgx_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 1.4142135623730951*fin[0]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgx_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 1.4142135623730951*fin[0]*subvol; 
  out[1] += 1.4142135623730951*fin[2]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 1.4142135623730951*fin[0]*subvol; 
  out[1] += 1.4142135623730951*fin[1]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgxy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 2.0*fin[0]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgx_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 1.4142135623730951*fin[0]*subvol; 
  out[1] += 1.4142135623730951*fin[2]*subvol; 
  out[2] += 1.4142135623730951*fin[3]*subvol; 
  out[3] += 1.4142135623730951*fin[6]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 1.4142135623730951*fin[0]*subvol; 
  out[1] += 1.4142135623730951*fin[1]*subvol; 
  out[2] += 1.4142135623730951*fin[3]*subvol; 
  out[3] += 1.4142135623730951*fin[5]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 1.4142135623730951*fin[0]*subvol; 
  out[1] += 1.4142135623730951*fin[1]*subvol; 
  out[2] += 1.4142135623730951*fin[2]*subvol; 
  out[3] += 1.4142135623730951*fin[4]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgxy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 2.0*fin[0]*subvol; 
  out[1] += 2.0*fin[3]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgxz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 2.0*fin[0]*subvol; 
  out[1] += 2.0*fin[2]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgyz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 2.0*fin[0]*subvol; 
  out[1] += 2.0*fin[1]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgxyz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 2.8284271247461907*fin[0]*subvol; 
} 
