#include <gkyl_array_average_kernels.h>

GKYL_CU_DH void gkyl_array_average_1x_ser_p1_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 1.4142135623730951*fin[0]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_y_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 1.4142135623730951*fin[0]*subvol; 
  out[1] += 1.4142135623730951*fin[2]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_x_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 1.4142135623730951*fin[0]*subvol; 
  out[1] += 1.4142135623730951*fin[1]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 2.0*fin[0]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_yz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 1.4142135623730951*fin[0]*subvol; 
  out[1] += 1.4142135623730951*fin[2]*subvol; 
  out[2] += 1.4142135623730951*fin[3]*subvol; 
  out[3] += 1.4142135623730951*fin[6]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_xz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 1.4142135623730951*fin[0]*subvol; 
  out[1] += 1.4142135623730951*fin[1]*subvol; 
  out[2] += 1.4142135623730951*fin[3]*subvol; 
  out[3] += 1.4142135623730951*fin[5]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_xy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 1.4142135623730951*fin[0]*subvol; 
  out[1] += 1.4142135623730951*fin[1]*subvol; 
  out[2] += 1.4142135623730951*fin[2]*subvol; 
  out[3] += 1.4142135623730951*fin[4]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_z_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 2.0*fin[0]*subvol; 
  out[1] += 2.0*fin[3]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_y_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 2.0*fin[0]*subvol; 
  out[1] += 2.0*fin[2]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_x_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 2.0*fin[0]*subvol; 
  out[1] += 2.0*fin[1]*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += 2.8284271247461907*fin[0]*subvol; 
} 
