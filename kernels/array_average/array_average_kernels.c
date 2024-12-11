#include <gkyl_array_average_kernels.h>

GKYL_CU_DH void gkyl_array_average_1x_ser_p1_avgx_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += (fin[1]*win[1]+fin[0]*win[0])*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgx_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += (0.7071067811865475*fin[3]*win[3]+0.7071067811865475*fin[2]*win[2]+0.7071067811865475*fin[1]*win[1]+0.7071067811865475*fin[0]*win[0])*subvol; 
  out[1] += (0.7071067811865475*fin[1]*win[3]+0.7071067811865475*win[1]*fin[3]+0.7071067811865475*fin[0]*win[2]+0.7071067811865475*win[0]*fin[2])*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += (0.7071067811865475*fin[3]*win[3]+0.7071067811865475*fin[2]*win[2]+0.7071067811865475*fin[1]*win[1]+0.7071067811865475*fin[0]*win[0])*subvol; 
  out[1] += (0.7071067811865475*fin[2]*win[3]+0.7071067811865475*win[2]*fin[3]+0.7071067811865475*fin[0]*win[1]+0.7071067811865475*win[0]*fin[1])*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_2x_ser_p1_avgxy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += (fin[3]*win[3]+fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgx_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += (0.5*fin[7]*win[7]+0.5*fin[6]*win[6]+0.5*fin[5]*win[5]+0.5*fin[4]*win[4]+0.5*fin[3]*win[3]+0.5*fin[2]*win[2]+0.5*fin[1]*win[1]+0.5*fin[0]*win[0])*subvol; 
  out[1] += (0.5*fin[5]*win[7]+0.5*win[5]*fin[7]+0.5*fin[3]*win[6]+0.5*win[3]*fin[6]+0.5*fin[1]*win[4]+0.5*win[1]*fin[4]+0.5*fin[0]*win[2]+0.5*win[0]*fin[2])*subvol; 
  out[2] += (0.5*fin[4]*win[7]+0.5*win[4]*fin[7]+0.5*fin[2]*win[6]+0.5*win[2]*fin[6]+0.5*fin[1]*win[5]+0.5*win[1]*fin[5]+0.5*fin[0]*win[3]+0.5*win[0]*fin[3])*subvol; 
  out[3] += (0.5*fin[1]*win[7]+0.5*win[1]*fin[7]+0.5*fin[0]*win[6]+0.5*win[0]*fin[6]+0.5*fin[4]*win[5]+0.5*win[4]*fin[5]+0.5*fin[2]*win[3]+0.5*win[2]*fin[3])*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += (0.5*fin[7]*win[7]+0.5*fin[6]*win[6]+0.5*fin[5]*win[5]+0.5*fin[4]*win[4]+0.5*fin[3]*win[3]+0.5*fin[2]*win[2]+0.5*fin[1]*win[1]+0.5*fin[0]*win[0])*subvol; 
  out[1] += (0.5*fin[6]*win[7]+0.5*win[6]*fin[7]+0.5*fin[3]*win[5]+0.5*win[3]*fin[5]+0.5*fin[2]*win[4]+0.5*win[2]*fin[4]+0.5*fin[0]*win[1]+0.5*win[0]*fin[1])*subvol; 
  out[2] += (0.5*fin[4]*win[7]+0.5*win[4]*fin[7]+0.5*fin[2]*win[6]+0.5*win[2]*fin[6]+0.5*fin[1]*win[5]+0.5*win[1]*fin[5]+0.5*fin[0]*win[3]+0.5*win[0]*fin[3])*subvol; 
  out[3] += (0.5*fin[2]*win[7]+0.5*win[2]*fin[7]+0.5*fin[4]*win[6]+0.5*win[4]*fin[6]+0.5*fin[0]*win[5]+0.5*win[0]*fin[5]+0.5*fin[1]*win[3]+0.5*win[1]*fin[3])*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += (0.5*fin[7]*win[7]+0.5*fin[6]*win[6]+0.5*fin[5]*win[5]+0.5*fin[4]*win[4]+0.5*fin[3]*win[3]+0.5*fin[2]*win[2]+0.5*fin[1]*win[1]+0.5*fin[0]*win[0])*subvol; 
  out[1] += (0.5*fin[6]*win[7]+0.5*win[6]*fin[7]+0.5*fin[3]*win[5]+0.5*win[3]*fin[5]+0.5*fin[2]*win[4]+0.5*win[2]*fin[4]+0.5*fin[0]*win[1]+0.5*win[0]*fin[1])*subvol; 
  out[2] += (0.5*fin[5]*win[7]+0.5*win[5]*fin[7]+0.5*fin[3]*win[6]+0.5*win[3]*fin[6]+0.5*fin[1]*win[4]+0.5*win[1]*fin[4]+0.5*fin[0]*win[2]+0.5*win[0]*fin[2])*subvol; 
  out[3] += (0.5*fin[3]*win[7]+0.5*win[3]*fin[7]+0.5*fin[5]*win[6]+0.5*win[5]*fin[6]+0.5*fin[0]*win[4]+0.5*win[0]*fin[4]+0.5*fin[1]*win[2]+0.5*win[1]*fin[2])*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgxy_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += (0.7071067811865475*fin[7]*win[7]+0.7071067811865475*fin[6]*win[6]+0.7071067811865475*fin[5]*win[5]+0.7071067811865475*fin[4]*win[4]+0.7071067811865475*fin[3]*win[3]+0.7071067811865475*fin[2]*win[2]+0.7071067811865475*fin[1]*win[1]+0.7071067811865475*fin[0]*win[0])*subvol; 
  out[1] += (0.7071067811865475*fin[4]*win[7]+0.7071067811865475*win[4]*fin[7]+0.7071067811865475*fin[2]*win[6]+0.7071067811865475*win[2]*fin[6]+0.7071067811865475*fin[1]*win[5]+0.7071067811865475*win[1]*fin[5]+0.7071067811865475*fin[0]*win[3]+0.7071067811865475*win[0]*fin[3])*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgxz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += (0.7071067811865475*fin[7]*win[7]+0.7071067811865475*fin[6]*win[6]+0.7071067811865475*fin[5]*win[5]+0.7071067811865475*fin[4]*win[4]+0.7071067811865475*fin[3]*win[3]+0.7071067811865475*fin[2]*win[2]+0.7071067811865475*fin[1]*win[1]+0.7071067811865475*fin[0]*win[0])*subvol; 
  out[1] += (0.7071067811865475*fin[5]*win[7]+0.7071067811865475*win[5]*fin[7]+0.7071067811865475*fin[3]*win[6]+0.7071067811865475*win[3]*fin[6]+0.7071067811865475*fin[1]*win[4]+0.7071067811865475*win[1]*fin[4]+0.7071067811865475*fin[0]*win[2]+0.7071067811865475*win[0]*fin[2])*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgyz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += (0.7071067811865475*fin[7]*win[7]+0.7071067811865475*fin[6]*win[6]+0.7071067811865475*fin[5]*win[5]+0.7071067811865475*fin[4]*win[4]+0.7071067811865475*fin[3]*win[3]+0.7071067811865475*fin[2]*win[2]+0.7071067811865475*fin[1]*win[1]+0.7071067811865475*fin[0]*win[0])*subvol; 
  out[1] += (0.7071067811865475*fin[6]*win[7]+0.7071067811865475*win[6]*fin[7]+0.7071067811865475*fin[3]*win[5]+0.7071067811865475*win[3]*fin[5]+0.7071067811865475*fin[2]*win[4]+0.7071067811865475*win[2]*fin[4]+0.7071067811865475*fin[0]*win[1]+0.7071067811865475*win[0]*fin[1])*subvol; 
} 
GKYL_CU_DH void gkyl_array_average_3x_ser_p1_avgxyz_ker(const double subvol, const double *win, const double *fin, double* GKYL_RESTRICT out) 
{ 
  out[0] += (fin[7]*win[7]+fin[6]*win[6]+fin[5]*win[5]+fin[4]*win[4]+fin[3]*win[3]+fin[2]*win[2]+fin[1]*win[1]+fin[0]*win[0])*subvol; 
} 
