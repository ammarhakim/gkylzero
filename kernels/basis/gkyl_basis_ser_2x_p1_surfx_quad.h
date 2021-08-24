GKYL_CU_DH static inline double 
ser_2x_p1_surfx_quad_0(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return (-0.8660254037844386*f[3])-0.5*f[2]+0.8660254037844386*f[1]+0.5*f[0]; 
  else 
    return 0.8660254037844386*f[3]-0.5*f[2]-0.8660254037844386*f[1]+0.5*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_2x_p1_surfx_quad_1(int side, const double* GKYL_RESTRICT f) { 
  if (side == 1) 
    return 0.8660254037844386*f[3]+0.5*f[2]+0.8660254037844386*f[1]+0.5*f[0]; 
  else 
    return (-0.8660254037844386*f[3])+0.5*f[2]-0.8660254037844386*f[1]+0.5*f[0]; 
} 
