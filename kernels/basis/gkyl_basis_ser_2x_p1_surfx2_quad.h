GKYL_CU_DH static inline double 
ser_2x_p1_surfx2_quad_0_r(const double* GKYL_RESTRICT f) { 
  return (-0.8660254037844386*f[3])+0.8660254037844386*f[2]-0.5*f[1]+0.5*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_2x_p1_surfx2_quad_0_l(const double* GKYL_RESTRICT f) { 
    return 0.8660254037844386*f[3]-0.8660254037844386*f[2]-0.5*f[1]+0.5*f[0]; 
} 
GKYL_CU_DH static inline double 
ser_2x_p1_surfx2_quad_1_r(const double* GKYL_RESTRICT f) { 
  return 0.8660254037844386*(f[3]+f[2])+0.5*(f[1]+f[0]); 
} 
GKYL_CU_DH static inline double 
ser_2x_p1_surfx2_quad_1_l(const double* GKYL_RESTRICT f) { 
    return 0.5*(f[1]+f[0])-0.8660254037844386*(f[3]+f[2]); 
} 
