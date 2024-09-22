#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_int_mom_2x2v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.39269908169872414*dxv[0]*dxv[1]*dxv[2]*dxv[3]/m_; 
 
  double tmp[8]; 
  tmp[0] = (2.8284271247461907*vmap[3]*f[4])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[8])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.8284271247461907*vmap[3]*f[9])/m_+(2.8284271247461907*f[2]*vmap[2])/m_; 
  tmp[3] = (2.8284271247461907*vmap[3]*f[16])/m_+(2.8284271247461907*vmap[2]*f[5])/m_; 
  tmp[4] = (2.828427124746191*vmap[3]*f[25])/m_+(2.8284271247461907*vmap[2]*f[11])/m_; 
  tmp[5] = (2.828427124746191*vmap[3]*f[26])/m_+(2.8284271247461907*vmap[2]*f[12])/m_; 
  tmp[6] = (2.828427124746191*vmap[3]*f[35])/m_+(2.8284271247461907*vmap[2]*f[19])/m_; 
  tmp[7] = (2.828427124746191*vmap[3]*f[36])/m_+(2.8284271247461907*vmap[2]*f[20])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 4.0*f[0]*volFact; 
  out[1] += (2.8284271247461907*vmap[1]*f[3]+2.8284271247461907*f[0]*vmap[0])*volFact; 
  out[2] += (1.7888543819998317*vmap1R2*f[13]+4.0*vmap[0]*vmap[1]*f[3]+2.0*f[0]*vmap1R2+2.0*f[0]*vmap0R2)*volFact; 
  out[3] += (bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 
