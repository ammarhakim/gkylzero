#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_int_mom_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[3]; 
  tmp[0] = (2.8284271247461907*f[3]*vmap[3])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.828427124746191*vmap[3]*f[13])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += (2.0*vmap[1]*f[2]+2.0*f[0]*vmap[0])*volFact; 
  out[2] += (1.264911064067352*vmap1R2*f[8]+2.8284271247461907*vmap[0]*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap1R2+1.4142135623730951*f[0]*vmap0R2)*volFact; 
  out[3] += (bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 
