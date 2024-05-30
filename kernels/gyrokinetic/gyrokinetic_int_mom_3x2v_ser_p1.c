#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_int_mom_3x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = M_PI/m_*0.0625*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]; 
 
  double tmp[8]; 
  tmp[0] = (2.828427124746191*vmap[3]*f[5])/m_+(2.828427124746191*f[0]*vmap[2])/m_; 
  tmp[1] = (2.828427124746191*vmap[3]*f[12])/m_+(2.828427124746191*f[1]*vmap[2])/m_; 
  tmp[2] = (2.828427124746191*vmap[3]*f[13])/m_+(2.828427124746191*f[2]*vmap[2])/m_; 
  tmp[3] = (2.828427124746191*vmap[3]*f[14])/m_+(2.828427124746191*vmap[2]*f[3])/m_; 
  tmp[4] = (2.828427124746191*vmap[3]*f[20])/m_+(2.828427124746191*vmap[2]*f[6])/m_; 
  tmp[5] = (2.828427124746191*vmap[3]*f[21])/m_+(2.828427124746191*vmap[2]*f[7])/m_; 
  tmp[6] = (2.828427124746191*vmap[3]*f[22])/m_+(2.828427124746191*vmap[2]*f[8])/m_; 
  tmp[7] = (2.828427124746191*vmap[3]*f[27])/m_+(2.828427124746191*vmap[2]*f[16])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += (4.0*vmap[1]*f[4]+4.0*f[0]*vmap[0])*volFact; 
  out[2] += (2.529822128134704*vmap1R2*f[32]+5.656854249492382*vmap[0]*vmap[1]*f[4]+2.828427124746191*f[0]*vmap1R2+2.828427124746191*f[0]*vmap0R2)*volFact; 
  out[3] += (bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 
