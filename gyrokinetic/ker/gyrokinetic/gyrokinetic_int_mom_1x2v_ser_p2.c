#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_int_M0_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 

  out[0] += 2.8284271247461907*f[0]*volFact; 
} 

GKYL_CU_DH void gyrokinetic_int_M1_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 

  out[0] += (2.0*vmap[1]*f[2]+2.0*f[0]*vmap[0])*volFact; 
} 

GKYL_CU_DH void gyrokinetic_int_M2_par_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.264911064067352*vmap1R2*f[8]+2.8284271247461907*vmap[0]*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap1R2+1.4142135623730951*f[0]*vmap0R2)*volFact; 
} 

GKYL_CU_DH void gyrokinetic_int_M2_perp_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[3]; 
  tmp[0] = (2.8284271247461907*f[3]*vmap[3])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.828427124746191*vmap[3]*f[13])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
 

  out[0] += (bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void gyrokinetic_int_M2_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[3]; 
  tmp[0] = (2.8284271247461907*f[3]*vmap[3])/m_+(2.8284271247461907*f[0]*vmap[2])/m_; 
  tmp[1] = (2.8284271247461907*vmap[3]*f[5])/m_+(2.8284271247461907*f[1]*vmap[2])/m_; 
  tmp[2] = (2.828427124746191*vmap[3]*f[13])/m_+(2.8284271247461907*vmap[2]*f[7])/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (1.264911064067352*vmap1R2*f[8]+bmag[2]*tmp[2]+2.8284271247461907*vmap[0]*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap1R2+bmag[1]*tmp[1]+1.4142135623730951*f[0]*vmap0R2+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void gyrokinetic_int_M3_par_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += (2.6832815729997477*vmap[0]*vmap1R2*f[8]+1.8*vmap1R3*f[2]+3.0*vmap0R2*vmap[1]*f[2]+3.0*f[0]*vmap[0]*vmap1R2+f[0]*vmap0R3)*volFact; 
} 

GKYL_CU_DH void gyrokinetic_int_M3_perp_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[3]; 
  tmp[0] = (2.0*vmap[1]*vmap[3]*f[6])/m_+(2.0*vmap[0]*f[3]*vmap[3])/m_+(2.0*vmap[1]*f[2]*vmap[2])/m_+(2.0*f[0]*vmap[0]*vmap[2])/m_; 
  tmp[1] = (2.0*vmap[1]*vmap[3]*f[10])/m_+(2.0*vmap[0]*vmap[3]*f[5])/m_+(2.0*vmap[1]*vmap[2]*f[4])/m_+(2.0*vmap[0]*f[1]*vmap[2])/m_; 
  tmp[2] = (2.0*vmap[1]*vmap[3]*f[17])/m_+(2.0000000000000004*vmap[0]*vmap[3]*f[13])/m_+(2.0000000000000004*vmap[1]*vmap[2]*f[11])/m_+(2.0*vmap[0]*vmap[2]*f[7])/m_; 
 

  out[0] += (bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void gyrokinetic_int_three_moments_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
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
  out[2] += (1.264911064067352*vmap1R2*f[8]+bmag[2]*tmp[2]+2.8284271247461907*vmap[0]*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap1R2+bmag[1]*tmp[1]+1.4142135623730951*f[0]*vmap0R2+bmag[0]*tmp[0])*volFact; 
} 

GKYL_CU_DH void gyrokinetic_int_four_moments_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
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

GKYL_CU_DH void gyrokinetic_int_hamiltonian_moments_1x2v_ser_p2(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.7853981633974483*dxv[0]*dxv[1]*dxv[2]/m_; 
 
  double tmp[3]; 
  tmp[0] = 1.4142135623730951*f[3]*vmap[3]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[5]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[13]+1.4142135623730951*vmap[2]*f[7]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += (2.0*vmap[1]*f[2]*m_+2.0*f[0]*vmap[0]*m_)*volFact; 
  out[2] += (2.0*phi[2]*f[7]*q_+2.0*f[1]*phi[1]*q_+2.0*f[0]*phi[0]*q_+0.6324555320336759*vmap1R2*f[8]*m_+1.4142135623730951*vmap[0]*vmap[1]*f[2]*m_+0.7071067811865475*f[0]*vmap1R2*m_+0.7071067811865475*f[0]*vmap0R2*m_+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 

