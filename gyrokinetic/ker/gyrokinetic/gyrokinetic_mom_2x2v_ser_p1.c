#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_M0_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
} 

GKYL_CU_DH void gyrokinetic_M1_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  out[0] += (1.4142135623730951*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*f[6]+1.4142135623730951*vmap[0]*f[1])*volFact; 
  out[2] += (1.4142135623730951*vmap[1]*f[7]+1.4142135623730951*vmap[0]*f[2])*volFact; 
  out[3] += (1.4142135623730951*vmap[1]*f[11]+1.4142135623730951*vmap[0]*f[5])*volFact; 
} 

GKYL_CU_DH void gyrokinetic_M2_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*f[16]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[1] += (0.8944271909999161*vmap1R2*f[17]+2.0*vmap[0]*vmap[1]*f[6]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  out[2] += (0.8944271909999161*vmap1R2*f[18]+2.0*vmap[0]*vmap[1]*f[7]+vmap1R2*f[2]+vmap0R2*f[2])*volFact; 
  out[3] += (0.8944271909999159*vmap1R2*f[20]+2.0*vmap[0]*vmap[1]*f[11]+vmap1R2*f[5]+vmap0R2*f[5])*volFact; 
  double tmp[4]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[12]+1.4142135623730951*vmap[2]*f[5]; 
  out[0] += (2.0*(0.5*bmag[3]*tmp[3]+0.5*bmag[2]*tmp[2]+0.5*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.5*bmag[2]*tmp[3]+0.5*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[1]+0.5*tmp[0]*bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.5*bmag[1]*tmp[3]+0.5*tmp[1]*bmag[3]+0.5*bmag[0]*tmp[2]+0.5*tmp[0]*bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.5*bmag[0]*tmp[3]+0.5*tmp[0]*bmag[3]+0.5*bmag[1]*tmp[2]+0.5*tmp[1]*bmag[2])*volFact)/m_; 
} 

GKYL_CU_DH void gyrokinetic_M2_par_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*f[16]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[1] += (0.8944271909999161*vmap1R2*f[17]+2.0*vmap[0]*vmap[1]*f[6]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  out[2] += (0.8944271909999161*vmap1R2*f[18]+2.0*vmap[0]*vmap[1]*f[7]+vmap1R2*f[2]+vmap0R2*f[2])*volFact; 
  out[3] += (0.8944271909999159*vmap1R2*f[20]+2.0*vmap[0]*vmap[1]*f[11]+vmap1R2*f[5]+vmap0R2*f[5])*volFact; 
} 

GKYL_CU_DH void gyrokinetic_M2_perp_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp[4]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[12]+1.4142135623730951*vmap[2]*f[5]; 
  out[0] += (2.0*(0.5*bmag[3]*tmp[3]+0.5*bmag[2]*tmp[2]+0.5*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.5*bmag[2]*tmp[3]+0.5*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[1]+0.5*tmp[0]*bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.5*bmag[1]*tmp[3]+0.5*tmp[1]*bmag[3]+0.5*bmag[0]*tmp[2]+0.5*tmp[0]*bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.5*bmag[0]*tmp[3]+0.5*tmp[0]*bmag[3]+0.5*bmag[1]*tmp[2]+0.5*tmp[1]*bmag[2])*volFact)/m_; 
} 

GKYL_CU_DH void gyrokinetic_M3_par_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += (1.8973665961010278*vmap[0]*vmap1R2*f[16]+1.2727922061357855*vmap1R3*f[3]+2.1213203435596424*vmap0R2*vmap[1]*f[3]+2.1213203435596424*f[0]*vmap[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R3)*volFact; 
  out[1] += (1.897366596101028*vmap[0]*vmap1R2*f[17]+1.2727922061357855*vmap1R3*f[6]+2.1213203435596424*vmap0R2*vmap[1]*f[6]+2.1213203435596424*vmap[0]*f[1]*vmap1R2+0.7071067811865475*vmap0R3*f[1])*volFact; 
  out[2] += (1.897366596101028*vmap[0]*vmap1R2*f[18]+1.2727922061357855*vmap1R3*f[7]+2.1213203435596424*vmap0R2*vmap[1]*f[7]+2.1213203435596424*vmap[0]*vmap1R2*f[2]+0.7071067811865475*vmap0R3*f[2])*volFact; 
  out[3] += (1.8973665961010278*vmap[0]*vmap1R2*f[20]+1.2727922061357855*vmap1R3*f[11]+2.1213203435596424*vmap0R2*vmap[1]*f[11]+2.1213203435596424*vmap[0]*vmap1R2*f[5]+0.7071067811865475*vmap0R3*f[5])*volFact; 
} 

GKYL_CU_DH void gyrokinetic_M3_perp_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 

  out[0] += ((vmap[1]*bmag[3]*vmap[3]*f[15]+vmap[1]*bmag[2]*vmap[3]*f[14]+bmag[1]*vmap[1]*vmap[3]*f[13]+vmap[0]*bmag[3]*vmap[3]*f[12]+vmap[1]*vmap[2]*bmag[3]*f[11]+bmag[0]*vmap[1]*vmap[3]*f[10]+vmap[0]*bmag[2]*vmap[3]*f[9]+vmap[0]*bmag[1]*vmap[3]*f[8]+vmap[1]*bmag[2]*vmap[2]*f[7]+bmag[1]*vmap[1]*vmap[2]*f[6]+vmap[0]*vmap[2]*bmag[3]*f[5]+bmag[0]*vmap[0]*vmap[3]*f[4]+bmag[0]*vmap[1]*vmap[2]*f[3]+vmap[0]*bmag[2]*f[2]*vmap[2]+vmap[0]*bmag[1]*f[1]*vmap[2]+bmag[0]*f[0]*vmap[0]*vmap[2])*volFact)/m_; 
  out[1] += ((vmap[1]*bmag[2]*vmap[3]*f[15]+vmap[1]*bmag[3]*vmap[3]*f[14]+bmag[0]*vmap[1]*vmap[3]*f[13]+vmap[0]*bmag[2]*vmap[3]*f[12]+vmap[1]*bmag[2]*vmap[2]*f[11]+bmag[1]*vmap[1]*vmap[3]*f[10]+vmap[0]*bmag[3]*vmap[3]*f[9]+bmag[0]*vmap[0]*vmap[3]*f[8]+vmap[1]*vmap[2]*bmag[3]*f[7]+bmag[0]*vmap[1]*vmap[2]*f[6]+vmap[0]*bmag[2]*vmap[2]*f[5]+vmap[0]*bmag[1]*vmap[3]*f[4]+bmag[1]*vmap[1]*vmap[2]*f[3]+vmap[0]*f[2]*vmap[2]*bmag[3]+bmag[0]*vmap[0]*f[1]*vmap[2]+f[0]*vmap[0]*bmag[1]*vmap[2])*volFact)/m_; 
  out[2] += ((bmag[1]*vmap[1]*vmap[3]*f[15]+bmag[0]*vmap[1]*vmap[3]*f[14]+vmap[1]*bmag[3]*vmap[3]*f[13]+vmap[0]*bmag[1]*vmap[3]*f[12]+bmag[1]*vmap[1]*vmap[2]*f[11]+vmap[1]*bmag[2]*vmap[3]*f[10]+bmag[0]*vmap[0]*vmap[3]*f[9]+vmap[0]*bmag[3]*vmap[3]*f[8]+bmag[0]*vmap[1]*vmap[2]*f[7]+vmap[1]*vmap[2]*bmag[3]*f[6]+vmap[0]*bmag[1]*vmap[2]*f[5]+vmap[0]*bmag[2]*vmap[3]*f[4]+vmap[1]*bmag[2]*vmap[2]*f[3]+vmap[0]*f[1]*vmap[2]*bmag[3]+bmag[0]*vmap[0]*f[2]*vmap[2]+f[0]*vmap[0]*bmag[2]*vmap[2])*volFact)/m_; 
  out[3] += ((bmag[0]*vmap[1]*vmap[3]*f[15]+bmag[1]*vmap[1]*vmap[3]*f[14]+vmap[1]*bmag[2]*vmap[3]*f[13]+bmag[0]*vmap[0]*vmap[3]*f[12]+bmag[0]*vmap[1]*vmap[2]*f[11]+vmap[1]*bmag[3]*vmap[3]*f[10]+vmap[0]*bmag[1]*vmap[3]*f[9]+vmap[0]*bmag[2]*vmap[3]*f[8]+bmag[1]*vmap[1]*vmap[2]*f[7]+vmap[1]*bmag[2]*vmap[2]*f[6]+bmag[0]*vmap[0]*vmap[2]*f[5]+vmap[0]*bmag[3]*vmap[3]*f[4]+vmap[1]*vmap[2]*bmag[3]*f[3]+f[0]*vmap[0]*vmap[2]*bmag[3]+vmap[0]*bmag[1]*f[2]*vmap[2]+vmap[0]*f[1]*bmag[2]*vmap[2])*volFact)/m_; 
} 

GKYL_CU_DH void gyrokinetic_three_moments_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp[4]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[12]+1.4142135623730951*vmap[2]*f[5]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += (1.4142135623730951*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[5] += (1.4142135623730951*vmap[1]*f[6]+1.4142135623730951*vmap[0]*f[1])*volFact; 
  out[6] += (1.4142135623730951*vmap[1]*f[7]+1.4142135623730951*vmap[0]*f[2])*volFact; 
  out[7] += (1.4142135623730951*vmap[1]*f[11]+1.4142135623730951*vmap[0]*f[5])*volFact; 
  out[8] += ((bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact)/m_+(0.8944271909999159*vmap1R2*f[16]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[9] += ((bmag[2]*tmp[3]+tmp[2]*bmag[3]+bmag[0]*tmp[1]+tmp[0]*bmag[1])*volFact)/m_+(0.8944271909999161*vmap1R2*f[17]+2.0*vmap[0]*vmap[1]*f[6]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  out[10] += ((bmag[1]*tmp[3]+tmp[1]*bmag[3]+bmag[0]*tmp[2]+tmp[0]*bmag[2])*volFact)/m_+(0.8944271909999161*vmap1R2*f[18]+2.0*vmap[0]*vmap[1]*f[7]+vmap1R2*f[2]+vmap0R2*f[2])*volFact; 
  out[11] += ((bmag[0]*tmp[3]+tmp[0]*bmag[3]+bmag[1]*tmp[2]+tmp[1]*bmag[2])*volFact)/m_+(0.8944271909999159*vmap1R2*f[20]+2.0*vmap[0]*vmap[1]*f[11]+vmap1R2*f[5]+vmap0R2*f[5])*volFact; 
} 

GKYL_CU_DH void gyrokinetic_four_moments_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp[4]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[12]+1.4142135623730951*vmap[2]*f[5]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += (1.4142135623730951*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[5] += (1.4142135623730951*vmap[1]*f[6]+1.4142135623730951*vmap[0]*f[1])*volFact; 
  out[6] += (1.4142135623730951*vmap[1]*f[7]+1.4142135623730951*vmap[0]*f[2])*volFact; 
  out[7] += (1.4142135623730951*vmap[1]*f[11]+1.4142135623730951*vmap[0]*f[5])*volFact; 
  out[8] += (0.8944271909999159*vmap1R2*f[16]+2.0*vmap[0]*vmap[1]*f[3]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[9] += (0.8944271909999161*vmap1R2*f[17]+2.0*vmap[0]*vmap[1]*f[6]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  out[10] += (0.8944271909999161*vmap1R2*f[18]+2.0*vmap[0]*vmap[1]*f[7]+vmap1R2*f[2]+vmap0R2*f[2])*volFact; 
  out[11] += (0.8944271909999159*vmap1R2*f[20]+2.0*vmap[0]*vmap[1]*f[11]+vmap1R2*f[5]+vmap0R2*f[5])*volFact; 
  out[12] += ((bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact)/m_; 
  out[13] += ((bmag[2]*tmp[3]+tmp[2]*bmag[3]+bmag[0]*tmp[1]+tmp[0]*bmag[1])*volFact)/m_; 
  out[14] += ((bmag[1]*tmp[3]+tmp[1]*bmag[3]+bmag[0]*tmp[2]+tmp[0]*bmag[2])*volFact)/m_; 
  out[15] += ((bmag[0]*tmp[3]+tmp[0]*bmag[3]+bmag[1]*tmp[2]+tmp[1]*bmag[2])*volFact)/m_; 
} 

GKYL_CU_DH void gyrokinetic_M0_step1_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  out[0] += 1.4142135623730951*f[0]*volFact; 
  out[1] += 1.4142135623730951*f[1]*volFact; 
  out[2] += 1.4142135623730951*f[2]*volFact; 
  out[3] += 1.4142135623730951*f[4]*volFact; 
  out[4] += 1.4142135623730951*f[5]*volFact; 
  out[5] += 1.4142135623730951*f[8]*volFact; 
  out[6] += 1.4142135623730951*f[9]*volFact; 
  out[7] += 1.4142135623730951*f[12]*volFact; 
} 

GKYL_CU_DH void gyrokinetic_M0_step2_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 3.141592653589793*dxv[3]/m_; 
  out[0] += 2.8284271247461907*f[0]*volFact; 
  out[1] += 2.8284271247461907*f[1]*volFact; 
  out[2] += 2.8284271247461907*f[2]*volFact; 
  out[3] += 2.8284271247461907*f[4]*volFact; 
} 

GKYL_CU_DH void gyrokinetic_hamiltonian_moments_2x2v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.5707963267948966*dxv[2]*dxv[3]/m_; 
  double tmp[4]; 
  tmp[0] = 1.4142135623730951*vmap[3]*f[4]+1.4142135623730951*f[0]*vmap[2]; 
  tmp[1] = 1.4142135623730951*vmap[3]*f[8]+1.4142135623730951*f[1]*vmap[2]; 
  tmp[2] = 1.4142135623730951*vmap[3]*f[9]+1.4142135623730951*f[2]*vmap[2]; 
  tmp[3] = 1.4142135623730951*vmap[3]*f[12]+1.4142135623730951*vmap[2]*f[5]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += (1.4142135623730951*vmap[1]*f[3]+1.4142135623730951*f[0]*vmap[0])*m_*volFact; 
  out[5] += (1.4142135623730951*vmap[1]*f[6]+1.4142135623730951*vmap[0]*f[1])*m_*volFact; 
  out[6] += (1.4142135623730951*vmap[1]*f[7]+1.4142135623730951*vmap[0]*f[2])*m_*volFact; 
  out[7] += (1.4142135623730951*vmap[1]*f[11]+1.4142135623730951*vmap[0]*f[5])*m_*volFact; 
  out[8] += (phi[3]*f[5]+f[2]*phi[2]+f[1]*phi[1]+f[0]*phi[0])*q_*volFact+(0.4472135954999579*vmap1R2*f[16]+vmap[0]*vmap[1]*f[3]+0.5*f[0]*vmap1R2+0.5*f[0]*vmap0R2)*m_*volFact+(0.5*bmag[3]*tmp[3]+0.5*bmag[2]*tmp[2]+0.5*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0])*volFact; 
  out[9] += (phi[2]*f[5]+f[2]*phi[3]+f[0]*phi[1]+phi[0]*f[1])*q_*volFact+(0.44721359549995804*vmap1R2*f[17]+vmap[0]*vmap[1]*f[6]+0.5*f[1]*vmap1R2+0.5*vmap0R2*f[1])*m_*volFact+(0.5*bmag[2]*tmp[3]+0.5*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[1]+0.5*tmp[0]*bmag[1])*volFact; 
  out[10] += (phi[1]*f[5]+f[1]*phi[3]+f[0]*phi[2]+phi[0]*f[2])*q_*volFact+(0.44721359549995804*vmap1R2*f[18]+vmap[0]*vmap[1]*f[7]+0.5*vmap1R2*f[2]+0.5*vmap0R2*f[2])*m_*volFact+(0.5*bmag[1]*tmp[3]+0.5*tmp[1]*bmag[3]+0.5*bmag[0]*tmp[2]+0.5*tmp[0]*bmag[2])*volFact; 
  out[11] += (phi[0]*f[5]+f[0]*phi[3]+f[1]*phi[2]+phi[1]*f[2])*q_*volFact+(0.4472135954999579*vmap1R2*f[20]+vmap[0]*vmap[1]*f[11]+0.5*vmap1R2*f[5]+0.5*vmap0R2*f[5])*m_*volFact+(0.5*bmag[0]*tmp[3]+0.5*tmp[0]*bmag[3]+0.5*bmag[1]*tmp[2]+0.5*tmp[1]*bmag[2])*volFact; 
} 

