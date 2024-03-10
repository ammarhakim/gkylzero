#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_M0_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.570796326794897*dxv[1]*dxv[2]/m_; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M1_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.570796326794897*dxv[1]*dxv[2]/m_; 
  out[0] += (1.414213562373095*vmap[1]*f[2]+1.414213562373095*f[0]*vmap[0])*volFact; 
  out[1] += (1.414213562373095*vmap[1]*f[4]+1.414213562373095*vmap[0]*f[1])*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M2_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.570796326794897*dxv[1]*dxv[2]/m_; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[1] += (0.8944271909999161*vmap1R2*f[9]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
  double tmp[2]; 
  tmp[0] = 1.414213562373095*f[3]*vmap[3]+1.414213562373095*f[0]*vmap[2]; 
  tmp[1] = 1.414213562373095*vmap[3]*f[5]+1.414213562373095*f[1]*vmap[2]; 
  out[0] += (2.0*(0.7071067811865475*bmag[1]*tmp[1]+0.7071067811865475*bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.7071067811865475*bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*bmag[1])*volFact)/m_; 
} 
GKYL_CU_DH void gyrokinetic_M2_par_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.570796326794897*dxv[1]*dxv[2]/m_; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[1] += (0.8944271909999161*vmap1R2*f[9]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M2_perp_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.570796326794897*dxv[1]*dxv[2]/m_; 
  double tmp[2]; 
  tmp[0] = 1.414213562373095*f[3]*vmap[3]+1.414213562373095*f[0]*vmap[2]; 
  tmp[1] = 1.414213562373095*vmap[3]*f[5]+1.414213562373095*f[1]*vmap[2]; 
  out[0] += (2.0*(0.7071067811865475*bmag[1]*tmp[1]+0.7071067811865475*bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.7071067811865475*bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*bmag[1])*volFact)/m_; 
} 
GKYL_CU_DH void gyrokinetic_M3_par_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.570796326794897*dxv[1]*dxv[2]/m_; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += (1.897366596101028*vmap[0]*vmap1R2*f[8]+1.272792206135785*vmap1R3*f[2]+2.121320343559642*vmap0R2*vmap[1]*f[2]+2.121320343559642*f[0]*vmap[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R3)*volFact; 
  out[1] += (1.897366596101028*vmap[0]*vmap1R2*f[9]+1.272792206135785*vmap1R3*f[4]+2.121320343559642*vmap0R2*vmap[1]*f[4]+2.121320343559642*vmap[0]*f[1]*vmap1R2+0.7071067811865475*vmap0R3*f[1])*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M3_perp_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.570796326794897*dxv[1]*dxv[2]/m_; 

  out[0] += ((1.414213562373095*bmag[1]*vmap[1]*vmap[3]*f[7]+1.414213562373095*bmag[0]*vmap[1]*vmap[3]*f[6]+1.414213562373095*vmap[0]*bmag[1]*vmap[3]*f[5]+1.414213562373095*bmag[1]*vmap[1]*vmap[2]*f[4]+1.414213562373095*bmag[0]*vmap[0]*f[3]*vmap[3]+1.414213562373095*bmag[0]*vmap[1]*f[2]*vmap[2]+1.414213562373095*vmap[0]*bmag[1]*f[1]*vmap[2]+1.414213562373095*bmag[0]*f[0]*vmap[0]*vmap[2])*volFact)/m_; 
  out[1] += ((1.414213562373095*bmag[0]*vmap[1]*vmap[3]*f[7]+1.414213562373095*bmag[1]*vmap[1]*vmap[3]*f[6]+1.414213562373095*bmag[0]*vmap[0]*vmap[3]*f[5]+1.414213562373095*bmag[0]*vmap[1]*vmap[2]*f[4]+1.414213562373095*vmap[0]*bmag[1]*f[3]*vmap[3]+1.414213562373095*bmag[1]*vmap[1]*f[2]*vmap[2]+1.414213562373095*bmag[0]*vmap[0]*f[1]*vmap[2]+1.414213562373095*f[0]*vmap[0]*bmag[1]*vmap[2])*volFact)/m_; 
} 
GKYL_CU_DH void gyrokinetic_three_moments_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 1.570796326794897*dxv[1]*dxv[2]/m_; 
  double tmp[2]; 
  tmp[0] = 1.414213562373095*f[3]*vmap[3]+1.414213562373095*f[0]*vmap[2]; 
  tmp[1] = 1.414213562373095*vmap[3]*f[5]+1.414213562373095*f[1]*vmap[2]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += (1.414213562373095*vmap[1]*f[2]+1.414213562373095*f[0]*vmap[0])*volFact; 
  out[3] += (1.414213562373095*vmap[1]*f[4]+1.414213562373095*vmap[0]*f[1])*volFact; 
  out[4] += ((1.414213562373095*bmag[1]*tmp[1]+1.414213562373095*bmag[0]*tmp[0])*volFact)/m_+(0.8944271909999159*vmap1R2*f[8]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
  out[5] += ((1.414213562373095*bmag[0]*tmp[1]+1.414213562373095*tmp[0]*bmag[1])*volFact)/m_+(0.8944271909999161*vmap1R2*f[9]+2.0*vmap[0]*vmap[1]*f[4]+f[1]*vmap1R2+vmap0R2*f[1])*volFact; 
} 

GKYL_CU_DH void gyrokinetic_M0_step1_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[3]*volFact; 
  out[3] += 1.414213562373095*f[5]*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M0_step2_1x2v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]/2; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
} 
