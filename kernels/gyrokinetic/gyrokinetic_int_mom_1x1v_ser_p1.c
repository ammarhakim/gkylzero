#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_int_M0_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 


  out[0] += 2.0*f[0]*volFact; 
} 

GKYL_CU_DH void gyrokinetic_int_M1_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 

  out[0] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*volFact; 
} 

GKYL_CU_DH void gyrokinetic_int_M2_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.8944271909999159*vmap1R2*f[4]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
} 

GKYL_CU_DH void gyrokinetic_int_three_moments_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[2] += (0.8944271909999159*vmap1R2*f[4]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
} 

GKYL_CU_DH void gyrokinetic_int_four_moments_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*f[2]+1.4142135623730951*f[0]*vmap[0])*volFact; 
  out[2] += (0.8944271909999159*vmap1R2*f[4]+2.0*vmap[0]*vmap[1]*f[2]+f[0]*vmap1R2+f[0]*vmap0R2)*volFact; 
} 

GKYL_CU_DH void gyrokinetic_int_hamiltonian_moments_1x1v_ser_p1(const double *dxv, const double *vmap, double m_, double q_, const double *bmag, const double *phi, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.25*dxv[0]*dxv[1]; 
 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 2.0*f[0]*volFact; 
  out[1] += (1.4142135623730951*vmap[1]*f[2]*m_+1.4142135623730951*f[0]*vmap[0]*m_)*volFact; 
  out[2] += (1.4142135623730951*f[1]*phi[1]*q_+1.4142135623730951*f[0]*phi[0]*q_+0.4472135954999579*vmap1R2*f[4]*m_+vmap[0]*vmap[1]*f[2]*m_+0.5*f[0]*vmap1R2*m_+0.5*f[0]*vmap0R2*m_)*volFact; 
} 

