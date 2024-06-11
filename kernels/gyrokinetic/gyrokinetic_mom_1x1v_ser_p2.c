#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_M0_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[4]*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M1_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  out[0] += (vmap[1]*f[2]+f[0]*vmap[0])*volFact; 
  out[1] += (vmap[1]*f[3]+vmap[0]*f[1])*volFact; 
  out[2] += (1.0*vmap[1]*f[6]+vmap[0]*f[4])*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M2_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.6324555320336759*vmap1R2*f[5]+1.414213562373095*vmap[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2)*volFact; 
  out[1] += (0.632455532033676*vmap1R2*f[7]+1.414213562373095*vmap[0]*vmap[1]*f[3]+0.7071067811865475*f[1]*vmap1R2+0.7071067811865475*vmap0R2*f[1])*volFact; 
  out[2] += (1.414213562373095*vmap[0]*vmap[1]*f[6]+0.7071067811865475*vmap1R2*f[4]+0.7071067811865475*vmap0R2*f[4])*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M2_par_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += (0.6324555320336759*vmap1R2*f[5]+1.414213562373095*vmap[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2)*volFact; 
  out[1] += (0.632455532033676*vmap1R2*f[7]+1.414213562373095*vmap[0]*vmap[1]*f[3]+0.7071067811865475*f[1]*vmap1R2+0.7071067811865475*vmap0R2*f[1])*volFact; 
  out[2] += (1.414213562373095*vmap[0]*vmap[1]*f[6]+0.7071067811865475*vmap1R2*f[4]+0.7071067811865475*vmap0R2*f[4])*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M3_par_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap0R3 = pow(vmap[0],3);
  const double vmap1R2 = pow(vmap[1],2);
  const double vmap1R3 = pow(vmap[1],3);

  out[0] += (1.341640786499874*vmap[0]*vmap1R2*f[5]+0.9*vmap1R3*f[2]+1.5*vmap0R2*vmap[1]*f[2]+1.5*f[0]*vmap[0]*vmap1R2+0.5*f[0]*vmap0R3)*volFact; 
  out[1] += (1.341640786499874*vmap[0]*vmap1R2*f[7]+0.9*vmap1R3*f[3]+1.5*vmap0R2*vmap[1]*f[3]+1.5*vmap[0]*f[1]*vmap1R2+0.5*vmap0R3*f[1])*volFact; 
  out[2] += (0.8999999999999998*vmap1R3*f[6]+1.5*vmap0R2*vmap[1]*f[6]+1.5*vmap[0]*vmap1R2*f[4]+0.5*vmap0R3*f[4])*volFact; 
} 
GKYL_CU_DH void gyrokinetic_three_moments_1x1v_ser_p2(const double *dxv, const double *vmap, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 0.5*dxv[1]; 
  const double vmap0R2 = pow(vmap[0],2);
  const double vmap1R2 = pow(vmap[1],2);

  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[4]*volFact; 
  out[3] += (vmap[1]*f[2]+f[0]*vmap[0])*volFact; 
  out[4] += (vmap[1]*f[3]+vmap[0]*f[1])*volFact; 
  out[5] += (1.0*vmap[1]*f[6]+vmap[0]*f[4])*volFact; 
  out[6] += (0.6324555320336759*vmap1R2*f[5]+1.414213562373095*vmap[0]*vmap[1]*f[2]+0.7071067811865475*f[0]*vmap1R2+0.7071067811865475*f[0]*vmap0R2)*volFact; 
  out[7] += (0.632455532033676*vmap1R2*f[7]+1.414213562373095*vmap[0]*vmap[1]*f[3]+0.7071067811865475*f[1]*vmap1R2+0.7071067811865475*vmap0R2*f[1])*volFact; 
  out[8] += (1.414213562373095*vmap[0]*vmap[1]*f[6]+0.7071067811865475*vmap1R2*f[4]+0.7071067811865475*vmap0R2*f[4])*volFact; 
} 

