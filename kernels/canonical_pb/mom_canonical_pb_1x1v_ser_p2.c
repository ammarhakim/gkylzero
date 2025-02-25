#include <gkyl_mom_canonical_pb_kernels.h> 
GKYL_CU_DH void canonical_pb_M1i_from_H_1x1v_ser_p2(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double rdvx2 = 2.0/dxv[1];

  out[0] += (2.738612787525831*f[3]*hamil[7]*rdvx2+1.224744871391589*f[4]*hamil[6]*rdvx2+2.738612787525831*f[2]*hamil[5]*rdvx2+1.224744871391589*f[1]*hamil[3]*rdvx2+1.224744871391589*f[0]*hamil[2]*rdvx2)*volFact; 
  out[1] += (2.449489742783178*f[6]*hamil[7]*rdvx2+2.738612787525831*f[2]*hamil[7]*rdvx2+1.095445115010332*f[1]*hamil[6]*rdvx2+2.738612787525831*f[3]*hamil[5]*rdvx2+1.095445115010332*hamil[3]*f[4]*rdvx2+1.224744871391589*f[0]*hamil[3]*rdvx2+1.224744871391589*f[1]*hamil[2]*rdvx2)*volFact; 
  out[2] += (2.449489742783178*f[3]*hamil[7]*rdvx2+0.7824607964359517*f[4]*hamil[6]*rdvx2+1.224744871391589*f[0]*hamil[6]*rdvx2+2.738612787525831*hamil[5]*f[6]*rdvx2+1.224744871391589*hamil[2]*f[4]*rdvx2+1.095445115010332*f[1]*hamil[3]*rdvx2)*volFact; 
} 
GKYL_CU_DH void canonical_pb_MEnergy_1x1v_ser_p2(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double dv10 = 2.0/dxv[1]; 

  out[0] += (0.7071067811865475*f[7]*hamil[7]+0.7071067811865475*f[6]*hamil[6]+0.7071067811865475*f[5]*hamil[5]+0.7071067811865475*f[4]*hamil[4]+0.7071067811865475*f[3]*hamil[3]+0.7071067811865475*f[2]*hamil[2]+0.7071067811865475*f[1]*hamil[1]+0.7071067811865475*f[0]*hamil[0])*volFact; 
  out[1] += (0.7071067811865475*f[5]*hamil[7]+0.7071067811865475*hamil[5]*f[7]+0.632455532033676*f[3]*hamil[6]+0.632455532033676*hamil[3]*f[6]+0.6324555320336759*f[1]*hamil[4]+0.6324555320336759*hamil[1]*f[4]+0.7071067811865475*f[2]*hamil[3]+0.7071067811865475*hamil[2]*f[3]+0.7071067811865475*f[0]*hamil[1]+0.7071067811865475*hamil[0]*f[1])*volFact; 
  out[2] += (0.6324555320336759*f[7]*hamil[7]+0.4517539514526256*f[6]*hamil[6]+0.7071067811865475*f[2]*hamil[6]+0.7071067811865475*hamil[2]*f[6]+0.4517539514526256*f[4]*hamil[4]+0.7071067811865475*f[0]*hamil[4]+0.7071067811865475*hamil[0]*f[4]+0.6324555320336759*f[3]*hamil[3]+0.6324555320336759*f[1]*hamil[1])*volFact; 
} 
GKYL_CU_DH void canonical_pb_int_mom_1x1v_ser_p2(const double *dxv, const double *hamil, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*0.25; 
  const double dv1 = dxv[1]; 
  const double rdvx2 = 2.0/dxv[1];

  out[0] += 2.0*f[0]*volFact; 
  out[1] += (3.872983346207417*f[3]*hamil[7]*rdvx2+1.732050807568877*f[4]*hamil[6]*rdvx2+3.872983346207417*f[2]*hamil[5]*rdvx2+1.732050807568877*f[1]*hamil[3]*rdvx2+1.732050807568877*f[0]*hamil[2]*rdvx2)*volFact; 
  out[2] += (f[7]*hamil[7]+f[6]*hamil[6]+f[5]*hamil[5]+f[4]*hamil[4]+f[3]*hamil[3]+f[2]*hamil[2]+f[1]*hamil[1]+f[0]*hamil[0])*volFact; 
} 
