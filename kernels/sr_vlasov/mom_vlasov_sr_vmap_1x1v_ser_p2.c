#include <gkyl_mom_vlasov_sr_kernels.h> 
GKYL_CU_DH void vlasov_sr_vmap_M1i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *jacob_vel_inv, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double dv10 = 2.0/dxv[1]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  double p0_over_gamma[3] = {0.0}; 
  p0_over_gamma[0] = 2.738612787525831*jacob_vel_inv0[1]*gamma[2]*dv10+1.224744871391589*jacob_vel_inv0[0]*gamma[1]*dv10; 
  p0_over_gamma[1] = 2.449489742783178*jacob_vel_inv0[2]*gamma[2]*dv10+2.738612787525831*jacob_vel_inv0[0]*gamma[2]*dv10+1.224744871391589*jacob_vel_inv0[1]*gamma[1]*dv10; 
  p0_over_gamma[2] = 2.449489742783178*jacob_vel_inv0[1]*gamma[2]*dv10+1.224744871391589*gamma[1]*jacob_vel_inv0[2]*dv10; 

  out[0] += (p0_over_gamma[2]*f[5]+p0_over_gamma[1]*f[2]+f[0]*p0_over_gamma[0])*volFact; 
  out[1] += (1.0*p0_over_gamma[2]*f[7]+p0_over_gamma[1]*f[3]+p0_over_gamma[0]*f[1])*volFact; 
  out[2] += (1.0*p0_over_gamma[1]*f[6]+p0_over_gamma[0]*f[4])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_vmap_M3i_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *vmap, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  const double *p0 = &vmap[0]; 
  out[0] += (p0[2]*f[5]+p0[1]*f[2]+f[0]*p0[0])*volFact; 
  out[1] += (1.0*p0[2]*f[7]+p0[1]*f[3]+p0[0]*f[1])*volFact; 
  out[2] += (1.0*p0[1]*f[6]+p0[0]*f[4])*volFact; 
} 
GKYL_CU_DH void vlasov_sr_vmap_int_mom_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, const double *vmap, const double *gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*0.25; 
  const double *p0 = &vmap[0]; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += (1.414213562373095*gamma[2]*f[5]+1.414213562373095*gamma[1]*f[2]+1.414213562373095*f[0]*gamma[0])*volFact; 
  out[2] += (1.414213562373095*p0[2]*f[5]+1.414213562373095*p0[1]*f[2]+1.414213562373095*f[0]*p0[0])*volFact; 
} 
