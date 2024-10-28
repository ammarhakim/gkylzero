#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_vmap_drag_vol_1x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_vel_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]: Cell-center coordinates. 
  // dxv[3]: Cell spacing. 
  // vmap: Velocity-space nonuniform mapping in each dimension. 
  // jacob_vel_inv: Inverse of velocity space Jacobian in each dimension. 
  // nuSum: collisionalities added (self and cross species collisionalities). 
  // nuPrimMomsSum: sum of bulk velocities and thermal speeds (squared) times their respective collisionalities. 
  // f: Input distribution function.
  // out: Incremented output 
  const double *nuUSum = nuPrimMomsSum;

  const double rdvx2 = 2.0/dxv[1]; 
  const double *v0 = &vmap[0]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 
  const double rdvy2 = 2.0/dxv[2]; 
  const double *v1 = &vmap[4]; 
  const double *jacob_vel_inv1 = &jacob_vel_inv[3]; 

  double alphaDrag[32]; 
  alphaDrag[0] = (jacob_vel_inv0[0]*((-1.0*nuSum[0]*v0[0])-1.414213562373095*nuUSum[0])-1.0*nuSum[0]*(jacob_vel_inv0[2]*v0[2]+jacob_vel_inv0[1]*v0[1]))*rdvx2; 
  alphaDrag[1] = (jacob_vel_inv0[0]*((-1.414213562373095*nuUSum[1])-1.0*v0[0]*nuSum[1])-1.0*nuSum[1]*(jacob_vel_inv0[2]*v0[2]+jacob_vel_inv0[1]*v0[1]))*rdvx2; 
  alphaDrag[2] = (nuSum[0]*((-0.8783100656536796*jacob_vel_inv0[2]*v0[3])-0.8944271909999159*jacob_vel_inv0[1]*v0[2]+v0[1]*((-0.8944271909999159*jacob_vel_inv0[2])-1.0*jacob_vel_inv0[0]))+((-1.0*nuSum[0]*v0[0])-1.414213562373095*nuUSum[0])*jacob_vel_inv0[1])*rdvx2; 
  alphaDrag[4] = (nuSum[1]*((-0.8783100656536796*jacob_vel_inv0[2]*v0[3])-0.8944271909999159*jacob_vel_inv0[1]*v0[2]+v0[1]*((-0.8944271909999159*jacob_vel_inv0[2])-1.0*jacob_vel_inv0[0]))+jacob_vel_inv0[1]*((-1.414213562373095*nuUSum[1])-1.0*v0[0]*nuSum[1]))*rdvx2; 
  alphaDrag[8] = (nuSum[0]*(((-0.6388765649999399*jacob_vel_inv0[2])-1.0*jacob_vel_inv0[0])*v0[2]-0.8783100656536796*jacob_vel_inv0[1]*v0[3])+((-1.0*nuSum[0]*v0[0])-1.414213562373095*nuUSum[0])*jacob_vel_inv0[2]-0.8944271909999159*nuSum[0]*jacob_vel_inv0[1]*v0[1])*rdvx2; 
  alphaDrag[9] = (nuSum[1]*(((-0.6388765649999399*jacob_vel_inv0[2])-1.0*jacob_vel_inv0[0])*v0[2]-0.8783100656536796*jacob_vel_inv0[1]*v0[3])-1.414213562373095*nuUSum[1]*jacob_vel_inv0[2]+nuSum[1]*((-1.0*v0[0]*jacob_vel_inv0[2])-0.8944271909999159*jacob_vel_inv0[1]*v0[1]))*rdvx2; 

  alphaDrag[16] = ((-1.0*nuSum[0]*jacob_vel_inv1[2]*v1[2])-1.414213562373095*jacob_vel_inv1[0]*nuUSum[2]-1.0*nuSum[0]*(jacob_vel_inv1[1]*v1[1]+jacob_vel_inv1[0]*v1[0]))*rdvy2; 
  alphaDrag[17] = ((-1.414213562373095*jacob_vel_inv1[0]*nuUSum[3])-1.0*nuSum[1]*(jacob_vel_inv1[2]*v1[2]+jacob_vel_inv1[1]*v1[1]+jacob_vel_inv1[0]*v1[0]))*rdvy2; 
  alphaDrag[19] = ((-0.8783100656536796*nuSum[0]*jacob_vel_inv1[2]*v1[3])+jacob_vel_inv1[1]*((-0.8944271909999159*nuSum[0]*v1[2])-1.414213562373095*nuUSum[2])+nuSum[0]*((-0.8944271909999159*v1[1]*jacob_vel_inv1[2])-1.0*(jacob_vel_inv1[0]*v1[1]+v1[0]*jacob_vel_inv1[1])))*rdvy2; 
  alphaDrag[21] = ((-0.8783100656536796*nuSum[1]*jacob_vel_inv1[2]*v1[3])-1.414213562373095*jacob_vel_inv1[1]*nuUSum[3]+nuSum[1]*((-0.8944271909999159*(jacob_vel_inv1[1]*v1[2]+v1[1]*jacob_vel_inv1[2]))-1.0*(jacob_vel_inv1[0]*v1[1]+v1[0]*jacob_vel_inv1[1])))*rdvy2; 
  alphaDrag[28] = (nuSum[0]*(((-0.6388765649999399*jacob_vel_inv1[2])-1.0*jacob_vel_inv1[0])*v1[2]-0.8783100656536796*jacob_vel_inv1[1]*v1[3])-1.414213562373095*jacob_vel_inv1[2]*nuUSum[2]+nuSum[0]*((-1.0*v1[0]*jacob_vel_inv1[2])-0.8944271909999159*jacob_vel_inv1[1]*v1[1]))*rdvy2; 
  alphaDrag[29] = ((-0.8783100656536796*jacob_vel_inv1[1]*nuSum[1]*v1[3])-1.414213562373095*jacob_vel_inv1[2]*nuUSum[3]+nuSum[1]*((-1.0*(jacob_vel_inv1[0]*v1[2]+v1[0]*jacob_vel_inv1[2]))-0.6388765649999399*jacob_vel_inv1[2]*v1[2]-0.8944271909999159*jacob_vel_inv1[1]*v1[1]))*rdvy2; 

  out[2] += 0.6123724356957944*(alphaDrag[9]*f[9]+alphaDrag[8]*f[8]+alphaDrag[4]*f[4]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.6123724356957944*(f[13]*alphaDrag[29]+f[12]*alphaDrag[28]+f[5]*alphaDrag[21]+f[3]*alphaDrag[19]+f[1]*alphaDrag[17]+f[0]*alphaDrag[16]); 
  out[4] += 0.6123724356957944*(alphaDrag[8]*f[9]+f[8]*alphaDrag[9]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 0.6123724356957944*(f[12]*alphaDrag[29]+f[13]*alphaDrag[28]+f[3]*alphaDrag[21]+f[5]*alphaDrag[19]+f[0]*alphaDrag[17]+f[1]*alphaDrag[16]); 
  out[6] += 0.6123724356957944*(f[15]*alphaDrag[29]+f[14]*alphaDrag[28]+f[7]*alphaDrag[21]+f[6]*alphaDrag[19]+f[4]*alphaDrag[17]+f[2]*alphaDrag[16]+alphaDrag[9]*f[11]+alphaDrag[8]*f[10]+alphaDrag[4]*f[7]+alphaDrag[2]*f[6]+alphaDrag[1]*f[5]+alphaDrag[0]*f[3]); 
  out[7] += 0.6123724356957944*(f[14]*alphaDrag[29]+f[15]*alphaDrag[28]+f[6]*alphaDrag[21]+f[7]*alphaDrag[19]+f[2]*alphaDrag[17]+f[4]*alphaDrag[16]+alphaDrag[8]*f[11]+alphaDrag[9]*f[10]+alphaDrag[2]*f[7]+alphaDrag[4]*f[6]+alphaDrag[0]*f[5]+alphaDrag[1]*f[3]); 
  out[8] += 1.224744871391589*(alphaDrag[4]*f[9]+f[4]*alphaDrag[9]+alphaDrag[2]*f[8]+f[2]*alphaDrag[8])+1.369306393762915*(alphaDrag[1]*f[4]+f[1]*alphaDrag[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[9] += 1.224744871391589*(alphaDrag[2]*f[9]+f[2]*alphaDrag[9]+alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+1.369306393762915*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[10] += 0.6123724356957944*(f[11]*alphaDrag[21]+f[10]*alphaDrag[19]+f[9]*alphaDrag[17]+f[8]*alphaDrag[16])+1.224744871391589*(alphaDrag[4]*f[11]+alphaDrag[2]*f[10]+f[7]*alphaDrag[9]+f[6]*alphaDrag[8])+1.369306393762915*(alphaDrag[1]*f[7]+alphaDrag[0]*f[6]+alphaDrag[4]*f[5]+alphaDrag[2]*f[3]); 
  out[11] += 0.6123724356957944*(f[10]*alphaDrag[21]+f[11]*alphaDrag[19]+f[8]*alphaDrag[17]+f[9]*alphaDrag[16])+1.224744871391589*(alphaDrag[2]*f[11]+alphaDrag[4]*f[10]+f[6]*alphaDrag[9]+f[7]*alphaDrag[8])+1.369306393762915*(alphaDrag[0]*f[7]+alphaDrag[1]*f[6]+alphaDrag[2]*f[5]+f[3]*alphaDrag[4]); 
  out[12] += 1.224744871391589*(f[5]*alphaDrag[29]+f[3]*alphaDrag[28])+(1.224744871391589*f[13]+1.369306393762915*f[1])*alphaDrag[21]+1.224744871391589*f[12]*alphaDrag[19]+1.369306393762915*(f[0]*alphaDrag[19]+f[5]*alphaDrag[17]+f[3]*alphaDrag[16]); 
  out[13] += 1.224744871391589*(f[3]*alphaDrag[29]+f[5]*alphaDrag[28])+(1.224744871391589*f[12]+1.369306393762915*f[0])*alphaDrag[21]+1.224744871391589*f[13]*alphaDrag[19]+1.369306393762915*(f[1]*alphaDrag[19]+f[3]*alphaDrag[17]+f[5]*alphaDrag[16]); 
  out[14] += 1.224744871391589*(f[7]*alphaDrag[29]+f[6]*alphaDrag[28])+(1.224744871391589*f[15]+1.369306393762915*f[4])*alphaDrag[21]+1.224744871391589*f[14]*alphaDrag[19]+1.369306393762915*(f[2]*alphaDrag[19]+f[7]*alphaDrag[17]+f[6]*alphaDrag[16])+0.6123724356957944*(alphaDrag[4]*f[15]+alphaDrag[2]*f[14]+alphaDrag[1]*f[13]+alphaDrag[0]*f[12]); 
  out[15] += 1.224744871391589*(f[6]*alphaDrag[29]+f[7]*alphaDrag[28])+(1.224744871391589*f[14]+1.369306393762915*f[2])*alphaDrag[21]+1.224744871391589*f[15]*alphaDrag[19]+1.369306393762915*(f[4]*alphaDrag[19]+f[6]*alphaDrag[17]+f[7]*alphaDrag[16])+0.6123724356957944*(alphaDrag[2]*f[15]+alphaDrag[4]*f[14]+alphaDrag[0]*f[13]+alphaDrag[1]*f[12]); 

  return fabs(0.8838834764831842*alphaDrag[0]-0.9882117688026182*alphaDrag[8])+fabs(0.8838834764831842*alphaDrag[16]-0.9882117688026182*alphaDrag[28]); 

} 
