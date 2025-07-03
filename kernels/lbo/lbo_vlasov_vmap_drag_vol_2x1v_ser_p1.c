#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_vmap_drag_vol_2x1v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_vel_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
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

  const double rdvx2 = 2.0/dxv[2]; 
  const double *v0 = &vmap[0]; 
  const double *jacob_vel_inv0 = &jacob_vel_inv[0]; 

  double alphaDrag[12]; 
  alphaDrag[0] = (jacob_vel_inv0[0]*((-0.7071067811865475*nuSum[0]*v0[0])-1.0*nuUSum[0])-0.7071067811865475*nuSum[0]*(jacob_vel_inv0[2]*v0[2]+jacob_vel_inv0[1]*v0[1]))*rdvx2; 
  alphaDrag[1] = (jacob_vel_inv0[0]*((-1.0*nuUSum[1])-0.7071067811865475*v0[0]*nuSum[1])-0.7071067811865475*nuSum[1]*(jacob_vel_inv0[2]*v0[2]+jacob_vel_inv0[1]*v0[1]))*rdvx2; 
  alphaDrag[2] = ((-0.7071067811865475*jacob_vel_inv0[2]*nuSum[2]*v0[2])-1.0*jacob_vel_inv0[0]*nuUSum[2]-0.7071067811865475*(jacob_vel_inv0[1]*v0[1]+jacob_vel_inv0[0]*v0[0])*nuSum[2])*rdvx2; 
  alphaDrag[3] = (nuSum[0]*((-0.6210590034081186*jacob_vel_inv0[2]*v0[3])-0.6324555320336759*jacob_vel_inv0[1]*v0[2]+v0[1]*((-0.6324555320336759*jacob_vel_inv0[2])-0.7071067811865475*jacob_vel_inv0[0]))+((-0.7071067811865475*nuSum[0]*v0[0])-1.0*nuUSum[0])*jacob_vel_inv0[1])*rdvx2; 
  alphaDrag[4] = ((-1.0*jacob_vel_inv0[0]*nuUSum[3])-0.7071067811865475*(jacob_vel_inv0[2]*v0[2]+jacob_vel_inv0[1]*v0[1]+jacob_vel_inv0[0]*v0[0])*nuSum[3])*rdvx2; 
  alphaDrag[5] = (nuSum[1]*((-0.6210590034081186*jacob_vel_inv0[2]*v0[3])-0.6324555320336759*jacob_vel_inv0[1]*v0[2]+v0[1]*((-0.6324555320336759*jacob_vel_inv0[2])-0.7071067811865475*jacob_vel_inv0[0]))+jacob_vel_inv0[1]*((-1.0*nuUSum[1])-0.7071067811865475*v0[0]*nuSum[1]))*rdvx2; 
  alphaDrag[6] = ((-0.6210590034081186*jacob_vel_inv0[2]*nuSum[2]*v0[3])+jacob_vel_inv0[1]*((-0.6324555320336759*nuSum[2]*v0[2])-1.0*nuUSum[2])+((-0.6324555320336759*v0[1]*jacob_vel_inv0[2])-0.7071067811865475*(jacob_vel_inv0[0]*v0[1]+v0[0]*jacob_vel_inv0[1]))*nuSum[2])*rdvx2; 
  alphaDrag[7] = (nuSum[0]*(((-0.4517539514526256*jacob_vel_inv0[2])-0.7071067811865475*jacob_vel_inv0[0])*v0[2]-0.6210590034081186*jacob_vel_inv0[1]*v0[3])+((-0.7071067811865475*nuSum[0]*v0[0])-1.0*nuUSum[0])*jacob_vel_inv0[2]-0.6324555320336759*nuSum[0]*jacob_vel_inv0[1]*v0[1])*rdvx2; 
  alphaDrag[8] = ((-0.6210590034081186*jacob_vel_inv0[2]*nuSum[3]*v0[3])-1.0*jacob_vel_inv0[1]*nuUSum[3]+((-0.6324555320336759*(jacob_vel_inv0[1]*v0[2]+v0[1]*jacob_vel_inv0[2]))-0.7071067811865475*(jacob_vel_inv0[0]*v0[1]+v0[0]*jacob_vel_inv0[1]))*nuSum[3])*rdvx2; 
  alphaDrag[9] = (nuSum[1]*(((-0.4517539514526256*jacob_vel_inv0[2])-0.7071067811865475*jacob_vel_inv0[0])*v0[2]-0.6210590034081186*jacob_vel_inv0[1]*v0[3])-1.0*nuUSum[1]*jacob_vel_inv0[2]+nuSum[1]*((-0.7071067811865475*v0[0]*jacob_vel_inv0[2])-0.6324555320336759*jacob_vel_inv0[1]*v0[1]))*rdvx2; 
  alphaDrag[10] = (nuSum[2]*(((-0.4517539514526256*jacob_vel_inv0[2])-0.7071067811865475*jacob_vel_inv0[0])*v0[2]-0.6210590034081186*jacob_vel_inv0[1]*v0[3])-1.0*jacob_vel_inv0[2]*nuUSum[2]+((-0.7071067811865475*v0[0]*jacob_vel_inv0[2])-0.6324555320336759*jacob_vel_inv0[1]*v0[1])*nuSum[2])*rdvx2; 
  alphaDrag[11] = ((-0.6210590034081186*jacob_vel_inv0[1]*nuSum[3]*v0[3])-1.0*jacob_vel_inv0[2]*nuUSum[3]+((-0.7071067811865475*(jacob_vel_inv0[0]*v0[2]+v0[0]*jacob_vel_inv0[2]))-0.4517539514526256*jacob_vel_inv0[2]*v0[2]-0.6324555320336759*jacob_vel_inv0[1]*v0[1])*nuSum[3])*rdvx2; 

  out[3] += 0.6123724356957944*(alphaDrag[11]*f[11]+alphaDrag[10]*f[10]+alphaDrag[9]*f[9]+alphaDrag[8]*f[8]+alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[5] += 0.6123724356957944*(alphaDrag[10]*f[11]+f[10]*alphaDrag[11]+alphaDrag[7]*f[9]+f[7]*alphaDrag[9]+alphaDrag[6]*f[8]+f[6]*alphaDrag[8]+alphaDrag[3]*f[5]+f[3]*alphaDrag[5]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[6] += 0.6123724356957944*(alphaDrag[9]*f[11]+f[9]*alphaDrag[11]+alphaDrag[7]*f[10]+f[7]*alphaDrag[10]+alphaDrag[5]*f[8]+f[5]*alphaDrag[8]+alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[7] += 1.224744871391589*(alphaDrag[8]*f[11]+f[8]*alphaDrag[11]+alphaDrag[6]*f[10]+f[6]*alphaDrag[10]+alphaDrag[5]*f[9]+f[5]*alphaDrag[9])+1.369306393762915*(alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+1.224744871391589*(alphaDrag[3]*f[7]+f[3]*alphaDrag[7])+1.369306393762915*(alphaDrag[2]*f[6]+f[2]*alphaDrag[6]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[8] += 0.6123724356957944*(alphaDrag[7]*f[11]+f[7]*alphaDrag[11]+alphaDrag[9]*f[10]+f[9]*alphaDrag[10]+alphaDrag[3]*f[8]+f[3]*alphaDrag[8]+alphaDrag[5]*f[6]+f[5]*alphaDrag[6]+alphaDrag[0]*f[4]+f[0]*alphaDrag[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[9] += 1.224744871391589*(alphaDrag[6]*f[11]+f[6]*alphaDrag[11]+alphaDrag[8]*f[10]+f[8]*alphaDrag[10]+alphaDrag[3]*f[9]+f[3]*alphaDrag[9])+1.369306393762915*(alphaDrag[2]*f[8]+f[2]*alphaDrag[8])+1.224744871391589*(alphaDrag[5]*f[7]+f[5]*alphaDrag[7])+1.369306393762915*(alphaDrag[4]*f[6]+f[4]*alphaDrag[6]+alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[10] += 1.224744871391589*(alphaDrag[5]*f[11]+f[5]*alphaDrag[11]+alphaDrag[3]*f[10]+f[3]*alphaDrag[10]+alphaDrag[8]*f[9]+f[8]*alphaDrag[9])+1.369306393762915*(alphaDrag[1]*f[8]+f[1]*alphaDrag[8])+1.224744871391589*(alphaDrag[6]*f[7]+f[6]*alphaDrag[7])+1.369306393762915*(alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[4]*f[5]+f[4]*alphaDrag[5]+alphaDrag[2]*f[3]+f[2]*alphaDrag[3]); 
  out[11] += 1.224744871391589*(alphaDrag[3]*f[11]+f[3]*alphaDrag[11]+alphaDrag[5]*f[10]+f[5]*alphaDrag[10]+alphaDrag[6]*f[9]+f[6]*alphaDrag[9])+(1.224744871391589*alphaDrag[7]+1.369306393762915*alphaDrag[0])*f[8]+1.224744871391589*f[7]*alphaDrag[8]+1.369306393762915*(f[0]*alphaDrag[8]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[3]*f[4]+f[3]*alphaDrag[4]); 

  return fabs(0.8838834764831842*alphaDrag[0]-0.9882117688026182*alphaDrag[7]); 

} 
