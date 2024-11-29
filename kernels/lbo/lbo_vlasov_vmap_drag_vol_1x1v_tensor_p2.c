#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_vmap_drag_vol_1x1v_tensor_p2(const double *w, const double *dxv, const double *vmap, const double *jacob_vel_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]: Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
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

  double alphaDrag[9]; 
  alphaDrag[0] = (jacob_vel_inv0[0]*((-0.7071067811865475*nuSum[0]*v0[0])-1.0*nuUSum[0])-0.7071067811865475*nuSum[0]*(jacob_vel_inv0[2]*v0[2]+jacob_vel_inv0[1]*v0[1]))*rdvx2; 
  alphaDrag[1] = (jacob_vel_inv0[0]*((-1.0*nuUSum[1])-0.7071067811865475*v0[0]*nuSum[1])-0.7071067811865475*nuSum[1]*(jacob_vel_inv0[2]*v0[2]+jacob_vel_inv0[1]*v0[1]))*rdvx2; 
  alphaDrag[2] = (nuSum[0]*((-0.6210590034081186*jacob_vel_inv0[2]*v0[3])-0.6324555320336759*jacob_vel_inv0[1]*v0[2]+v0[1]*((-0.6324555320336759*jacob_vel_inv0[2])-0.7071067811865475*jacob_vel_inv0[0]))+((-0.7071067811865475*nuSum[0]*v0[0])-1.0*nuUSum[0])*jacob_vel_inv0[1])*rdvx2; 
  alphaDrag[3] = (nuSum[1]*((-0.6210590034081186*jacob_vel_inv0[2]*v0[3])-0.6324555320336759*jacob_vel_inv0[1]*v0[2]+v0[1]*((-0.6324555320336759*jacob_vel_inv0[2])-0.7071067811865475*jacob_vel_inv0[0]))+jacob_vel_inv0[1]*((-1.0*nuUSum[1])-0.7071067811865475*v0[0]*nuSum[1]))*rdvx2; 
  alphaDrag[4] = ((-0.7071067811865475*jacob_vel_inv0[2]*nuSum[2]*v0[2])-1.0*jacob_vel_inv0[0]*nuUSum[2]-0.7071067811865475*(jacob_vel_inv0[1]*v0[1]+jacob_vel_inv0[0]*v0[0])*nuSum[2])*rdvx2; 
  alphaDrag[5] = (nuSum[0]*(((-0.4517539514526256*jacob_vel_inv0[2])-0.7071067811865475*jacob_vel_inv0[0])*v0[2]-0.6210590034081186*jacob_vel_inv0[1]*v0[3])+((-0.7071067811865475*nuSum[0]*v0[0])-1.0*nuUSum[0])*jacob_vel_inv0[2]-0.6324555320336759*nuSum[0]*jacob_vel_inv0[1]*v0[1])*rdvx2; 
  alphaDrag[6] = ((-0.6210590034081186*jacob_vel_inv0[2]*nuSum[2]*v0[3])+jacob_vel_inv0[1]*((-0.6324555320336759*nuSum[2]*v0[2])-1.0*nuUSum[2])+((-0.6324555320336759*v0[1]*jacob_vel_inv0[2])-0.7071067811865475*(jacob_vel_inv0[0]*v0[1]+v0[0]*jacob_vel_inv0[1]))*nuSum[2])*rdvx2; 
  alphaDrag[7] = (nuSum[1]*(((-0.4517539514526256*jacob_vel_inv0[2])-0.7071067811865475*jacob_vel_inv0[0])*v0[2]-0.6210590034081186*jacob_vel_inv0[1]*v0[3])-1.0*nuUSum[1]*jacob_vel_inv0[2]+nuSum[1]*((-0.7071067811865475*v0[0]*jacob_vel_inv0[2])-0.6324555320336759*jacob_vel_inv0[1]*v0[1]))*rdvx2; 
  alphaDrag[8] = (nuSum[2]*(((-0.4517539514526256*jacob_vel_inv0[2])-0.7071067811865475*jacob_vel_inv0[0])*v0[2]-0.6210590034081186*jacob_vel_inv0[1]*v0[3])-1.0*jacob_vel_inv0[2]*nuUSum[2]+((-0.7071067811865475*v0[0]*jacob_vel_inv0[2])-0.6324555320336759*jacob_vel_inv0[1]*v0[1])*nuSum[2])*rdvx2; 

  out[2] += 0.8660254037844386*(alphaDrag[8]*f[8]+alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.7745966692414833*(alphaDrag[7]*f[8]+f[7]*alphaDrag[8])+0.8660254037844386*(alphaDrag[5]*f[7]+f[5]*alphaDrag[7])+0.7745966692414833*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4])+0.8660254037844386*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 1.732050807568877*(alphaDrag[6]*f[8]+f[6]*alphaDrag[8]+alphaDrag[3]*f[7]+f[3]*alphaDrag[7])+1.936491673103709*(alphaDrag[4]*f[6]+f[4]*alphaDrag[6])+1.732050807568877*(alphaDrag[2]*f[5]+f[2]*alphaDrag[5])+1.936491673103709*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[6] += 0.5532833351724881*alphaDrag[8]*f[8]+0.8660254037844386*(alphaDrag[5]*f[8]+f[5]*alphaDrag[8])+0.7745966692414833*alphaDrag[7]*f[7]+0.5532833351724881*alphaDrag[6]*f[6]+0.8660254037844386*(alphaDrag[2]*f[6]+f[2]*alphaDrag[6])+0.5532833351724881*alphaDrag[4]*f[4]+0.8660254037844386*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4])+0.7745966692414833*(alphaDrag[3]*f[3]+alphaDrag[1]*f[1]); 
  out[7] += 1.549193338482967*(alphaDrag[3]*f[8]+f[3]*alphaDrag[8])+(1.549193338482967*alphaDrag[6]+1.732050807568877*alphaDrag[2])*f[7]+1.549193338482967*f[6]*alphaDrag[7]+1.732050807568877*(f[2]*alphaDrag[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[3]*f[5]+f[3]*alphaDrag[5]+alphaDrag[3]*f[4]+f[3]*alphaDrag[4])+1.936491673103709*(alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[8] += (1.106566670344976*alphaDrag[6]+1.732050807568877*alphaDrag[2])*f[8]+(1.106566670344976*f[6]+1.732050807568877*f[2])*alphaDrag[8]+1.549193338482967*(alphaDrag[3]*f[7]+f[3]*alphaDrag[7])+(1.732050807568877*alphaDrag[5]+1.237179148263484*alphaDrag[4]+1.936491673103709*alphaDrag[0])*f[6]+(1.732050807568877*f[5]+1.237179148263484*f[4])*alphaDrag[6]+1.936491673103709*(f[0]*alphaDrag[6]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4])+1.732050807568877*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 

  return fabs(1.5625*alphaDrag[8]-1.397542485937369*(alphaDrag[5]+alphaDrag[4])+1.25*alphaDrag[0]); 

} 
