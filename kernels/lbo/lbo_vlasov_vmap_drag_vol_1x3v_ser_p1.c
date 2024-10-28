#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_vmap_drag_vol_1x3v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *jacob_vel_inv, const double *nuSum, const double *nuPrimMomsSum, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]: Cell-center coordinates. 
  // dxv[4]: Cell spacing. 
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
  const double rdvz2 = 2.0/dxv[3]; 
  const double *v2 = &vmap[8]; 
  const double *jacob_vel_inv2 = &jacob_vel_inv[6]; 

  double alphaDrag[120]; 
  alphaDrag[0] = (jacob_vel_inv0[0]*((-1.414213562373095*nuSum[0]*v0[0])-2.0*nuUSum[0])-1.414213562373095*nuSum[0]*(jacob_vel_inv0[2]*v0[2]+jacob_vel_inv0[1]*v0[1]))*rdvx2; 
  alphaDrag[1] = (jacob_vel_inv0[0]*((-2.0*nuUSum[1])-1.414213562373095*v0[0]*nuSum[1])-1.414213562373095*nuSum[1]*(jacob_vel_inv0[2]*v0[2]+jacob_vel_inv0[1]*v0[1]))*rdvx2; 
  alphaDrag[2] = (nuSum[0]*((-1.242118006816237*jacob_vel_inv0[2]*v0[3])-1.264911064067352*jacob_vel_inv0[1]*v0[2]+v0[1]*((-1.264911064067352*jacob_vel_inv0[2])-1.414213562373095*jacob_vel_inv0[0]))+((-1.414213562373095*nuSum[0]*v0[0])-2.0*nuUSum[0])*jacob_vel_inv0[1])*rdvx2; 
  alphaDrag[5] = (nuSum[1]*((-1.242118006816237*jacob_vel_inv0[2]*v0[3])-1.264911064067352*jacob_vel_inv0[1]*v0[2]+v0[1]*((-1.264911064067352*jacob_vel_inv0[2])-1.414213562373095*jacob_vel_inv0[0]))+jacob_vel_inv0[1]*((-2.0*nuUSum[1])-1.414213562373095*v0[0]*nuSum[1]))*rdvx2; 
  alphaDrag[16] = (nuSum[0]*(((-0.9035079029052515*jacob_vel_inv0[2])-1.414213562373095*jacob_vel_inv0[0])*v0[2]-1.242118006816237*jacob_vel_inv0[1]*v0[3])+((-1.414213562373095*nuSum[0]*v0[0])-2.0*nuUSum[0])*jacob_vel_inv0[2]-1.264911064067352*nuSum[0]*jacob_vel_inv0[1]*v0[1])*rdvx2; 
  alphaDrag[17] = (nuSum[1]*(((-0.9035079029052515*jacob_vel_inv0[2])-1.414213562373095*jacob_vel_inv0[0])*v0[2]-1.242118006816237*jacob_vel_inv0[1]*v0[3])-2.0*nuUSum[1]*jacob_vel_inv0[2]+nuSum[1]*((-1.414213562373095*v0[0]*jacob_vel_inv0[2])-1.264911064067352*jacob_vel_inv0[1]*v0[1]))*rdvx2; 

  alphaDrag[40] = ((-1.414213562373095*nuSum[0]*jacob_vel_inv1[2]*v1[2])-2.0*jacob_vel_inv1[0]*nuUSum[2]-1.414213562373095*nuSum[0]*(jacob_vel_inv1[1]*v1[1]+jacob_vel_inv1[0]*v1[0]))*rdvy2; 
  alphaDrag[41] = ((-2.0*jacob_vel_inv1[0]*nuUSum[3])-1.414213562373095*nuSum[1]*(jacob_vel_inv1[2]*v1[2]+jacob_vel_inv1[1]*v1[1]+jacob_vel_inv1[0]*v1[0]))*rdvy2; 
  alphaDrag[43] = ((-1.242118006816237*nuSum[0]*jacob_vel_inv1[2]*v1[3])+jacob_vel_inv1[1]*((-1.264911064067352*nuSum[0]*v1[2])-2.0*nuUSum[2])+nuSum[0]*((-1.264911064067352*v1[1]*jacob_vel_inv1[2])-1.414213562373095*(jacob_vel_inv1[0]*v1[1]+v1[0]*jacob_vel_inv1[1])))*rdvy2; 
  alphaDrag[46] = ((-1.242118006816237*nuSum[1]*jacob_vel_inv1[2]*v1[3])-2.0*jacob_vel_inv1[1]*nuUSum[3]+nuSum[1]*((-1.264911064067352*(jacob_vel_inv1[1]*v1[2]+v1[1]*jacob_vel_inv1[2]))-1.414213562373095*(jacob_vel_inv1[0]*v1[1]+v1[0]*jacob_vel_inv1[1])))*rdvy2; 
  alphaDrag[64] = (nuSum[0]*(((-0.9035079029052515*jacob_vel_inv1[2])-1.414213562373095*jacob_vel_inv1[0])*v1[2]-1.242118006816237*jacob_vel_inv1[1]*v1[3])-2.0*jacob_vel_inv1[2]*nuUSum[2]+nuSum[0]*((-1.414213562373095*v1[0]*jacob_vel_inv1[2])-1.264911064067352*jacob_vel_inv1[1]*v1[1]))*rdvy2; 
  alphaDrag[65] = ((-1.242118006816237*jacob_vel_inv1[1]*nuSum[1]*v1[3])-2.0*jacob_vel_inv1[2]*nuUSum[3]+nuSum[1]*((-1.414213562373095*(jacob_vel_inv1[0]*v1[2]+v1[0]*jacob_vel_inv1[2]))-0.9035079029052515*jacob_vel_inv1[2]*v1[2]-1.264911064067352*jacob_vel_inv1[1]*v1[1]))*rdvy2; 

  alphaDrag[80] = ((-2.0*jacob_vel_inv2[0]*nuUSum[4])-1.414213562373095*nuSum[0]*(jacob_vel_inv2[2]*v2[2]+jacob_vel_inv2[1]*v2[1]+jacob_vel_inv2[0]*v2[0]))*rdvz2; 
  alphaDrag[81] = ((-2.0*jacob_vel_inv2[0]*nuUSum[5])-1.414213562373095*nuSum[1]*(jacob_vel_inv2[2]*v2[2]+jacob_vel_inv2[1]*v2[1]+jacob_vel_inv2[0]*v2[0]))*rdvz2; 
  alphaDrag[84] = (nuSum[0]*((-1.242118006816237*jacob_vel_inv2[2]*v2[3])-1.264911064067352*(jacob_vel_inv2[1]*v2[2]+v2[1]*jacob_vel_inv2[2])-1.414213562373095*(jacob_vel_inv2[0]*v2[1]+v2[0]*jacob_vel_inv2[1]))-2.0*jacob_vel_inv2[1]*nuUSum[4])*rdvz2; 
  alphaDrag[88] = (nuSum[1]*((-1.242118006816237*jacob_vel_inv2[2]*v2[3])-1.264911064067352*(jacob_vel_inv2[1]*v2[2]+v2[1]*jacob_vel_inv2[2])-1.414213562373095*(jacob_vel_inv2[0]*v2[1]+v2[0]*jacob_vel_inv2[1]))-2.0*jacob_vel_inv2[1]*nuUSum[5])*rdvz2; 
  alphaDrag[112] = (nuSum[0]*((-1.242118006816237*jacob_vel_inv2[1]*v2[3])-1.414213562373095*(jacob_vel_inv2[0]*v2[2]+v2[0]*jacob_vel_inv2[2])-0.9035079029052515*jacob_vel_inv2[2]*v2[2]-1.264911064067352*jacob_vel_inv2[1]*v2[1])-2.0*jacob_vel_inv2[2]*nuUSum[4])*rdvz2; 
  alphaDrag[113] = (nuSum[1]*((-1.242118006816237*jacob_vel_inv2[1]*v2[3])-1.414213562373095*(jacob_vel_inv2[0]*v2[2]+v2[0]*jacob_vel_inv2[2])-0.9035079029052515*jacob_vel_inv2[2]*v2[2]-1.264911064067352*jacob_vel_inv2[1]*v2[1])-2.0*jacob_vel_inv2[2]*nuUSum[5])*rdvz2; 

  out[2] += 0.4330127018922193*(alphaDrag[17]*f[17]+alphaDrag[16]*f[16]+alphaDrag[5]*f[5]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[25]*alphaDrag[65]+f[24]*alphaDrag[64]+f[6]*alphaDrag[46]+f[3]*alphaDrag[43]+f[1]*alphaDrag[41]+f[0]*alphaDrag[40]); 
  out[4] += 0.4330127018922193*(f[33]*alphaDrag[113]+f[32]*alphaDrag[112]+f[8]*alphaDrag[88]+f[4]*alphaDrag[84]+f[1]*alphaDrag[81]+f[0]*alphaDrag[80]); 
  out[5] += 0.4330127018922193*(alphaDrag[16]*f[17]+f[16]*alphaDrag[17]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[6] += 0.4330127018922193*(f[24]*alphaDrag[65]+f[25]*alphaDrag[64]+f[3]*alphaDrag[46]+f[6]*alphaDrag[43]+f[0]*alphaDrag[41]+f[1]*alphaDrag[40]); 
  out[7] += 0.4330127018922193*(f[28]*alphaDrag[65]+f[26]*alphaDrag[64]+f[11]*alphaDrag[46]+f[7]*alphaDrag[43]+f[5]*alphaDrag[41]+f[2]*alphaDrag[40]+alphaDrag[17]*f[20]+alphaDrag[16]*f[18]+alphaDrag[5]*f[11]+alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+alphaDrag[0]*f[3]); 
  out[8] += 0.4330127018922193*(f[32]*alphaDrag[113]+f[33]*alphaDrag[112]+f[4]*alphaDrag[88]+f[8]*alphaDrag[84]+f[0]*alphaDrag[81]+f[1]*alphaDrag[80]); 
  out[9] += 0.4330127018922193*(f[36]*alphaDrag[113]+f[34]*alphaDrag[112]+f[12]*alphaDrag[88]+f[9]*alphaDrag[84]+f[5]*alphaDrag[81]+f[2]*alphaDrag[80]+alphaDrag[17]*f[21]+alphaDrag[16]*f[19]+alphaDrag[5]*f[12]+alphaDrag[2]*f[9]+alphaDrag[1]*f[8]+alphaDrag[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[37]*alphaDrag[113]+f[35]*alphaDrag[112]+f[13]*alphaDrag[88]+f[10]*alphaDrag[84]+f[6]*alphaDrag[81]+f[3]*alphaDrag[80]+f[29]*alphaDrag[65]+f[27]*alphaDrag[64]+f[13]*alphaDrag[46]+f[10]*alphaDrag[43]+f[8]*alphaDrag[41]+f[4]*alphaDrag[40]); 
  out[11] += 0.4330127018922193*(f[26]*alphaDrag[65]+f[28]*alphaDrag[64]+f[7]*alphaDrag[46]+f[11]*alphaDrag[43]+f[2]*alphaDrag[41]+f[5]*alphaDrag[40]+alphaDrag[16]*f[20]+alphaDrag[17]*f[18]+alphaDrag[2]*f[11]+alphaDrag[5]*f[7]+alphaDrag[0]*f[6]+alphaDrag[1]*f[3]); 
  out[12] += 0.4330127018922193*(f[34]*alphaDrag[113]+f[36]*alphaDrag[112]+f[9]*alphaDrag[88]+f[12]*alphaDrag[84]+f[2]*alphaDrag[81]+f[5]*alphaDrag[80]+alphaDrag[16]*f[21]+alphaDrag[17]*f[19]+alphaDrag[2]*f[12]+alphaDrag[5]*f[9]+alphaDrag[0]*f[8]+alphaDrag[1]*f[4]); 
  out[13] += 0.4330127018922193*(f[35]*alphaDrag[113]+f[37]*alphaDrag[112]+f[10]*alphaDrag[88]+f[13]*alphaDrag[84]+f[3]*alphaDrag[81]+f[6]*alphaDrag[80]+f[27]*alphaDrag[65]+f[29]*alphaDrag[64]+f[10]*alphaDrag[46]+f[13]*alphaDrag[43]+f[4]*alphaDrag[41]+f[8]*alphaDrag[40]); 
  out[14] += 0.4330127018922193*(f[39]*alphaDrag[113]+f[38]*alphaDrag[112]+f[15]*alphaDrag[88]+f[14]*alphaDrag[84]+f[11]*alphaDrag[81]+f[7]*alphaDrag[80]+f[31]*alphaDrag[65]+f[30]*alphaDrag[64]+f[15]*alphaDrag[46]+f[14]*alphaDrag[43]+f[12]*alphaDrag[41]+f[9]*alphaDrag[40]+alphaDrag[17]*f[23]+alphaDrag[16]*f[22]+alphaDrag[5]*f[15]+alphaDrag[2]*f[14]+alphaDrag[1]*f[13]+alphaDrag[0]*f[10]); 
  out[15] += 0.4330127018922193*(f[38]*alphaDrag[113]+f[39]*alphaDrag[112]+f[14]*alphaDrag[88]+f[15]*alphaDrag[84]+f[7]*alphaDrag[81]+f[11]*alphaDrag[80]+f[30]*alphaDrag[65]+f[31]*alphaDrag[64]+f[14]*alphaDrag[46]+f[15]*alphaDrag[43]+f[9]*alphaDrag[41]+f[12]*alphaDrag[40]+alphaDrag[16]*f[23]+alphaDrag[17]*f[22]+alphaDrag[2]*f[15]+alphaDrag[5]*f[14]+alphaDrag[0]*f[13]+alphaDrag[1]*f[10]); 
  out[16] += 0.8660254037844386*(alphaDrag[5]*f[17]+f[5]*alphaDrag[17]+alphaDrag[2]*f[16]+f[2]*alphaDrag[16])+0.9682458365518543*(alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[17] += 0.8660254037844386*(alphaDrag[2]*f[17]+f[2]*alphaDrag[17]+alphaDrag[5]*f[16]+f[5]*alphaDrag[16])+0.9682458365518543*(alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[18] += 0.4330127018922193*(f[20]*alphaDrag[46]+f[18]*alphaDrag[43]+f[17]*alphaDrag[41]+f[16]*alphaDrag[40])+0.8660254037844386*(alphaDrag[5]*f[20]+alphaDrag[2]*f[18]+f[11]*alphaDrag[17]+f[7]*alphaDrag[16])+0.9682458365518543*(alphaDrag[1]*f[11]+alphaDrag[0]*f[7]+alphaDrag[5]*f[6]+alphaDrag[2]*f[3]); 
  out[19] += 0.4330127018922193*(f[21]*alphaDrag[88]+f[19]*alphaDrag[84]+f[17]*alphaDrag[81]+f[16]*alphaDrag[80])+0.8660254037844386*(alphaDrag[5]*f[21]+alphaDrag[2]*f[19]+f[12]*alphaDrag[17]+f[9]*alphaDrag[16])+0.9682458365518543*(alphaDrag[1]*f[12]+alphaDrag[0]*f[9]+alphaDrag[5]*f[8]+alphaDrag[2]*f[4]); 
  out[20] += 0.4330127018922193*(f[18]*alphaDrag[46]+f[20]*alphaDrag[43]+f[16]*alphaDrag[41]+f[17]*alphaDrag[40])+0.8660254037844386*(alphaDrag[2]*f[20]+alphaDrag[5]*f[18]+f[7]*alphaDrag[17]+f[11]*alphaDrag[16])+0.9682458365518543*(alphaDrag[0]*f[11]+alphaDrag[1]*f[7]+alphaDrag[2]*f[6]+f[3]*alphaDrag[5]); 
  out[21] += 0.4330127018922193*(f[19]*alphaDrag[88]+f[21]*alphaDrag[84]+f[16]*alphaDrag[81]+f[17]*alphaDrag[80])+0.8660254037844386*(alphaDrag[2]*f[21]+alphaDrag[5]*f[19]+f[9]*alphaDrag[17]+f[12]*alphaDrag[16])+0.9682458365518543*(alphaDrag[0]*f[12]+alphaDrag[1]*f[9]+alphaDrag[2]*f[8]+f[4]*alphaDrag[5]); 
  out[22] += 0.4330127018922193*(f[23]*alphaDrag[88]+f[22]*alphaDrag[84]+f[20]*alphaDrag[81]+f[18]*alphaDrag[80]+f[23]*alphaDrag[46]+f[22]*alphaDrag[43]+f[21]*alphaDrag[41]+f[19]*alphaDrag[40])+0.8660254037844386*(alphaDrag[5]*f[23]+alphaDrag[2]*f[22]+f[15]*alphaDrag[17]+f[14]*alphaDrag[16])+0.9682458365518543*(alphaDrag[1]*f[15]+alphaDrag[0]*f[14]+alphaDrag[5]*f[13]+alphaDrag[2]*f[10]); 
  out[23] += 0.4330127018922193*(f[22]*alphaDrag[88]+f[23]*alphaDrag[84]+f[18]*alphaDrag[81]+f[20]*alphaDrag[80]+f[22]*alphaDrag[46]+f[23]*alphaDrag[43]+f[19]*alphaDrag[41]+f[21]*alphaDrag[40])+0.8660254037844386*(alphaDrag[2]*f[23]+alphaDrag[5]*f[22]+f[14]*alphaDrag[17]+f[15]*alphaDrag[16])+0.9682458365518543*(alphaDrag[0]*f[15]+alphaDrag[1]*f[14]+alphaDrag[2]*f[13]+alphaDrag[5]*f[10]); 
  out[24] += 0.8660254037844386*(f[6]*alphaDrag[65]+f[3]*alphaDrag[64])+(0.8660254037844386*f[25]+0.9682458365518543*f[1])*alphaDrag[46]+0.8660254037844386*f[24]*alphaDrag[43]+0.9682458365518543*(f[0]*alphaDrag[43]+f[6]*alphaDrag[41]+f[3]*alphaDrag[40]); 
  out[25] += 0.8660254037844386*(f[3]*alphaDrag[65]+f[6]*alphaDrag[64])+(0.8660254037844386*f[24]+0.9682458365518543*f[0])*alphaDrag[46]+0.8660254037844386*f[25]*alphaDrag[43]+0.9682458365518543*(f[1]*alphaDrag[43]+f[3]*alphaDrag[41]+f[6]*alphaDrag[40]); 
  out[26] += 0.8660254037844386*(f[11]*alphaDrag[65]+f[7]*alphaDrag[64])+(0.8660254037844386*f[28]+0.9682458365518543*f[5])*alphaDrag[46]+0.8660254037844386*f[26]*alphaDrag[43]+0.9682458365518543*(f[2]*alphaDrag[43]+f[11]*alphaDrag[41]+f[7]*alphaDrag[40])+0.4330127018922193*(alphaDrag[5]*f[28]+alphaDrag[2]*f[26]+alphaDrag[1]*f[25]+alphaDrag[0]*f[24]); 
  out[27] += 0.4330127018922193*(f[29]*alphaDrag[88]+f[27]*alphaDrag[84]+f[25]*alphaDrag[81]+f[24]*alphaDrag[80])+0.8660254037844386*(f[13]*alphaDrag[65]+f[10]*alphaDrag[64])+(0.8660254037844386*f[29]+0.9682458365518543*f[8])*alphaDrag[46]+0.8660254037844386*f[27]*alphaDrag[43]+0.9682458365518543*(f[4]*alphaDrag[43]+f[13]*alphaDrag[41]+f[10]*alphaDrag[40]); 
  out[28] += 0.8660254037844386*(f[7]*alphaDrag[65]+f[11]*alphaDrag[64])+(0.8660254037844386*f[26]+0.9682458365518543*f[2])*alphaDrag[46]+0.8660254037844386*f[28]*alphaDrag[43]+0.9682458365518543*(f[5]*alphaDrag[43]+f[7]*alphaDrag[41]+f[11]*alphaDrag[40])+0.4330127018922193*(alphaDrag[2]*f[28]+alphaDrag[5]*f[26]+alphaDrag[0]*f[25]+alphaDrag[1]*f[24]); 
  out[29] += 0.4330127018922193*(f[27]*alphaDrag[88]+f[29]*alphaDrag[84]+f[24]*alphaDrag[81]+f[25]*alphaDrag[80])+0.8660254037844386*(f[10]*alphaDrag[65]+f[13]*alphaDrag[64])+(0.8660254037844386*f[27]+0.9682458365518543*f[4])*alphaDrag[46]+0.8660254037844386*f[29]*alphaDrag[43]+0.9682458365518543*(f[8]*alphaDrag[43]+f[10]*alphaDrag[41]+f[13]*alphaDrag[40]); 
  out[30] += 0.4330127018922193*(f[31]*alphaDrag[88]+f[30]*alphaDrag[84]+f[28]*alphaDrag[81]+f[26]*alphaDrag[80])+0.8660254037844386*(f[15]*alphaDrag[65]+f[14]*alphaDrag[64])+(0.8660254037844386*f[31]+0.9682458365518543*f[12])*alphaDrag[46]+0.8660254037844386*f[30]*alphaDrag[43]+0.9682458365518543*(f[9]*alphaDrag[43]+f[15]*alphaDrag[41]+f[14]*alphaDrag[40])+0.4330127018922193*(alphaDrag[5]*f[31]+alphaDrag[2]*f[30]+alphaDrag[1]*f[29]+alphaDrag[0]*f[27]); 
  out[31] += 0.4330127018922193*(f[30]*alphaDrag[88]+f[31]*alphaDrag[84]+f[26]*alphaDrag[81]+f[28]*alphaDrag[80])+0.8660254037844386*(f[14]*alphaDrag[65]+f[15]*alphaDrag[64])+(0.8660254037844386*f[30]+0.9682458365518543*f[9])*alphaDrag[46]+0.8660254037844386*f[31]*alphaDrag[43]+0.9682458365518543*(f[12]*alphaDrag[43]+f[14]*alphaDrag[41]+f[15]*alphaDrag[40])+0.4330127018922193*(alphaDrag[2]*f[31]+alphaDrag[5]*f[30]+alphaDrag[0]*f[29]+alphaDrag[1]*f[27]); 
  out[32] += 0.8660254037844386*(f[8]*alphaDrag[113]+f[4]*alphaDrag[112])+(0.8660254037844386*f[33]+0.9682458365518543*f[1])*alphaDrag[88]+0.8660254037844386*f[32]*alphaDrag[84]+0.9682458365518543*(f[0]*alphaDrag[84]+f[8]*alphaDrag[81]+f[4]*alphaDrag[80]); 
  out[33] += 0.8660254037844386*(f[4]*alphaDrag[113]+f[8]*alphaDrag[112])+(0.8660254037844386*f[32]+0.9682458365518543*f[0])*alphaDrag[88]+0.8660254037844386*f[33]*alphaDrag[84]+0.9682458365518543*(f[1]*alphaDrag[84]+f[4]*alphaDrag[81]+f[8]*alphaDrag[80]); 
  out[34] += 0.8660254037844386*(f[12]*alphaDrag[113]+f[9]*alphaDrag[112])+(0.8660254037844386*f[36]+0.9682458365518543*f[5])*alphaDrag[88]+0.8660254037844386*f[34]*alphaDrag[84]+0.9682458365518543*(f[2]*alphaDrag[84]+f[12]*alphaDrag[81]+f[9]*alphaDrag[80])+0.4330127018922193*(alphaDrag[5]*f[36]+alphaDrag[2]*f[34]+alphaDrag[1]*f[33]+alphaDrag[0]*f[32]); 
  out[35] += 0.8660254037844386*(f[13]*alphaDrag[113]+f[10]*alphaDrag[112])+(0.8660254037844386*f[37]+0.9682458365518543*f[6])*alphaDrag[88]+0.8660254037844386*f[35]*alphaDrag[84]+0.9682458365518543*(f[3]*alphaDrag[84]+f[13]*alphaDrag[81]+f[10]*alphaDrag[80])+0.4330127018922193*(f[37]*alphaDrag[46]+f[35]*alphaDrag[43]+f[33]*alphaDrag[41]+f[32]*alphaDrag[40]); 
  out[36] += 0.8660254037844386*(f[9]*alphaDrag[113]+f[12]*alphaDrag[112])+(0.8660254037844386*f[34]+0.9682458365518543*f[2])*alphaDrag[88]+0.8660254037844386*f[36]*alphaDrag[84]+0.9682458365518543*(f[5]*alphaDrag[84]+f[9]*alphaDrag[81]+f[12]*alphaDrag[80])+0.4330127018922193*(alphaDrag[2]*f[36]+alphaDrag[5]*f[34]+alphaDrag[0]*f[33]+alphaDrag[1]*f[32]); 
  out[37] += 0.8660254037844386*(f[10]*alphaDrag[113]+f[13]*alphaDrag[112])+(0.8660254037844386*f[35]+0.9682458365518543*f[3])*alphaDrag[88]+0.8660254037844386*f[37]*alphaDrag[84]+0.9682458365518543*(f[6]*alphaDrag[84]+f[10]*alphaDrag[81]+f[13]*alphaDrag[80])+0.4330127018922193*(f[35]*alphaDrag[46]+f[37]*alphaDrag[43]+f[32]*alphaDrag[41]+f[33]*alphaDrag[40]); 
  out[38] += 0.8660254037844386*(f[15]*alphaDrag[113]+f[14]*alphaDrag[112])+(0.8660254037844386*f[39]+0.9682458365518543*f[11])*alphaDrag[88]+0.8660254037844386*f[38]*alphaDrag[84]+0.9682458365518543*(f[7]*alphaDrag[84]+f[15]*alphaDrag[81]+f[14]*alphaDrag[80])+0.4330127018922193*(f[39]*alphaDrag[46]+f[38]*alphaDrag[43]+f[36]*alphaDrag[41]+f[34]*alphaDrag[40]+alphaDrag[5]*f[39]+alphaDrag[2]*f[38]+alphaDrag[1]*f[37]+alphaDrag[0]*f[35]); 
  out[39] += 0.8660254037844386*(f[14]*alphaDrag[113]+f[15]*alphaDrag[112])+(0.8660254037844386*f[38]+0.9682458365518543*f[7])*alphaDrag[88]+0.8660254037844386*f[39]*alphaDrag[84]+0.9682458365518543*(f[11]*alphaDrag[84]+f[14]*alphaDrag[81]+f[15]*alphaDrag[80])+0.4330127018922193*(f[38]*alphaDrag[46]+f[39]*alphaDrag[43]+f[34]*alphaDrag[41]+f[36]*alphaDrag[40]+alphaDrag[2]*f[39]+alphaDrag[5]*f[38]+alphaDrag[0]*f[37]+alphaDrag[1]*f[35]); 

  return fabs(0.625*alphaDrag[0]-0.6987712429686843*alphaDrag[16])+fabs(0.625*alphaDrag[40]-0.6987712429686843*alphaDrag[64])+fabs(0.625*alphaDrag[80]-0.6987712429686843*alphaDrag[112]); 

} 
