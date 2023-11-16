#include <gkyl_rad_gyrokinetic_kernels.h> 
GKYL_CU_DH double rad_gyrokinetic_drag_vol_2x2v_ser_p1(const double *w, const double *dxv, const double *nI, const double *vnu, const double *vsqnu, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]: cell-center coordinates. 
  // dxv[4]: cell spacing. 
  // nI: atomic density 
  // vnu: 2/pi*v*nu(v) dg field representation (v'(v||,mu) in notes) 
  // vsqnu: sqrt(mu*me/2B)*v^2*nu(v) dg field representation (v''(v||,mu) in notes) 
  // f: input distribution function.
  // out: incremented output 

  double rdv2[2]; 
  rdv2[0] = 2.0/dxv[2]; 
  rdv2[1] = 2.0/dxv[3]; 

  double alphaDrag[24] = {0,0}; 
  alphaDrag[0] = 0.5*nI[3]*vnu[5]+0.5*nI[2]*vnu[2]+0.5*nI[1]*vnu[1]+0.5*nI[0]*vnu[0]; 
  alphaDrag[1] = 0.5*nI[2]*vnu[5]+0.5*vnu[2]*nI[3]+0.5*nI[0]*vnu[1]+0.5*vnu[0]*nI[1]; 
  alphaDrag[2] = 0.5*nI[1]*vnu[5]+0.5*vnu[1]*nI[3]+0.5*nI[0]*vnu[2]+0.5*vnu[0]*nI[2]; 
  alphaDrag[3] = 0.5*nI[3]*vnu[11]+0.5*nI[2]*vnu[7]+0.5*nI[1]*vnu[6]+0.5*nI[0]*vnu[3]; 
  alphaDrag[4] = 0.5*nI[3]*vnu[12]+0.5*nI[2]*vnu[9]+0.5*nI[1]*vnu[8]+0.5*nI[0]*vnu[4]; 
  alphaDrag[5] = 0.5*nI[0]*vnu[5]+0.5*vnu[0]*nI[3]+0.5*nI[1]*vnu[2]+0.5*vnu[1]*nI[2]; 
  alphaDrag[6] = 0.5*nI[2]*vnu[11]+0.5*nI[3]*vnu[7]+0.5*nI[0]*vnu[6]+0.5*nI[1]*vnu[3]; 
  alphaDrag[7] = 0.5*nI[1]*vnu[11]+0.5*nI[0]*vnu[7]+0.5*nI[3]*vnu[6]+0.5*nI[2]*vnu[3]; 
  alphaDrag[8] = 0.5*nI[2]*vnu[12]+0.5*nI[3]*vnu[9]+0.5*nI[0]*vnu[8]+0.5*nI[1]*vnu[4]; 
  alphaDrag[9] = 0.5*nI[1]*vnu[12]+0.5*nI[0]*vnu[9]+0.5*nI[3]*vnu[8]+0.5*nI[2]*vnu[4]; 
  alphaDrag[10] = 0.5*nI[3]*vnu[15]+0.5*nI[2]*vnu[14]+0.5*nI[1]*vnu[13]+0.5*nI[0]*vnu[10]; 
  alphaDrag[11] = 0.5*nI[0]*vnu[11]+0.5*nI[1]*vnu[7]+0.5*nI[2]*vnu[6]+0.5*nI[3]*vnu[3]; 
  alphaDrag[12] = 0.5*nI[0]*vnu[12]+0.5*nI[1]*vnu[9]+0.5*nI[2]*vnu[8]+0.5*nI[3]*vnu[4]; 
  alphaDrag[13] = 0.5*nI[2]*vnu[15]+0.5*nI[3]*vnu[14]+0.5*nI[0]*vnu[13]+0.5*nI[1]*vnu[10]; 
  alphaDrag[14] = 0.5*nI[1]*vnu[15]+0.5*nI[0]*vnu[14]+0.5*nI[3]*vnu[13]+0.5*nI[2]*vnu[10]; 
  alphaDrag[15] = 0.5*nI[0]*vnu[15]+0.5*nI[1]*vnu[14]+0.5*nI[2]*vnu[13]+0.5*nI[3]*vnu[10]; 
  alphaDrag[16] = 0.5*nI[3]*vnu[20]+0.5000000000000001*nI[2]*vnu[18]+0.5000000000000001*nI[1]*vnu[17]+0.5*nI[0]*vnu[16]; 
  alphaDrag[17] = 0.5000000000000001*nI[2]*vnu[20]+0.5*nI[3]*vnu[18]+0.5*nI[0]*vnu[17]+0.5000000000000001*nI[1]*vnu[16]; 
  alphaDrag[18] = 0.5000000000000001*nI[1]*vnu[20]+0.5*nI[0]*vnu[18]+0.5*nI[3]*vnu[17]+0.5000000000000001*nI[2]*vnu[16]; 
  alphaDrag[19] = 0.5*nI[3]*vnu[23]+0.5000000000000001*nI[2]*vnu[22]+0.5000000000000001*nI[1]*vnu[21]+0.5*nI[0]*vnu[19]; 
  alphaDrag[20] = 0.5*nI[0]*vnu[20]+0.5000000000000001*nI[1]*vnu[18]+0.5000000000000001*nI[2]*vnu[17]+0.5*nI[3]*vnu[16]; 
  alphaDrag[21] = 0.5000000000000001*nI[2]*vnu[23]+0.5*nI[3]*vnu[22]+0.5*nI[0]*vnu[21]+0.5000000000000001*nI[1]*vnu[19]; 
  alphaDrag[22] = 0.5000000000000001*nI[1]*vnu[23]+0.5*nI[0]*vnu[22]+0.5*nI[3]*vnu[21]+0.5000000000000001*nI[2]*vnu[19]; 
  alphaDrag[23] = 0.5*nI[0]*vnu[23]+0.5000000000000001*nI[1]*vnu[22]+0.5000000000000001*nI[2]*vnu[21]+0.5*nI[3]*vnu[19]; 

  alphaDrag[0] += 0.5*nI[3]*vsqnu[5]+0.5*nI[2]*vsqnu[2]+0.5*nI[1]*vsqnu[1]+0.5*nI[0]*vsqnu[0]; 
  alphaDrag[1] += 0.5*nI[2]*vsqnu[5]+0.5*vsqnu[2]*nI[3]+0.5*nI[0]*vsqnu[1]+0.5*vsqnu[0]*nI[1]; 
  alphaDrag[2] += 0.5*nI[1]*vsqnu[5]+0.5*vsqnu[1]*nI[3]+0.5*nI[0]*vsqnu[2]+0.5*vsqnu[0]*nI[2]; 
  alphaDrag[3] += 0.5*nI[3]*vsqnu[11]+0.5*nI[2]*vsqnu[7]+0.5*nI[1]*vsqnu[6]+0.5*nI[0]*vsqnu[3]; 
  alphaDrag[4] += 0.5*nI[3]*vsqnu[12]+0.5*nI[2]*vsqnu[9]+0.5*nI[1]*vsqnu[8]+0.5*nI[0]*vsqnu[4]; 
  alphaDrag[5] += 0.5*nI[0]*vsqnu[5]+0.5*vsqnu[0]*nI[3]+0.5*nI[1]*vsqnu[2]+0.5*vsqnu[1]*nI[2]; 
  alphaDrag[6] += 0.5*nI[2]*vsqnu[11]+0.5*nI[3]*vsqnu[7]+0.5*nI[0]*vsqnu[6]+0.5*nI[1]*vsqnu[3]; 
  alphaDrag[7] += 0.5*nI[1]*vsqnu[11]+0.5*nI[0]*vsqnu[7]+0.5*nI[3]*vsqnu[6]+0.5*nI[2]*vsqnu[3]; 
  alphaDrag[8] += 0.5*nI[2]*vsqnu[12]+0.5*nI[3]*vsqnu[9]+0.5*nI[0]*vsqnu[8]+0.5*nI[1]*vsqnu[4]; 
  alphaDrag[9] += 0.5*nI[1]*vsqnu[12]+0.5*nI[0]*vsqnu[9]+0.5*nI[3]*vsqnu[8]+0.5*nI[2]*vsqnu[4]; 
  alphaDrag[10] += 0.5*nI[3]*vsqnu[15]+0.5*nI[2]*vsqnu[14]+0.5*nI[1]*vsqnu[13]+0.5*nI[0]*vsqnu[10]; 
  alphaDrag[11] += 0.5*nI[0]*vsqnu[11]+0.5*nI[1]*vsqnu[7]+0.5*nI[2]*vsqnu[6]+0.5*nI[3]*vsqnu[3]; 
  alphaDrag[12] += 0.5*nI[0]*vsqnu[12]+0.5*nI[1]*vsqnu[9]+0.5*nI[2]*vsqnu[8]+0.5*nI[3]*vsqnu[4]; 
  alphaDrag[13] += 0.5*nI[2]*vsqnu[15]+0.5*nI[3]*vsqnu[14]+0.5*nI[0]*vsqnu[13]+0.5*nI[1]*vsqnu[10]; 
  alphaDrag[14] += 0.5*nI[1]*vsqnu[15]+0.5*nI[0]*vsqnu[14]+0.5*nI[3]*vsqnu[13]+0.5*nI[2]*vsqnu[10]; 
  alphaDrag[15] += 0.5*nI[0]*vsqnu[15]+0.5*nI[1]*vsqnu[14]+0.5*nI[2]*vsqnu[13]+0.5*nI[3]*vsqnu[10]; 
  alphaDrag[16] += 0.5*nI[3]*vsqnu[20]+0.5000000000000001*nI[2]*vsqnu[18]+0.5000000000000001*nI[1]*vsqnu[17]+0.5*nI[0]*vsqnu[16]; 
  alphaDrag[17] += 0.5000000000000001*nI[2]*vsqnu[20]+0.5*nI[3]*vsqnu[18]+0.5*nI[0]*vsqnu[17]+0.5000000000000001*nI[1]*vsqnu[16]; 
  alphaDrag[18] += 0.5000000000000001*nI[1]*vsqnu[20]+0.5*nI[0]*vsqnu[18]+0.5*nI[3]*vsqnu[17]+0.5000000000000001*nI[2]*vsqnu[16]; 
  alphaDrag[19] += 0.5*nI[3]*vsqnu[23]+0.5000000000000001*nI[2]*vsqnu[22]+0.5000000000000001*nI[1]*vsqnu[21]+0.5*nI[0]*vsqnu[19]; 
  alphaDrag[20] += 0.5*nI[0]*vsqnu[20]+0.5000000000000001*nI[1]*vsqnu[18]+0.5000000000000001*nI[2]*vsqnu[17]+0.5*nI[3]*vsqnu[16]; 
  alphaDrag[21] += 0.5000000000000001*nI[2]*vsqnu[23]+0.5*nI[3]*vsqnu[22]+0.5*nI[0]*vsqnu[21]+0.5000000000000001*nI[1]*vsqnu[19]; 
  alphaDrag[22] += 0.5000000000000001*nI[1]*vsqnu[23]+0.5*nI[0]*vsqnu[22]+0.5*nI[3]*vsqnu[21]+0.5000000000000001*nI[2]*vsqnu[19]; 
  alphaDrag[23] += 0.5*nI[0]*vsqnu[23]+0.5000000000000001*nI[1]*vsqnu[22]+0.5000000000000001*nI[2]*vsqnu[21]+0.5*nI[3]*vsqnu[19]; 

  out[3] += 0.8660254037844386*(alphaDrag[23]*f[23]+alphaDrag[22]*f[22]+alphaDrag[21]*f[21]+alphaDrag[20]*f[20]+alphaDrag[19]*f[19]+alphaDrag[18]*f[18]+alphaDrag[17]*f[17]+alphaDrag[16]*f[16]+alphaDrag[15]*f[15]+alphaDrag[14]*f[14]+alphaDrag[13]*f[13]+alphaDrag[12]*f[12]+alphaDrag[11]*f[11]+alphaDrag[10]*f[10]+alphaDrag[9]*f[9]+alphaDrag[8]*f[8]+alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[6] += 0.8660254037844386*(alphaDrag[22]*f[23]+f[22]*alphaDrag[23]+alphaDrag[19]*f[21]+f[19]*alphaDrag[21]+alphaDrag[18]*f[20]+f[18]*alphaDrag[20]+alphaDrag[16]*f[17]+f[16]*alphaDrag[17]+alphaDrag[14]*f[15]+f[14]*alphaDrag[15]+alphaDrag[10]*f[13]+f[10]*alphaDrag[13]+alphaDrag[9]*f[12]+f[9]*alphaDrag[12]+alphaDrag[7]*f[11]+f[7]*alphaDrag[11]+alphaDrag[4]*f[8]+f[4]*alphaDrag[8]+alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[7] += 0.8660254037844386*(alphaDrag[21]*f[23]+f[21]*alphaDrag[23]+alphaDrag[19]*f[22]+f[19]*alphaDrag[22]+alphaDrag[17]*f[20]+f[17]*alphaDrag[20]+alphaDrag[16]*f[18]+f[16]*alphaDrag[18]+alphaDrag[13]*f[15]+f[13]*alphaDrag[15]+alphaDrag[10]*f[14]+f[10]*alphaDrag[14]+alphaDrag[8]*f[12]+f[8]*alphaDrag[12]+alphaDrag[6]*f[11]+f[6]*alphaDrag[11]+alphaDrag[4]*f[9]+f[4]*alphaDrag[9]+alphaDrag[3]*f[7]+f[3]*alphaDrag[7]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[10] += 0.8660254037844386*(alphaDrag[20]*f[23]+f[20]*alphaDrag[23]+alphaDrag[18]*f[22]+f[18]*alphaDrag[22]+alphaDrag[17]*f[21]+f[17]*alphaDrag[21]+alphaDrag[16]*f[19]+f[16]*alphaDrag[19]+alphaDrag[11]*f[15]+f[11]*alphaDrag[15]+alphaDrag[7]*f[14]+f[7]*alphaDrag[14]+alphaDrag[6]*f[13]+f[6]*alphaDrag[13]+alphaDrag[5]*f[12]+f[5]*alphaDrag[12]+alphaDrag[3]*f[10]+f[3]*alphaDrag[10]+alphaDrag[2]*f[9]+f[2]*alphaDrag[9]+alphaDrag[1]*f[8]+f[1]*alphaDrag[8]+alphaDrag[0]*f[4]+f[0]*alphaDrag[4]); 
  out[11] += 0.8660254037844386*(alphaDrag[19]*f[23]+f[19]*alphaDrag[23]+alphaDrag[21]*f[22]+f[21]*alphaDrag[22]+alphaDrag[16]*f[20]+f[16]*alphaDrag[20]+alphaDrag[17]*f[18]+f[17]*alphaDrag[18]+alphaDrag[10]*f[15]+f[10]*alphaDrag[15]+alphaDrag[13]*f[14]+f[13]*alphaDrag[14]+alphaDrag[4]*f[12]+f[4]*alphaDrag[12]+alphaDrag[3]*f[11]+f[3]*alphaDrag[11]+alphaDrag[8]*f[9]+f[8]*alphaDrag[9]+alphaDrag[6]*f[7]+f[6]*alphaDrag[7]+alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[13] += 0.8660254037844386*(alphaDrag[18]*f[23]+f[18]*alphaDrag[23]+alphaDrag[20]*f[22]+f[20]*alphaDrag[22]+alphaDrag[16]*f[21]+f[16]*alphaDrag[21]+alphaDrag[17]*f[19]+f[17]*alphaDrag[19]+alphaDrag[7]*f[15]+f[7]*alphaDrag[15]+alphaDrag[11]*f[14]+f[11]*alphaDrag[14]+alphaDrag[3]*f[13]+f[3]*alphaDrag[13]+alphaDrag[2]*f[12]+f[2]*alphaDrag[12]+alphaDrag[6]*f[10]+f[6]*alphaDrag[10]+alphaDrag[5]*f[9]+f[5]*alphaDrag[9]+alphaDrag[0]*f[8]+f[0]*alphaDrag[8]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4]); 
  out[14] += 0.8660254037844386*(alphaDrag[17]*f[23]+f[17]*alphaDrag[23]+alphaDrag[16]*f[22]+f[16]*alphaDrag[22]+alphaDrag[20]*f[21]+f[20]*alphaDrag[21]+alphaDrag[18]*f[19]+f[18]*alphaDrag[19]+alphaDrag[6]*f[15]+f[6]*alphaDrag[15]+alphaDrag[3]*f[14]+f[3]*alphaDrag[14]+alphaDrag[11]*f[13]+f[11]*alphaDrag[13]+alphaDrag[1]*f[12]+f[1]*alphaDrag[12]+alphaDrag[7]*f[10]+f[7]*alphaDrag[10]+alphaDrag[0]*f[9]+f[0]*alphaDrag[9]+alphaDrag[5]*f[8]+f[5]*alphaDrag[8]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4]); 
  out[15] += 0.8660254037844386*(alphaDrag[16]*f[23]+f[16]*alphaDrag[23]+alphaDrag[17]*f[22]+f[17]*alphaDrag[22]+alphaDrag[18]*f[21]+f[18]*alphaDrag[21]+alphaDrag[19]*f[20]+f[19]*alphaDrag[20]+alphaDrag[3]*f[15]+f[3]*alphaDrag[15]+alphaDrag[6]*f[14]+f[6]*alphaDrag[14]+alphaDrag[7]*f[13]+f[7]*alphaDrag[13]+alphaDrag[0]*f[12]+f[0]*alphaDrag[12]+alphaDrag[10]*f[11]+f[10]*alphaDrag[11]+alphaDrag[1]*f[9]+f[1]*alphaDrag[9]+alphaDrag[2]*f[8]+f[2]*alphaDrag[8]+alphaDrag[4]*f[5]+f[4]*alphaDrag[5]); 
  out[16] += 1.732050807568877*(alphaDrag[15]*f[23]+f[15]*alphaDrag[23]+alphaDrag[14]*f[22]+f[14]*alphaDrag[22]+alphaDrag[13]*f[21]+f[13]*alphaDrag[21]+alphaDrag[11]*f[20]+f[11]*alphaDrag[20]+alphaDrag[10]*f[19]+f[10]*alphaDrag[19]+alphaDrag[7]*f[18]+f[7]*alphaDrag[18]+alphaDrag[6]*f[17]+f[6]*alphaDrag[17]+alphaDrag[3]*f[16]+f[3]*alphaDrag[16])+1.936491673103709*(alphaDrag[12]*f[15]+f[12]*alphaDrag[15]+alphaDrag[9]*f[14]+f[9]*alphaDrag[14]+alphaDrag[8]*f[13]+f[8]*alphaDrag[13]+alphaDrag[5]*f[11]+f[5]*alphaDrag[11]+alphaDrag[4]*f[10]+f[4]*alphaDrag[10]+alphaDrag[2]*f[7]+f[2]*alphaDrag[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[17] += 1.732050807568877*(alphaDrag[14]*f[23]+f[14]*alphaDrag[23]+alphaDrag[15]*f[22]+f[15]*alphaDrag[22]+alphaDrag[10]*f[21]+f[10]*alphaDrag[21]+alphaDrag[7]*f[20]+f[7]*alphaDrag[20]+alphaDrag[13]*f[19]+f[13]*alphaDrag[19]+alphaDrag[11]*f[18]+f[11]*alphaDrag[18]+alphaDrag[3]*f[17]+f[3]*alphaDrag[17]+alphaDrag[6]*f[16]+f[6]*alphaDrag[16])+1.936491673103709*(alphaDrag[9]*f[15]+f[9]*alphaDrag[15]+alphaDrag[12]*f[14]+f[12]*alphaDrag[14]+alphaDrag[4]*f[13]+f[4]*alphaDrag[13]+alphaDrag[2]*f[11]+f[2]*alphaDrag[11]+alphaDrag[8]*f[10]+f[8]*alphaDrag[10]+alphaDrag[5]*f[7]+f[5]*alphaDrag[7]+alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[18] += 1.732050807568877*(alphaDrag[13]*f[23]+f[13]*alphaDrag[23]+alphaDrag[10]*f[22]+f[10]*alphaDrag[22]+alphaDrag[15]*f[21]+f[15]*alphaDrag[21]+alphaDrag[6]*f[20]+f[6]*alphaDrag[20]+alphaDrag[14]*f[19]+f[14]*alphaDrag[19]+alphaDrag[3]*f[18]+f[3]*alphaDrag[18]+alphaDrag[11]*f[17]+f[11]*alphaDrag[17]+alphaDrag[7]*f[16]+f[7]*alphaDrag[16])+1.936491673103709*(alphaDrag[8]*f[15]+f[8]*alphaDrag[15]+alphaDrag[4]*f[14]+f[4]*alphaDrag[14]+alphaDrag[12]*f[13]+f[12]*alphaDrag[13]+alphaDrag[1]*f[11]+f[1]*alphaDrag[11]+alphaDrag[9]*f[10]+f[9]*alphaDrag[10]+alphaDrag[0]*f[7]+f[0]*alphaDrag[7]+alphaDrag[5]*f[6]+f[5]*alphaDrag[6]+alphaDrag[2]*f[3]+f[2]*alphaDrag[3]); 
  out[19] += 1.732050807568877*(alphaDrag[11]*f[23]+f[11]*alphaDrag[23]+alphaDrag[7]*f[22]+f[7]*alphaDrag[22]+alphaDrag[6]*f[21]+f[6]*alphaDrag[21]+alphaDrag[15]*f[20]+f[15]*alphaDrag[20]+alphaDrag[3]*f[19]+f[3]*alphaDrag[19]+alphaDrag[14]*f[18]+f[14]*alphaDrag[18]+alphaDrag[13]*f[17]+f[13]*alphaDrag[17]+alphaDrag[10]*f[16]+f[10]*alphaDrag[16])+1.936491673103709*(alphaDrag[5]*f[15]+f[5]*alphaDrag[15]+alphaDrag[2]*f[14]+f[2]*alphaDrag[14]+alphaDrag[1]*f[13]+f[1]*alphaDrag[13]+alphaDrag[11]*f[12]+f[11]*alphaDrag[12]+alphaDrag[0]*f[10]+f[0]*alphaDrag[10]+alphaDrag[7]*f[9]+f[7]*alphaDrag[9]+alphaDrag[6]*f[8]+f[6]*alphaDrag[8]+alphaDrag[3]*f[4]+f[3]*alphaDrag[4]); 
  out[20] += 1.732050807568877*(alphaDrag[10]*f[23]+f[10]*alphaDrag[23]+alphaDrag[13]*f[22]+f[13]*alphaDrag[22]+alphaDrag[14]*f[21]+f[14]*alphaDrag[21]+alphaDrag[3]*f[20]+f[3]*alphaDrag[20]+alphaDrag[15]*f[19]+f[15]*alphaDrag[19]+alphaDrag[6]*f[18]+f[6]*alphaDrag[18]+alphaDrag[7]*f[17]+f[7]*alphaDrag[17]+alphaDrag[11]*f[16]+f[11]*alphaDrag[16])+1.936491673103709*(alphaDrag[4]*f[15]+f[4]*alphaDrag[15]+alphaDrag[8]*f[14]+f[8]*alphaDrag[14]+alphaDrag[9]*f[13]+f[9]*alphaDrag[13]+alphaDrag[10]*f[12]+f[10]*alphaDrag[12]+alphaDrag[0]*f[11]+f[0]*alphaDrag[11]+alphaDrag[1]*f[7]+f[1]*alphaDrag[7]+alphaDrag[2]*f[6]+f[2]*alphaDrag[6]+alphaDrag[3]*f[5]+f[3]*alphaDrag[5]); 
  out[21] += 1.732050807568877*(alphaDrag[7]*f[23]+f[7]*alphaDrag[23]+alphaDrag[11]*f[22]+f[11]*alphaDrag[22]+alphaDrag[3]*f[21]+f[3]*alphaDrag[21]+alphaDrag[14]*f[20]+f[14]*alphaDrag[20]+alphaDrag[6]*f[19]+f[6]*alphaDrag[19]+alphaDrag[15]*f[18]+f[15]*alphaDrag[18]+alphaDrag[10]*f[17]+f[10]*alphaDrag[17]+alphaDrag[13]*f[16]+f[13]*alphaDrag[16])+1.936491673103709*(alphaDrag[2]*f[15]+f[2]*alphaDrag[15]+alphaDrag[5]*f[14]+f[5]*alphaDrag[14]+alphaDrag[0]*f[13]+f[0]*alphaDrag[13]+alphaDrag[7]*f[12]+f[7]*alphaDrag[12]+alphaDrag[9]*f[11]+f[9]*alphaDrag[11]+alphaDrag[1]*f[10]+f[1]*alphaDrag[10]+alphaDrag[3]*f[8]+f[3]*alphaDrag[8]+alphaDrag[4]*f[6]+f[4]*alphaDrag[6]); 
  out[22] += 1.732050807568877*(alphaDrag[6]*f[23]+f[6]*alphaDrag[23]+alphaDrag[3]*f[22]+f[3]*alphaDrag[22]+alphaDrag[11]*f[21]+f[11]*alphaDrag[21]+alphaDrag[13]*f[20]+f[13]*alphaDrag[20]+alphaDrag[7]*f[19]+f[7]*alphaDrag[19]+alphaDrag[10]*f[18]+f[10]*alphaDrag[18]+alphaDrag[15]*f[17]+f[15]*alphaDrag[17]+alphaDrag[14]*f[16]+f[14]*alphaDrag[16])+1.936491673103709*(alphaDrag[1]*f[15]+f[1]*alphaDrag[15]+alphaDrag[0]*f[14]+f[0]*alphaDrag[14]+alphaDrag[5]*f[13]+f[5]*alphaDrag[13]+alphaDrag[6]*f[12]+f[6]*alphaDrag[12]+alphaDrag[8]*f[11]+f[8]*alphaDrag[11]+alphaDrag[2]*f[10]+f[2]*alphaDrag[10]+alphaDrag[3]*f[9]+f[3]*alphaDrag[9]+alphaDrag[4]*f[7]+f[4]*alphaDrag[7]); 
  out[23] += 1.732050807568877*(alphaDrag[3]*f[23]+f[3]*alphaDrag[23]+alphaDrag[6]*f[22]+f[6]*alphaDrag[22]+alphaDrag[7]*f[21]+f[7]*alphaDrag[21]+alphaDrag[10]*f[20]+f[10]*alphaDrag[20]+alphaDrag[11]*f[19]+f[11]*alphaDrag[19]+alphaDrag[13]*f[18]+f[13]*alphaDrag[18]+alphaDrag[14]*f[17]+f[14]*alphaDrag[17]+alphaDrag[15]*f[16]+f[15]*alphaDrag[16])+1.936491673103709*(alphaDrag[0]*f[15]+f[0]*alphaDrag[15]+alphaDrag[1]*f[14]+f[1]*alphaDrag[14]+alphaDrag[2]*f[13]+f[2]*alphaDrag[13]+alphaDrag[3]*f[12]+f[3]*alphaDrag[12]+alphaDrag[4]*f[11]+f[4]*alphaDrag[11]+alphaDrag[5]*f[10]+f[5]*alphaDrag[10]+alphaDrag[6]*f[9]+f[6]*alphaDrag[9]+alphaDrag[7]*f[8]+f[7]*alphaDrag[8]); 

  return fabs(0.625*alphaDrag[0]-0.6987712429686843*alphaDrag[16])+fabs(nI[0]*(0.1875*vsqnu[0]-0.2096313728906053*vsqnu[16])); 

} 
