#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_pkpm_drag_vol_3x1v_ser_p1(const double *w, const double *dxv, const double *nu, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[4]:      Cell-center coordinates. 
  // dxv[4]:    Cell spacing. 
  // nu:     collisionality. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvpar = 2.0/dxv[3]; 
  const double dvpar = dxv[3], wvpar = w[3]; 
  double alphaDrag[24]; 
  // Expand rdvpar*(nu*vx) in phase basis.
  alphaDrag[0] = -1.414213562373095*nu[0]*rdvpar*wvpar; 
  alphaDrag[1] = -1.414213562373095*nu[1]*rdvpar*wvpar; 
  alphaDrag[2] = -1.414213562373095*nu[2]*rdvpar*wvpar; 
  alphaDrag[3] = -1.414213562373095*nu[3]*rdvpar*wvpar; 
  alphaDrag[4] = -0.408248290463863*nu[0]*dvpar*rdvpar; 
  alphaDrag[5] = -1.414213562373095*nu[4]*rdvpar*wvpar; 
  alphaDrag[6] = -1.414213562373095*nu[5]*rdvpar*wvpar; 
  alphaDrag[7] = -1.414213562373095*nu[6]*rdvpar*wvpar; 
  alphaDrag[8] = -0.408248290463863*nu[1]*dvpar*rdvpar; 
  alphaDrag[9] = -0.408248290463863*nu[2]*dvpar*rdvpar; 
  alphaDrag[10] = -0.408248290463863*nu[3]*dvpar*rdvpar; 
  alphaDrag[11] = -1.414213562373095*nu[7]*rdvpar*wvpar; 
  alphaDrag[12] = -0.408248290463863*nu[4]*dvpar*rdvpar; 
  alphaDrag[13] = -0.408248290463863*nu[5]*dvpar*rdvpar; 
  alphaDrag[14] = -0.408248290463863*nu[6]*dvpar*rdvpar; 
  alphaDrag[15] = -0.408248290463863*nu[7]*dvpar*rdvpar; 

  out[4] += 0.4330127018922193*(alphaDrag[15]*f[15]+alphaDrag[14]*f[14]+alphaDrag[13]*f[13]+alphaDrag[12]*f[12]+alphaDrag[11]*f[11]+alphaDrag[10]*f[10]+alphaDrag[9]*f[9]+alphaDrag[8]*f[8]+alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[8] += 0.4330127018922193*(alphaDrag[14]*f[15]+f[14]*alphaDrag[15]+alphaDrag[10]*f[13]+f[10]*alphaDrag[13]+alphaDrag[9]*f[12]+f[9]*alphaDrag[12]+alphaDrag[7]*f[11]+f[7]*alphaDrag[11]+alphaDrag[4]*f[8]+f[4]*alphaDrag[8]+alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[9] += 0.4330127018922193*(alphaDrag[13]*f[15]+f[13]*alphaDrag[15]+alphaDrag[10]*f[14]+f[10]*alphaDrag[14]+alphaDrag[8]*f[12]+f[8]*alphaDrag[12]+alphaDrag[6]*f[11]+f[6]*alphaDrag[11]+alphaDrag[4]*f[9]+f[4]*alphaDrag[9]+alphaDrag[3]*f[7]+f[3]*alphaDrag[7]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[10] += 0.4330127018922193*(alphaDrag[12]*f[15]+f[12]*alphaDrag[15]+alphaDrag[9]*f[14]+f[9]*alphaDrag[14]+alphaDrag[8]*f[13]+f[8]*alphaDrag[13]+alphaDrag[5]*f[11]+f[5]*alphaDrag[11]+alphaDrag[4]*f[10]+f[4]*alphaDrag[10]+alphaDrag[2]*f[7]+f[2]*alphaDrag[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[12] += 0.4330127018922193*(alphaDrag[10]*f[15]+f[10]*alphaDrag[15]+alphaDrag[13]*f[14]+f[13]*alphaDrag[14]+alphaDrag[4]*f[12]+f[4]*alphaDrag[12]+alphaDrag[3]*f[11]+f[3]*alphaDrag[11]+alphaDrag[8]*f[9]+f[8]*alphaDrag[9]+alphaDrag[6]*f[7]+f[6]*alphaDrag[7]+alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[13] += 0.4330127018922193*(alphaDrag[9]*f[15]+f[9]*alphaDrag[15]+alphaDrag[12]*f[14]+f[12]*alphaDrag[14]+alphaDrag[4]*f[13]+f[4]*alphaDrag[13]+alphaDrag[2]*f[11]+f[2]*alphaDrag[11]+alphaDrag[8]*f[10]+f[8]*alphaDrag[10]+alphaDrag[5]*f[7]+f[5]*alphaDrag[7]+alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[14] += 0.4330127018922193*(alphaDrag[8]*f[15]+f[8]*alphaDrag[15]+alphaDrag[4]*f[14]+f[4]*alphaDrag[14]+alphaDrag[12]*f[13]+f[12]*alphaDrag[13]+alphaDrag[1]*f[11]+f[1]*alphaDrag[11]+alphaDrag[9]*f[10]+f[9]*alphaDrag[10]+alphaDrag[0]*f[7]+f[0]*alphaDrag[7]+alphaDrag[5]*f[6]+f[5]*alphaDrag[6]+alphaDrag[2]*f[3]+f[2]*alphaDrag[3]); 
  out[15] += 0.4330127018922193*(alphaDrag[4]*f[15]+f[4]*alphaDrag[15]+alphaDrag[8]*f[14]+f[8]*alphaDrag[14]+alphaDrag[9]*f[13]+f[9]*alphaDrag[13]+alphaDrag[10]*f[12]+f[10]*alphaDrag[12]+alphaDrag[0]*f[11]+f[0]*alphaDrag[11]+alphaDrag[1]*f[7]+f[1]*alphaDrag[7]+alphaDrag[2]*f[6]+f[2]*alphaDrag[6]+alphaDrag[3]*f[5]+f[3]*alphaDrag[5]); 
  out[16] += 0.8660254037844386*(alphaDrag[15]*f[23]+alphaDrag[14]*f[22]+alphaDrag[13]*f[21]+alphaDrag[12]*f[20]+alphaDrag[10]*f[19]+alphaDrag[9]*f[18]+alphaDrag[8]*f[17]+alphaDrag[4]*f[16])+0.9682458365518543*(alphaDrag[11]*f[15]+f[11]*alphaDrag[15]+alphaDrag[7]*f[14]+f[7]*alphaDrag[14]+alphaDrag[6]*f[13]+f[6]*alphaDrag[13]+alphaDrag[5]*f[12]+f[5]*alphaDrag[12]+alphaDrag[3]*f[10]+f[3]*alphaDrag[10]+alphaDrag[2]*f[9]+f[2]*alphaDrag[9]+alphaDrag[1]*f[8]+f[1]*alphaDrag[8]+alphaDrag[0]*f[4]+f[0]*alphaDrag[4]); 
  out[17] += 0.8660254037844386*(alphaDrag[14]*f[23]+alphaDrag[15]*f[22]+alphaDrag[10]*f[21]+alphaDrag[9]*f[20]+alphaDrag[13]*f[19]+alphaDrag[12]*f[18]+alphaDrag[4]*f[17]+alphaDrag[8]*f[16])+0.9682458365518543*(alphaDrag[7]*f[15]+f[7]*alphaDrag[15]+alphaDrag[11]*f[14]+f[11]*alphaDrag[14]+alphaDrag[3]*f[13]+f[3]*alphaDrag[13]+alphaDrag[2]*f[12]+f[2]*alphaDrag[12]+alphaDrag[6]*f[10]+f[6]*alphaDrag[10]+alphaDrag[5]*f[9]+f[5]*alphaDrag[9]+alphaDrag[0]*f[8]+f[0]*alphaDrag[8]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4]); 
  out[18] += 0.8660254037844386*(alphaDrag[13]*f[23]+alphaDrag[10]*f[22]+alphaDrag[15]*f[21]+alphaDrag[8]*f[20]+alphaDrag[14]*f[19]+alphaDrag[4]*f[18]+alphaDrag[12]*f[17]+alphaDrag[9]*f[16])+0.9682458365518543*(alphaDrag[6]*f[15]+f[6]*alphaDrag[15]+alphaDrag[3]*f[14]+f[3]*alphaDrag[14]+alphaDrag[11]*f[13]+f[11]*alphaDrag[13]+alphaDrag[1]*f[12]+f[1]*alphaDrag[12]+alphaDrag[7]*f[10]+f[7]*alphaDrag[10]+alphaDrag[0]*f[9]+f[0]*alphaDrag[9]+alphaDrag[5]*f[8]+f[5]*alphaDrag[8]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4]); 
  out[19] += 0.8660254037844386*(alphaDrag[12]*f[23]+alphaDrag[9]*f[22]+alphaDrag[8]*f[21]+alphaDrag[15]*f[20]+alphaDrag[4]*f[19]+alphaDrag[14]*f[18]+alphaDrag[13]*f[17]+alphaDrag[10]*f[16])+0.9682458365518543*(alphaDrag[5]*f[15]+f[5]*alphaDrag[15]+alphaDrag[2]*f[14]+f[2]*alphaDrag[14]+alphaDrag[1]*f[13]+f[1]*alphaDrag[13]+alphaDrag[11]*f[12]+f[11]*alphaDrag[12]+alphaDrag[0]*f[10]+f[0]*alphaDrag[10]+alphaDrag[7]*f[9]+f[7]*alphaDrag[9]+alphaDrag[6]*f[8]+f[6]*alphaDrag[8]+alphaDrag[3]*f[4]+f[3]*alphaDrag[4]); 
  out[20] += 0.8660254037844386*(alphaDrag[10]*f[23]+alphaDrag[13]*f[22]+alphaDrag[14]*f[21]+alphaDrag[4]*f[20]+alphaDrag[15]*f[19]+alphaDrag[8]*f[18]+alphaDrag[9]*f[17]+alphaDrag[12]*f[16])+0.9682458365518543*(alphaDrag[3]*f[15]+f[3]*alphaDrag[15]+alphaDrag[6]*f[14]+f[6]*alphaDrag[14]+alphaDrag[7]*f[13]+f[7]*alphaDrag[13]+alphaDrag[0]*f[12]+f[0]*alphaDrag[12]+alphaDrag[10]*f[11]+f[10]*alphaDrag[11]+alphaDrag[1]*f[9]+f[1]*alphaDrag[9]+alphaDrag[2]*f[8]+f[2]*alphaDrag[8]+alphaDrag[4]*f[5]+f[4]*alphaDrag[5]); 
  out[21] += 0.8660254037844386*(alphaDrag[9]*f[23]+alphaDrag[12]*f[22]+alphaDrag[4]*f[21]+alphaDrag[14]*f[20]+alphaDrag[8]*f[19]+alphaDrag[15]*f[18]+alphaDrag[10]*f[17]+alphaDrag[13]*f[16])+0.9682458365518543*(alphaDrag[2]*f[15]+f[2]*alphaDrag[15]+alphaDrag[5]*f[14]+f[5]*alphaDrag[14]+alphaDrag[0]*f[13]+f[0]*alphaDrag[13]+alphaDrag[7]*f[12]+f[7]*alphaDrag[12]+alphaDrag[9]*f[11]+f[9]*alphaDrag[11]+alphaDrag[1]*f[10]+f[1]*alphaDrag[10]+alphaDrag[3]*f[8]+f[3]*alphaDrag[8]+alphaDrag[4]*f[6]+f[4]*alphaDrag[6]); 
  out[22] += 0.8660254037844386*(alphaDrag[8]*f[23]+alphaDrag[4]*f[22]+alphaDrag[12]*f[21]+alphaDrag[13]*f[20]+alphaDrag[9]*f[19]+alphaDrag[10]*f[18]+alphaDrag[15]*f[17]+alphaDrag[14]*f[16])+0.9682458365518543*(alphaDrag[1]*f[15]+f[1]*alphaDrag[15]+alphaDrag[0]*f[14]+f[0]*alphaDrag[14]+alphaDrag[5]*f[13]+f[5]*alphaDrag[13]+alphaDrag[6]*f[12]+f[6]*alphaDrag[12]+alphaDrag[8]*f[11]+f[8]*alphaDrag[11]+alphaDrag[2]*f[10]+f[2]*alphaDrag[10]+alphaDrag[3]*f[9]+f[3]*alphaDrag[9]+alphaDrag[4]*f[7]+f[4]*alphaDrag[7]); 
  out[23] += 0.8660254037844386*(alphaDrag[4]*f[23]+alphaDrag[8]*f[22]+alphaDrag[9]*f[21]+alphaDrag[10]*f[20]+alphaDrag[12]*f[19]+alphaDrag[13]*f[18]+alphaDrag[14]*f[17]+alphaDrag[15]*f[16])+0.9682458365518543*(alphaDrag[0]*f[15]+f[0]*alphaDrag[15]+alphaDrag[1]*f[14]+f[1]*alphaDrag[14]+alphaDrag[2]*f[13]+f[2]*alphaDrag[13]+alphaDrag[3]*f[12]+f[3]*alphaDrag[12]+alphaDrag[4]*f[11]+f[4]*alphaDrag[11]+alphaDrag[5]*f[10]+f[5]*alphaDrag[10]+alphaDrag[6]*f[9]+f[6]*alphaDrag[9]+alphaDrag[7]*f[8]+f[7]*alphaDrag[8]); 

  return fabs(0.625*alphaDrag[0]); 

} 
