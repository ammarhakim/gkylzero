#include <gkyl_vlasov_lbo_kernels.h> 
GKYL_CU_DH double vlasov_lbo_vol_1x3v_ser_p2(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const double *f, double GKYL_RESTRICT *out) 
{ 
  // w[4]:      Cell-center coordinates. 
  // dxv[4]:    Cell spacing. 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvx2 = 2.0/dxv[1]; 
  const double rdvxSq4 = 4.0/(dxv[1]*dxv[1]); 
  const double rdvy2 = 2.0/dxv[2]; 
  const double rdvySq4 = 4.0/(dxv[2]*dxv[2]); 
  const double rdvz2 = 2.0/dxv[3]; 
  const double rdvzSq4 = 4.0/(dxv[3]*dxv[3]); 

  double alphaDrag[144]; 
  // Expand rdv2*(nu*vx-nuUSumx) in phase basis.
  alphaDrag[0] = (2.828427124746191*nuUSum[0]-4.0*w[1]*nuSum)*rdvx2; 
  alphaDrag[1] = 2.828427124746191*nuUSum[1]*rdvx2; 
  alphaDrag[2] = -1.154700538379252*dxv[1]*nuSum*rdvx2; 
  alphaDrag[11] = 2.828427124746191*nuUSum[2]*rdvx2; 

  // Expand rdv2*(nu*vy-nuUSumy) in phase basis.
  alphaDrag[48] = (2.828427124746191*nuUSum[3]-4.0*w[2]*nuSum)*rdvy2; 
  alphaDrag[49] = 2.828427124746191*nuUSum[4]*rdvy2; 
  alphaDrag[51] = -1.154700538379252*dxv[2]*nuSum*rdvy2; 
  alphaDrag[59] = 2.828427124746191*nuUSum[5]*rdvy2; 

  // Expand rdv2*(nu*vz-nuUSumz) in phase basis.
  alphaDrag[96] = (2.828427124746191*nuUSum[6]-4.0*w[3]*nuSum)*rdvz2; 
  alphaDrag[97] = 2.828427124746191*nuUSum[7]*rdvz2; 
  alphaDrag[100] = -1.154700538379252*dxv[3]*nuSum*rdvz2; 
  alphaDrag[107] = 2.828427124746191*nuUSum[8]*rdvz2; 

  double facDiff[3]; 
  // Expand nuVtSqSum in phase basis.
  facDiff[0] = nuVtSqSum[0]; 
  facDiff[1] = nuVtSqSum[1]; 
  facDiff[2] = nuVtSqSum[2]; 

  // Put together updates due to drag and diffusion terms.
  out[2] += 0.4330127018922193*(alphaDrag[11]*f[11]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.4330127018922193*(f[11]*alphaDrag[59]+f[3]*alphaDrag[51]+f[1]*alphaDrag[49]+f[0]*alphaDrag[48]); 
  out[4] += 0.4330127018922193*(f[11]*alphaDrag[107]+f[4]*alphaDrag[100]+f[1]*alphaDrag[97]+f[0]*alphaDrag[96]); 
  out[5] += 0.3872983346207416*(alphaDrag[1]*f[11]+f[1]*alphaDrag[11])+0.4330127018922193*(alphaDrag[2]*f[5]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[6] += 0.3872983346207416*f[1]*alphaDrag[59]+0.4330127018922193*f[6]*alphaDrag[51]+0.3872983346207416*f[11]*alphaDrag[49]+0.4330127018922193*(f[0]*alphaDrag[49]+f[1]*alphaDrag[48]); 
  out[7] += 0.4330127018922193*(f[19]*alphaDrag[59]+f[7]*alphaDrag[51]+f[5]*alphaDrag[49]+f[2]*alphaDrag[48]+alphaDrag[11]*f[21]+alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+alphaDrag[0]*f[3]); 
  out[8] += 0.3872983346207416*f[1]*alphaDrag[107]+0.4330127018922193*f[8]*alphaDrag[100]+0.3872983346207416*f[11]*alphaDrag[97]+0.4330127018922193*(f[0]*alphaDrag[97]+f[1]*alphaDrag[96]); 
  out[9] += 0.4330127018922193*(f[19]*alphaDrag[107]+f[9]*alphaDrag[100]+f[5]*alphaDrag[97]+f[2]*alphaDrag[96]+alphaDrag[11]*f[25]+alphaDrag[2]*f[9]+alphaDrag[1]*f[8]+alphaDrag[0]*f[4]); 
  out[10] += 0.4330127018922193*(f[21]*alphaDrag[107]+f[10]*alphaDrag[100]+f[6]*alphaDrag[97]+f[3]*alphaDrag[96]+f[25]*alphaDrag[59]+f[10]*alphaDrag[51]+f[8]*alphaDrag[49]+f[4]*alphaDrag[48]); 
  out[12] += 4.743416490252569*(facDiff[2]*f[11]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvxSq4+0.9682458365518543*alphaDrag[11]*f[19]+0.8660254037844386*alphaDrag[2]*f[12]+0.9682458365518543*(alphaDrag[1]*f[5]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[13] += 4.743416490252569*(facDiff[2]*f[11]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvySq4+0.9682458365518543*f[21]*alphaDrag[59]+0.8660254037844386*f[13]*alphaDrag[51]+0.9682458365518543*(f[0]*alphaDrag[51]+f[6]*alphaDrag[49]+f[3]*alphaDrag[48]); 
  out[14] += 4.743416490252569*(facDiff[2]*f[11]+f[1]*facDiff[1]+f[0]*facDiff[0])*rdvzSq4+0.9682458365518543*f[25]*alphaDrag[107]+0.8660254037844386*f[14]*alphaDrag[100]+0.9682458365518543*(f[0]*alphaDrag[100]+f[8]*alphaDrag[97]+f[4]*alphaDrag[96]); 
  out[15] += 0.3872983346207416*f[5]*alphaDrag[59]+0.4330127018922193*f[15]*alphaDrag[51]+0.3872983346207416*f[19]*alphaDrag[49]+0.4330127018922193*(f[2]*alphaDrag[49]+f[5]*alphaDrag[48])+0.3872983346207416*alphaDrag[1]*f[21]+0.4330127018922193*alphaDrag[2]*f[15]+0.3872983346207416*f[6]*alphaDrag[11]+0.4330127018922193*(alphaDrag[0]*f[6]+alphaDrag[1]*f[3]); 
  out[16] += 0.3872983346207416*f[5]*alphaDrag[107]+0.4330127018922193*f[16]*alphaDrag[100]+0.3872983346207416*f[19]*alphaDrag[97]+0.4330127018922193*(f[2]*alphaDrag[97]+f[5]*alphaDrag[96])+0.3872983346207416*alphaDrag[1]*f[25]+0.4330127018922193*alphaDrag[2]*f[16]+0.3872983346207416*f[8]*alphaDrag[11]+0.4330127018922193*(alphaDrag[0]*f[8]+alphaDrag[1]*f[4]); 
  out[17] += 0.3872983346207416*f[6]*alphaDrag[107]+0.4330127018922193*f[17]*alphaDrag[100]+0.3872983346207416*f[21]*alphaDrag[97]+0.4330127018922193*(f[3]*alphaDrag[97]+f[6]*alphaDrag[96])+0.3872983346207416*f[8]*alphaDrag[59]+0.4330127018922193*f[17]*alphaDrag[51]+0.3872983346207416*f[25]*alphaDrag[49]+0.4330127018922193*(f[4]*alphaDrag[49]+f[8]*alphaDrag[48]); 
  out[18] += 0.4330127018922193*(f[32]*alphaDrag[107]+f[18]*alphaDrag[100]+f[15]*alphaDrag[97]+f[7]*alphaDrag[96]+f[35]*alphaDrag[59]+f[18]*alphaDrag[51]+f[16]*alphaDrag[49]+f[9]*alphaDrag[48]+alphaDrag[11]*f[37]+alphaDrag[2]*f[18]+alphaDrag[1]*f[17]+alphaDrag[0]*f[10]); 
  out[19] += 0.4330127018922193*alphaDrag[2]*f[19]+0.276641667586244*alphaDrag[11]*f[11]+0.4330127018922193*(alphaDrag[0]*f[11]+f[0]*alphaDrag[11])+0.3872983346207416*alphaDrag[1]*f[1]; 
  out[20] += (4.242640687119286*(facDiff[1]*f[11]+f[1]*facDiff[2])+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvxSq4+0.8660254037844386*(alphaDrag[2]*f[20]+alphaDrag[1]*f[19]+f[5]*alphaDrag[11])+0.9682458365518543*(alphaDrag[0]*f[5]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[21] += 0.276641667586244*f[11]*alphaDrag[59]+0.4330127018922193*(f[0]*alphaDrag[59]+f[21]*alphaDrag[51])+0.3872983346207416*f[1]*alphaDrag[49]+0.4330127018922193*f[11]*alphaDrag[48]; 
  out[22] += 4.743416490252569*(facDiff[2]*f[21]+facDiff[1]*f[6]+facDiff[0]*f[3])*rdvxSq4+0.4330127018922193*(f[22]*alphaDrag[51]+f[20]*alphaDrag[49]+f[12]*alphaDrag[48])+0.9682458365518543*alphaDrag[11]*f[32]+0.8660254037844386*alphaDrag[2]*f[22]+0.9682458365518543*(alphaDrag[1]*f[15]+alphaDrag[0]*f[7]+alphaDrag[2]*f[3]); 
  out[23] += (4.242640687119286*(facDiff[1]*f[11]+f[1]*facDiff[2])+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvySq4+0.8660254037844386*f[6]*alphaDrag[59]+(0.8660254037844386*f[23]+0.9682458365518543*f[1])*alphaDrag[51]+0.8660254037844386*f[21]*alphaDrag[49]+0.9682458365518543*(f[3]*alphaDrag[49]+f[6]*alphaDrag[48]); 
  out[24] += 4.743416490252569*(facDiff[2]*f[19]+facDiff[1]*f[5]+facDiff[0]*f[2])*rdvySq4+0.9682458365518543*f[32]*alphaDrag[59]+0.8660254037844386*f[24]*alphaDrag[51]+0.9682458365518543*(f[2]*alphaDrag[51]+f[15]*alphaDrag[49]+f[7]*alphaDrag[48])+0.4330127018922193*(alphaDrag[2]*f[24]+alphaDrag[1]*f[23]+alphaDrag[0]*f[13]); 
  out[25] += 0.276641667586244*f[11]*alphaDrag[107]+0.4330127018922193*(f[0]*alphaDrag[107]+f[25]*alphaDrag[100])+0.3872983346207416*f[1]*alphaDrag[97]+0.4330127018922193*f[11]*alphaDrag[96]; 
  out[26] += 4.743416490252569*(facDiff[2]*f[25]+facDiff[1]*f[8]+facDiff[0]*f[4])*rdvxSq4+0.4330127018922193*(f[26]*alphaDrag[100]+f[20]*alphaDrag[97]+f[12]*alphaDrag[96])+0.9682458365518543*alphaDrag[11]*f[35]+0.8660254037844386*alphaDrag[2]*f[26]+0.9682458365518543*(alphaDrag[1]*f[16]+alphaDrag[0]*f[9]+alphaDrag[2]*f[4]); 
  out[27] += 4.743416490252569*(facDiff[2]*f[25]+facDiff[1]*f[8]+facDiff[0]*f[4])*rdvySq4+0.4330127018922193*(f[27]*alphaDrag[100]+f[23]*alphaDrag[97]+f[13]*alphaDrag[96])+0.9682458365518543*f[37]*alphaDrag[59]+0.8660254037844386*f[27]*alphaDrag[51]+0.9682458365518543*(f[4]*alphaDrag[51]+f[17]*alphaDrag[49]+f[10]*alphaDrag[48]); 
  out[28] += (4.242640687119286*(facDiff[1]*f[11]+f[1]*facDiff[2])+4.743416490252569*(f[0]*facDiff[1]+facDiff[0]*f[1]))*rdvzSq4+0.8660254037844386*f[8]*alphaDrag[107]+(0.8660254037844386*f[28]+0.9682458365518543*f[1])*alphaDrag[100]+0.8660254037844386*f[25]*alphaDrag[97]+0.9682458365518543*(f[4]*alphaDrag[97]+f[8]*alphaDrag[96]); 
  out[29] += 4.743416490252569*(facDiff[2]*f[19]+facDiff[1]*f[5]+facDiff[0]*f[2])*rdvzSq4+0.9682458365518543*f[35]*alphaDrag[107]+0.8660254037844386*f[29]*alphaDrag[100]+0.9682458365518543*(f[2]*alphaDrag[100]+f[16]*alphaDrag[97]+f[9]*alphaDrag[96])+0.4330127018922193*(alphaDrag[2]*f[29]+alphaDrag[1]*f[28]+alphaDrag[0]*f[14]); 
  out[30] += 4.743416490252569*(facDiff[2]*f[21]+facDiff[1]*f[6]+facDiff[0]*f[3])*rdvzSq4+0.9682458365518543*f[37]*alphaDrag[107]+0.8660254037844386*f[30]*alphaDrag[100]+0.9682458365518543*(f[3]*alphaDrag[100]+f[17]*alphaDrag[97]+f[10]*alphaDrag[96])+0.4330127018922193*(f[30]*alphaDrag[51]+f[28]*alphaDrag[49]+f[14]*alphaDrag[48]); 
  out[31] += 0.3872983346207416*f[15]*alphaDrag[107]+0.4330127018922193*f[31]*alphaDrag[100]+0.3872983346207416*f[32]*alphaDrag[97]+0.4330127018922193*(f[7]*alphaDrag[97]+f[15]*alphaDrag[96])+0.3872983346207416*f[16]*alphaDrag[59]+0.4330127018922193*f[31]*alphaDrag[51]+0.3872983346207416*f[35]*alphaDrag[49]+0.4330127018922193*(f[9]*alphaDrag[49]+f[16]*alphaDrag[48])+0.3872983346207416*alphaDrag[1]*f[37]+0.4330127018922193*alphaDrag[2]*f[31]+0.3872983346207416*alphaDrag[11]*f[17]+0.4330127018922193*(alphaDrag[0]*f[17]+alphaDrag[1]*f[10]); 
  out[32] += 0.276641667586244*f[19]*alphaDrag[59]+0.4330127018922193*(f[2]*alphaDrag[59]+f[32]*alphaDrag[51])+0.3872983346207416*f[5]*alphaDrag[49]+0.4330127018922193*(f[19]*alphaDrag[48]+alphaDrag[2]*f[32])+0.276641667586244*alphaDrag[11]*f[21]+0.4330127018922193*(alphaDrag[0]*f[21]+f[3]*alphaDrag[11])+0.3872983346207416*alphaDrag[1]*f[6]; 
  out[33] += (4.242640687119286*facDiff[1]*f[21]+4.743416490252569*(facDiff[0]*f[6]+facDiff[1]*f[3])+4.242640687119286*facDiff[2]*f[6])*rdvxSq4+0.3872983346207416*f[20]*alphaDrag[59]+0.4330127018922193*(f[33]*alphaDrag[51]+f[12]*alphaDrag[49]+f[20]*alphaDrag[48])+0.8660254037844386*(alphaDrag[2]*f[33]+alphaDrag[1]*f[32]+alphaDrag[11]*f[15])+0.9682458365518543*(alphaDrag[0]*f[15]+alphaDrag[1]*f[7]+alphaDrag[2]*f[6]); 
  out[34] += (4.242640687119286*facDiff[1]*f[19]+4.743416490252569*(facDiff[0]*f[5]+facDiff[1]*f[2])+4.242640687119286*facDiff[2]*f[5])*rdvySq4+0.8660254037844386*f[15]*alphaDrag[59]+(0.8660254037844386*f[34]+0.9682458365518543*f[5])*alphaDrag[51]+0.8660254037844386*f[32]*alphaDrag[49]+0.9682458365518543*(f[7]*alphaDrag[49]+f[15]*alphaDrag[48])+0.4330127018922193*alphaDrag[2]*f[34]+0.3872983346207416*alphaDrag[11]*f[23]+0.4330127018922193*(alphaDrag[0]*f[23]+alphaDrag[1]*f[13]); 
  out[35] += 0.276641667586244*f[19]*alphaDrag[107]+0.4330127018922193*(f[2]*alphaDrag[107]+f[35]*alphaDrag[100])+0.3872983346207416*f[5]*alphaDrag[97]+0.4330127018922193*(f[19]*alphaDrag[96]+alphaDrag[2]*f[35])+0.276641667586244*alphaDrag[11]*f[25]+0.4330127018922193*(alphaDrag[0]*f[25]+f[4]*alphaDrag[11])+0.3872983346207416*alphaDrag[1]*f[8]; 
  out[36] += (4.242640687119286*facDiff[1]*f[25]+4.743416490252569*(facDiff[0]*f[8]+facDiff[1]*f[4])+4.242640687119286*facDiff[2]*f[8])*rdvxSq4+0.3872983346207416*f[20]*alphaDrag[107]+0.4330127018922193*(f[36]*alphaDrag[100]+f[12]*alphaDrag[97]+f[20]*alphaDrag[96])+0.8660254037844386*(alphaDrag[2]*f[36]+alphaDrag[1]*f[35]+alphaDrag[11]*f[16])+0.9682458365518543*(alphaDrag[0]*f[16]+alphaDrag[1]*f[9]+alphaDrag[2]*f[8]); 
  out[37] += 0.276641667586244*f[21]*alphaDrag[107]+0.4330127018922193*(f[3]*alphaDrag[107]+f[37]*alphaDrag[100])+0.3872983346207416*f[6]*alphaDrag[97]+0.4330127018922193*f[21]*alphaDrag[96]+0.276641667586244*f[25]*alphaDrag[59]+0.4330127018922193*(f[4]*alphaDrag[59]+f[37]*alphaDrag[51])+0.3872983346207416*f[8]*alphaDrag[49]+0.4330127018922193*f[25]*alphaDrag[48]; 
  out[38] += 4.743416490252569*(facDiff[2]*f[37]+facDiff[1]*f[17]+facDiff[0]*f[10])*rdvxSq4+0.4330127018922193*(f[38]*alphaDrag[100]+f[33]*alphaDrag[97]+f[22]*alphaDrag[96]+f[38]*alphaDrag[51]+f[36]*alphaDrag[49]+f[26]*alphaDrag[48])+0.9682458365518543*alphaDrag[11]*f[44]+0.8660254037844386*alphaDrag[2]*f[38]+0.9682458365518543*(alphaDrag[1]*f[31]+alphaDrag[0]*f[18]+alphaDrag[2]*f[10]); 
  out[39] += (4.242640687119286*facDiff[1]*f[25]+4.743416490252569*(facDiff[0]*f[8]+facDiff[1]*f[4])+4.242640687119286*facDiff[2]*f[8])*rdvySq4+0.3872983346207416*f[23]*alphaDrag[107]+0.4330127018922193*(f[39]*alphaDrag[100]+f[13]*alphaDrag[97]+f[23]*alphaDrag[96])+0.8660254037844386*f[17]*alphaDrag[59]+(0.8660254037844386*f[39]+0.9682458365518543*f[8])*alphaDrag[51]+0.8660254037844386*f[37]*alphaDrag[49]+0.9682458365518543*(f[10]*alphaDrag[49]+f[17]*alphaDrag[48]); 
  out[40] += 4.743416490252569*(facDiff[2]*f[35]+facDiff[1]*f[16]+facDiff[0]*f[9])*rdvySq4+0.4330127018922193*(f[40]*alphaDrag[100]+f[34]*alphaDrag[97]+f[24]*alphaDrag[96])+0.9682458365518543*f[44]*alphaDrag[59]+0.8660254037844386*f[40]*alphaDrag[51]+0.9682458365518543*(f[9]*alphaDrag[51]+f[31]*alphaDrag[49]+f[18]*alphaDrag[48])+0.4330127018922193*(alphaDrag[2]*f[40]+alphaDrag[1]*f[39]+alphaDrag[0]*f[27]); 
  out[41] += (4.242640687119286*facDiff[1]*f[19]+4.743416490252569*(facDiff[0]*f[5]+facDiff[1]*f[2])+4.242640687119286*facDiff[2]*f[5])*rdvzSq4+0.8660254037844386*f[16]*alphaDrag[107]+(0.8660254037844386*f[41]+0.9682458365518543*f[5])*alphaDrag[100]+0.8660254037844386*f[35]*alphaDrag[97]+0.9682458365518543*(f[9]*alphaDrag[97]+f[16]*alphaDrag[96])+0.4330127018922193*alphaDrag[2]*f[41]+0.3872983346207416*alphaDrag[11]*f[28]+0.4330127018922193*(alphaDrag[0]*f[28]+alphaDrag[1]*f[14]); 
  out[42] += (4.242640687119286*facDiff[1]*f[21]+4.743416490252569*(facDiff[0]*f[6]+facDiff[1]*f[3])+4.242640687119286*facDiff[2]*f[6])*rdvzSq4+0.8660254037844386*f[17]*alphaDrag[107]+(0.8660254037844386*f[42]+0.9682458365518543*f[6])*alphaDrag[100]+0.8660254037844386*f[37]*alphaDrag[97]+0.9682458365518543*(f[10]*alphaDrag[97]+f[17]*alphaDrag[96])+0.3872983346207416*f[28]*alphaDrag[59]+0.4330127018922193*(f[42]*alphaDrag[51]+f[14]*alphaDrag[49]+f[28]*alphaDrag[48]); 
  out[43] += 4.743416490252569*(facDiff[2]*f[32]+facDiff[1]*f[15]+facDiff[0]*f[7])*rdvzSq4+0.9682458365518543*f[44]*alphaDrag[107]+0.8660254037844386*f[43]*alphaDrag[100]+0.9682458365518543*(f[7]*alphaDrag[100]+f[31]*alphaDrag[97]+f[18]*alphaDrag[96])+0.4330127018922193*(f[43]*alphaDrag[51]+f[41]*alphaDrag[49]+f[29]*alphaDrag[48]+alphaDrag[2]*f[43]+alphaDrag[1]*f[42]+alphaDrag[0]*f[30]); 
  out[44] += 0.276641667586244*f[32]*alphaDrag[107]+0.4330127018922193*(f[7]*alphaDrag[107]+f[44]*alphaDrag[100])+0.3872983346207416*f[15]*alphaDrag[97]+0.4330127018922193*f[32]*alphaDrag[96]+0.276641667586244*f[35]*alphaDrag[59]+0.4330127018922193*(f[9]*alphaDrag[59]+f[44]*alphaDrag[51])+0.3872983346207416*f[16]*alphaDrag[49]+0.4330127018922193*(f[35]*alphaDrag[48]+alphaDrag[2]*f[44])+(0.276641667586244*alphaDrag[11]+0.4330127018922193*alphaDrag[0])*f[37]+0.3872983346207416*alphaDrag[1]*f[17]+0.4330127018922193*f[10]*alphaDrag[11]; 
  out[45] += (4.242640687119286*facDiff[1]*f[37]+4.743416490252569*(facDiff[0]*f[17]+facDiff[1]*f[10])+4.242640687119286*facDiff[2]*f[17])*rdvxSq4+0.3872983346207416*f[33]*alphaDrag[107]+0.4330127018922193*(f[45]*alphaDrag[100]+f[22]*alphaDrag[97]+f[33]*alphaDrag[96])+0.3872983346207416*f[36]*alphaDrag[59]+0.4330127018922193*(f[45]*alphaDrag[51]+f[26]*alphaDrag[49]+f[36]*alphaDrag[48])+0.8660254037844386*(alphaDrag[2]*f[45]+alphaDrag[1]*f[44]+alphaDrag[11]*f[31])+0.9682458365518543*(alphaDrag[0]*f[31]+alphaDrag[1]*f[18]+alphaDrag[2]*f[17]); 
  out[46] += (4.242640687119286*facDiff[1]*f[35]+4.743416490252569*(facDiff[0]*f[16]+facDiff[1]*f[9])+4.242640687119286*facDiff[2]*f[16])*rdvySq4+0.3872983346207416*f[34]*alphaDrag[107]+0.4330127018922193*(f[46]*alphaDrag[100]+f[24]*alphaDrag[97]+f[34]*alphaDrag[96])+0.8660254037844386*f[31]*alphaDrag[59]+(0.8660254037844386*f[46]+0.9682458365518543*f[16])*alphaDrag[51]+0.8660254037844386*f[44]*alphaDrag[49]+0.9682458365518543*(f[18]*alphaDrag[49]+f[31]*alphaDrag[48])+0.4330127018922193*alphaDrag[2]*f[46]+0.3872983346207416*alphaDrag[11]*f[39]+0.4330127018922193*(alphaDrag[0]*f[39]+alphaDrag[1]*f[27]); 
  out[47] += (4.242640687119286*facDiff[1]*f[32]+4.743416490252569*(facDiff[0]*f[15]+facDiff[1]*f[7])+4.242640687119286*facDiff[2]*f[15])*rdvzSq4+0.8660254037844386*f[31]*alphaDrag[107]+(0.8660254037844386*f[47]+0.9682458365518543*f[15])*alphaDrag[100]+0.8660254037844386*f[44]*alphaDrag[97]+0.9682458365518543*(f[18]*alphaDrag[97]+f[31]*alphaDrag[96])+0.3872983346207416*f[41]*alphaDrag[59]+0.4330127018922193*(f[47]*alphaDrag[51]+f[29]*alphaDrag[49]+f[41]*alphaDrag[48]+alphaDrag[2]*f[47])+0.3872983346207416*alphaDrag[11]*f[42]+0.4330127018922193*(alphaDrag[0]*f[42]+alphaDrag[1]*f[30]); 

  return fabs(0.125*alphaDrag[0]-0.1397542485937369*alphaDrag[11])+fabs(0.125*alphaDrag[48]-0.1397542485937369*alphaDrag[59])+fabs(0.125*alphaDrag[96]-0.1397542485937369*alphaDrag[107])+fabs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvxSq4)+fabs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvySq4)+fabs((1.272792206135785*facDiff[0]-1.42302494707577*facDiff[2])*rdvzSq4); 

} 
