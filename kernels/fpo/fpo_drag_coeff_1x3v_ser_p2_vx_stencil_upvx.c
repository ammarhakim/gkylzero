#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_drag_coeff_1x3v_vx_ser_p2_upvx(const double *dxv, const double gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_h_stencil[3]: 3 cell stencil of Rosenbluth potential H. 
  // fpo_dhdv_surf: Surface projection of dH/dv in center cell. 
  // drag_coeff: Output array for drag coefficient. 

  double dv1 = 2.0/dxv[1]; 

  const double* H_L = fpo_h_stencil[0]; 
  const double* H_C = fpo_h_stencil[1]; 
  
  const double *fpo_dHdv_surf_C_vx = &fpo_dhdv_surf[0]; 
  const double *fpo_dHdv_surf_C_vy = &fpo_dhdv_surf[20]; 
  const double *fpo_dHdv_surf_C_vz = &fpo_dhdv_surf[40]; 
  
  const double* dHdv_surf_C = fpo_dHdv_surf_C_vx; 
  
  double *drag_coeff_vx = &drag_coeff[0]; 
  double *drag_coeff_vy = &drag_coeff[48]; 
  double *drag_coeff_vz = &drag_coeff[96]; 
  
  double *out = drag_coeff_vx; 
  
  out[0] = -(0.8441156615061707*H_L[12]*dv1*gamma)+2.01805134969356*H_C[12]*dv1*gamma-1.1691342951089918*H_L[2]*dv1*gamma+1.6454482671904334*H_C[2]*dv1*gamma-0.8125*H_L[0]*dv1*gamma+0.8125*H_C[0]*dv1*gamma+0.1414213562373096*dHdv_surf_C[0]*gamma; 
  out[1] = -(0.8441156615061707*H_L[20]*dv1*gamma)+2.0180513496935606*H_C[20]*dv1*gamma-1.1691342951089916*H_L[5]*dv1*gamma+1.6454482671904334*H_C[5]*dv1*gamma-0.8125*H_L[1]*dv1*gamma+0.8125*H_C[1]*dv1*gamma+0.1414213562373096*dHdv_surf_C[1]*gamma; 
  out[2] = 0.232379000772445*H_L[12]*dv1*gamma+5.189797683917939*H_C[12]*dv1*gamma+0.41250000000000003*H_L[2]*dv1*gamma+0.41250000000000003*H_C[2]*dv1*gamma+0.32475952641916445*H_L[0]*dv1*gamma-0.32475952641916445*H_C[0]*dv1*gamma+0.24494897427831794*dHdv_surf_C[0]*gamma; 
  out[3] = -(0.8441156615061707*H_L[22]*dv1*gamma)+2.0180513496935606*H_C[22]*dv1*gamma-1.1691342951089916*H_L[7]*dv1*gamma+1.6454482671904334*H_C[7]*dv1*gamma-0.8125*H_L[3]*dv1*gamma+0.8125*H_C[3]*dv1*gamma+0.1414213562373096*dHdv_surf_C[2]*gamma; 
  out[4] = -(0.8441156615061707*H_L[26]*dv1*gamma)+2.0180513496935606*H_C[26]*dv1*gamma-1.1691342951089916*H_L[9]*dv1*gamma+1.6454482671904334*H_C[9]*dv1*gamma-0.8125*H_L[4]*dv1*gamma+0.8125*H_C[4]*dv1*gamma+0.1414213562373096*dHdv_surf_C[3]*gamma; 
  out[5] = 0.232379000772445*H_L[20]*dv1*gamma+5.189797683917939*H_C[20]*dv1*gamma+0.41250000000000003*H_L[5]*dv1*gamma+0.41250000000000003*H_C[5]*dv1*gamma+0.32475952641916445*H_L[1]*dv1*gamma-0.32475952641916445*H_C[1]*dv1*gamma+0.24494897427831794*dHdv_surf_C[1]*gamma; 
  out[6] = -(0.8441156615061707*H_L[33]*dv1*gamma)+2.01805134969356*H_C[33]*dv1*gamma-1.1691342951089918*H_L[15]*dv1*gamma+1.6454482671904334*H_C[15]*dv1*gamma-0.8125*H_L[6]*dv1*gamma+0.8125*H_C[6]*dv1*gamma+0.1414213562373096*dHdv_surf_C[4]*gamma; 
  out[7] = 0.232379000772445*H_L[22]*dv1*gamma+5.189797683917939*H_C[22]*dv1*gamma+0.41250000000000003*H_L[7]*dv1*gamma+0.41250000000000003*H_C[7]*dv1*gamma+0.32475952641916445*H_L[3]*dv1*gamma-0.32475952641916445*H_C[3]*dv1*gamma+0.24494897427831794*dHdv_surf_C[2]*gamma; 
  out[8] = -(0.8441156615061707*H_L[36]*dv1*gamma)+2.01805134969356*H_C[36]*dv1*gamma-1.1691342951089918*H_L[16]*dv1*gamma+1.6454482671904334*H_C[16]*dv1*gamma-0.8125*H_L[8]*dv1*gamma+0.8125*H_C[8]*dv1*gamma+0.1414213562373096*dHdv_surf_C[5]*gamma; 
  out[9] = 0.232379000772445*H_L[26]*dv1*gamma+5.189797683917939*H_C[26]*dv1*gamma+0.41250000000000003*H_L[9]*dv1*gamma+0.41250000000000003*H_C[9]*dv1*gamma+0.32475952641916445*H_L[4]*dv1*gamma-0.32475952641916445*H_C[4]*dv1*gamma+0.24494897427831794*dHdv_surf_C[3]*gamma; 
  out[10] = -(0.8441156615061707*H_L[38]*dv1*gamma)+2.01805134969356*H_C[38]*dv1*gamma-1.1691342951089918*H_L[18]*dv1*gamma+1.6454482671904334*H_C[18]*dv1*gamma-0.8125*H_L[10]*dv1*gamma+0.8125*H_C[10]*dv1*gamma+0.1414213562373096*dHdv_surf_C[6]*gamma; 
  out[11] = -(1.1691342951089922*H_L[19]*dv1*gamma)+1.6454482671904336*H_C[19]*dv1*gamma-0.8125*H_L[11]*dv1*gamma+0.8125*H_C[11]*dv1*gamma+0.1414213562373096*dHdv_surf_C[7]*gamma; 
  out[12] = -(1.8875000000000002*H_L[12]*dv1*gamma)+4.5125*H_C[12]*dv1*gamma-2.6142637586900057*H_L[2]*dv1*gamma-4.066632513517788*H_C[2]*dv1*gamma-1.8168052317185797*H_L[0]*dv1*gamma+1.8168052317185797*H_C[0]*dv1*gamma+0.3162277660168381*dHdv_surf_C[0]*gamma; 
  out[13] = -(1.1691342951089922*H_L[24]*dv1*gamma)+1.6454482671904336*H_C[24]*dv1*gamma-0.8125*H_L[13]*dv1*gamma+0.8125*H_C[13]*dv1*gamma+0.1414213562373096*dHdv_surf_C[8]*gamma; 
  out[14] = -(1.1691342951089922*H_L[29]*dv1*gamma)+1.6454482671904336*H_C[29]*dv1*gamma-0.8125*H_L[14]*dv1*gamma+0.8125*H_C[14]*dv1*gamma+0.1414213562373096*dHdv_surf_C[9]*gamma; 
  out[15] = 0.232379000772445*H_L[33]*dv1*gamma+5.189797683917939*H_C[33]*dv1*gamma+0.41250000000000003*H_L[15]*dv1*gamma+0.41250000000000003*H_C[15]*dv1*gamma+0.32475952641916445*H_L[6]*dv1*gamma-0.32475952641916445*H_C[6]*dv1*gamma+0.24494897427831794*dHdv_surf_C[4]*gamma; 
  out[16] = 0.232379000772445*H_L[36]*dv1*gamma+5.189797683917939*H_C[36]*dv1*gamma+0.41250000000000003*H_L[16]*dv1*gamma+0.41250000000000003*H_C[16]*dv1*gamma+0.32475952641916445*H_L[8]*dv1*gamma-0.32475952641916445*H_C[8]*dv1*gamma+0.24494897427831794*dHdv_surf_C[5]*gamma; 
  out[17] = -(0.8441156615061707*H_L[45]*dv1*gamma)+2.0180513496935606*H_C[45]*dv1*gamma-1.1691342951089916*H_L[31]*dv1*gamma+1.6454482671904334*H_C[31]*dv1*gamma-0.8125*H_L[17]*dv1*gamma+0.8125*H_C[17]*dv1*gamma+0.1414213562373096*dHdv_surf_C[10]*gamma; 
  out[18] = 0.232379000772445*H_L[38]*dv1*gamma+5.189797683917939*H_C[38]*dv1*gamma+0.41250000000000003*H_L[18]*dv1*gamma+0.41250000000000003*H_C[18]*dv1*gamma+0.32475952641916445*H_L[10]*dv1*gamma-0.32475952641916445*H_C[10]*dv1*gamma+0.24494897427831794*dHdv_surf_C[6]*gamma; 
  out[19] = 0.41250000000000003*H_L[19]*dv1*gamma+0.41250000000000003*H_C[19]*dv1*gamma+0.32475952641916456*H_L[11]*dv1*gamma-0.32475952641916456*H_C[11]*dv1*gamma+0.24494897427831797*dHdv_surf_C[7]*gamma; 
  out[20] = -(1.8875*H_L[20]*dv1*gamma)+4.5125*H_C[20]*dv1*gamma-2.6142637586900066*H_L[5]*dv1*gamma-4.066632513517788*H_C[5]*dv1*gamma-1.816805231718579*H_L[1]*dv1*gamma+1.816805231718579*H_C[1]*dv1*gamma+0.3162277660168381*dHdv_surf_C[1]*gamma; 
  out[21] = -(1.1691342951089922*H_L[32]*dv1*gamma)+1.6454482671904336*H_C[32]*dv1*gamma-0.8125*H_L[21]*dv1*gamma+0.8125*H_C[21]*dv1*gamma+0.1414213562373096*dHdv_surf_C[11]*gamma; 
  out[22] = -(1.8875*H_L[22]*dv1*gamma)+4.5125*H_C[22]*dv1*gamma-2.6142637586900066*H_L[7]*dv1*gamma-4.066632513517788*H_C[7]*dv1*gamma-1.816805231718579*H_L[3]*dv1*gamma+1.816805231718579*H_C[3]*dv1*gamma+0.3162277660168381*dHdv_surf_C[2]*gamma; 
  out[23] = -(1.1691342951089922*H_L[34]*dv1*gamma)+1.6454482671904336*H_C[34]*dv1*gamma-0.8125*H_L[23]*dv1*gamma+0.8125*H_C[23]*dv1*gamma+0.1414213562373096*dHdv_surf_C[12]*gamma; 
  out[24] = 0.41250000000000003*H_L[24]*dv1*gamma+0.41250000000000003*H_C[24]*dv1*gamma+0.32475952641916456*H_L[13]*dv1*gamma-0.32475952641916456*H_C[13]*dv1*gamma+0.24494897427831797*dHdv_surf_C[8]*gamma; 
  out[25] = -(1.1691342951089922*H_L[35]*dv1*gamma)+1.6454482671904336*H_C[35]*dv1*gamma-0.8125*H_L[25]*dv1*gamma+0.8125*H_C[25]*dv1*gamma+0.1414213562373096*dHdv_surf_C[13]*gamma; 
  out[26] = -(1.8875*H_L[26]*dv1*gamma)+4.5125*H_C[26]*dv1*gamma-2.6142637586900066*H_L[9]*dv1*gamma-4.066632513517788*H_C[9]*dv1*gamma-1.816805231718579*H_L[4]*dv1*gamma+1.816805231718579*H_C[4]*dv1*gamma+0.3162277660168381*dHdv_surf_C[3]*gamma; 
  out[27] = -(1.1691342951089922*H_L[40]*dv1*gamma)+1.6454482671904336*H_C[40]*dv1*gamma-0.8125*H_L[27]*dv1*gamma+0.8125*H_C[27]*dv1*gamma+0.1414213562373096*dHdv_surf_C[14]*gamma; 
  out[28] = -(1.1691342951089922*H_L[41]*dv1*gamma)+1.6454482671904336*H_C[41]*dv1*gamma-0.8125*H_L[28]*dv1*gamma+0.8125*H_C[28]*dv1*gamma+0.1414213562373096*dHdv_surf_C[15]*gamma; 
  out[29] = 0.41250000000000003*H_L[29]*dv1*gamma+0.41250000000000003*H_C[29]*dv1*gamma+0.32475952641916456*H_L[14]*dv1*gamma-0.32475952641916456*H_C[14]*dv1*gamma+0.24494897427831797*dHdv_surf_C[9]*gamma; 
  out[30] = -(1.1691342951089922*H_L[43]*dv1*gamma)+1.6454482671904336*H_C[43]*dv1*gamma-0.8125*H_L[30]*dv1*gamma+0.8125*H_C[30]*dv1*gamma+0.1414213562373096*dHdv_surf_C[16]*gamma; 
  out[31] = 0.232379000772445*H_L[45]*dv1*gamma+5.189797683917939*H_C[45]*dv1*gamma+0.41250000000000003*H_L[31]*dv1*gamma+0.41250000000000003*H_C[31]*dv1*gamma+0.32475952641916445*H_L[17]*dv1*gamma-0.32475952641916445*H_C[17]*dv1*gamma+0.24494897427831794*dHdv_surf_C[10]*gamma; 
  out[32] = 0.41250000000000003*H_L[32]*dv1*gamma+0.41250000000000003*H_C[32]*dv1*gamma+0.32475952641916456*H_L[21]*dv1*gamma-0.32475952641916456*H_C[21]*dv1*gamma+0.24494897427831797*dHdv_surf_C[11]*gamma; 
  out[33] = -(1.8875000000000002*H_L[33]*dv1*gamma)+4.5125*H_C[33]*dv1*gamma-2.6142637586900057*H_L[15]*dv1*gamma-4.066632513517788*H_C[15]*dv1*gamma-1.8168052317185797*H_L[6]*dv1*gamma+1.8168052317185797*H_C[6]*dv1*gamma+0.3162277660168381*dHdv_surf_C[4]*gamma; 
  out[34] = 0.41250000000000003*H_L[34]*dv1*gamma+0.41250000000000003*H_C[34]*dv1*gamma+0.32475952641916456*H_L[23]*dv1*gamma-0.32475952641916456*H_C[23]*dv1*gamma+0.24494897427831797*dHdv_surf_C[12]*gamma; 
  out[35] = 0.41250000000000003*H_L[35]*dv1*gamma+0.41250000000000003*H_C[35]*dv1*gamma+0.32475952641916456*H_L[25]*dv1*gamma-0.32475952641916456*H_C[25]*dv1*gamma+0.24494897427831797*dHdv_surf_C[13]*gamma; 
  out[36] = -(1.8875000000000002*H_L[36]*dv1*gamma)+4.5125*H_C[36]*dv1*gamma-2.6142637586900057*H_L[16]*dv1*gamma-4.066632513517788*H_C[16]*dv1*gamma-1.8168052317185797*H_L[8]*dv1*gamma+1.8168052317185797*H_C[8]*dv1*gamma+0.3162277660168381*dHdv_surf_C[5]*gamma; 
  out[37] = -(1.1691342951089922*H_L[44]*dv1*gamma)+1.6454482671904336*H_C[44]*dv1*gamma-0.8125*H_L[37]*dv1*gamma+0.8125*H_C[37]*dv1*gamma+0.1414213562373096*dHdv_surf_C[17]*gamma; 
  out[38] = -(1.8875000000000002*H_L[38]*dv1*gamma)+4.5125*H_C[38]*dv1*gamma-2.6142637586900057*H_L[18]*dv1*gamma-4.066632513517788*H_C[18]*dv1*gamma-1.8168052317185797*H_L[10]*dv1*gamma+1.8168052317185797*H_C[10]*dv1*gamma+0.3162277660168381*dHdv_surf_C[6]*gamma; 
  out[39] = -(1.1691342951089922*H_L[46]*dv1*gamma)+1.6454482671904336*H_C[46]*dv1*gamma-0.8125*H_L[39]*dv1*gamma+0.8125*H_C[39]*dv1*gamma+0.1414213562373096*dHdv_surf_C[18]*gamma; 
  out[40] = 0.41250000000000003*H_L[40]*dv1*gamma+0.41250000000000003*H_C[40]*dv1*gamma+0.32475952641916456*H_L[27]*dv1*gamma-0.32475952641916456*H_C[27]*dv1*gamma+0.24494897427831797*dHdv_surf_C[14]*gamma; 
  out[41] = 0.41250000000000003*H_L[41]*dv1*gamma+0.41250000000000003*H_C[41]*dv1*gamma+0.32475952641916456*H_L[28]*dv1*gamma-0.32475952641916456*H_C[28]*dv1*gamma+0.24494897427831797*dHdv_surf_C[15]*gamma; 
  out[42] = -(1.1691342951089922*H_L[47]*dv1*gamma)+1.6454482671904336*H_C[47]*dv1*gamma-0.8125*H_L[42]*dv1*gamma+0.8125*H_C[42]*dv1*gamma+0.1414213562373096*dHdv_surf_C[19]*gamma; 
  out[43] = 0.41250000000000003*H_L[43]*dv1*gamma+0.41250000000000003*H_C[43]*dv1*gamma+0.32475952641916456*H_L[30]*dv1*gamma-0.32475952641916456*H_C[30]*dv1*gamma+0.24494897427831797*dHdv_surf_C[16]*gamma; 
  out[44] = 0.41250000000000003*H_L[44]*dv1*gamma+0.41250000000000003*H_C[44]*dv1*gamma+0.32475952641916456*H_L[37]*dv1*gamma-0.32475952641916456*H_C[37]*dv1*gamma+0.24494897427831797*dHdv_surf_C[17]*gamma; 
  out[45] = -(1.8875*H_L[45]*dv1*gamma)+4.5125*H_C[45]*dv1*gamma-2.6142637586900066*H_L[31]*dv1*gamma-4.066632513517788*H_C[31]*dv1*gamma-1.816805231718579*H_L[17]*dv1*gamma+1.816805231718579*H_C[17]*dv1*gamma+0.3162277660168381*dHdv_surf_C[10]*gamma; 
  out[46] = 0.41250000000000003*H_L[46]*dv1*gamma+0.41250000000000003*H_C[46]*dv1*gamma+0.32475952641916456*H_L[39]*dv1*gamma-0.32475952641916456*H_C[39]*dv1*gamma+0.24494897427831797*dHdv_surf_C[18]*gamma; 
  out[47] = 0.41250000000000003*H_L[47]*dv1*gamma+0.41250000000000003*H_C[47]*dv1*gamma+0.32475952641916456*H_L[42]*dv1*gamma-0.32475952641916456*H_C[42]*dv1*gamma+0.24494897427831797*dHdv_surf_C[19]*gamma; 
} 

