#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_drag_coeff_1x3v_vz_ser_p2_lovz(const double *dxv, const double gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_h_stencil[3]: 3 cell stencil of Rosenbluth potential H. 
  // fpo_dhdv_surf: Surface projection of dH/dv in center cell. 
  // drag_coeff: Output array for drag coefficient. 

  double dv1 = 2.0/dxv[3]; 

  const double* H_C = fpo_h_stencil[0]; 
  const double* H_R = fpo_h_stencil[1]; 
  
  const double *fpo_dHdv_surf_C_vx = &fpo_dhdv_surf[0]; 
  const double *fpo_dHdv_surf_C_vy = &fpo_dhdv_surf[20]; 
  const double *fpo_dHdv_surf_C_vz = &fpo_dhdv_surf[40]; 
  
  const double* dHdv_surf_C = fpo_dHdv_surf_C_vz; 
  
  double *drag_coeff_vx = &drag_coeff[0]; 
  double *drag_coeff_vy = &drag_coeff[48]; 
  double *drag_coeff_vz = &drag_coeff[96]; 
  
  double *out = drag_coeff_vz; 
  
  out[0] = 0.8441156615061707*H_R[14]*dv1*gamma-2.01805134969356*H_C[14]*dv1*gamma-1.1691342951089918*H_R[4]*dv1*gamma+1.6454482671904334*H_C[4]*dv1*gamma+0.8125*H_R[0]*dv1*gamma-0.8125*H_C[0]*dv1*gamma+0.1414213562373096*dHdv_surf_C[0]*gamma; 
  out[1] = 0.8441156615061707*H_R[28]*dv1*gamma-2.0180513496935606*H_C[28]*dv1*gamma-1.1691342951089916*H_R[8]*dv1*gamma+1.6454482671904334*H_C[8]*dv1*gamma+0.8125*H_R[1]*dv1*gamma-0.8125*H_C[1]*dv1*gamma+0.1414213562373096*dHdv_surf_C[1]*gamma; 
  out[2] = 0.8441156615061707*H_R[29]*dv1*gamma-2.0180513496935606*H_C[29]*dv1*gamma-1.1691342951089916*H_R[9]*dv1*gamma+1.6454482671904334*H_C[9]*dv1*gamma+0.8125*H_R[2]*dv1*gamma-0.8125*H_C[2]*dv1*gamma+0.1414213562373096*dHdv_surf_C[2]*gamma; 
  out[3] = 0.8441156615061707*H_R[30]*dv1*gamma-2.0180513496935606*H_C[30]*dv1*gamma-1.1691342951089916*H_R[10]*dv1*gamma+1.6454482671904334*H_C[10]*dv1*gamma+0.8125*H_R[3]*dv1*gamma-0.8125*H_C[3]*dv1*gamma+0.1414213562373096*dHdv_surf_C[3]*gamma; 
  out[4] = 0.232379000772445*H_R[14]*dv1*gamma+5.189797683917939*H_C[14]*dv1*gamma-0.41250000000000003*H_R[4]*dv1*gamma-0.41250000000000003*H_C[4]*dv1*gamma+0.32475952641916445*H_R[0]*dv1*gamma-0.32475952641916445*H_C[0]*dv1*gamma-0.24494897427831794*dHdv_surf_C[0]*gamma; 
  out[5] = 0.8441156615061707*H_R[41]*dv1*gamma-2.01805134969356*H_C[41]*dv1*gamma-1.1691342951089918*H_R[16]*dv1*gamma+1.6454482671904334*H_C[16]*dv1*gamma+0.8125*H_R[5]*dv1*gamma-0.8125*H_C[5]*dv1*gamma+0.1414213562373096*dHdv_surf_C[4]*gamma; 
  out[6] = 0.8441156615061707*H_R[42]*dv1*gamma-2.01805134969356*H_C[42]*dv1*gamma-1.1691342951089918*H_R[17]*dv1*gamma+1.6454482671904334*H_C[17]*dv1*gamma+0.8125*H_R[6]*dv1*gamma-0.8125*H_C[6]*dv1*gamma+0.1414213562373096*dHdv_surf_C[5]*gamma; 
  out[7] = 0.8441156615061707*H_R[43]*dv1*gamma-2.01805134969356*H_C[43]*dv1*gamma-1.1691342951089918*H_R[18]*dv1*gamma+1.6454482671904334*H_C[18]*dv1*gamma+0.8125*H_R[7]*dv1*gamma-0.8125*H_C[7]*dv1*gamma+0.1414213562373096*dHdv_surf_C[6]*gamma; 
  out[8] = 0.232379000772445*H_R[28]*dv1*gamma+5.189797683917939*H_C[28]*dv1*gamma-0.41250000000000003*H_R[8]*dv1*gamma-0.41250000000000003*H_C[8]*dv1*gamma+0.32475952641916445*H_R[1]*dv1*gamma-0.32475952641916445*H_C[1]*dv1*gamma-0.24494897427831794*dHdv_surf_C[1]*gamma; 
  out[9] = 0.232379000772445*H_R[29]*dv1*gamma+5.189797683917939*H_C[29]*dv1*gamma-0.41250000000000003*H_R[9]*dv1*gamma-0.41250000000000003*H_C[9]*dv1*gamma+0.32475952641916445*H_R[2]*dv1*gamma-0.32475952641916445*H_C[2]*dv1*gamma-0.24494897427831794*dHdv_surf_C[2]*gamma; 
  out[10] = 0.232379000772445*H_R[30]*dv1*gamma+5.189797683917939*H_C[30]*dv1*gamma-0.41250000000000003*H_R[10]*dv1*gamma-0.41250000000000003*H_C[10]*dv1*gamma+0.32475952641916445*H_R[3]*dv1*gamma-0.32475952641916445*H_C[3]*dv1*gamma-0.24494897427831794*dHdv_surf_C[3]*gamma; 
  out[11] = -(1.1691342951089922*H_R[25]*dv1*gamma)+1.6454482671904336*H_C[25]*dv1*gamma+0.8125*H_R[11]*dv1*gamma-0.8125*H_C[11]*dv1*gamma+0.1414213562373096*dHdv_surf_C[7]*gamma; 
  out[12] = -(1.1691342951089922*H_R[26]*dv1*gamma)+1.6454482671904336*H_C[26]*dv1*gamma+0.8125*H_R[12]*dv1*gamma-0.8125*H_C[12]*dv1*gamma+0.1414213562373096*dHdv_surf_C[8]*gamma; 
  out[13] = -(1.1691342951089922*H_R[27]*dv1*gamma)+1.6454482671904336*H_C[27]*dv1*gamma+0.8125*H_R[13]*dv1*gamma-0.8125*H_C[13]*dv1*gamma+0.1414213562373096*dHdv_surf_C[9]*gamma; 
  out[14] = 1.8875000000000002*H_R[14]*dv1*gamma-4.5125*H_C[14]*dv1*gamma-2.6142637586900057*H_R[4]*dv1*gamma-4.066632513517788*H_C[4]*dv1*gamma+1.8168052317185797*H_R[0]*dv1*gamma-1.8168052317185797*H_C[0]*dv1*gamma+0.3162277660168381*dHdv_surf_C[0]*gamma; 
  out[15] = 0.8441156615061707*H_R[47]*dv1*gamma-2.0180513496935606*H_C[47]*dv1*gamma-1.1691342951089916*H_R[31]*dv1*gamma+1.6454482671904334*H_C[31]*dv1*gamma+0.8125*H_R[15]*dv1*gamma-0.8125*H_C[15]*dv1*gamma+0.1414213562373096*dHdv_surf_C[10]*gamma; 
  out[16] = 0.232379000772445*H_R[41]*dv1*gamma+5.189797683917939*H_C[41]*dv1*gamma-0.41250000000000003*H_R[16]*dv1*gamma-0.41250000000000003*H_C[16]*dv1*gamma+0.32475952641916445*H_R[5]*dv1*gamma-0.32475952641916445*H_C[5]*dv1*gamma-0.24494897427831794*dHdv_surf_C[4]*gamma; 
  out[17] = 0.232379000772445*H_R[42]*dv1*gamma+5.189797683917939*H_C[42]*dv1*gamma-0.41250000000000003*H_R[17]*dv1*gamma-0.41250000000000003*H_C[17]*dv1*gamma+0.32475952641916445*H_R[6]*dv1*gamma-0.32475952641916445*H_C[6]*dv1*gamma-0.24494897427831794*dHdv_surf_C[5]*gamma; 
  out[18] = 0.232379000772445*H_R[43]*dv1*gamma+5.189797683917939*H_C[43]*dv1*gamma-0.41250000000000003*H_R[18]*dv1*gamma-0.41250000000000003*H_C[18]*dv1*gamma+0.32475952641916445*H_R[7]*dv1*gamma-0.32475952641916445*H_C[7]*dv1*gamma-0.24494897427831794*dHdv_surf_C[6]*gamma; 
  out[19] = -(1.1691342951089922*H_R[35]*dv1*gamma)+1.6454482671904336*H_C[35]*dv1*gamma+0.8125*H_R[19]*dv1*gamma-0.8125*H_C[19]*dv1*gamma+0.1414213562373096*dHdv_surf_C[11]*gamma; 
  out[20] = -(1.1691342951089922*H_R[36]*dv1*gamma)+1.6454482671904336*H_C[36]*dv1*gamma+0.8125*H_R[20]*dv1*gamma-0.8125*H_C[20]*dv1*gamma+0.1414213562373096*dHdv_surf_C[12]*gamma; 
  out[21] = -(1.1691342951089922*H_R[37]*dv1*gamma)+1.6454482671904336*H_C[37]*dv1*gamma+0.8125*H_R[21]*dv1*gamma-0.8125*H_C[21]*dv1*gamma+0.1414213562373096*dHdv_surf_C[13]*gamma; 
  out[22] = -(1.1691342951089922*H_R[38]*dv1*gamma)+1.6454482671904336*H_C[38]*dv1*gamma+0.8125*H_R[22]*dv1*gamma-0.8125*H_C[22]*dv1*gamma+0.1414213562373096*dHdv_surf_C[14]*gamma; 
  out[23] = -(1.1691342951089922*H_R[39]*dv1*gamma)+1.6454482671904336*H_C[39]*dv1*gamma+0.8125*H_R[23]*dv1*gamma-0.8125*H_C[23]*dv1*gamma+0.1414213562373096*dHdv_surf_C[15]*gamma; 
  out[24] = -(1.1691342951089922*H_R[40]*dv1*gamma)+1.6454482671904336*H_C[40]*dv1*gamma+0.8125*H_R[24]*dv1*gamma-0.8125*H_C[24]*dv1*gamma+0.1414213562373096*dHdv_surf_C[16]*gamma; 
  out[25] = -(0.41250000000000003*H_R[25]*dv1*gamma)-0.41250000000000003*H_C[25]*dv1*gamma+0.32475952641916456*H_R[11]*dv1*gamma-0.32475952641916456*H_C[11]*dv1*gamma-0.24494897427831797*dHdv_surf_C[7]*gamma; 
  out[26] = -(0.41250000000000003*H_R[26]*dv1*gamma)-0.41250000000000003*H_C[26]*dv1*gamma+0.32475952641916456*H_R[12]*dv1*gamma-0.32475952641916456*H_C[12]*dv1*gamma-0.24494897427831797*dHdv_surf_C[8]*gamma; 
  out[27] = -(0.41250000000000003*H_R[27]*dv1*gamma)-0.41250000000000003*H_C[27]*dv1*gamma+0.32475952641916456*H_R[13]*dv1*gamma-0.32475952641916456*H_C[13]*dv1*gamma-0.24494897427831797*dHdv_surf_C[9]*gamma; 
  out[28] = 1.8875*H_R[28]*dv1*gamma-4.5125*H_C[28]*dv1*gamma-2.6142637586900066*H_R[8]*dv1*gamma-4.066632513517788*H_C[8]*dv1*gamma+1.816805231718579*H_R[1]*dv1*gamma-1.816805231718579*H_C[1]*dv1*gamma+0.3162277660168381*dHdv_surf_C[1]*gamma; 
  out[29] = 1.8875*H_R[29]*dv1*gamma-4.5125*H_C[29]*dv1*gamma-2.6142637586900066*H_R[9]*dv1*gamma-4.066632513517788*H_C[9]*dv1*gamma+1.816805231718579*H_R[2]*dv1*gamma-1.816805231718579*H_C[2]*dv1*gamma+0.3162277660168381*dHdv_surf_C[2]*gamma; 
  out[30] = 1.8875*H_R[30]*dv1*gamma-4.5125*H_C[30]*dv1*gamma-2.6142637586900066*H_R[10]*dv1*gamma-4.066632513517788*H_C[10]*dv1*gamma+1.816805231718579*H_R[3]*dv1*gamma-1.816805231718579*H_C[3]*dv1*gamma+0.3162277660168381*dHdv_surf_C[3]*gamma; 
  out[31] = 0.232379000772445*H_R[47]*dv1*gamma+5.189797683917939*H_C[47]*dv1*gamma-0.41250000000000003*H_R[31]*dv1*gamma-0.41250000000000003*H_C[31]*dv1*gamma+0.32475952641916445*H_R[15]*dv1*gamma-0.32475952641916445*H_C[15]*dv1*gamma-0.24494897427831794*dHdv_surf_C[10]*gamma; 
  out[32] = -(1.1691342951089922*H_R[44]*dv1*gamma)+1.6454482671904336*H_C[44]*dv1*gamma+0.8125*H_R[32]*dv1*gamma-0.8125*H_C[32]*dv1*gamma+0.1414213562373096*dHdv_surf_C[17]*gamma; 
  out[33] = -(1.1691342951089922*H_R[45]*dv1*gamma)+1.6454482671904336*H_C[45]*dv1*gamma+0.8125*H_R[33]*dv1*gamma-0.8125*H_C[33]*dv1*gamma+0.1414213562373096*dHdv_surf_C[18]*gamma; 
  out[34] = -(1.1691342951089922*H_R[46]*dv1*gamma)+1.6454482671904336*H_C[46]*dv1*gamma+0.8125*H_R[34]*dv1*gamma-0.8125*H_C[34]*dv1*gamma+0.1414213562373096*dHdv_surf_C[19]*gamma; 
  out[35] = -(0.41250000000000003*H_R[35]*dv1*gamma)-0.41250000000000003*H_C[35]*dv1*gamma+0.32475952641916456*H_R[19]*dv1*gamma-0.32475952641916456*H_C[19]*dv1*gamma-0.24494897427831797*dHdv_surf_C[11]*gamma; 
  out[36] = -(0.41250000000000003*H_R[36]*dv1*gamma)-0.41250000000000003*H_C[36]*dv1*gamma+0.32475952641916456*H_R[20]*dv1*gamma-0.32475952641916456*H_C[20]*dv1*gamma-0.24494897427831797*dHdv_surf_C[12]*gamma; 
  out[37] = -(0.41250000000000003*H_R[37]*dv1*gamma)-0.41250000000000003*H_C[37]*dv1*gamma+0.32475952641916456*H_R[21]*dv1*gamma-0.32475952641916456*H_C[21]*dv1*gamma-0.24494897427831797*dHdv_surf_C[13]*gamma; 
  out[38] = -(0.41250000000000003*H_R[38]*dv1*gamma)-0.41250000000000003*H_C[38]*dv1*gamma+0.32475952641916456*H_R[22]*dv1*gamma-0.32475952641916456*H_C[22]*dv1*gamma-0.24494897427831797*dHdv_surf_C[14]*gamma; 
  out[39] = -(0.41250000000000003*H_R[39]*dv1*gamma)-0.41250000000000003*H_C[39]*dv1*gamma+0.32475952641916456*H_R[23]*dv1*gamma-0.32475952641916456*H_C[23]*dv1*gamma-0.24494897427831797*dHdv_surf_C[15]*gamma; 
  out[40] = -(0.41250000000000003*H_R[40]*dv1*gamma)-0.41250000000000003*H_C[40]*dv1*gamma+0.32475952641916456*H_R[24]*dv1*gamma-0.32475952641916456*H_C[24]*dv1*gamma-0.24494897427831797*dHdv_surf_C[16]*gamma; 
  out[41] = 1.8875000000000002*H_R[41]*dv1*gamma-4.5125*H_C[41]*dv1*gamma-2.6142637586900057*H_R[16]*dv1*gamma-4.066632513517788*H_C[16]*dv1*gamma+1.8168052317185797*H_R[5]*dv1*gamma-1.8168052317185797*H_C[5]*dv1*gamma+0.3162277660168381*dHdv_surf_C[4]*gamma; 
  out[42] = 1.8875000000000002*H_R[42]*dv1*gamma-4.5125*H_C[42]*dv1*gamma-2.6142637586900057*H_R[17]*dv1*gamma-4.066632513517788*H_C[17]*dv1*gamma+1.8168052317185797*H_R[6]*dv1*gamma-1.8168052317185797*H_C[6]*dv1*gamma+0.3162277660168381*dHdv_surf_C[5]*gamma; 
  out[43] = 1.8875000000000002*H_R[43]*dv1*gamma-4.5125*H_C[43]*dv1*gamma-2.6142637586900057*H_R[18]*dv1*gamma-4.066632513517788*H_C[18]*dv1*gamma+1.8168052317185797*H_R[7]*dv1*gamma-1.8168052317185797*H_C[7]*dv1*gamma+0.3162277660168381*dHdv_surf_C[6]*gamma; 
  out[44] = -(0.41250000000000003*H_R[44]*dv1*gamma)-0.41250000000000003*H_C[44]*dv1*gamma+0.32475952641916456*H_R[32]*dv1*gamma-0.32475952641916456*H_C[32]*dv1*gamma-0.24494897427831797*dHdv_surf_C[17]*gamma; 
  out[45] = -(0.41250000000000003*H_R[45]*dv1*gamma)-0.41250000000000003*H_C[45]*dv1*gamma+0.32475952641916456*H_R[33]*dv1*gamma-0.32475952641916456*H_C[33]*dv1*gamma-0.24494897427831797*dHdv_surf_C[18]*gamma; 
  out[46] = -(0.41250000000000003*H_R[46]*dv1*gamma)-0.41250000000000003*H_C[46]*dv1*gamma+0.32475952641916456*H_R[34]*dv1*gamma-0.32475952641916456*H_C[34]*dv1*gamma-0.24494897427831797*dHdv_surf_C[19]*gamma; 
  out[47] = 1.8875*H_R[47]*dv1*gamma-4.5125*H_C[47]*dv1*gamma-2.6142637586900066*H_R[31]*dv1*gamma-4.066632513517788*H_C[31]*dv1*gamma+1.816805231718579*H_R[15]*dv1*gamma-1.816805231718579*H_C[15]*dv1*gamma+0.3162277660168381*dHdv_surf_C[10]*gamma; 
} 
