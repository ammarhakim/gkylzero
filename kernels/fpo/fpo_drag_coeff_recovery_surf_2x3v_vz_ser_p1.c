#include <gkyl_fpo_vlasov_kernels.h> 
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vz_2x3v_ser_p1(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff) {
  // dxv[NDIM]: Cell spacing. 
  // H_skin/edge:   Input potential in skin/edge cells in recovery direction.
  
  const double dv1 = 2.0/dxv[4]; 
  double *drag_coeff_x = &drag_coeff[0]; 
  double *drag_coeff_y = &drag_coeff[32]; 
  double *drag_coeff_z = &drag_coeff[64]; 
  if (edge == 1) {
  drag_coeff_z[0] = (-0.2886751345948129*H_skin[5]*dv1)+1.154700538379252*H_edge[5]*dv1-0.25*H_skin[0]*dv1+0.25*H_edge[0]*dv1; 
  drag_coeff_z[1] = (-0.2886751345948129*H_skin[12]*dv1)+1.154700538379252*H_edge[12]*dv1-0.25*H_skin[1]*dv1+0.25*H_edge[1]*dv1; 
  drag_coeff_z[2] = (-0.2886751345948129*H_skin[13]*dv1)+1.154700538379252*H_edge[13]*dv1-0.25*H_skin[2]*dv1+0.25*H_edge[2]*dv1; 
  drag_coeff_z[3] = (-0.2886751345948129*H_skin[14]*dv1)+1.154700538379252*H_edge[14]*dv1-0.25*H_skin[3]*dv1+0.25*H_edge[3]*dv1; 
  drag_coeff_z[4] = (-0.2886751345948129*H_skin[15]*dv1)+1.154700538379252*H_edge[15]*dv1-0.25*H_skin[4]*dv1+0.25*H_edge[4]*dv1; 
  drag_coeff_z[5] = 0.5*H_skin[5]*dv1+H_edge[5]*dv1+0.4330127018922193*H_skin[0]*dv1-0.4330127018922193*H_edge[0]*dv1; 
  drag_coeff_z[6] = (-0.2886751345948129*H_skin[20]*dv1)+1.154700538379252*H_edge[20]*dv1-0.25*H_skin[6]*dv1+0.25*H_edge[6]*dv1; 
  drag_coeff_z[7] = (-0.2886751345948129*H_skin[21]*dv1)+1.154700538379252*H_edge[21]*dv1-0.25*H_skin[7]*dv1+0.25*H_edge[7]*dv1; 
  drag_coeff_z[8] = (-0.2886751345948129*H_skin[22]*dv1)+1.154700538379252*H_edge[22]*dv1-0.25*H_skin[8]*dv1+0.25*H_edge[8]*dv1; 
  drag_coeff_z[9] = (-0.2886751345948129*H_skin[23]*dv1)+1.154700538379252*H_edge[23]*dv1-0.25*H_skin[9]*dv1+0.25*H_edge[9]*dv1; 
  drag_coeff_z[10] = (-0.2886751345948129*H_skin[24]*dv1)+1.154700538379252*H_edge[24]*dv1-0.25*H_skin[10]*dv1+0.25*H_edge[10]*dv1; 
  drag_coeff_z[11] = (-0.2886751345948129*H_skin[25]*dv1)+1.154700538379252*H_edge[25]*dv1-0.25*H_skin[11]*dv1+0.25*H_edge[11]*dv1; 
  drag_coeff_z[12] = 0.5*H_skin[12]*dv1+H_edge[12]*dv1+0.4330127018922193*H_skin[1]*dv1-0.4330127018922193*H_edge[1]*dv1; 
  drag_coeff_z[13] = 0.5*H_skin[13]*dv1+H_edge[13]*dv1+0.4330127018922193*H_skin[2]*dv1-0.4330127018922193*H_edge[2]*dv1; 
  drag_coeff_z[14] = 0.5*H_skin[14]*dv1+H_edge[14]*dv1+0.4330127018922193*H_skin[3]*dv1-0.4330127018922193*H_edge[3]*dv1; 
  drag_coeff_z[15] = 0.5*H_skin[15]*dv1+H_edge[15]*dv1+0.4330127018922193*H_skin[4]*dv1-0.4330127018922193*H_edge[4]*dv1; 
  drag_coeff_z[16] = (-0.2886751345948129*H_skin[27]*dv1)+1.154700538379252*H_edge[27]*dv1-0.25*H_skin[16]*dv1+0.25*H_edge[16]*dv1; 
  drag_coeff_z[17] = (-0.2886751345948129*H_skin[28]*dv1)+1.154700538379252*H_edge[28]*dv1-0.25*H_skin[17]*dv1+0.25*H_edge[17]*dv1; 
  drag_coeff_z[18] = (-0.2886751345948129*H_skin[29]*dv1)+1.154700538379252*H_edge[29]*dv1-0.25*H_skin[18]*dv1+0.25*H_edge[18]*dv1; 
  drag_coeff_z[19] = (-0.2886751345948129*H_skin[30]*dv1)+1.154700538379252*H_edge[30]*dv1-0.25*H_skin[19]*dv1+0.25*H_edge[19]*dv1; 
  drag_coeff_z[20] = 0.5*H_skin[20]*dv1+H_edge[20]*dv1+0.4330127018922193*H_skin[6]*dv1-0.4330127018922193*H_edge[6]*dv1; 
  drag_coeff_z[21] = 0.5*H_skin[21]*dv1+H_edge[21]*dv1+0.4330127018922193*H_skin[7]*dv1-0.4330127018922193*H_edge[7]*dv1; 
  drag_coeff_z[22] = 0.5*H_skin[22]*dv1+H_edge[22]*dv1+0.4330127018922193*H_skin[8]*dv1-0.4330127018922193*H_edge[8]*dv1; 
  drag_coeff_z[23] = 0.5*H_skin[23]*dv1+H_edge[23]*dv1+0.4330127018922193*H_skin[9]*dv1-0.4330127018922193*H_edge[9]*dv1; 
  drag_coeff_z[24] = 0.5*H_skin[24]*dv1+H_edge[24]*dv1+0.4330127018922193*H_skin[10]*dv1-0.4330127018922193*H_edge[10]*dv1; 
  drag_coeff_z[25] = 0.5*H_skin[25]*dv1+H_edge[25]*dv1+0.4330127018922193*H_skin[11]*dv1-0.4330127018922193*H_edge[11]*dv1; 
  drag_coeff_z[26] = (-0.2886751345948129*H_skin[31]*dv1)+1.154700538379252*H_edge[31]*dv1-0.25*H_skin[26]*dv1+0.25*H_edge[26]*dv1; 
  drag_coeff_z[27] = 0.5*H_skin[27]*dv1+H_edge[27]*dv1+0.4330127018922193*H_skin[16]*dv1-0.4330127018922193*H_edge[16]*dv1; 
  drag_coeff_z[28] = 0.5*H_skin[28]*dv1+H_edge[28]*dv1+0.4330127018922193*H_skin[17]*dv1-0.4330127018922193*H_edge[17]*dv1; 
  drag_coeff_z[29] = 0.5*H_skin[29]*dv1+H_edge[29]*dv1+0.4330127018922193*H_skin[18]*dv1-0.4330127018922193*H_edge[18]*dv1; 
  drag_coeff_z[30] = 0.5*H_skin[30]*dv1+H_edge[30]*dv1+0.4330127018922193*H_skin[19]*dv1-0.4330127018922193*H_edge[19]*dv1; 
  drag_coeff_z[31] = 0.5*H_skin[31]*dv1+H_edge[31]*dv1+0.4330127018922193*H_skin[26]*dv1-0.4330127018922193*H_edge[26]*dv1; 
  } else {
  drag_coeff_z[0] = (-0.2886751345948129*H_skin[5]*dv1)+1.154700538379252*H_edge[5]*dv1+0.25*H_skin[0]*dv1-0.25*H_edge[0]*dv1; 
  drag_coeff_z[1] = (-0.2886751345948129*H_skin[12]*dv1)+1.154700538379252*H_edge[12]*dv1+0.25*H_skin[1]*dv1-0.25*H_edge[1]*dv1; 
  drag_coeff_z[2] = (-0.2886751345948129*H_skin[13]*dv1)+1.154700538379252*H_edge[13]*dv1+0.25*H_skin[2]*dv1-0.25*H_edge[2]*dv1; 
  drag_coeff_z[3] = (-0.2886751345948129*H_skin[14]*dv1)+1.154700538379252*H_edge[14]*dv1+0.25*H_skin[3]*dv1-0.25*H_edge[3]*dv1; 
  drag_coeff_z[4] = (-0.2886751345948129*H_skin[15]*dv1)+1.154700538379252*H_edge[15]*dv1+0.25*H_skin[4]*dv1-0.25*H_edge[4]*dv1; 
  drag_coeff_z[5] = (-0.5*H_skin[5]*dv1)-1.0*H_edge[5]*dv1+0.4330127018922193*H_skin[0]*dv1-0.4330127018922193*H_edge[0]*dv1; 
  drag_coeff_z[6] = (-0.2886751345948129*H_skin[20]*dv1)+1.154700538379252*H_edge[20]*dv1+0.25*H_skin[6]*dv1-0.25*H_edge[6]*dv1; 
  drag_coeff_z[7] = (-0.2886751345948129*H_skin[21]*dv1)+1.154700538379252*H_edge[21]*dv1+0.25*H_skin[7]*dv1-0.25*H_edge[7]*dv1; 
  drag_coeff_z[8] = (-0.2886751345948129*H_skin[22]*dv1)+1.154700538379252*H_edge[22]*dv1+0.25*H_skin[8]*dv1-0.25*H_edge[8]*dv1; 
  drag_coeff_z[9] = (-0.2886751345948129*H_skin[23]*dv1)+1.154700538379252*H_edge[23]*dv1+0.25*H_skin[9]*dv1-0.25*H_edge[9]*dv1; 
  drag_coeff_z[10] = (-0.2886751345948129*H_skin[24]*dv1)+1.154700538379252*H_edge[24]*dv1+0.25*H_skin[10]*dv1-0.25*H_edge[10]*dv1; 
  drag_coeff_z[11] = (-0.2886751345948129*H_skin[25]*dv1)+1.154700538379252*H_edge[25]*dv1+0.25*H_skin[11]*dv1-0.25*H_edge[11]*dv1; 
  drag_coeff_z[12] = (-0.5*H_skin[12]*dv1)-1.0*H_edge[12]*dv1+0.4330127018922193*H_skin[1]*dv1-0.4330127018922193*H_edge[1]*dv1; 
  drag_coeff_z[13] = (-0.5*H_skin[13]*dv1)-1.0*H_edge[13]*dv1+0.4330127018922193*H_skin[2]*dv1-0.4330127018922193*H_edge[2]*dv1; 
  drag_coeff_z[14] = (-0.5*H_skin[14]*dv1)-1.0*H_edge[14]*dv1+0.4330127018922193*H_skin[3]*dv1-0.4330127018922193*H_edge[3]*dv1; 
  drag_coeff_z[15] = (-0.5*H_skin[15]*dv1)-1.0*H_edge[15]*dv1+0.4330127018922193*H_skin[4]*dv1-0.4330127018922193*H_edge[4]*dv1; 
  drag_coeff_z[16] = (-0.2886751345948129*H_skin[27]*dv1)+1.154700538379252*H_edge[27]*dv1+0.25*H_skin[16]*dv1-0.25*H_edge[16]*dv1; 
  drag_coeff_z[17] = (-0.2886751345948129*H_skin[28]*dv1)+1.154700538379252*H_edge[28]*dv1+0.25*H_skin[17]*dv1-0.25*H_edge[17]*dv1; 
  drag_coeff_z[18] = (-0.2886751345948129*H_skin[29]*dv1)+1.154700538379252*H_edge[29]*dv1+0.25*H_skin[18]*dv1-0.25*H_edge[18]*dv1; 
  drag_coeff_z[19] = (-0.2886751345948129*H_skin[30]*dv1)+1.154700538379252*H_edge[30]*dv1+0.25*H_skin[19]*dv1-0.25*H_edge[19]*dv1; 
  drag_coeff_z[20] = (-0.5*H_skin[20]*dv1)-1.0*H_edge[20]*dv1+0.4330127018922193*H_skin[6]*dv1-0.4330127018922193*H_edge[6]*dv1; 
  drag_coeff_z[21] = (-0.5*H_skin[21]*dv1)-1.0*H_edge[21]*dv1+0.4330127018922193*H_skin[7]*dv1-0.4330127018922193*H_edge[7]*dv1; 
  drag_coeff_z[22] = (-0.5*H_skin[22]*dv1)-1.0*H_edge[22]*dv1+0.4330127018922193*H_skin[8]*dv1-0.4330127018922193*H_edge[8]*dv1; 
  drag_coeff_z[23] = (-0.5*H_skin[23]*dv1)-1.0*H_edge[23]*dv1+0.4330127018922193*H_skin[9]*dv1-0.4330127018922193*H_edge[9]*dv1; 
  drag_coeff_z[24] = (-0.5*H_skin[24]*dv1)-1.0*H_edge[24]*dv1+0.4330127018922193*H_skin[10]*dv1-0.4330127018922193*H_edge[10]*dv1; 
  drag_coeff_z[25] = (-0.5*H_skin[25]*dv1)-1.0*H_edge[25]*dv1+0.4330127018922193*H_skin[11]*dv1-0.4330127018922193*H_edge[11]*dv1; 
  drag_coeff_z[26] = (-0.2886751345948129*H_skin[31]*dv1)+1.154700538379252*H_edge[31]*dv1+0.25*H_skin[26]*dv1-0.25*H_edge[26]*dv1; 
  drag_coeff_z[27] = (-0.5*H_skin[27]*dv1)-1.0*H_edge[27]*dv1+0.4330127018922193*H_skin[16]*dv1-0.4330127018922193*H_edge[16]*dv1; 
  drag_coeff_z[28] = (-0.5*H_skin[28]*dv1)-1.0*H_edge[28]*dv1+0.4330127018922193*H_skin[17]*dv1-0.4330127018922193*H_edge[17]*dv1; 
  drag_coeff_z[29] = (-0.5*H_skin[29]*dv1)-1.0*H_edge[29]*dv1+0.4330127018922193*H_skin[18]*dv1-0.4330127018922193*H_edge[18]*dv1; 
  drag_coeff_z[30] = (-0.5*H_skin[30]*dv1)-1.0*H_edge[30]*dv1+0.4330127018922193*H_skin[19]*dv1-0.4330127018922193*H_edge[19]*dv1; 
  drag_coeff_z[31] = (-0.5*H_skin[31]*dv1)-1.0*H_edge[31]*dv1+0.4330127018922193*H_skin[26]*dv1-0.4330127018922193*H_edge[26]*dv1; 
  } 
}
