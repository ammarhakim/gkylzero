#include <gkyl_fpo_vlasov_kernels.h> 
GKYL_CU_DH void fpo_drag_coeff_recov_surf_vx_1x3v_ser_p1(const int edge, const double *dxv, const double *H_skin, const double *H_edge, double *drag_coeff) {
  // dxv[NDIM]: Cell spacing. 
  // H_skin/edge:   Input potential in skin/edge cells in recovery direction.
  
  const double dv1 = 2.0/dxv[1]; 
  double *drag_coeff_x = &drag_coeff[0]; 
  double *drag_coeff_y = &drag_coeff[16]; 
  double *drag_coeff_z = &drag_coeff[32]; 
  if (edge == 1) {
  drag_coeff_x[0] = (-0.2886751345948129*H_skin[2]*dv1)+1.154700538379252*H_edge[2]*dv1-0.25*H_skin[0]*dv1+0.25*H_edge[0]*dv1; 
  drag_coeff_x[1] = (-0.2886751345948129*H_skin[5]*dv1)+1.154700538379252*H_edge[5]*dv1-0.25*H_skin[1]*dv1+0.25*H_edge[1]*dv1; 
  drag_coeff_x[2] = 0.5*H_skin[2]*dv1+H_edge[2]*dv1+0.4330127018922193*H_skin[0]*dv1-0.4330127018922193*H_edge[0]*dv1; 
  drag_coeff_x[3] = (-0.2886751345948129*H_skin[7]*dv1)+1.154700538379252*H_edge[7]*dv1-0.25*H_skin[3]*dv1+0.25*H_edge[3]*dv1; 
  drag_coeff_x[4] = (-0.2886751345948129*H_skin[9]*dv1)+1.154700538379252*H_edge[9]*dv1-0.25*H_skin[4]*dv1+0.25*H_edge[4]*dv1; 
  drag_coeff_x[5] = 0.5*H_skin[5]*dv1+H_edge[5]*dv1+0.4330127018922193*H_skin[1]*dv1-0.4330127018922193*H_edge[1]*dv1; 
  drag_coeff_x[6] = (-0.2886751345948129*H_skin[11]*dv1)+1.154700538379252*H_edge[11]*dv1-0.25*H_skin[6]*dv1+0.25*H_edge[6]*dv1; 
  drag_coeff_x[7] = 0.5*H_skin[7]*dv1+H_edge[7]*dv1+0.4330127018922193*H_skin[3]*dv1-0.4330127018922193*H_edge[3]*dv1; 
  drag_coeff_x[8] = (-0.2886751345948129*H_skin[12]*dv1)+1.154700538379252*H_edge[12]*dv1-0.25*H_skin[8]*dv1+0.25*H_edge[8]*dv1; 
  drag_coeff_x[9] = 0.5*H_skin[9]*dv1+H_edge[9]*dv1+0.4330127018922193*H_skin[4]*dv1-0.4330127018922193*H_edge[4]*dv1; 
  drag_coeff_x[10] = (-0.2886751345948129*H_skin[14]*dv1)+1.154700538379252*H_edge[14]*dv1-0.25*H_skin[10]*dv1+0.25*H_edge[10]*dv1; 
  drag_coeff_x[11] = 0.5*H_skin[11]*dv1+H_edge[11]*dv1+0.4330127018922193*H_skin[6]*dv1-0.4330127018922193*H_edge[6]*dv1; 
  drag_coeff_x[12] = 0.5*H_skin[12]*dv1+H_edge[12]*dv1+0.4330127018922193*H_skin[8]*dv1-0.4330127018922193*H_edge[8]*dv1; 
  drag_coeff_x[13] = (-0.2886751345948129*H_skin[15]*dv1)+1.154700538379252*H_edge[15]*dv1-0.25*H_skin[13]*dv1+0.25*H_edge[13]*dv1; 
  drag_coeff_x[14] = 0.5*H_skin[14]*dv1+H_edge[14]*dv1+0.4330127018922193*H_skin[10]*dv1-0.4330127018922193*H_edge[10]*dv1; 
  drag_coeff_x[15] = 0.5*H_skin[15]*dv1+H_edge[15]*dv1+0.4330127018922193*H_skin[13]*dv1-0.4330127018922193*H_edge[13]*dv1; 
  } else {
  drag_coeff_x[0] = (-0.2886751345948129*H_skin[2]*dv1)+1.154700538379252*H_edge[2]*dv1+0.25*H_skin[0]*dv1-0.25*H_edge[0]*dv1; 
  drag_coeff_x[1] = (-0.2886751345948129*H_skin[5]*dv1)+1.154700538379252*H_edge[5]*dv1+0.25*H_skin[1]*dv1-0.25*H_edge[1]*dv1; 
  drag_coeff_x[2] = (-0.5*H_skin[2]*dv1)-1.0*H_edge[2]*dv1+0.4330127018922193*H_skin[0]*dv1-0.4330127018922193*H_edge[0]*dv1; 
  drag_coeff_x[3] = (-0.2886751345948129*H_skin[7]*dv1)+1.154700538379252*H_edge[7]*dv1+0.25*H_skin[3]*dv1-0.25*H_edge[3]*dv1; 
  drag_coeff_x[4] = (-0.2886751345948129*H_skin[9]*dv1)+1.154700538379252*H_edge[9]*dv1+0.25*H_skin[4]*dv1-0.25*H_edge[4]*dv1; 
  drag_coeff_x[5] = (-0.5*H_skin[5]*dv1)-1.0*H_edge[5]*dv1+0.4330127018922193*H_skin[1]*dv1-0.4330127018922193*H_edge[1]*dv1; 
  drag_coeff_x[6] = (-0.2886751345948129*H_skin[11]*dv1)+1.154700538379252*H_edge[11]*dv1+0.25*H_skin[6]*dv1-0.25*H_edge[6]*dv1; 
  drag_coeff_x[7] = (-0.5*H_skin[7]*dv1)-1.0*H_edge[7]*dv1+0.4330127018922193*H_skin[3]*dv1-0.4330127018922193*H_edge[3]*dv1; 
  drag_coeff_x[8] = (-0.2886751345948129*H_skin[12]*dv1)+1.154700538379252*H_edge[12]*dv1+0.25*H_skin[8]*dv1-0.25*H_edge[8]*dv1; 
  drag_coeff_x[9] = (-0.5*H_skin[9]*dv1)-1.0*H_edge[9]*dv1+0.4330127018922193*H_skin[4]*dv1-0.4330127018922193*H_edge[4]*dv1; 
  drag_coeff_x[10] = (-0.2886751345948129*H_skin[14]*dv1)+1.154700538379252*H_edge[14]*dv1+0.25*H_skin[10]*dv1-0.25*H_edge[10]*dv1; 
  drag_coeff_x[11] = (-0.5*H_skin[11]*dv1)-1.0*H_edge[11]*dv1+0.4330127018922193*H_skin[6]*dv1-0.4330127018922193*H_edge[6]*dv1; 
  drag_coeff_x[12] = (-0.5*H_skin[12]*dv1)-1.0*H_edge[12]*dv1+0.4330127018922193*H_skin[8]*dv1-0.4330127018922193*H_edge[8]*dv1; 
  drag_coeff_x[13] = (-0.2886751345948129*H_skin[15]*dv1)+1.154700538379252*H_edge[15]*dv1+0.25*H_skin[13]*dv1-0.25*H_edge[13]*dv1; 
  drag_coeff_x[14] = (-0.5*H_skin[14]*dv1)-1.0*H_edge[14]*dv1+0.4330127018922193*H_skin[10]*dv1-0.4330127018922193*H_edge[10]*dv1; 
  drag_coeff_x[15] = (-0.5*H_skin[15]*dv1)-1.0*H_edge[15]*dv1+0.4330127018922193*H_skin[13]*dv1-0.4330127018922193*H_edge[13]*dv1; 
  } 
}
