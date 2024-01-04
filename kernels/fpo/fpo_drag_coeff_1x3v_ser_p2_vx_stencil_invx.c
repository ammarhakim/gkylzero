#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_drag_coeff_1x3v_vx_ser_p2_invx(const double *dxv, const double gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_h_stencil[3]: 3 cell stencil of Rosenbluth potential H. 
  // fpo_dhdv_surf: Surface projection of dH/dv in center cell. 
  // drag_coeff: Output array for drag coefficient. 

  double dv1 = 2.0/dxv[1]; 

  const double* H_L = fpo_h_stencil[0]; 
  const double* H_C = fpo_h_stencil[1]; 
  const double* H_R = fpo_h_stencil[2]; 
  
  const double *fpo_dHdv_surf_C_vx = &fpo_dhdv_surf[0]; 
  const double *fpo_dHdv_surf_C_vy = &fpo_dhdv_surf[20]; 
  const double *fpo_dHdv_surf_C_vz = &fpo_dhdv_surf[40]; 
  
  const double* dHdv_surf_C = fpo_dHdv_surf_C_vx; 
  
  double *drag_coeff_vx = &drag_coeff[0]; 
  double *drag_coeff_vy = &drag_coeff[48]; 
  double *drag_coeff_vz = &drag_coeff[96]; 
  
  double *out = drag_coeff_vx; 
  
  out[0] = 0.489139870078079*H_R[12]*dv1*gamma-0.489139870078079*H_L[12]*dv1*gamma-0.7036456405748563*H_R[2]*dv1*gamma-0.7036456405748563*H_L[2]*dv1*gamma+1.4072912811497127*H_C[2]*dv1*gamma+0.5*H_R[0]*dv1*gamma-0.5*H_L[0]*dv1*gamma; 
  out[1] = 0.4891398700780789*H_R[20]*dv1*gamma-0.4891398700780789*H_L[20]*dv1*gamma-0.7036456405748562*H_R[5]*dv1*gamma-0.7036456405748562*H_L[5]*dv1*gamma+1.4072912811497125*H_C[5]*dv1*gamma+0.5*H_R[1]*dv1*gamma-0.5*H_L[1]*dv1*gamma; 
  out[2] = 0.8472151069828725*H_R[12]*dv1*gamma+0.8472151069828725*H_L[12]*dv1*gamma+1.694430213965745*H_C[12]*dv1*gamma-1.21875*H_R[2]*dv1*gamma+1.21875*H_L[2]*dv1*gamma+0.8660254037844386*H_R[0]*dv1*gamma+0.8660254037844386*H_L[0]*dv1*gamma-1.7320508075688772*H_C[0]*dv1*gamma; 
  out[3] = 0.4891398700780789*H_R[22]*dv1*gamma-0.4891398700780789*H_L[22]*dv1*gamma-0.7036456405748562*H_R[7]*dv1*gamma-0.7036456405748562*H_L[7]*dv1*gamma+1.4072912811497125*H_C[7]*dv1*gamma+0.5*H_R[3]*dv1*gamma-0.5*H_L[3]*dv1*gamma; 
  out[4] = 0.4891398700780789*H_R[26]*dv1*gamma-0.4891398700780789*H_L[26]*dv1*gamma-0.7036456405748562*H_R[9]*dv1*gamma-0.7036456405748562*H_L[9]*dv1*gamma+1.4072912811497125*H_C[9]*dv1*gamma+0.5*H_R[4]*dv1*gamma-0.5*H_L[4]*dv1*gamma; 
  out[5] = 0.8472151069828725*H_R[20]*dv1*gamma+0.8472151069828725*H_L[20]*dv1*gamma+1.694430213965745*H_C[20]*dv1*gamma-1.21875*H_R[5]*dv1*gamma+1.21875*H_L[5]*dv1*gamma+0.8660254037844386*H_R[1]*dv1*gamma+0.8660254037844386*H_L[1]*dv1*gamma-1.7320508075688772*H_C[1]*dv1*gamma; 
  out[6] = 0.489139870078079*H_R[33]*dv1*gamma-0.489139870078079*H_L[33]*dv1*gamma-0.7036456405748563*H_R[15]*dv1*gamma-0.7036456405748563*H_L[15]*dv1*gamma+1.4072912811497127*H_C[15]*dv1*gamma+0.5*H_R[6]*dv1*gamma-0.5*H_L[6]*dv1*gamma; 
  out[7] = 0.8472151069828725*H_R[22]*dv1*gamma+0.8472151069828725*H_L[22]*dv1*gamma+1.694430213965745*H_C[22]*dv1*gamma-1.21875*H_R[7]*dv1*gamma+1.21875*H_L[7]*dv1*gamma+0.8660254037844386*H_R[3]*dv1*gamma+0.8660254037844386*H_L[3]*dv1*gamma-1.7320508075688772*H_C[3]*dv1*gamma; 
  out[8] = 0.489139870078079*H_R[36]*dv1*gamma-0.489139870078079*H_L[36]*dv1*gamma-0.7036456405748563*H_R[16]*dv1*gamma-0.7036456405748563*H_L[16]*dv1*gamma+1.4072912811497127*H_C[16]*dv1*gamma+0.5*H_R[8]*dv1*gamma-0.5*H_L[8]*dv1*gamma; 
  out[9] = 0.8472151069828725*H_R[26]*dv1*gamma+0.8472151069828725*H_L[26]*dv1*gamma+1.694430213965745*H_C[26]*dv1*gamma-1.21875*H_R[9]*dv1*gamma+1.21875*H_L[9]*dv1*gamma+0.8660254037844386*H_R[4]*dv1*gamma+0.8660254037844386*H_L[4]*dv1*gamma-1.7320508075688772*H_C[4]*dv1*gamma; 
  out[10] = 0.489139870078079*H_R[38]*dv1*gamma-0.489139870078079*H_L[38]*dv1*gamma-0.7036456405748563*H_R[18]*dv1*gamma-0.7036456405748563*H_L[18]*dv1*gamma+1.4072912811497127*H_C[18]*dv1*gamma+0.5*H_R[10]*dv1*gamma-0.5*H_L[10]*dv1*gamma; 
  out[11] = -(0.7036456405748565*H_R[19]*dv1*gamma)-0.7036456405748565*H_L[19]*dv1*gamma+1.407291281149713*H_C[19]*dv1*gamma+0.5*H_R[11]*dv1*gamma-0.5*H_L[11]*dv1*gamma; 
  out[12] = 1.09375*H_R[12]*dv1*gamma-1.09375*H_L[12]*dv1*gamma-1.5733994843967631*H_R[2]*dv1*gamma-1.5733994843967631*H_L[2]*dv1*gamma-4.599167723621307*H_C[2]*dv1*gamma+1.118033988749895*H_R[0]*dv1*gamma-1.118033988749895*H_L[0]*dv1*gamma; 
  out[13] = -(0.7036456405748565*H_R[24]*dv1*gamma)-0.7036456405748565*H_L[24]*dv1*gamma+1.407291281149713*H_C[24]*dv1*gamma+0.5*H_R[13]*dv1*gamma-0.5*H_L[13]*dv1*gamma; 
  out[14] = -(0.7036456405748565*H_R[29]*dv1*gamma)-0.7036456405748565*H_L[29]*dv1*gamma+1.407291281149713*H_C[29]*dv1*gamma+0.5*H_R[14]*dv1*gamma-0.5*H_L[14]*dv1*gamma; 
  out[15] = 0.8472151069828725*H_R[33]*dv1*gamma+0.8472151069828725*H_L[33]*dv1*gamma+1.694430213965745*H_C[33]*dv1*gamma-1.21875*H_R[15]*dv1*gamma+1.21875*H_L[15]*dv1*gamma+0.8660254037844386*H_R[6]*dv1*gamma+0.8660254037844386*H_L[6]*dv1*gamma-1.7320508075688772*H_C[6]*dv1*gamma; 
  out[16] = 0.8472151069828725*H_R[36]*dv1*gamma+0.8472151069828725*H_L[36]*dv1*gamma+1.694430213965745*H_C[36]*dv1*gamma-1.21875*H_R[16]*dv1*gamma+1.21875*H_L[16]*dv1*gamma+0.8660254037844386*H_R[8]*dv1*gamma+0.8660254037844386*H_L[8]*dv1*gamma-1.7320508075688772*H_C[8]*dv1*gamma; 
  out[17] = 0.4891398700780789*H_R[45]*dv1*gamma-0.4891398700780789*H_L[45]*dv1*gamma-0.7036456405748562*H_R[31]*dv1*gamma-0.7036456405748562*H_L[31]*dv1*gamma+1.4072912811497125*H_C[31]*dv1*gamma+0.5*H_R[17]*dv1*gamma-0.5*H_L[17]*dv1*gamma; 
  out[18] = 0.8472151069828725*H_R[38]*dv1*gamma+0.8472151069828725*H_L[38]*dv1*gamma+1.694430213965745*H_C[38]*dv1*gamma-1.21875*H_R[18]*dv1*gamma+1.21875*H_L[18]*dv1*gamma+0.8660254037844386*H_R[10]*dv1*gamma+0.8660254037844386*H_L[10]*dv1*gamma-1.7320508075688772*H_C[10]*dv1*gamma; 
  out[19] = -(1.21875*H_R[19]*dv1*gamma)+1.21875*H_L[19]*dv1*gamma+0.8660254037844387*H_R[11]*dv1*gamma+0.8660254037844387*H_L[11]*dv1*gamma-1.7320508075688774*H_C[11]*dv1*gamma; 
  out[20] = 1.09375*H_R[20]*dv1*gamma-1.09375*H_L[20]*dv1*gamma-1.5733994843967631*H_R[5]*dv1*gamma-1.5733994843967631*H_L[5]*dv1*gamma-4.599167723621307*H_C[5]*dv1*gamma+1.118033988749895*H_R[1]*dv1*gamma-1.118033988749895*H_L[1]*dv1*gamma; 
  out[21] = -(0.7036456405748565*H_R[32]*dv1*gamma)-0.7036456405748565*H_L[32]*dv1*gamma+1.407291281149713*H_C[32]*dv1*gamma+0.5*H_R[21]*dv1*gamma-0.5*H_L[21]*dv1*gamma; 
  out[22] = 1.09375*H_R[22]*dv1*gamma-1.09375*H_L[22]*dv1*gamma-1.5733994843967631*H_R[7]*dv1*gamma-1.5733994843967631*H_L[7]*dv1*gamma-4.599167723621307*H_C[7]*dv1*gamma+1.118033988749895*H_R[3]*dv1*gamma-1.118033988749895*H_L[3]*dv1*gamma; 
  out[23] = -(0.7036456405748565*H_R[34]*dv1*gamma)-0.7036456405748565*H_L[34]*dv1*gamma+1.407291281149713*H_C[34]*dv1*gamma+0.5*H_R[23]*dv1*gamma-0.5*H_L[23]*dv1*gamma; 
  out[24] = -(1.21875*H_R[24]*dv1*gamma)+1.21875*H_L[24]*dv1*gamma+0.8660254037844387*H_R[13]*dv1*gamma+0.8660254037844387*H_L[13]*dv1*gamma-1.7320508075688774*H_C[13]*dv1*gamma; 
  out[25] = -(0.7036456405748565*H_R[35]*dv1*gamma)-0.7036456405748565*H_L[35]*dv1*gamma+1.407291281149713*H_C[35]*dv1*gamma+0.5*H_R[25]*dv1*gamma-0.5*H_L[25]*dv1*gamma; 
  out[26] = 1.09375*H_R[26]*dv1*gamma-1.09375*H_L[26]*dv1*gamma-1.5733994843967631*H_R[9]*dv1*gamma-1.5733994843967631*H_L[9]*dv1*gamma-4.599167723621307*H_C[9]*dv1*gamma+1.118033988749895*H_R[4]*dv1*gamma-1.118033988749895*H_L[4]*dv1*gamma; 
  out[27] = -(0.7036456405748565*H_R[40]*dv1*gamma)-0.7036456405748565*H_L[40]*dv1*gamma+1.407291281149713*H_C[40]*dv1*gamma+0.5*H_R[27]*dv1*gamma-0.5*H_L[27]*dv1*gamma; 
  out[28] = -(0.7036456405748565*H_R[41]*dv1*gamma)-0.7036456405748565*H_L[41]*dv1*gamma+1.407291281149713*H_C[41]*dv1*gamma+0.5*H_R[28]*dv1*gamma-0.5*H_L[28]*dv1*gamma; 
  out[29] = -(1.21875*H_R[29]*dv1*gamma)+1.21875*H_L[29]*dv1*gamma+0.8660254037844387*H_R[14]*dv1*gamma+0.8660254037844387*H_L[14]*dv1*gamma-1.7320508075688774*H_C[14]*dv1*gamma; 
  out[30] = -(0.7036456405748565*H_R[43]*dv1*gamma)-0.7036456405748565*H_L[43]*dv1*gamma+1.407291281149713*H_C[43]*dv1*gamma+0.5*H_R[30]*dv1*gamma-0.5*H_L[30]*dv1*gamma; 
  out[31] = 0.8472151069828725*H_R[45]*dv1*gamma+0.8472151069828725*H_L[45]*dv1*gamma+1.694430213965745*H_C[45]*dv1*gamma-1.21875*H_R[31]*dv1*gamma+1.21875*H_L[31]*dv1*gamma+0.8660254037844386*H_R[17]*dv1*gamma+0.8660254037844386*H_L[17]*dv1*gamma-1.7320508075688772*H_C[17]*dv1*gamma; 
  out[32] = -(1.21875*H_R[32]*dv1*gamma)+1.21875*H_L[32]*dv1*gamma+0.8660254037844387*H_R[21]*dv1*gamma+0.8660254037844387*H_L[21]*dv1*gamma-1.7320508075688774*H_C[21]*dv1*gamma; 
  out[33] = 1.09375*H_R[33]*dv1*gamma-1.09375*H_L[33]*dv1*gamma-1.5733994843967631*H_R[15]*dv1*gamma-1.5733994843967631*H_L[15]*dv1*gamma-4.599167723621307*H_C[15]*dv1*gamma+1.118033988749895*H_R[6]*dv1*gamma-1.118033988749895*H_L[6]*dv1*gamma; 
  out[34] = -(1.21875*H_R[34]*dv1*gamma)+1.21875*H_L[34]*dv1*gamma+0.8660254037844387*H_R[23]*dv1*gamma+0.8660254037844387*H_L[23]*dv1*gamma-1.7320508075688774*H_C[23]*dv1*gamma; 
  out[35] = -(1.21875*H_R[35]*dv1*gamma)+1.21875*H_L[35]*dv1*gamma+0.8660254037844387*H_R[25]*dv1*gamma+0.8660254037844387*H_L[25]*dv1*gamma-1.7320508075688774*H_C[25]*dv1*gamma; 
  out[36] = 1.09375*H_R[36]*dv1*gamma-1.09375*H_L[36]*dv1*gamma-1.5733994843967631*H_R[16]*dv1*gamma-1.5733994843967631*H_L[16]*dv1*gamma-4.599167723621307*H_C[16]*dv1*gamma+1.118033988749895*H_R[8]*dv1*gamma-1.118033988749895*H_L[8]*dv1*gamma; 
  out[37] = -(0.7036456405748565*H_R[44]*dv1*gamma)-0.7036456405748565*H_L[44]*dv1*gamma+1.407291281149713*H_C[44]*dv1*gamma+0.5*H_R[37]*dv1*gamma-0.5*H_L[37]*dv1*gamma; 
  out[38] = 1.09375*H_R[38]*dv1*gamma-1.09375*H_L[38]*dv1*gamma-1.5733994843967631*H_R[18]*dv1*gamma-1.5733994843967631*H_L[18]*dv1*gamma-4.599167723621307*H_C[18]*dv1*gamma+1.118033988749895*H_R[10]*dv1*gamma-1.118033988749895*H_L[10]*dv1*gamma; 
  out[39] = -(0.7036456405748565*H_R[46]*dv1*gamma)-0.7036456405748565*H_L[46]*dv1*gamma+1.407291281149713*H_C[46]*dv1*gamma+0.5*H_R[39]*dv1*gamma-0.5*H_L[39]*dv1*gamma; 
  out[40] = -(1.21875*H_R[40]*dv1*gamma)+1.21875*H_L[40]*dv1*gamma+0.8660254037844387*H_R[27]*dv1*gamma+0.8660254037844387*H_L[27]*dv1*gamma-1.7320508075688774*H_C[27]*dv1*gamma; 
  out[41] = -(1.21875*H_R[41]*dv1*gamma)+1.21875*H_L[41]*dv1*gamma+0.8660254037844387*H_R[28]*dv1*gamma+0.8660254037844387*H_L[28]*dv1*gamma-1.7320508075688774*H_C[28]*dv1*gamma; 
  out[42] = -(0.7036456405748565*H_R[47]*dv1*gamma)-0.7036456405748565*H_L[47]*dv1*gamma+1.407291281149713*H_C[47]*dv1*gamma+0.5*H_R[42]*dv1*gamma-0.5*H_L[42]*dv1*gamma; 
  out[43] = -(1.21875*H_R[43]*dv1*gamma)+1.21875*H_L[43]*dv1*gamma+0.8660254037844387*H_R[30]*dv1*gamma+0.8660254037844387*H_L[30]*dv1*gamma-1.7320508075688774*H_C[30]*dv1*gamma; 
  out[44] = -(1.21875*H_R[44]*dv1*gamma)+1.21875*H_L[44]*dv1*gamma+0.8660254037844387*H_R[37]*dv1*gamma+0.8660254037844387*H_L[37]*dv1*gamma-1.7320508075688774*H_C[37]*dv1*gamma; 
  out[45] = 1.09375*H_R[45]*dv1*gamma-1.09375*H_L[45]*dv1*gamma-1.5733994843967631*H_R[31]*dv1*gamma-1.5733994843967631*H_L[31]*dv1*gamma-4.599167723621307*H_C[31]*dv1*gamma+1.118033988749895*H_R[17]*dv1*gamma-1.118033988749895*H_L[17]*dv1*gamma; 
  out[46] = -(1.21875*H_R[46]*dv1*gamma)+1.21875*H_L[46]*dv1*gamma+0.8660254037844387*H_R[39]*dv1*gamma+0.8660254037844387*H_L[39]*dv1*gamma-1.7320508075688774*H_C[39]*dv1*gamma; 
  out[47] = -(1.21875*H_R[47]*dv1*gamma)+1.21875*H_L[47]*dv1*gamma+0.8660254037844387*H_R[42]*dv1*gamma+0.8660254037844387*H_L[42]*dv1*gamma-1.7320508075688774*H_C[42]*dv1*gamma; 
} 

