#include <gkyl_fpo_vlasov_kernels.h> 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vx_2x3v_ser_p1(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing. 
  // G_skin/edge:   Input potential skin/edge cells in recovery direction.

  const double dv1 = 2.0/dxv[2]; 
  const double dv1_sq = dv1*dv1; 
  double *diff_coeff_xx = &diff_coeff[0]; 
  double *diff_coeff_yy = &diff_coeff[32]; 
  double *diff_coeff_zz = &diff_coeff[64]; 
  if (edge == 1) {
  diff_coeff_xx[0] = (-0.6495190528383289*G_skin[3]*dv1_sq)-0.6495190528383289*G_edge[3]*dv1_sq-0.375*G_skin[0]*dv1_sq+0.375*G_edge[0]*dv1_sq; 
  diff_coeff_xx[1] = (-0.6495190528383289*G_skin[7]*dv1_sq)-0.6495190528383289*G_edge[7]*dv1_sq-0.375*G_skin[1]*dv1_sq+0.375*G_edge[1]*dv1_sq; 
  diff_coeff_xx[2] = (-0.6495190528383289*G_skin[8]*dv1_sq)-0.6495190528383289*G_edge[8]*dv1_sq-0.375*G_skin[2]*dv1_sq+0.375*G_edge[2]*dv1_sq; 
  diff_coeff_xx[3] = (-2.5*G_skin[3]*dv1_sq)-5.0*G_edge[3]*dv1_sq-2.165063509461096*G_skin[0]*dv1_sq+2.165063509461096*G_edge[0]*dv1_sq; 
  diff_coeff_xx[4] = (-0.6495190528383289*G_skin[11]*dv1_sq)-0.6495190528383289*G_edge[11]*dv1_sq-0.375*G_skin[4]*dv1_sq+0.375*G_edge[4]*dv1_sq; 
  diff_coeff_xx[5] = (-0.6495190528383289*G_skin[14]*dv1_sq)-0.6495190528383289*G_edge[14]*dv1_sq-0.375*G_skin[5]*dv1_sq+0.375*G_edge[5]*dv1_sq; 
  diff_coeff_xx[6] = (-0.6495190528383289*G_skin[16]*dv1_sq)-0.6495190528383289*G_edge[16]*dv1_sq-0.375*G_skin[6]*dv1_sq+0.375*G_edge[6]*dv1_sq; 
  diff_coeff_xx[7] = (-2.5*G_skin[7]*dv1_sq)-5.0*G_edge[7]*dv1_sq-2.165063509461096*G_skin[1]*dv1_sq+2.165063509461096*G_edge[1]*dv1_sq; 
  diff_coeff_xx[8] = (-2.5*G_skin[8]*dv1_sq)-5.0*G_edge[8]*dv1_sq-2.165063509461096*G_skin[2]*dv1_sq+2.165063509461096*G_edge[2]*dv1_sq; 
  diff_coeff_xx[9] = (-0.6495190528383289*G_skin[18]*dv1_sq)-0.6495190528383289*G_edge[18]*dv1_sq-0.375*G_skin[9]*dv1_sq+0.375*G_edge[9]*dv1_sq; 
  diff_coeff_xx[10] = (-0.6495190528383289*G_skin[19]*dv1_sq)-0.6495190528383289*G_edge[19]*dv1_sq-0.375*G_skin[10]*dv1_sq+0.375*G_edge[10]*dv1_sq; 
  diff_coeff_xx[11] = (-2.5*G_skin[11]*dv1_sq)-5.0*G_edge[11]*dv1_sq-2.165063509461096*G_skin[4]*dv1_sq+2.165063509461096*G_edge[4]*dv1_sq; 
  diff_coeff_xx[12] = (-0.6495190528383289*G_skin[21]*dv1_sq)-0.6495190528383289*G_edge[21]*dv1_sq-0.375*G_skin[12]*dv1_sq+0.375*G_edge[12]*dv1_sq; 
  diff_coeff_xx[13] = (-0.6495190528383289*G_skin[22]*dv1_sq)-0.6495190528383289*G_edge[22]*dv1_sq-0.375*G_skin[13]*dv1_sq+0.375*G_edge[13]*dv1_sq; 
  diff_coeff_xx[14] = (-2.5*G_skin[14]*dv1_sq)-5.0*G_edge[14]*dv1_sq-2.165063509461096*G_skin[5]*dv1_sq+2.165063509461096*G_edge[5]*dv1_sq; 
  diff_coeff_xx[15] = (-0.6495190528383289*G_skin[25]*dv1_sq)-0.6495190528383289*G_edge[25]*dv1_sq-0.375*G_skin[15]*dv1_sq+0.375*G_edge[15]*dv1_sq; 
  diff_coeff_xx[16] = (-2.5*G_skin[16]*dv1_sq)-5.0*G_edge[16]*dv1_sq-2.165063509461096*G_skin[6]*dv1_sq+2.165063509461096*G_edge[6]*dv1_sq; 
  diff_coeff_xx[17] = (-0.6495190528383289*G_skin[26]*dv1_sq)-0.6495190528383289*G_edge[26]*dv1_sq-0.375*G_skin[17]*dv1_sq+0.375*G_edge[17]*dv1_sq; 
  diff_coeff_xx[18] = (-2.5*G_skin[18]*dv1_sq)-5.0*G_edge[18]*dv1_sq-2.165063509461096*G_skin[9]*dv1_sq+2.165063509461096*G_edge[9]*dv1_sq; 
  diff_coeff_xx[19] = (-2.5*G_skin[19]*dv1_sq)-5.0*G_edge[19]*dv1_sq-2.165063509461096*G_skin[10]*dv1_sq+2.165063509461096*G_edge[10]*dv1_sq; 
  diff_coeff_xx[20] = (-0.6495190528383289*G_skin[27]*dv1_sq)-0.6495190528383289*G_edge[27]*dv1_sq-0.375*G_skin[20]*dv1_sq+0.375*G_edge[20]*dv1_sq; 
  diff_coeff_xx[21] = (-2.5*G_skin[21]*dv1_sq)-5.0*G_edge[21]*dv1_sq-2.165063509461096*G_skin[12]*dv1_sq+2.165063509461096*G_edge[12]*dv1_sq; 
  diff_coeff_xx[22] = (-2.5*G_skin[22]*dv1_sq)-5.0*G_edge[22]*dv1_sq-2.165063509461096*G_skin[13]*dv1_sq+2.165063509461096*G_edge[13]*dv1_sq; 
  diff_coeff_xx[23] = (-0.6495190528383289*G_skin[29]*dv1_sq)-0.6495190528383289*G_edge[29]*dv1_sq-0.375*G_skin[23]*dv1_sq+0.375*G_edge[23]*dv1_sq; 
  diff_coeff_xx[24] = (-0.6495190528383289*G_skin[30]*dv1_sq)-0.6495190528383289*G_edge[30]*dv1_sq-0.375*G_skin[24]*dv1_sq+0.375*G_edge[24]*dv1_sq; 
  diff_coeff_xx[25] = (-2.5*G_skin[25]*dv1_sq)-5.0*G_edge[25]*dv1_sq-2.165063509461096*G_skin[15]*dv1_sq+2.165063509461096*G_edge[15]*dv1_sq; 
  diff_coeff_xx[26] = (-2.5*G_skin[26]*dv1_sq)-5.0*G_edge[26]*dv1_sq-2.165063509461096*G_skin[17]*dv1_sq+2.165063509461096*G_edge[17]*dv1_sq; 
  diff_coeff_xx[27] = (-2.5*G_skin[27]*dv1_sq)-5.0*G_edge[27]*dv1_sq-2.165063509461096*G_skin[20]*dv1_sq+2.165063509461096*G_edge[20]*dv1_sq; 
  diff_coeff_xx[28] = (-0.6495190528383289*G_skin[31]*dv1_sq)-0.6495190528383289*G_edge[31]*dv1_sq-0.375*G_skin[28]*dv1_sq+0.375*G_edge[28]*dv1_sq; 
  diff_coeff_xx[29] = (-2.5*G_skin[29]*dv1_sq)-5.0*G_edge[29]*dv1_sq-2.165063509461096*G_skin[23]*dv1_sq+2.165063509461096*G_edge[23]*dv1_sq; 
  diff_coeff_xx[30] = (-2.5*G_skin[30]*dv1_sq)-5.0*G_edge[30]*dv1_sq-2.165063509461096*G_skin[24]*dv1_sq+2.165063509461096*G_edge[24]*dv1_sq; 
  diff_coeff_xx[31] = (-2.5*G_skin[31]*dv1_sq)-5.0*G_edge[31]*dv1_sq-2.165063509461096*G_skin[28]*dv1_sq+2.165063509461096*G_edge[28]*dv1_sq; 
 } else {
  diff_coeff_xx[0] = 0.6495190528383289*G_skin[3]*dv1_sq+0.6495190528383289*G_edge[3]*dv1_sq-0.375*G_skin[0]*dv1_sq+0.375*G_edge[0]*dv1_sq; 
  diff_coeff_xx[1] = 0.6495190528383289*G_skin[7]*dv1_sq+0.6495190528383289*G_edge[7]*dv1_sq-0.375*G_skin[1]*dv1_sq+0.375*G_edge[1]*dv1_sq; 
  diff_coeff_xx[2] = 0.6495190528383289*G_skin[8]*dv1_sq+0.6495190528383289*G_edge[8]*dv1_sq-0.375*G_skin[2]*dv1_sq+0.375*G_edge[2]*dv1_sq; 
  diff_coeff_xx[3] = (-2.5*G_skin[3]*dv1_sq)-5.0*G_edge[3]*dv1_sq+2.165063509461096*G_skin[0]*dv1_sq-2.165063509461096*G_edge[0]*dv1_sq; 
  diff_coeff_xx[4] = 0.6495190528383289*G_skin[11]*dv1_sq+0.6495190528383289*G_edge[11]*dv1_sq-0.375*G_skin[4]*dv1_sq+0.375*G_edge[4]*dv1_sq; 
  diff_coeff_xx[5] = 0.6495190528383289*G_skin[14]*dv1_sq+0.6495190528383289*G_edge[14]*dv1_sq-0.375*G_skin[5]*dv1_sq+0.375*G_edge[5]*dv1_sq; 
  diff_coeff_xx[6] = 0.6495190528383289*G_skin[16]*dv1_sq+0.6495190528383289*G_edge[16]*dv1_sq-0.375*G_skin[6]*dv1_sq+0.375*G_edge[6]*dv1_sq; 
  diff_coeff_xx[7] = (-2.5*G_skin[7]*dv1_sq)-5.0*G_edge[7]*dv1_sq+2.165063509461096*G_skin[1]*dv1_sq-2.165063509461096*G_edge[1]*dv1_sq; 
  diff_coeff_xx[8] = (-2.5*G_skin[8]*dv1_sq)-5.0*G_edge[8]*dv1_sq+2.165063509461096*G_skin[2]*dv1_sq-2.165063509461096*G_edge[2]*dv1_sq; 
  diff_coeff_xx[9] = 0.6495190528383289*G_skin[18]*dv1_sq+0.6495190528383289*G_edge[18]*dv1_sq-0.375*G_skin[9]*dv1_sq+0.375*G_edge[9]*dv1_sq; 
  diff_coeff_xx[10] = 0.6495190528383289*G_skin[19]*dv1_sq+0.6495190528383289*G_edge[19]*dv1_sq-0.375*G_skin[10]*dv1_sq+0.375*G_edge[10]*dv1_sq; 
  diff_coeff_xx[11] = (-2.5*G_skin[11]*dv1_sq)-5.0*G_edge[11]*dv1_sq+2.165063509461096*G_skin[4]*dv1_sq-2.165063509461096*G_edge[4]*dv1_sq; 
  diff_coeff_xx[12] = 0.6495190528383289*G_skin[21]*dv1_sq+0.6495190528383289*G_edge[21]*dv1_sq-0.375*G_skin[12]*dv1_sq+0.375*G_edge[12]*dv1_sq; 
  diff_coeff_xx[13] = 0.6495190528383289*G_skin[22]*dv1_sq+0.6495190528383289*G_edge[22]*dv1_sq-0.375*G_skin[13]*dv1_sq+0.375*G_edge[13]*dv1_sq; 
  diff_coeff_xx[14] = (-2.5*G_skin[14]*dv1_sq)-5.0*G_edge[14]*dv1_sq+2.165063509461096*G_skin[5]*dv1_sq-2.165063509461096*G_edge[5]*dv1_sq; 
  diff_coeff_xx[15] = 0.6495190528383289*G_skin[25]*dv1_sq+0.6495190528383289*G_edge[25]*dv1_sq-0.375*G_skin[15]*dv1_sq+0.375*G_edge[15]*dv1_sq; 
  diff_coeff_xx[16] = (-2.5*G_skin[16]*dv1_sq)-5.0*G_edge[16]*dv1_sq+2.165063509461096*G_skin[6]*dv1_sq-2.165063509461096*G_edge[6]*dv1_sq; 
  diff_coeff_xx[17] = 0.6495190528383289*G_skin[26]*dv1_sq+0.6495190528383289*G_edge[26]*dv1_sq-0.375*G_skin[17]*dv1_sq+0.375*G_edge[17]*dv1_sq; 
  diff_coeff_xx[18] = (-2.5*G_skin[18]*dv1_sq)-5.0*G_edge[18]*dv1_sq+2.165063509461096*G_skin[9]*dv1_sq-2.165063509461096*G_edge[9]*dv1_sq; 
  diff_coeff_xx[19] = (-2.5*G_skin[19]*dv1_sq)-5.0*G_edge[19]*dv1_sq+2.165063509461096*G_skin[10]*dv1_sq-2.165063509461096*G_edge[10]*dv1_sq; 
  diff_coeff_xx[20] = 0.6495190528383289*G_skin[27]*dv1_sq+0.6495190528383289*G_edge[27]*dv1_sq-0.375*G_skin[20]*dv1_sq+0.375*G_edge[20]*dv1_sq; 
  diff_coeff_xx[21] = (-2.5*G_skin[21]*dv1_sq)-5.0*G_edge[21]*dv1_sq+2.165063509461096*G_skin[12]*dv1_sq-2.165063509461096*G_edge[12]*dv1_sq; 
  diff_coeff_xx[22] = (-2.5*G_skin[22]*dv1_sq)-5.0*G_edge[22]*dv1_sq+2.165063509461096*G_skin[13]*dv1_sq-2.165063509461096*G_edge[13]*dv1_sq; 
  diff_coeff_xx[23] = 0.6495190528383289*G_skin[29]*dv1_sq+0.6495190528383289*G_edge[29]*dv1_sq-0.375*G_skin[23]*dv1_sq+0.375*G_edge[23]*dv1_sq; 
  diff_coeff_xx[24] = 0.6495190528383289*G_skin[30]*dv1_sq+0.6495190528383289*G_edge[30]*dv1_sq-0.375*G_skin[24]*dv1_sq+0.375*G_edge[24]*dv1_sq; 
  diff_coeff_xx[25] = (-2.5*G_skin[25]*dv1_sq)-5.0*G_edge[25]*dv1_sq+2.165063509461096*G_skin[15]*dv1_sq-2.165063509461096*G_edge[15]*dv1_sq; 
  diff_coeff_xx[26] = (-2.5*G_skin[26]*dv1_sq)-5.0*G_edge[26]*dv1_sq+2.165063509461096*G_skin[17]*dv1_sq-2.165063509461096*G_edge[17]*dv1_sq; 
  diff_coeff_xx[27] = (-2.5*G_skin[27]*dv1_sq)-5.0*G_edge[27]*dv1_sq+2.165063509461096*G_skin[20]*dv1_sq-2.165063509461096*G_edge[20]*dv1_sq; 
  diff_coeff_xx[28] = 0.6495190528383289*G_skin[31]*dv1_sq+0.6495190528383289*G_edge[31]*dv1_sq-0.375*G_skin[28]*dv1_sq+0.375*G_edge[28]*dv1_sq; 
  diff_coeff_xx[29] = (-2.5*G_skin[29]*dv1_sq)-5.0*G_edge[29]*dv1_sq+2.165063509461096*G_skin[23]*dv1_sq-2.165063509461096*G_edge[23]*dv1_sq; 
  diff_coeff_xx[30] = (-2.5*G_skin[30]*dv1_sq)-5.0*G_edge[30]*dv1_sq+2.165063509461096*G_skin[24]*dv1_sq-2.165063509461096*G_edge[24]*dv1_sq; 
  diff_coeff_xx[31] = (-2.5*G_skin[31]*dv1_sq)-5.0*G_edge[31]*dv1_sq+2.165063509461096*G_skin[28]*dv1_sq-2.165063509461096*G_edge[28]*dv1_sq; 
  } 
}
