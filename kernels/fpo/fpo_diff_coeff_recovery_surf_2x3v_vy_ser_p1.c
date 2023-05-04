#include <gkyl_fpo_vlasov_kernels.h> 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vy_2x3v_ser_p1(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing. 
  // G_skin/edge:   Input potential skin/edge cells in recovery direction.

  const double dv1 = 2.0/dxv[3]; 
  const double dv1_sq = dv1*dv1; 
  double *diff_coeff_xx = &diff_coeff[0]; 
  double *diff_coeff_yy = &diff_coeff[32]; 
  double *diff_coeff_zz = &diff_coeff[64]; 
  if (edge == 1) {
  diff_coeff_yy[0] = (-0.6495190528383289*G_skin[4]*dv1_sq)-0.6495190528383289*G_edge[4]*dv1_sq-0.375*G_skin[0]*dv1_sq+0.375*G_edge[0]*dv1_sq; 
  diff_coeff_yy[1] = (-0.6495190528383289*G_skin[9]*dv1_sq)-0.6495190528383289*G_edge[9]*dv1_sq-0.375*G_skin[1]*dv1_sq+0.375*G_edge[1]*dv1_sq; 
  diff_coeff_yy[2] = (-0.6495190528383289*G_skin[10]*dv1_sq)-0.6495190528383289*G_edge[10]*dv1_sq-0.375*G_skin[2]*dv1_sq+0.375*G_edge[2]*dv1_sq; 
  diff_coeff_yy[3] = (-0.6495190528383289*G_skin[11]*dv1_sq)-0.6495190528383289*G_edge[11]*dv1_sq-0.375*G_skin[3]*dv1_sq+0.375*G_edge[3]*dv1_sq; 
  diff_coeff_yy[4] = (-2.5*G_skin[4]*dv1_sq)-5.0*G_edge[4]*dv1_sq-2.165063509461096*G_skin[0]*dv1_sq+2.165063509461096*G_edge[0]*dv1_sq; 
  diff_coeff_yy[5] = (-0.6495190528383289*G_skin[15]*dv1_sq)-0.6495190528383289*G_edge[15]*dv1_sq-0.375*G_skin[5]*dv1_sq+0.375*G_edge[5]*dv1_sq; 
  diff_coeff_yy[6] = (-0.6495190528383289*G_skin[17]*dv1_sq)-0.6495190528383289*G_edge[17]*dv1_sq-0.375*G_skin[6]*dv1_sq+0.375*G_edge[6]*dv1_sq; 
  diff_coeff_yy[7] = (-0.6495190528383289*G_skin[18]*dv1_sq)-0.6495190528383289*G_edge[18]*dv1_sq-0.375*G_skin[7]*dv1_sq+0.375*G_edge[7]*dv1_sq; 
  diff_coeff_yy[8] = (-0.6495190528383289*G_skin[19]*dv1_sq)-0.6495190528383289*G_edge[19]*dv1_sq-0.375*G_skin[8]*dv1_sq+0.375*G_edge[8]*dv1_sq; 
  diff_coeff_yy[9] = (-2.5*G_skin[9]*dv1_sq)-5.0*G_edge[9]*dv1_sq-2.165063509461096*G_skin[1]*dv1_sq+2.165063509461096*G_edge[1]*dv1_sq; 
  diff_coeff_yy[10] = (-2.5*G_skin[10]*dv1_sq)-5.0*G_edge[10]*dv1_sq-2.165063509461096*G_skin[2]*dv1_sq+2.165063509461096*G_edge[2]*dv1_sq; 
  diff_coeff_yy[11] = (-2.5*G_skin[11]*dv1_sq)-5.0*G_edge[11]*dv1_sq-2.165063509461096*G_skin[3]*dv1_sq+2.165063509461096*G_edge[3]*dv1_sq; 
  diff_coeff_yy[12] = (-0.6495190528383289*G_skin[23]*dv1_sq)-0.6495190528383289*G_edge[23]*dv1_sq-0.375*G_skin[12]*dv1_sq+0.375*G_edge[12]*dv1_sq; 
  diff_coeff_yy[13] = (-0.6495190528383289*G_skin[24]*dv1_sq)-0.6495190528383289*G_edge[24]*dv1_sq-0.375*G_skin[13]*dv1_sq+0.375*G_edge[13]*dv1_sq; 
  diff_coeff_yy[14] = (-0.6495190528383289*G_skin[25]*dv1_sq)-0.6495190528383289*G_edge[25]*dv1_sq-0.375*G_skin[14]*dv1_sq+0.375*G_edge[14]*dv1_sq; 
  diff_coeff_yy[15] = (-2.5*G_skin[15]*dv1_sq)-5.0*G_edge[15]*dv1_sq-2.165063509461096*G_skin[5]*dv1_sq+2.165063509461096*G_edge[5]*dv1_sq; 
  diff_coeff_yy[16] = (-0.6495190528383289*G_skin[26]*dv1_sq)-0.6495190528383289*G_edge[26]*dv1_sq-0.375*G_skin[16]*dv1_sq+0.375*G_edge[16]*dv1_sq; 
  diff_coeff_yy[17] = (-2.5*G_skin[17]*dv1_sq)-5.0*G_edge[17]*dv1_sq-2.165063509461096*G_skin[6]*dv1_sq+2.165063509461096*G_edge[6]*dv1_sq; 
  diff_coeff_yy[18] = (-2.5*G_skin[18]*dv1_sq)-5.0*G_edge[18]*dv1_sq-2.165063509461096*G_skin[7]*dv1_sq+2.165063509461096*G_edge[7]*dv1_sq; 
  diff_coeff_yy[19] = (-2.5*G_skin[19]*dv1_sq)-5.0*G_edge[19]*dv1_sq-2.165063509461096*G_skin[8]*dv1_sq+2.165063509461096*G_edge[8]*dv1_sq; 
  diff_coeff_yy[20] = (-0.6495190528383289*G_skin[28]*dv1_sq)-0.6495190528383289*G_edge[28]*dv1_sq-0.375*G_skin[20]*dv1_sq+0.375*G_edge[20]*dv1_sq; 
  diff_coeff_yy[21] = (-0.6495190528383289*G_skin[29]*dv1_sq)-0.6495190528383289*G_edge[29]*dv1_sq-0.375*G_skin[21]*dv1_sq+0.375*G_edge[21]*dv1_sq; 
  diff_coeff_yy[22] = (-0.6495190528383289*G_skin[30]*dv1_sq)-0.6495190528383289*G_edge[30]*dv1_sq-0.375*G_skin[22]*dv1_sq+0.375*G_edge[22]*dv1_sq; 
  diff_coeff_yy[23] = (-2.5*G_skin[23]*dv1_sq)-5.0*G_edge[23]*dv1_sq-2.165063509461096*G_skin[12]*dv1_sq+2.165063509461096*G_edge[12]*dv1_sq; 
  diff_coeff_yy[24] = (-2.5*G_skin[24]*dv1_sq)-5.0*G_edge[24]*dv1_sq-2.165063509461096*G_skin[13]*dv1_sq+2.165063509461096*G_edge[13]*dv1_sq; 
  diff_coeff_yy[25] = (-2.5*G_skin[25]*dv1_sq)-5.0*G_edge[25]*dv1_sq-2.165063509461096*G_skin[14]*dv1_sq+2.165063509461096*G_edge[14]*dv1_sq; 
  diff_coeff_yy[26] = (-2.5*G_skin[26]*dv1_sq)-5.0*G_edge[26]*dv1_sq-2.165063509461096*G_skin[16]*dv1_sq+2.165063509461096*G_edge[16]*dv1_sq; 
  diff_coeff_yy[27] = (-0.6495190528383289*G_skin[31]*dv1_sq)-0.6495190528383289*G_edge[31]*dv1_sq-0.375*G_skin[27]*dv1_sq+0.375*G_edge[27]*dv1_sq; 
  diff_coeff_yy[28] = (-2.5*G_skin[28]*dv1_sq)-5.0*G_edge[28]*dv1_sq-2.165063509461096*G_skin[20]*dv1_sq+2.165063509461096*G_edge[20]*dv1_sq; 
  diff_coeff_yy[29] = (-2.5*G_skin[29]*dv1_sq)-5.0*G_edge[29]*dv1_sq-2.165063509461096*G_skin[21]*dv1_sq+2.165063509461096*G_edge[21]*dv1_sq; 
  diff_coeff_yy[30] = (-2.5*G_skin[30]*dv1_sq)-5.0*G_edge[30]*dv1_sq-2.165063509461096*G_skin[22]*dv1_sq+2.165063509461096*G_edge[22]*dv1_sq; 
  diff_coeff_yy[31] = (-2.5*G_skin[31]*dv1_sq)-5.0*G_edge[31]*dv1_sq-2.165063509461096*G_skin[27]*dv1_sq+2.165063509461096*G_edge[27]*dv1_sq; 
 } else {
  diff_coeff_yy[0] = 0.6495190528383289*G_skin[4]*dv1_sq+0.6495190528383289*G_edge[4]*dv1_sq-0.375*G_skin[0]*dv1_sq+0.375*G_edge[0]*dv1_sq; 
  diff_coeff_yy[1] = 0.6495190528383289*G_skin[9]*dv1_sq+0.6495190528383289*G_edge[9]*dv1_sq-0.375*G_skin[1]*dv1_sq+0.375*G_edge[1]*dv1_sq; 
  diff_coeff_yy[2] = 0.6495190528383289*G_skin[10]*dv1_sq+0.6495190528383289*G_edge[10]*dv1_sq-0.375*G_skin[2]*dv1_sq+0.375*G_edge[2]*dv1_sq; 
  diff_coeff_yy[3] = 0.6495190528383289*G_skin[11]*dv1_sq+0.6495190528383289*G_edge[11]*dv1_sq-0.375*G_skin[3]*dv1_sq+0.375*G_edge[3]*dv1_sq; 
  diff_coeff_yy[4] = (-2.5*G_skin[4]*dv1_sq)-5.0*G_edge[4]*dv1_sq+2.165063509461096*G_skin[0]*dv1_sq-2.165063509461096*G_edge[0]*dv1_sq; 
  diff_coeff_yy[5] = 0.6495190528383289*G_skin[15]*dv1_sq+0.6495190528383289*G_edge[15]*dv1_sq-0.375*G_skin[5]*dv1_sq+0.375*G_edge[5]*dv1_sq; 
  diff_coeff_yy[6] = 0.6495190528383289*G_skin[17]*dv1_sq+0.6495190528383289*G_edge[17]*dv1_sq-0.375*G_skin[6]*dv1_sq+0.375*G_edge[6]*dv1_sq; 
  diff_coeff_yy[7] = 0.6495190528383289*G_skin[18]*dv1_sq+0.6495190528383289*G_edge[18]*dv1_sq-0.375*G_skin[7]*dv1_sq+0.375*G_edge[7]*dv1_sq; 
  diff_coeff_yy[8] = 0.6495190528383289*G_skin[19]*dv1_sq+0.6495190528383289*G_edge[19]*dv1_sq-0.375*G_skin[8]*dv1_sq+0.375*G_edge[8]*dv1_sq; 
  diff_coeff_yy[9] = (-2.5*G_skin[9]*dv1_sq)-5.0*G_edge[9]*dv1_sq+2.165063509461096*G_skin[1]*dv1_sq-2.165063509461096*G_edge[1]*dv1_sq; 
  diff_coeff_yy[10] = (-2.5*G_skin[10]*dv1_sq)-5.0*G_edge[10]*dv1_sq+2.165063509461096*G_skin[2]*dv1_sq-2.165063509461096*G_edge[2]*dv1_sq; 
  diff_coeff_yy[11] = (-2.5*G_skin[11]*dv1_sq)-5.0*G_edge[11]*dv1_sq+2.165063509461096*G_skin[3]*dv1_sq-2.165063509461096*G_edge[3]*dv1_sq; 
  diff_coeff_yy[12] = 0.6495190528383289*G_skin[23]*dv1_sq+0.6495190528383289*G_edge[23]*dv1_sq-0.375*G_skin[12]*dv1_sq+0.375*G_edge[12]*dv1_sq; 
  diff_coeff_yy[13] = 0.6495190528383289*G_skin[24]*dv1_sq+0.6495190528383289*G_edge[24]*dv1_sq-0.375*G_skin[13]*dv1_sq+0.375*G_edge[13]*dv1_sq; 
  diff_coeff_yy[14] = 0.6495190528383289*G_skin[25]*dv1_sq+0.6495190528383289*G_edge[25]*dv1_sq-0.375*G_skin[14]*dv1_sq+0.375*G_edge[14]*dv1_sq; 
  diff_coeff_yy[15] = (-2.5*G_skin[15]*dv1_sq)-5.0*G_edge[15]*dv1_sq+2.165063509461096*G_skin[5]*dv1_sq-2.165063509461096*G_edge[5]*dv1_sq; 
  diff_coeff_yy[16] = 0.6495190528383289*G_skin[26]*dv1_sq+0.6495190528383289*G_edge[26]*dv1_sq-0.375*G_skin[16]*dv1_sq+0.375*G_edge[16]*dv1_sq; 
  diff_coeff_yy[17] = (-2.5*G_skin[17]*dv1_sq)-5.0*G_edge[17]*dv1_sq+2.165063509461096*G_skin[6]*dv1_sq-2.165063509461096*G_edge[6]*dv1_sq; 
  diff_coeff_yy[18] = (-2.5*G_skin[18]*dv1_sq)-5.0*G_edge[18]*dv1_sq+2.165063509461096*G_skin[7]*dv1_sq-2.165063509461096*G_edge[7]*dv1_sq; 
  diff_coeff_yy[19] = (-2.5*G_skin[19]*dv1_sq)-5.0*G_edge[19]*dv1_sq+2.165063509461096*G_skin[8]*dv1_sq-2.165063509461096*G_edge[8]*dv1_sq; 
  diff_coeff_yy[20] = 0.6495190528383289*G_skin[28]*dv1_sq+0.6495190528383289*G_edge[28]*dv1_sq-0.375*G_skin[20]*dv1_sq+0.375*G_edge[20]*dv1_sq; 
  diff_coeff_yy[21] = 0.6495190528383289*G_skin[29]*dv1_sq+0.6495190528383289*G_edge[29]*dv1_sq-0.375*G_skin[21]*dv1_sq+0.375*G_edge[21]*dv1_sq; 
  diff_coeff_yy[22] = 0.6495190528383289*G_skin[30]*dv1_sq+0.6495190528383289*G_edge[30]*dv1_sq-0.375*G_skin[22]*dv1_sq+0.375*G_edge[22]*dv1_sq; 
  diff_coeff_yy[23] = (-2.5*G_skin[23]*dv1_sq)-5.0*G_edge[23]*dv1_sq+2.165063509461096*G_skin[12]*dv1_sq-2.165063509461096*G_edge[12]*dv1_sq; 
  diff_coeff_yy[24] = (-2.5*G_skin[24]*dv1_sq)-5.0*G_edge[24]*dv1_sq+2.165063509461096*G_skin[13]*dv1_sq-2.165063509461096*G_edge[13]*dv1_sq; 
  diff_coeff_yy[25] = (-2.5*G_skin[25]*dv1_sq)-5.0*G_edge[25]*dv1_sq+2.165063509461096*G_skin[14]*dv1_sq-2.165063509461096*G_edge[14]*dv1_sq; 
  diff_coeff_yy[26] = (-2.5*G_skin[26]*dv1_sq)-5.0*G_edge[26]*dv1_sq+2.165063509461096*G_skin[16]*dv1_sq-2.165063509461096*G_edge[16]*dv1_sq; 
  diff_coeff_yy[27] = 0.6495190528383289*G_skin[31]*dv1_sq+0.6495190528383289*G_edge[31]*dv1_sq-0.375*G_skin[27]*dv1_sq+0.375*G_edge[27]*dv1_sq; 
  diff_coeff_yy[28] = (-2.5*G_skin[28]*dv1_sq)-5.0*G_edge[28]*dv1_sq+2.165063509461096*G_skin[20]*dv1_sq-2.165063509461096*G_edge[20]*dv1_sq; 
  diff_coeff_yy[29] = (-2.5*G_skin[29]*dv1_sq)-5.0*G_edge[29]*dv1_sq+2.165063509461096*G_skin[21]*dv1_sq-2.165063509461096*G_edge[21]*dv1_sq; 
  diff_coeff_yy[30] = (-2.5*G_skin[30]*dv1_sq)-5.0*G_edge[30]*dv1_sq+2.165063509461096*G_skin[22]*dv1_sq-2.165063509461096*G_edge[22]*dv1_sq; 
  diff_coeff_yy[31] = (-2.5*G_skin[31]*dv1_sq)-5.0*G_edge[31]*dv1_sq+2.165063509461096*G_skin[27]*dv1_sq-2.165063509461096*G_edge[27]*dv1_sq; 
  } 
}
