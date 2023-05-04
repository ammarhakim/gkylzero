#include <gkyl_fpo_vlasov_kernels.h> 
GKYL_CU_DH void fpo_diff_coeff_recov_surf_vy_1x3v_ser_p1(const int edge, const double *dxv, const double *G_skin, const double *G_edge, double *diff_coeff) { 
  // dxv[NDIM]: Cell spacing. 
  // G_skin/edge:   Input potential skin/edge cells in recovery direction.

  const double dv1 = 2.0/dxv[2]; 
  const double dv1_sq = dv1*dv1; 
  double *diff_coeff_xx = &diff_coeff[0]; 
  double *diff_coeff_yy = &diff_coeff[16]; 
  double *diff_coeff_zz = &diff_coeff[32]; 
  if (edge == 1) {
  diff_coeff_yy[0] = (-0.6495190528383289*G_skin[3]*dv1_sq)-0.6495190528383289*G_edge[3]*dv1_sq-0.375*G_skin[0]*dv1_sq+0.375*G_edge[0]*dv1_sq; 
  diff_coeff_yy[1] = (-0.6495190528383289*G_skin[6]*dv1_sq)-0.6495190528383289*G_edge[6]*dv1_sq-0.375*G_skin[1]*dv1_sq+0.375*G_edge[1]*dv1_sq; 
  diff_coeff_yy[2] = (-0.6495190528383289*G_skin[7]*dv1_sq)-0.6495190528383289*G_edge[7]*dv1_sq-0.375*G_skin[2]*dv1_sq+0.375*G_edge[2]*dv1_sq; 
  diff_coeff_yy[3] = (-2.5*G_skin[3]*dv1_sq)-5.0*G_edge[3]*dv1_sq-2.165063509461096*G_skin[0]*dv1_sq+2.165063509461096*G_edge[0]*dv1_sq; 
  diff_coeff_yy[4] = (-0.6495190528383289*G_skin[10]*dv1_sq)-0.6495190528383289*G_edge[10]*dv1_sq-0.375*G_skin[4]*dv1_sq+0.375*G_edge[4]*dv1_sq; 
  diff_coeff_yy[5] = (-0.6495190528383289*G_skin[11]*dv1_sq)-0.6495190528383289*G_edge[11]*dv1_sq-0.375*G_skin[5]*dv1_sq+0.375*G_edge[5]*dv1_sq; 
  diff_coeff_yy[6] = (-2.5*G_skin[6]*dv1_sq)-5.0*G_edge[6]*dv1_sq-2.165063509461096*G_skin[1]*dv1_sq+2.165063509461096*G_edge[1]*dv1_sq; 
  diff_coeff_yy[7] = (-2.5*G_skin[7]*dv1_sq)-5.0*G_edge[7]*dv1_sq-2.165063509461096*G_skin[2]*dv1_sq+2.165063509461096*G_edge[2]*dv1_sq; 
  diff_coeff_yy[8] = (-0.6495190528383289*G_skin[13]*dv1_sq)-0.6495190528383289*G_edge[13]*dv1_sq-0.375*G_skin[8]*dv1_sq+0.375*G_edge[8]*dv1_sq; 
  diff_coeff_yy[9] = (-0.6495190528383289*G_skin[14]*dv1_sq)-0.6495190528383289*G_edge[14]*dv1_sq-0.375*G_skin[9]*dv1_sq+0.375*G_edge[9]*dv1_sq; 
  diff_coeff_yy[10] = (-2.5*G_skin[10]*dv1_sq)-5.0*G_edge[10]*dv1_sq-2.165063509461096*G_skin[4]*dv1_sq+2.165063509461096*G_edge[4]*dv1_sq; 
  diff_coeff_yy[11] = (-2.5*G_skin[11]*dv1_sq)-5.0*G_edge[11]*dv1_sq-2.165063509461096*G_skin[5]*dv1_sq+2.165063509461096*G_edge[5]*dv1_sq; 
  diff_coeff_yy[12] = (-0.6495190528383289*G_skin[15]*dv1_sq)-0.6495190528383289*G_edge[15]*dv1_sq-0.375*G_skin[12]*dv1_sq+0.375*G_edge[12]*dv1_sq; 
  diff_coeff_yy[13] = (-2.5*G_skin[13]*dv1_sq)-5.0*G_edge[13]*dv1_sq-2.165063509461096*G_skin[8]*dv1_sq+2.165063509461096*G_edge[8]*dv1_sq; 
  diff_coeff_yy[14] = (-2.5*G_skin[14]*dv1_sq)-5.0*G_edge[14]*dv1_sq-2.165063509461096*G_skin[9]*dv1_sq+2.165063509461096*G_edge[9]*dv1_sq; 
  diff_coeff_yy[15] = (-2.5*G_skin[15]*dv1_sq)-5.0*G_edge[15]*dv1_sq-2.165063509461096*G_skin[12]*dv1_sq+2.165063509461096*G_edge[12]*dv1_sq; 
 } else {
  diff_coeff_yy[0] = 0.6495190528383289*G_skin[3]*dv1_sq+0.6495190528383289*G_edge[3]*dv1_sq-0.375*G_skin[0]*dv1_sq+0.375*G_edge[0]*dv1_sq; 
  diff_coeff_yy[1] = 0.6495190528383289*G_skin[6]*dv1_sq+0.6495190528383289*G_edge[6]*dv1_sq-0.375*G_skin[1]*dv1_sq+0.375*G_edge[1]*dv1_sq; 
  diff_coeff_yy[2] = 0.6495190528383289*G_skin[7]*dv1_sq+0.6495190528383289*G_edge[7]*dv1_sq-0.375*G_skin[2]*dv1_sq+0.375*G_edge[2]*dv1_sq; 
  diff_coeff_yy[3] = (-2.5*G_skin[3]*dv1_sq)-5.0*G_edge[3]*dv1_sq+2.165063509461096*G_skin[0]*dv1_sq-2.165063509461096*G_edge[0]*dv1_sq; 
  diff_coeff_yy[4] = 0.6495190528383289*G_skin[10]*dv1_sq+0.6495190528383289*G_edge[10]*dv1_sq-0.375*G_skin[4]*dv1_sq+0.375*G_edge[4]*dv1_sq; 
  diff_coeff_yy[5] = 0.6495190528383289*G_skin[11]*dv1_sq+0.6495190528383289*G_edge[11]*dv1_sq-0.375*G_skin[5]*dv1_sq+0.375*G_edge[5]*dv1_sq; 
  diff_coeff_yy[6] = (-2.5*G_skin[6]*dv1_sq)-5.0*G_edge[6]*dv1_sq+2.165063509461096*G_skin[1]*dv1_sq-2.165063509461096*G_edge[1]*dv1_sq; 
  diff_coeff_yy[7] = (-2.5*G_skin[7]*dv1_sq)-5.0*G_edge[7]*dv1_sq+2.165063509461096*G_skin[2]*dv1_sq-2.165063509461096*G_edge[2]*dv1_sq; 
  diff_coeff_yy[8] = 0.6495190528383289*G_skin[13]*dv1_sq+0.6495190528383289*G_edge[13]*dv1_sq-0.375*G_skin[8]*dv1_sq+0.375*G_edge[8]*dv1_sq; 
  diff_coeff_yy[9] = 0.6495190528383289*G_skin[14]*dv1_sq+0.6495190528383289*G_edge[14]*dv1_sq-0.375*G_skin[9]*dv1_sq+0.375*G_edge[9]*dv1_sq; 
  diff_coeff_yy[10] = (-2.5*G_skin[10]*dv1_sq)-5.0*G_edge[10]*dv1_sq+2.165063509461096*G_skin[4]*dv1_sq-2.165063509461096*G_edge[4]*dv1_sq; 
  diff_coeff_yy[11] = (-2.5*G_skin[11]*dv1_sq)-5.0*G_edge[11]*dv1_sq+2.165063509461096*G_skin[5]*dv1_sq-2.165063509461096*G_edge[5]*dv1_sq; 
  diff_coeff_yy[12] = 0.6495190528383289*G_skin[15]*dv1_sq+0.6495190528383289*G_edge[15]*dv1_sq-0.375*G_skin[12]*dv1_sq+0.375*G_edge[12]*dv1_sq; 
  diff_coeff_yy[13] = (-2.5*G_skin[13]*dv1_sq)-5.0*G_edge[13]*dv1_sq+2.165063509461096*G_skin[8]*dv1_sq-2.165063509461096*G_edge[8]*dv1_sq; 
  diff_coeff_yy[14] = (-2.5*G_skin[14]*dv1_sq)-5.0*G_edge[14]*dv1_sq+2.165063509461096*G_skin[9]*dv1_sq-2.165063509461096*G_edge[9]*dv1_sq; 
  diff_coeff_yy[15] = (-2.5*G_skin[15]*dv1_sq)-5.0*G_edge[15]*dv1_sq+2.165063509461096*G_skin[12]*dv1_sq-2.165063509461096*G_edge[12]*dv1_sq; 
  } 
}
