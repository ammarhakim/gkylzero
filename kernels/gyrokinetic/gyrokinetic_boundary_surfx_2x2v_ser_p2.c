#include <gkyl_gyrokinetic_kernels.h>
#include <gkyl_basis_ser_4x_p2_surfx1_eval_quad.h> 
#include <gkyl_basis_ser_4x_p2_upwind_quad_to_modal.h> 
GKYL_CU_DH void gyrokinetic_boundary_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *apar, const double *apardot, const int edge, const double *fedge, const double *fskin, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // jacobtot_inv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // cmag: coefficient multiplying parallel gradient.
  // b_i: covariant components of the field aligned unit vector.
  // phi: electrostatic potential .
  // apar: parallel component of magnetic vector potential.
  // apardot: time derivative of Apar.
  // edge: determines if the update is for the left edge (-1) or right edge (+1).
  // fskin,fedge: distribution function in skin cell/last edge cell.
  // out: output increment in center cell.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wy = w[1];
  double rdy2 = 2.0/dxv[1];
  double wvpar = w[2];
  double rdvpar2 = 2.0/dxv[2];
  double wmu = w[3];
  double rdmu2 = 2.0/dxv[3];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wySq = w[1]*w[1];
  double rdy2Sq = rdy2*rdy2;
  double wvparSq = w[2]*w[2];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[3]*w[3];
  double rdmu2Sq = rdmu2*rdmu2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[8];
  const double *b_z = &b_i[16];

  double hamil[48] = {0.}; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[11] = 2.0*(bmag[4]*wmu+phi[4]*q_); 
  hamil[12] = 2.0*phi[5]*q_; 
  hamil[13] = (0.5962847939999438*m_)/rdvpar2Sq; 
  hamil[19] = 2.0*phi[6]*q_; 
  hamil[20] = 2.0*phi[7]*q_; 
  hamil[25] = (1.154700538379251*bmag[4])/rdmu2; 

  double BstarXdBmag[48] = {0.}; 
  BstarXdBmag[0] = 0.01428571428571429*(3.872983346207417*((10.0*b_z[4]+15.65247584249853*b_z[0])*jacobtot_inv[4]+7.0*(2.23606797749979*jacobtot_inv[0]*b_z[4]+2.0*b_z[1]*jacobtot_inv[1]))*apar[6]+12.12435565298214*((5.0*apar[2]*b_z[4]+4.47213595499958*b_z[1]*apar[3])*jacobtot_inv[4]+4.47213595499958*jacobtot_inv[1]*apar[3]*b_z[4]+5.0*((b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*apar[3]+(b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*apar[2])))*rdy2; 
  BstarXdBmag[1] = 0.01428571428571429*(3.872983346207417*(24.59674775249769*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+14.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*apar[6]+1.732050807568877*((55.0*apar[3]*b_z[4]+31.30495168499706*(b_z[0]*apar[3]+b_z[1]*apar[2]))*jacobtot_inv[4]+7.0*(4.47213595499958*(jacobtot_inv[0]*apar[3]+jacobtot_inv[1]*apar[2])*b_z[4]+(9.0*b_z[1]*jacobtot_inv[1]+5.0*b_z[0]*jacobtot_inv[0])*apar[3]+5.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*apar[2])))*rdy2; 
  BstarXdBmag[2] = 0.1*(3.872983346207417*(4.47213595499958*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+5.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*apar[7]+19.36491673103709*(b_z[4]*jacobtot_inv[4]+b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*apar[5])*rdy2; 
  BstarXdBmag[5] = 0.01428571428571429*(3.872983346207417*((55.0*b_z[4]+31.30495168499706*b_z[0])*jacobtot_inv[4]+7.0*(4.47213595499958*jacobtot_inv[0]*b_z[4]+9.0*b_z[1]*jacobtot_inv[1]+5.0*b_z[0]*jacobtot_inv[0]))*apar[7]+12.12435565298214*(10.0*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+11.18033988749895*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*apar[5])*rdy2; 
  BstarXdBmag[11] = 0.01428571428571429*(3.872983346207417*((33.54101966249685*b_z[4]+10.0*b_z[0])*jacobtot_inv[4]+10.0*jacobtot_inv[0]*b_z[4]+2.23606797749979*(11.0*b_z[1]*jacobtot_inv[1]+7.0*b_z[0]*jacobtot_inv[0]))*apar[6]+1.732050807568877*((22.3606797749979*apar[2]*b_z[4]+5.0*(11.0*b_z[1]*apar[3]+7.0*b_z[0]*apar[2]))*jacobtot_inv[4]+5.0*(11.0*jacobtot_inv[1]*apar[3]+7.0*jacobtot_inv[0]*apar[2])*b_z[4]+31.30495168499706*((b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])*apar[3]+b_z[1]*jacobtot_inv[1]*apar[2])))*rdy2; 
  BstarXdBmag[19] = 0.01428571428571429*(1.732050807568877*(122.9837387624885*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+70.0*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*apar[7]+3.872983346207417*((22.3606797749979*b_z[4]+35.0*b_z[0])*jacobtot_inv[4]+7.0*(5.0*jacobtot_inv[0]*b_z[4]+4.47213595499958*b_z[1]*jacobtot_inv[1]))*apar[5])*rdy2; 

  if (edge == -1) { 

  double alphaR[20] = {0.}; 
  alphaR[0] = -(0.1767766952966368*((((19.36491673103708*b_z[4]+15.0*b_z[1]+8.660254037844387*b_z[0])*jacobtot_inv[4]+(15.0*jacobtot_inv[1]+8.660254037844387*jacobtot_inv[0])*b_z[4]+(11.61895003862225*b_z[1]+6.708203932499369*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(6.708203932499369*b_z[1]+3.872983346207417*b_z[0]))*hamil[19]+((15.0*b_z[4]+11.61895003862225*b_z[1]+6.708203932499369*b_z[0])*jacobtot_inv[4]+(11.61895003862225*jacobtot_inv[1]+6.708203932499369*jacobtot_inv[0])*b_z[4]+(9.0*b_z[1]+5.196152422706631*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(5.196152422706631*b_z[1]+3.0*b_z[0]))*hamil[5]+hamil[2]*((8.660254037844386*b_z[4]+6.708203932499369*b_z[1]+3.872983346207417*b_z[0])*jacobtot_inv[4]+(6.708203932499369*jacobtot_inv[1]+3.872983346207417*jacobtot_inv[0])*b_z[4]+(5.196152422706631*b_z[1]+3.0*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(3.0*b_z[1]+1.732050807568877*b_z[0])))*m_*rdy2+hamil[3]*((-3.872983346207417*BstarXdBmag[11])-3.0*BstarXdBmag[1]-1.732050807568877*BstarXdBmag[0])*q_*rdvpar2))/(m_*q_); 
  alphaR[1] = -(0.1767766952966368*((((33.54101966249684*b_z[4]+25.98076211353316*b_z[1]+15.0*b_z[0])*jacobtot_inv[4]+(25.98076211353316*jacobtot_inv[1]+15.0*jacobtot_inv[0])*b_z[4]+(20.1246117974981*b_z[1]+11.61895003862225*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(11.61895003862225*b_z[1]+6.708203932499369*b_z[0]))*hamil[20]+((19.36491673103709*b_z[4]+15.0*b_z[1]+8.660254037844386*b_z[0])*jacobtot_inv[4]+(15.0*jacobtot_inv[1]+8.660254037844386*jacobtot_inv[0])*b_z[4]+(11.61895003862225*b_z[1]+6.708203932499369*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(6.708203932499369*b_z[1]+3.872983346207417*b_z[0]))*hamil[12])*m_*rdy2+hamil[3]*((-3.872983346207417*BstarXdBmag[19])-3.0*BstarXdBmag[5]-1.732050807568877*BstarXdBmag[2])*q_*rdvpar2))/(m_*q_); 
  alphaR[2] = (0.1767766952966368*(8.660254037844386*BstarXdBmag[11]+6.708203932499369*BstarXdBmag[1]+3.872983346207417*BstarXdBmag[0])*hamil[13]*rdvpar2)/m_; 
  alphaR[4] = (0.1767766952966368*hamil[13]*(8.660254037844387*BstarXdBmag[19]+6.708203932499369*BstarXdBmag[5]+3.872983346207417*BstarXdBmag[2])*rdvpar2)/m_; 

  double fUpOrdR[27] = {0.};
  if (alphaR[4]-0.7453559924999297*alphaR[2]-0.7453559924999297*alphaR[1]+0.555555555555555*alphaR[0] > 0.) {
    fUpOrdR[0] = ser_4x_p2_surfx1_eval_quad_node_0_r(fskin); 
  } else { 
    fUpOrdR[0] = ser_4x_p2_surfx1_eval_quad_node_0_l(fedge); 
  } 
  if (alphaR[4]-0.7453559924999297*alphaR[2]-0.7453559924999297*alphaR[1]+0.555555555555555*alphaR[0] > 0.) {
    fUpOrdR[1] = ser_4x_p2_surfx1_eval_quad_node_1_r(fskin); 
  } else { 
    fUpOrdR[1] = ser_4x_p2_surfx1_eval_quad_node_1_l(fedge); 
  } 
  if (alphaR[4]-0.7453559924999297*alphaR[2]-0.7453559924999297*alphaR[1]+0.555555555555555*alphaR[0] > 0.) {
    fUpOrdR[2] = ser_4x_p2_surfx1_eval_quad_node_2_r(fskin); 
  } else { 
    fUpOrdR[2] = ser_4x_p2_surfx1_eval_quad_node_2_l(fedge); 
  } 
  if (alphaR[0]-1.341640786499875*alphaR[1] > 0.) {
    fUpOrdR[3] = ser_4x_p2_surfx1_eval_quad_node_3_r(fskin); 
  } else { 
    fUpOrdR[3] = ser_4x_p2_surfx1_eval_quad_node_3_l(fedge); 
  } 
  if (alphaR[0]-1.341640786499875*alphaR[1] > 0.) {
    fUpOrdR[4] = ser_4x_p2_surfx1_eval_quad_node_4_r(fskin); 
  } else { 
    fUpOrdR[4] = ser_4x_p2_surfx1_eval_quad_node_4_l(fedge); 
  } 
  if (alphaR[0]-1.341640786499875*alphaR[1] > 0.) {
    fUpOrdR[5] = ser_4x_p2_surfx1_eval_quad_node_5_r(fskin); 
  } else { 
    fUpOrdR[5] = ser_4x_p2_surfx1_eval_quad_node_5_l(fedge); 
  } 
  if ((-1.0*alphaR[4])+0.7453559924999297*alphaR[2]-0.7453559924999297*alphaR[1]+0.555555555555555*alphaR[0] > 0.) {
    fUpOrdR[6] = ser_4x_p2_surfx1_eval_quad_node_6_r(fskin); 
  } else { 
    fUpOrdR[6] = ser_4x_p2_surfx1_eval_quad_node_6_l(fedge); 
  } 
  if ((-1.0*alphaR[4])+0.7453559924999297*alphaR[2]-0.7453559924999297*alphaR[1]+0.555555555555555*alphaR[0] > 0.) {
    fUpOrdR[7] = ser_4x_p2_surfx1_eval_quad_node_7_r(fskin); 
  } else { 
    fUpOrdR[7] = ser_4x_p2_surfx1_eval_quad_node_7_l(fedge); 
  } 
  if ((-1.0*alphaR[4])+0.7453559924999297*alphaR[2]-0.7453559924999297*alphaR[1]+0.555555555555555*alphaR[0] > 0.) {
    fUpOrdR[8] = ser_4x_p2_surfx1_eval_quad_node_8_r(fskin); 
  } else { 
    fUpOrdR[8] = ser_4x_p2_surfx1_eval_quad_node_8_l(fedge); 
  } 
  if (alphaR[0]-1.341640786499875*alphaR[2] > 0.) {
    fUpOrdR[9] = ser_4x_p2_surfx1_eval_quad_node_9_r(fskin); 
  } else { 
    fUpOrdR[9] = ser_4x_p2_surfx1_eval_quad_node_9_l(fedge); 
  } 
  if (alphaR[0]-1.341640786499875*alphaR[2] > 0.) {
    fUpOrdR[10] = ser_4x_p2_surfx1_eval_quad_node_10_r(fskin); 
  } else { 
    fUpOrdR[10] = ser_4x_p2_surfx1_eval_quad_node_10_l(fedge); 
  } 
  if (alphaR[0]-1.341640786499875*alphaR[2] > 0.) {
    fUpOrdR[11] = ser_4x_p2_surfx1_eval_quad_node_11_r(fskin); 
  } else { 
    fUpOrdR[11] = ser_4x_p2_surfx1_eval_quad_node_11_l(fedge); 
  } 
  if (alphaR[0] > 0.) {
    fUpOrdR[12] = ser_4x_p2_surfx1_eval_quad_node_12_r(fskin); 
  } else { 
    fUpOrdR[12] = ser_4x_p2_surfx1_eval_quad_node_12_l(fedge); 
  } 
  if (alphaR[0] > 0.) {
    fUpOrdR[13] = ser_4x_p2_surfx1_eval_quad_node_13_r(fskin); 
  } else { 
    fUpOrdR[13] = ser_4x_p2_surfx1_eval_quad_node_13_l(fedge); 
  } 
  if (alphaR[0] > 0.) {
    fUpOrdR[14] = ser_4x_p2_surfx1_eval_quad_node_14_r(fskin); 
  } else { 
    fUpOrdR[14] = ser_4x_p2_surfx1_eval_quad_node_14_l(fedge); 
  } 
  if (alphaR[2]+0.7453559924999295*alphaR[0] > 0.) {
    fUpOrdR[15] = ser_4x_p2_surfx1_eval_quad_node_15_r(fskin); 
  } else { 
    fUpOrdR[15] = ser_4x_p2_surfx1_eval_quad_node_15_l(fedge); 
  } 
  if (alphaR[2]+0.7453559924999295*alphaR[0] > 0.) {
    fUpOrdR[16] = ser_4x_p2_surfx1_eval_quad_node_16_r(fskin); 
  } else { 
    fUpOrdR[16] = ser_4x_p2_surfx1_eval_quad_node_16_l(fedge); 
  } 
  if (alphaR[2]+0.7453559924999295*alphaR[0] > 0.) {
    fUpOrdR[17] = ser_4x_p2_surfx1_eval_quad_node_17_r(fskin); 
  } else { 
    fUpOrdR[17] = ser_4x_p2_surfx1_eval_quad_node_17_l(fedge); 
  } 
  if ((-1.0*alphaR[4])-0.7453559924999297*alphaR[2]+0.7453559924999297*alphaR[1]+0.555555555555555*alphaR[0] > 0.) {
    fUpOrdR[18] = ser_4x_p2_surfx1_eval_quad_node_18_r(fskin); 
  } else { 
    fUpOrdR[18] = ser_4x_p2_surfx1_eval_quad_node_18_l(fedge); 
  } 
  if ((-1.0*alphaR[4])-0.7453559924999297*alphaR[2]+0.7453559924999297*alphaR[1]+0.555555555555555*alphaR[0] > 0.) {
    fUpOrdR[19] = ser_4x_p2_surfx1_eval_quad_node_19_r(fskin); 
  } else { 
    fUpOrdR[19] = ser_4x_p2_surfx1_eval_quad_node_19_l(fedge); 
  } 
  if ((-1.0*alphaR[4])-0.7453559924999297*alphaR[2]+0.7453559924999297*alphaR[1]+0.555555555555555*alphaR[0] > 0.) {
    fUpOrdR[20] = ser_4x_p2_surfx1_eval_quad_node_20_r(fskin); 
  } else { 
    fUpOrdR[20] = ser_4x_p2_surfx1_eval_quad_node_20_l(fedge); 
  } 
  if (alphaR[1]+0.7453559924999295*alphaR[0] > 0.) {
    fUpOrdR[21] = ser_4x_p2_surfx1_eval_quad_node_21_r(fskin); 
  } else { 
    fUpOrdR[21] = ser_4x_p2_surfx1_eval_quad_node_21_l(fedge); 
  } 
  if (alphaR[1]+0.7453559924999295*alphaR[0] > 0.) {
    fUpOrdR[22] = ser_4x_p2_surfx1_eval_quad_node_22_r(fskin); 
  } else { 
    fUpOrdR[22] = ser_4x_p2_surfx1_eval_quad_node_22_l(fedge); 
  } 
  if (alphaR[1]+0.7453559924999295*alphaR[0] > 0.) {
    fUpOrdR[23] = ser_4x_p2_surfx1_eval_quad_node_23_r(fskin); 
  } else { 
    fUpOrdR[23] = ser_4x_p2_surfx1_eval_quad_node_23_l(fedge); 
  } 
  if (alphaR[4]+0.7453559924999297*alphaR[2]+0.7453559924999297*alphaR[1]+0.555555555555555*alphaR[0] > 0.) {
    fUpOrdR[24] = ser_4x_p2_surfx1_eval_quad_node_24_r(fskin); 
  } else { 
    fUpOrdR[24] = ser_4x_p2_surfx1_eval_quad_node_24_l(fedge); 
  } 
  if (alphaR[4]+0.7453559924999297*alphaR[2]+0.7453559924999297*alphaR[1]+0.555555555555555*alphaR[0] > 0.) {
    fUpOrdR[25] = ser_4x_p2_surfx1_eval_quad_node_25_r(fskin); 
  } else { 
    fUpOrdR[25] = ser_4x_p2_surfx1_eval_quad_node_25_l(fedge); 
  } 
  if (alphaR[4]+0.7453559924999297*alphaR[2]+0.7453559924999297*alphaR[1]+0.555555555555555*alphaR[0] > 0.) {
    fUpOrdR[26] = ser_4x_p2_surfx1_eval_quad_node_26_r(fskin); 
  } else { 
    fUpOrdR[26] = ser_4x_p2_surfx1_eval_quad_node_26_l(fedge); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpR[20] = {0.};
  ser_4x_p2_upwind_quad_to_modal(fUpOrdR, fUpR); 

  double GhatR[48] = {0.}; 
  GhatR[0] = 0.3535533905932737*alphaR[4]*fUpR[4]+0.3535533905932737*alphaR[2]*fUpR[2]+0.3535533905932737*alphaR[1]*fUpR[1]+0.3535533905932737*alphaR[0]*fUpR[0]; 
  GhatR[1] = 0.3162277660168379*alphaR[4]*fUpR[11]+0.3162277660168379*alphaR[1]*fUpR[7]+0.3535533905932737*alphaR[2]*fUpR[4]+0.3535533905932737*fUpR[2]*alphaR[4]+0.3535533905932737*alphaR[0]*fUpR[1]+0.3535533905932737*fUpR[0]*alphaR[1]; 
  GhatR[2] = 0.3162277660168379*alphaR[4]*fUpR[12]+0.3162277660168379*alphaR[2]*fUpR[8]+0.3535533905932737*alphaR[1]*fUpR[4]+0.3535533905932737*fUpR[1]*alphaR[4]+0.3535533905932737*alphaR[0]*fUpR[2]+0.3535533905932737*fUpR[0]*alphaR[2]; 
  GhatR[3] = 0.3535533905932737*alphaR[4]*fUpR[10]+0.3535533905932737*alphaR[2]*fUpR[6]+0.3535533905932737*alphaR[1]*fUpR[5]+0.3535533905932737*alphaR[0]*fUpR[3]; 
  GhatR[4] = 0.3162277660168379*alphaR[2]*fUpR[12]+0.3162277660168379*alphaR[1]*fUpR[11]+0.3162277660168379*alphaR[4]*fUpR[8]+0.3162277660168379*alphaR[4]*fUpR[7]+0.3535533905932737*alphaR[0]*fUpR[4]+0.3535533905932737*fUpR[0]*alphaR[4]+0.3535533905932737*alphaR[1]*fUpR[2]+0.3535533905932737*fUpR[1]*alphaR[2]; 
  GhatR[5] = 0.3162277660168379*alphaR[4]*fUpR[17]+0.3162277660168379*alphaR[1]*fUpR[13]+0.3535533905932737*alphaR[2]*fUpR[10]+0.3535533905932737*alphaR[4]*fUpR[6]+0.3535533905932737*alphaR[0]*fUpR[5]+0.3535533905932737*alphaR[1]*fUpR[3]; 
  GhatR[6] = 0.3162277660168379*alphaR[4]*fUpR[18]+0.3162277660168379*alphaR[2]*fUpR[14]+0.3535533905932737*alphaR[1]*fUpR[10]+0.3535533905932737*alphaR[0]*fUpR[6]+0.3535533905932737*alphaR[4]*fUpR[5]+0.3535533905932737*alphaR[2]*fUpR[3]; 
  GhatR[7] = 0.3535533905932737*alphaR[2]*fUpR[11]+0.3535533905932737*alphaR[0]*fUpR[7]+0.3162277660168379*alphaR[4]*fUpR[4]+0.3162277660168379*alphaR[1]*fUpR[1]; 
  GhatR[8] = 0.3535533905932737*alphaR[1]*fUpR[12]+0.3535533905932737*alphaR[0]*fUpR[8]+0.3162277660168379*alphaR[4]*fUpR[4]+0.3162277660168379*alphaR[2]*fUpR[2]; 
  GhatR[9] = 0.3535533905932737*alphaR[4]*fUpR[19]+0.3535533905932737*alphaR[2]*fUpR[16]+0.3535533905932737*alphaR[1]*fUpR[15]+0.3535533905932737*alphaR[0]*fUpR[9]; 
  GhatR[10] = 0.3162277660168379*alphaR[2]*fUpR[18]+0.3162277660168379*alphaR[1]*fUpR[17]+0.3162277660168379*alphaR[4]*fUpR[14]+0.3162277660168379*alphaR[4]*fUpR[13]+0.3535533905932737*alphaR[0]*fUpR[10]+0.3535533905932737*alphaR[1]*fUpR[6]+0.3535533905932737*alphaR[2]*fUpR[5]+0.3535533905932737*fUpR[3]*alphaR[4]; 
  GhatR[11] = 0.2828427124746191*alphaR[4]*fUpR[12]+0.3535533905932737*alphaR[0]*fUpR[11]+0.3535533905932737*alphaR[2]*fUpR[7]+0.3162277660168379*alphaR[1]*fUpR[4]+0.3162277660168379*fUpR[1]*alphaR[4]; 
  GhatR[12] = 0.3535533905932737*alphaR[0]*fUpR[12]+0.2828427124746191*alphaR[4]*fUpR[11]+0.3535533905932737*alphaR[1]*fUpR[8]+0.3162277660168379*alphaR[2]*fUpR[4]+0.3162277660168379*fUpR[2]*alphaR[4]; 
  GhatR[13] = 0.3535533905932737*alphaR[2]*fUpR[17]+0.3535533905932737*alphaR[0]*fUpR[13]+0.3162277660168379*alphaR[4]*fUpR[10]+0.3162277660168379*alphaR[1]*fUpR[5]; 
  GhatR[14] = 0.3535533905932737*alphaR[1]*fUpR[18]+0.3535533905932737*alphaR[0]*fUpR[14]+0.3162277660168379*alphaR[4]*fUpR[10]+0.3162277660168379*alphaR[2]*fUpR[6]; 
  GhatR[15] = 0.3535533905932737*alphaR[2]*fUpR[19]+0.3535533905932737*alphaR[4]*fUpR[16]+0.3535533905932737*alphaR[0]*fUpR[15]+0.3535533905932737*alphaR[1]*fUpR[9]; 
  GhatR[16] = 0.3535533905932737*alphaR[1]*fUpR[19]+0.3535533905932737*alphaR[0]*fUpR[16]+0.3535533905932737*alphaR[4]*fUpR[15]+0.3535533905932737*alphaR[2]*fUpR[9]; 
  GhatR[17] = 0.2828427124746191*alphaR[4]*fUpR[18]+0.3535533905932737*alphaR[0]*fUpR[17]+0.3535533905932737*alphaR[2]*fUpR[13]+0.3162277660168379*alphaR[1]*fUpR[10]+0.3162277660168379*alphaR[4]*fUpR[5]; 
  GhatR[18] = 0.3535533905932737*alphaR[0]*fUpR[18]+0.2828427124746191*alphaR[4]*fUpR[17]+0.3535533905932737*alphaR[1]*fUpR[14]+0.3162277660168379*alphaR[2]*fUpR[10]+0.3162277660168379*alphaR[4]*fUpR[6]; 
  GhatR[19] = 0.3535533905932737*alphaR[0]*fUpR[19]+0.3535533905932737*alphaR[1]*fUpR[16]+0.3535533905932737*alphaR[2]*fUpR[15]+0.3535533905932737*alphaR[4]*fUpR[9]; 

  out[0] += -0.7071067811865475*GhatR[0]*rdx2; 
  out[1] += -1.224744871391589*GhatR[0]*rdx2; 
  out[2] += -0.7071067811865475*GhatR[1]*rdx2; 
  out[3] += -0.7071067811865475*GhatR[2]*rdx2; 
  out[4] += -0.7071067811865475*GhatR[3]*rdx2; 
  out[5] += -1.224744871391589*GhatR[1]*rdx2; 
  out[6] += -1.224744871391589*GhatR[2]*rdx2; 
  out[7] += -0.7071067811865475*GhatR[4]*rdx2; 
  out[8] += -1.224744871391589*GhatR[3]*rdx2; 
  out[9] += -0.7071067811865475*GhatR[5]*rdx2; 
  out[10] += -0.7071067811865475*GhatR[6]*rdx2; 
  out[11] += -1.58113883008419*GhatR[0]*rdx2; 
  out[12] += -0.7071067811865475*GhatR[7]*rdx2; 
  out[13] += -0.7071067811865475*GhatR[8]*rdx2; 
  out[14] += -0.7071067811865475*GhatR[9]*rdx2; 
  out[15] += -1.224744871391589*GhatR[4]*rdx2; 
  out[16] += -1.224744871391589*GhatR[5]*rdx2; 
  out[17] += -1.224744871391589*GhatR[6]*rdx2; 
  out[18] += -0.7071067811865475*GhatR[10]*rdx2; 
  out[19] += -1.58113883008419*GhatR[1]*rdx2; 
  out[20] += -1.224744871391589*GhatR[7]*rdx2; 
  out[21] += -1.58113883008419*GhatR[2]*rdx2; 
  out[22] += -0.7071067811865475*GhatR[11]*rdx2; 
  out[23] += -1.224744871391589*GhatR[8]*rdx2; 
  out[24] += -0.7071067811865475*GhatR[12]*rdx2; 
  out[25] += -1.58113883008419*GhatR[3]*rdx2; 
  out[26] += -0.7071067811865475*GhatR[13]*rdx2; 
  out[27] += -0.7071067811865475*GhatR[14]*rdx2; 
  out[28] += -1.224744871391589*GhatR[9]*rdx2; 
  out[29] += -0.7071067811865475*GhatR[15]*rdx2; 
  out[30] += -0.7071067811865475*GhatR[16]*rdx2; 
  out[31] += -1.224744871391589*GhatR[10]*rdx2; 
  out[32] += -1.58113883008419*GhatR[4]*rdx2; 
  out[33] += -1.224744871391589*GhatR[11]*rdx2; 
  out[34] += -1.224744871391589*GhatR[12]*rdx2; 
  out[35] += -1.58113883008419*GhatR[5]*rdx2; 
  out[36] += -1.224744871391589*GhatR[13]*rdx2; 
  out[37] += -1.58113883008419*GhatR[6]*rdx2; 
  out[38] += -0.7071067811865475*GhatR[17]*rdx2; 
  out[39] += -1.224744871391589*GhatR[14]*rdx2; 
  out[40] += -0.7071067811865475*GhatR[18]*rdx2; 
  out[41] += -1.224744871391589*GhatR[15]*rdx2; 
  out[42] += -1.224744871391589*GhatR[16]*rdx2; 
  out[43] += -0.7071067811865475*GhatR[19]*rdx2; 
  out[44] += -1.58113883008419*GhatR[10]*rdx2; 
  out[45] += -1.224744871391589*GhatR[17]*rdx2; 
  out[46] += -1.224744871391589*GhatR[18]*rdx2; 
  out[47] += -1.224744871391589*GhatR[19]*rdx2; 

  } else { 

  double alphaL[20] = {0.}; 
  alphaL[0] = -(0.1767766952966368*((((19.36491673103708*b_z[4]-15.0*b_z[1]+8.660254037844387*b_z[0])*jacobtot_inv[4]+(8.660254037844387*jacobtot_inv[0]-15.0*jacobtot_inv[1])*b_z[4]+(11.61895003862225*b_z[1]-6.708203932499369*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(3.872983346207417*b_z[0]-6.708203932499369*b_z[1]))*hamil[19]+(((-15.0*b_z[4])+11.61895003862225*b_z[1]-6.708203932499369*b_z[0])*jacobtot_inv[4]+(11.61895003862225*jacobtot_inv[1]-6.708203932499369*jacobtot_inv[0])*b_z[4]+(5.196152422706631*b_z[0]-9.0*b_z[1])*jacobtot_inv[1]+jacobtot_inv[0]*(5.196152422706631*b_z[1]-3.0*b_z[0]))*hamil[5]+hamil[2]*((8.660254037844386*b_z[4]-6.708203932499369*b_z[1]+3.872983346207417*b_z[0])*jacobtot_inv[4]+(3.872983346207417*jacobtot_inv[0]-6.708203932499369*jacobtot_inv[1])*b_z[4]+(5.196152422706631*b_z[1]-3.0*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(1.732050807568877*b_z[0]-3.0*b_z[1])))*m_*rdy2+hamil[3]*((-3.872983346207417*BstarXdBmag[11])+3.0*BstarXdBmag[1]-1.732050807568877*BstarXdBmag[0])*q_*rdvpar2))/(m_*q_); 
  alphaL[1] = (0.1767766952966368*((((33.54101966249684*b_z[4]-25.98076211353316*b_z[1]+15.0*b_z[0])*jacobtot_inv[4]+(15.0*jacobtot_inv[0]-25.98076211353316*jacobtot_inv[1])*b_z[4]+(20.1246117974981*b_z[1]-11.61895003862225*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(6.708203932499369*b_z[0]-11.61895003862225*b_z[1]))*hamil[20]+(((-19.36491673103709*b_z[4])+15.0*b_z[1]-8.660254037844386*b_z[0])*jacobtot_inv[4]+(15.0*jacobtot_inv[1]-8.660254037844386*jacobtot_inv[0])*b_z[4]+(6.708203932499369*b_z[0]-11.61895003862225*b_z[1])*jacobtot_inv[1]+jacobtot_inv[0]*(6.708203932499369*b_z[1]-3.872983346207417*b_z[0]))*hamil[12])*m_*rdy2+hamil[3]*(3.872983346207417*BstarXdBmag[19]-3.0*BstarXdBmag[5]+1.732050807568877*BstarXdBmag[2])*q_*rdvpar2))/(m_*q_); 
  alphaL[2] = (0.1767766952966368*(8.660254037844386*BstarXdBmag[11]-6.708203932499369*BstarXdBmag[1]+3.872983346207417*BstarXdBmag[0])*hamil[13]*rdvpar2)/m_; 
  alphaL[4] = (0.1767766952966368*hamil[13]*(8.660254037844387*BstarXdBmag[19]-6.708203932499369*BstarXdBmag[5]+3.872983346207417*BstarXdBmag[2])*rdvpar2)/m_; 

  double fUpOrdL[27] = {0.};
  if (alphaL[4]-0.7453559924999297*alphaL[2]-0.7453559924999297*alphaL[1]+0.555555555555555*alphaL[0] > 0.) {
    fUpOrdL[0] = ser_4x_p2_surfx1_eval_quad_node_0_r(fedge); 
  } else { 
    fUpOrdL[0] = ser_4x_p2_surfx1_eval_quad_node_0_l(fskin); 
  } 
  if (alphaL[4]-0.7453559924999297*alphaL[2]-0.7453559924999297*alphaL[1]+0.555555555555555*alphaL[0] > 0.) {
    fUpOrdL[1] = ser_4x_p2_surfx1_eval_quad_node_1_r(fedge); 
  } else { 
    fUpOrdL[1] = ser_4x_p2_surfx1_eval_quad_node_1_l(fskin); 
  } 
  if (alphaL[4]-0.7453559924999297*alphaL[2]-0.7453559924999297*alphaL[1]+0.555555555555555*alphaL[0] > 0.) {
    fUpOrdL[2] = ser_4x_p2_surfx1_eval_quad_node_2_r(fedge); 
  } else { 
    fUpOrdL[2] = ser_4x_p2_surfx1_eval_quad_node_2_l(fskin); 
  } 
  if (alphaL[0]-1.341640786499875*alphaL[1] > 0.) {
    fUpOrdL[3] = ser_4x_p2_surfx1_eval_quad_node_3_r(fedge); 
  } else { 
    fUpOrdL[3] = ser_4x_p2_surfx1_eval_quad_node_3_l(fskin); 
  } 
  if (alphaL[0]-1.341640786499875*alphaL[1] > 0.) {
    fUpOrdL[4] = ser_4x_p2_surfx1_eval_quad_node_4_r(fedge); 
  } else { 
    fUpOrdL[4] = ser_4x_p2_surfx1_eval_quad_node_4_l(fskin); 
  } 
  if (alphaL[0]-1.341640786499875*alphaL[1] > 0.) {
    fUpOrdL[5] = ser_4x_p2_surfx1_eval_quad_node_5_r(fedge); 
  } else { 
    fUpOrdL[5] = ser_4x_p2_surfx1_eval_quad_node_5_l(fskin); 
  } 
  if ((-1.0*alphaL[4])+0.7453559924999297*alphaL[2]-0.7453559924999297*alphaL[1]+0.555555555555555*alphaL[0] > 0.) {
    fUpOrdL[6] = ser_4x_p2_surfx1_eval_quad_node_6_r(fedge); 
  } else { 
    fUpOrdL[6] = ser_4x_p2_surfx1_eval_quad_node_6_l(fskin); 
  } 
  if ((-1.0*alphaL[4])+0.7453559924999297*alphaL[2]-0.7453559924999297*alphaL[1]+0.555555555555555*alphaL[0] > 0.) {
    fUpOrdL[7] = ser_4x_p2_surfx1_eval_quad_node_7_r(fedge); 
  } else { 
    fUpOrdL[7] = ser_4x_p2_surfx1_eval_quad_node_7_l(fskin); 
  } 
  if ((-1.0*alphaL[4])+0.7453559924999297*alphaL[2]-0.7453559924999297*alphaL[1]+0.555555555555555*alphaL[0] > 0.) {
    fUpOrdL[8] = ser_4x_p2_surfx1_eval_quad_node_8_r(fedge); 
  } else { 
    fUpOrdL[8] = ser_4x_p2_surfx1_eval_quad_node_8_l(fskin); 
  } 
  if (alphaL[0]-1.341640786499875*alphaL[2] > 0.) {
    fUpOrdL[9] = ser_4x_p2_surfx1_eval_quad_node_9_r(fedge); 
  } else { 
    fUpOrdL[9] = ser_4x_p2_surfx1_eval_quad_node_9_l(fskin); 
  } 
  if (alphaL[0]-1.341640786499875*alphaL[2] > 0.) {
    fUpOrdL[10] = ser_4x_p2_surfx1_eval_quad_node_10_r(fedge); 
  } else { 
    fUpOrdL[10] = ser_4x_p2_surfx1_eval_quad_node_10_l(fskin); 
  } 
  if (alphaL[0]-1.341640786499875*alphaL[2] > 0.) {
    fUpOrdL[11] = ser_4x_p2_surfx1_eval_quad_node_11_r(fedge); 
  } else { 
    fUpOrdL[11] = ser_4x_p2_surfx1_eval_quad_node_11_l(fskin); 
  } 
  if (alphaL[0] > 0.) {
    fUpOrdL[12] = ser_4x_p2_surfx1_eval_quad_node_12_r(fedge); 
  } else { 
    fUpOrdL[12] = ser_4x_p2_surfx1_eval_quad_node_12_l(fskin); 
  } 
  if (alphaL[0] > 0.) {
    fUpOrdL[13] = ser_4x_p2_surfx1_eval_quad_node_13_r(fedge); 
  } else { 
    fUpOrdL[13] = ser_4x_p2_surfx1_eval_quad_node_13_l(fskin); 
  } 
  if (alphaL[0] > 0.) {
    fUpOrdL[14] = ser_4x_p2_surfx1_eval_quad_node_14_r(fedge); 
  } else { 
    fUpOrdL[14] = ser_4x_p2_surfx1_eval_quad_node_14_l(fskin); 
  } 
  if (alphaL[2]+0.7453559924999295*alphaL[0] > 0.) {
    fUpOrdL[15] = ser_4x_p2_surfx1_eval_quad_node_15_r(fedge); 
  } else { 
    fUpOrdL[15] = ser_4x_p2_surfx1_eval_quad_node_15_l(fskin); 
  } 
  if (alphaL[2]+0.7453559924999295*alphaL[0] > 0.) {
    fUpOrdL[16] = ser_4x_p2_surfx1_eval_quad_node_16_r(fedge); 
  } else { 
    fUpOrdL[16] = ser_4x_p2_surfx1_eval_quad_node_16_l(fskin); 
  } 
  if (alphaL[2]+0.7453559924999295*alphaL[0] > 0.) {
    fUpOrdL[17] = ser_4x_p2_surfx1_eval_quad_node_17_r(fedge); 
  } else { 
    fUpOrdL[17] = ser_4x_p2_surfx1_eval_quad_node_17_l(fskin); 
  } 
  if ((-1.0*alphaL[4])-0.7453559924999297*alphaL[2]+0.7453559924999297*alphaL[1]+0.555555555555555*alphaL[0] > 0.) {
    fUpOrdL[18] = ser_4x_p2_surfx1_eval_quad_node_18_r(fedge); 
  } else { 
    fUpOrdL[18] = ser_4x_p2_surfx1_eval_quad_node_18_l(fskin); 
  } 
  if ((-1.0*alphaL[4])-0.7453559924999297*alphaL[2]+0.7453559924999297*alphaL[1]+0.555555555555555*alphaL[0] > 0.) {
    fUpOrdL[19] = ser_4x_p2_surfx1_eval_quad_node_19_r(fedge); 
  } else { 
    fUpOrdL[19] = ser_4x_p2_surfx1_eval_quad_node_19_l(fskin); 
  } 
  if ((-1.0*alphaL[4])-0.7453559924999297*alphaL[2]+0.7453559924999297*alphaL[1]+0.555555555555555*alphaL[0] > 0.) {
    fUpOrdL[20] = ser_4x_p2_surfx1_eval_quad_node_20_r(fedge); 
  } else { 
    fUpOrdL[20] = ser_4x_p2_surfx1_eval_quad_node_20_l(fskin); 
  } 
  if (alphaL[1]+0.7453559924999295*alphaL[0] > 0.) {
    fUpOrdL[21] = ser_4x_p2_surfx1_eval_quad_node_21_r(fedge); 
  } else { 
    fUpOrdL[21] = ser_4x_p2_surfx1_eval_quad_node_21_l(fskin); 
  } 
  if (alphaL[1]+0.7453559924999295*alphaL[0] > 0.) {
    fUpOrdL[22] = ser_4x_p2_surfx1_eval_quad_node_22_r(fedge); 
  } else { 
    fUpOrdL[22] = ser_4x_p2_surfx1_eval_quad_node_22_l(fskin); 
  } 
  if (alphaL[1]+0.7453559924999295*alphaL[0] > 0.) {
    fUpOrdL[23] = ser_4x_p2_surfx1_eval_quad_node_23_r(fedge); 
  } else { 
    fUpOrdL[23] = ser_4x_p2_surfx1_eval_quad_node_23_l(fskin); 
  } 
  if (alphaL[4]+0.7453559924999297*alphaL[2]+0.7453559924999297*alphaL[1]+0.555555555555555*alphaL[0] > 0.) {
    fUpOrdL[24] = ser_4x_p2_surfx1_eval_quad_node_24_r(fedge); 
  } else { 
    fUpOrdL[24] = ser_4x_p2_surfx1_eval_quad_node_24_l(fskin); 
  } 
  if (alphaL[4]+0.7453559924999297*alphaL[2]+0.7453559924999297*alphaL[1]+0.555555555555555*alphaL[0] > 0.) {
    fUpOrdL[25] = ser_4x_p2_surfx1_eval_quad_node_25_r(fedge); 
  } else { 
    fUpOrdL[25] = ser_4x_p2_surfx1_eval_quad_node_25_l(fskin); 
  } 
  if (alphaL[4]+0.7453559924999297*alphaL[2]+0.7453559924999297*alphaL[1]+0.555555555555555*alphaL[0] > 0.) {
    fUpOrdL[26] = ser_4x_p2_surfx1_eval_quad_node_26_r(fedge); 
  } else { 
    fUpOrdL[26] = ser_4x_p2_surfx1_eval_quad_node_26_l(fskin); 
  } 

  // Project tensor nodal quadrature basis back onto modal basis. 
  double fUpL[20] = {0.};
  ser_4x_p2_upwind_quad_to_modal(fUpOrdL, fUpL); 

  double GhatL[48] = {0.}; 
  GhatL[0] = 0.3535533905932737*alphaL[4]*fUpL[4]+0.3535533905932737*alphaL[2]*fUpL[2]+0.3535533905932737*alphaL[1]*fUpL[1]+0.3535533905932737*alphaL[0]*fUpL[0]; 
  GhatL[1] = 0.3162277660168379*alphaL[4]*fUpL[11]+0.3162277660168379*alphaL[1]*fUpL[7]+0.3535533905932737*alphaL[2]*fUpL[4]+0.3535533905932737*fUpL[2]*alphaL[4]+0.3535533905932737*alphaL[0]*fUpL[1]+0.3535533905932737*fUpL[0]*alphaL[1]; 
  GhatL[2] = 0.3162277660168379*alphaL[4]*fUpL[12]+0.3162277660168379*alphaL[2]*fUpL[8]+0.3535533905932737*alphaL[1]*fUpL[4]+0.3535533905932737*fUpL[1]*alphaL[4]+0.3535533905932737*alphaL[0]*fUpL[2]+0.3535533905932737*fUpL[0]*alphaL[2]; 
  GhatL[3] = 0.3535533905932737*alphaL[4]*fUpL[10]+0.3535533905932737*alphaL[2]*fUpL[6]+0.3535533905932737*alphaL[1]*fUpL[5]+0.3535533905932737*alphaL[0]*fUpL[3]; 
  GhatL[4] = 0.3162277660168379*alphaL[2]*fUpL[12]+0.3162277660168379*alphaL[1]*fUpL[11]+0.3162277660168379*alphaL[4]*fUpL[8]+0.3162277660168379*alphaL[4]*fUpL[7]+0.3535533905932737*alphaL[0]*fUpL[4]+0.3535533905932737*fUpL[0]*alphaL[4]+0.3535533905932737*alphaL[1]*fUpL[2]+0.3535533905932737*fUpL[1]*alphaL[2]; 
  GhatL[5] = 0.3162277660168379*alphaL[4]*fUpL[17]+0.3162277660168379*alphaL[1]*fUpL[13]+0.3535533905932737*alphaL[2]*fUpL[10]+0.3535533905932737*alphaL[4]*fUpL[6]+0.3535533905932737*alphaL[0]*fUpL[5]+0.3535533905932737*alphaL[1]*fUpL[3]; 
  GhatL[6] = 0.3162277660168379*alphaL[4]*fUpL[18]+0.3162277660168379*alphaL[2]*fUpL[14]+0.3535533905932737*alphaL[1]*fUpL[10]+0.3535533905932737*alphaL[0]*fUpL[6]+0.3535533905932737*alphaL[4]*fUpL[5]+0.3535533905932737*alphaL[2]*fUpL[3]; 
  GhatL[7] = 0.3535533905932737*alphaL[2]*fUpL[11]+0.3535533905932737*alphaL[0]*fUpL[7]+0.3162277660168379*alphaL[4]*fUpL[4]+0.3162277660168379*alphaL[1]*fUpL[1]; 
  GhatL[8] = 0.3535533905932737*alphaL[1]*fUpL[12]+0.3535533905932737*alphaL[0]*fUpL[8]+0.3162277660168379*alphaL[4]*fUpL[4]+0.3162277660168379*alphaL[2]*fUpL[2]; 
  GhatL[9] = 0.3535533905932737*alphaL[4]*fUpL[19]+0.3535533905932737*alphaL[2]*fUpL[16]+0.3535533905932737*alphaL[1]*fUpL[15]+0.3535533905932737*alphaL[0]*fUpL[9]; 
  GhatL[10] = 0.3162277660168379*alphaL[2]*fUpL[18]+0.3162277660168379*alphaL[1]*fUpL[17]+0.3162277660168379*alphaL[4]*fUpL[14]+0.3162277660168379*alphaL[4]*fUpL[13]+0.3535533905932737*alphaL[0]*fUpL[10]+0.3535533905932737*alphaL[1]*fUpL[6]+0.3535533905932737*alphaL[2]*fUpL[5]+0.3535533905932737*fUpL[3]*alphaL[4]; 
  GhatL[11] = 0.2828427124746191*alphaL[4]*fUpL[12]+0.3535533905932737*alphaL[0]*fUpL[11]+0.3535533905932737*alphaL[2]*fUpL[7]+0.3162277660168379*alphaL[1]*fUpL[4]+0.3162277660168379*fUpL[1]*alphaL[4]; 
  GhatL[12] = 0.3535533905932737*alphaL[0]*fUpL[12]+0.2828427124746191*alphaL[4]*fUpL[11]+0.3535533905932737*alphaL[1]*fUpL[8]+0.3162277660168379*alphaL[2]*fUpL[4]+0.3162277660168379*fUpL[2]*alphaL[4]; 
  GhatL[13] = 0.3535533905932737*alphaL[2]*fUpL[17]+0.3535533905932737*alphaL[0]*fUpL[13]+0.3162277660168379*alphaL[4]*fUpL[10]+0.3162277660168379*alphaL[1]*fUpL[5]; 
  GhatL[14] = 0.3535533905932737*alphaL[1]*fUpL[18]+0.3535533905932737*alphaL[0]*fUpL[14]+0.3162277660168379*alphaL[4]*fUpL[10]+0.3162277660168379*alphaL[2]*fUpL[6]; 
  GhatL[15] = 0.3535533905932737*alphaL[2]*fUpL[19]+0.3535533905932737*alphaL[4]*fUpL[16]+0.3535533905932737*alphaL[0]*fUpL[15]+0.3535533905932737*alphaL[1]*fUpL[9]; 
  GhatL[16] = 0.3535533905932737*alphaL[1]*fUpL[19]+0.3535533905932737*alphaL[0]*fUpL[16]+0.3535533905932737*alphaL[4]*fUpL[15]+0.3535533905932737*alphaL[2]*fUpL[9]; 
  GhatL[17] = 0.2828427124746191*alphaL[4]*fUpL[18]+0.3535533905932737*alphaL[0]*fUpL[17]+0.3535533905932737*alphaL[2]*fUpL[13]+0.3162277660168379*alphaL[1]*fUpL[10]+0.3162277660168379*alphaL[4]*fUpL[5]; 
  GhatL[18] = 0.3535533905932737*alphaL[0]*fUpL[18]+0.2828427124746191*alphaL[4]*fUpL[17]+0.3535533905932737*alphaL[1]*fUpL[14]+0.3162277660168379*alphaL[2]*fUpL[10]+0.3162277660168379*alphaL[4]*fUpL[6]; 
  GhatL[19] = 0.3535533905932737*alphaL[0]*fUpL[19]+0.3535533905932737*alphaL[1]*fUpL[16]+0.3535533905932737*alphaL[2]*fUpL[15]+0.3535533905932737*alphaL[4]*fUpL[9]; 

  out[0] += 0.7071067811865475*GhatL[0]*rdx2; 
  out[1] += -1.224744871391589*GhatL[0]*rdx2; 
  out[2] += 0.7071067811865475*GhatL[1]*rdx2; 
  out[3] += 0.7071067811865475*GhatL[2]*rdx2; 
  out[4] += 0.7071067811865475*GhatL[3]*rdx2; 
  out[5] += -1.224744871391589*GhatL[1]*rdx2; 
  out[6] += -1.224744871391589*GhatL[2]*rdx2; 
  out[7] += 0.7071067811865475*GhatL[4]*rdx2; 
  out[8] += -1.224744871391589*GhatL[3]*rdx2; 
  out[9] += 0.7071067811865475*GhatL[5]*rdx2; 
  out[10] += 0.7071067811865475*GhatL[6]*rdx2; 
  out[11] += 1.58113883008419*GhatL[0]*rdx2; 
  out[12] += 0.7071067811865475*GhatL[7]*rdx2; 
  out[13] += 0.7071067811865475*GhatL[8]*rdx2; 
  out[14] += 0.7071067811865475*GhatL[9]*rdx2; 
  out[15] += -1.224744871391589*GhatL[4]*rdx2; 
  out[16] += -1.224744871391589*GhatL[5]*rdx2; 
  out[17] += -1.224744871391589*GhatL[6]*rdx2; 
  out[18] += 0.7071067811865475*GhatL[10]*rdx2; 
  out[19] += 1.58113883008419*GhatL[1]*rdx2; 
  out[20] += -1.224744871391589*GhatL[7]*rdx2; 
  out[21] += 1.58113883008419*GhatL[2]*rdx2; 
  out[22] += 0.7071067811865475*GhatL[11]*rdx2; 
  out[23] += -1.224744871391589*GhatL[8]*rdx2; 
  out[24] += 0.7071067811865475*GhatL[12]*rdx2; 
  out[25] += 1.58113883008419*GhatL[3]*rdx2; 
  out[26] += 0.7071067811865475*GhatL[13]*rdx2; 
  out[27] += 0.7071067811865475*GhatL[14]*rdx2; 
  out[28] += -1.224744871391589*GhatL[9]*rdx2; 
  out[29] += 0.7071067811865475*GhatL[15]*rdx2; 
  out[30] += 0.7071067811865475*GhatL[16]*rdx2; 
  out[31] += -1.224744871391589*GhatL[10]*rdx2; 
  out[32] += 1.58113883008419*GhatL[4]*rdx2; 
  out[33] += -1.224744871391589*GhatL[11]*rdx2; 
  out[34] += -1.224744871391589*GhatL[12]*rdx2; 
  out[35] += 1.58113883008419*GhatL[5]*rdx2; 
  out[36] += -1.224744871391589*GhatL[13]*rdx2; 
  out[37] += 1.58113883008419*GhatL[6]*rdx2; 
  out[38] += 0.7071067811865475*GhatL[17]*rdx2; 
  out[39] += -1.224744871391589*GhatL[14]*rdx2; 
  out[40] += 0.7071067811865475*GhatL[18]*rdx2; 
  out[41] += -1.224744871391589*GhatL[15]*rdx2; 
  out[42] += -1.224744871391589*GhatL[16]*rdx2; 
  out[43] += 0.7071067811865475*GhatL[19]*rdx2; 
  out[44] += 1.58113883008419*GhatL[10]*rdx2; 
  out[45] += -1.224744871391589*GhatL[17]*rdx2; 
  out[46] += -1.224744871391589*GhatL[18]*rdx2; 
  out[47] += -1.224744871391589*GhatL[19]*rdx2; 

  } 

} 
