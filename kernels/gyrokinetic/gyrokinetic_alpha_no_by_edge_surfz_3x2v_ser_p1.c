#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfz_3x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, 
  const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
  const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // jacobtot_inv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // cmag: coefficient multiplying parallel gradient.
  // b_i: covariant components of the field aligned unit vector.
  // phi: electrostatic potential.
  // alpha_surf: output surface phase space flux in each direction (cdim + 1 components).
  //             Note: Each cell owns their *lower* edge surface evaluation.
  // sgn_alpha_surf: output sign(alpha_surf) in each direction at quadrature points (cdim + 1 components).
  //                 Note: Each cell owns their *lower* edge sign(alpha_surf).
  // returns int const_sgn_alpha (true if sign(alpha_surf) is only one sign, either +1 or -1).

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wy = w[1];
  double rdy2 = 2.0/dxv[1];
  double wz = w[2];
  double rdz2 = 2.0/dxv[2];
  double wvpar = w[3];
  double rdvpar2 = 2.0/dxv[3];
  double wmu = w[4];
  double rdmu2 = 2.0/dxv[4];

  double wvparSq = wvpar*wvpar;
  double rdvpar2Sq = rdvpar2*rdvpar2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[8];
  const double *b_z = &b_i[16];

  double hamil[48] = {0.}; 
  hamil[0] = 2.828427124746191*m_*wvparSq+2.0*bmag[0]*wmu+(0.9428090415820636*m_)/rdvpar2Sq+2.0*phi[0]*q_; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = 2.0*(bmag[3]*wmu+phi[3]*q_); 
  hamil[4] = (3.265986323710906*m_*wvpar)/rdvpar2; 
  hamil[5] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[6] = 2.0*phi[4]*q_; 
  hamil[7] = 2.0*(bmag[5]*wmu+phi[5]*q_); 
  hamil[8] = 2.0*phi[6]*q_; 
  hamil[12] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[14] = (1.154700538379252*bmag[3])/rdmu2; 
  hamil[16] = 2.0*phi[7]*q_; 
  hamil[21] = (1.154700538379252*bmag[5])/rdmu2; 
  hamil[32] = (0.8432740427115681*m_)/rdvpar2Sq; 

  double *alphaR = &alpha_surf[48];
  double *sgn_alpha_surfR = &sgn_alpha_surf[48];
  alphaR[0] = ((0.7954951288348656*b_x[3]*jacobtot_inv[5]*hamil[16]+0.4592793267718456*b_x[0]*jacobtot_inv[5]*hamil[16]+0.7954951288348656*jacobtot_inv[3]*b_x[5]*hamil[16]+0.4592793267718456*jacobtot_inv[0]*b_x[5]*hamil[16]+0.4592793267718456*b_x[1]*jacobtot_inv[3]*hamil[16]+0.4592793267718456*jacobtot_inv[1]*b_x[3]*hamil[16]+0.2651650429449552*b_x[0]*jacobtot_inv[1]*hamil[16]+0.2651650429449552*jacobtot_inv[0]*b_x[1]*hamil[16]+0.7954951288348656*b_x[5]*jacobtot_inv[5]*hamil[8]+0.4592793267718456*b_x[1]*jacobtot_inv[5]*hamil[8]+0.4592793267718456*jacobtot_inv[1]*b_x[5]*hamil[8]+0.7954951288348656*b_x[3]*jacobtot_inv[3]*hamil[8]+0.4592793267718456*b_x[0]*jacobtot_inv[3]*hamil[8]+0.4592793267718456*jacobtot_inv[0]*b_x[3]*hamil[8]+0.2651650429449552*b_x[1]*jacobtot_inv[1]*hamil[8]+0.2651650429449552*b_x[0]*jacobtot_inv[0]*hamil[8]+0.4592793267718456*b_x[3]*jacobtot_inv[5]*hamil[6]+0.2651650429449552*b_x[0]*jacobtot_inv[5]*hamil[6]+0.4592793267718456*jacobtot_inv[3]*b_x[5]*hamil[6]+0.2651650429449552*jacobtot_inv[0]*b_x[5]*hamil[6]+0.2651650429449552*b_x[1]*jacobtot_inv[3]*hamil[6]+0.2651650429449552*jacobtot_inv[1]*b_x[3]*hamil[6]+0.1530931089239486*b_x[0]*jacobtot_inv[1]*hamil[6]+0.1530931089239486*jacobtot_inv[0]*b_x[1]*hamil[6]+0.4592793267718456*hamil[2]*b_x[5]*jacobtot_inv[5]+0.2651650429449552*b_x[1]*hamil[2]*jacobtot_inv[5]+0.2651650429449552*jacobtot_inv[1]*hamil[2]*b_x[5]+0.4592793267718456*hamil[2]*b_x[3]*jacobtot_inv[3]+0.2651650429449552*b_x[0]*hamil[2]*jacobtot_inv[3]+0.2651650429449552*jacobtot_inv[0]*hamil[2]*b_x[3]+0.1530931089239486*b_x[1]*jacobtot_inv[1]*hamil[2]+0.1530931089239486*b_x[0]*jacobtot_inv[0]*hamil[2])*rdy2)/q_+((0.4592793267718456*hamil[4]*cmag[5]*jacobtot_inv[5]+0.2651650429449552*cmag[1]*hamil[4]*jacobtot_inv[5]+0.2651650429449552*jacobtot_inv[1]*hamil[4]*cmag[5]+0.4592793267718456*cmag[3]*jacobtot_inv[3]*hamil[4]+0.2651650429449552*cmag[0]*jacobtot_inv[3]*hamil[4]+0.2651650429449552*jacobtot_inv[0]*cmag[3]*hamil[4]+0.1530931089239486*cmag[1]*jacobtot_inv[1]*hamil[4]+0.1530931089239486*cmag[0]*jacobtot_inv[0]*hamil[4])*rdvpar2)/m_; 
  alphaR[1] = ((1.431891231902758*b_x[5]*jacobtot_inv[5]*hamil[16]+0.826702788189322*b_x[1]*jacobtot_inv[5]*hamil[16]+0.826702788189322*jacobtot_inv[1]*b_x[5]*hamil[16]+0.7954951288348656*b_x[3]*jacobtot_inv[3]*hamil[16]+0.4592793267718456*b_x[0]*jacobtot_inv[3]*hamil[16]+0.4592793267718456*jacobtot_inv[0]*b_x[3]*hamil[16]+0.4772970773009194*b_x[1]*jacobtot_inv[1]*hamil[16]+0.2651650429449552*b_x[0]*jacobtot_inv[0]*hamil[16]+0.7954951288348656*b_x[3]*jacobtot_inv[5]*hamil[8]+0.4592793267718456*b_x[0]*jacobtot_inv[5]*hamil[8]+0.7954951288348656*jacobtot_inv[3]*b_x[5]*hamil[8]+0.4592793267718456*jacobtot_inv[0]*b_x[5]*hamil[8]+0.4592793267718456*b_x[1]*jacobtot_inv[3]*hamil[8]+0.4592793267718456*jacobtot_inv[1]*b_x[3]*hamil[8]+0.2651650429449552*b_x[0]*jacobtot_inv[1]*hamil[8]+0.2651650429449552*jacobtot_inv[0]*b_x[1]*hamil[8]+0.826702788189322*b_x[5]*jacobtot_inv[5]*hamil[6]+0.4772970773009194*b_x[1]*jacobtot_inv[5]*hamil[6]+0.4772970773009194*jacobtot_inv[1]*b_x[5]*hamil[6]+0.4592793267718456*b_x[3]*jacobtot_inv[3]*hamil[6]+0.2651650429449552*b_x[0]*jacobtot_inv[3]*hamil[6]+0.2651650429449552*jacobtot_inv[0]*b_x[3]*hamil[6]+0.2755675960631073*b_x[1]*jacobtot_inv[1]*hamil[6]+0.1530931089239486*b_x[0]*jacobtot_inv[0]*hamil[6]+0.4592793267718456*hamil[2]*b_x[3]*jacobtot_inv[5]+0.2651650429449552*b_x[0]*hamil[2]*jacobtot_inv[5]+0.4592793267718456*hamil[2]*jacobtot_inv[3]*b_x[5]+0.2651650429449552*jacobtot_inv[0]*hamil[2]*b_x[5]+0.2651650429449552*b_x[1]*hamil[2]*jacobtot_inv[3]+0.2651650429449552*jacobtot_inv[1]*hamil[2]*b_x[3]+0.1530931089239486*b_x[0]*jacobtot_inv[1]*hamil[2]+0.1530931089239486*jacobtot_inv[0]*b_x[1]*hamil[2])*rdy2)/q_+((0.4592793267718456*cmag[3]*hamil[4]*jacobtot_inv[5]+0.2651650429449552*cmag[0]*hamil[4]*jacobtot_inv[5]+0.4592793267718456*jacobtot_inv[3]*hamil[4]*cmag[5]+0.2651650429449552*jacobtot_inv[0]*hamil[4]*cmag[5]+0.2651650429449552*cmag[1]*jacobtot_inv[3]*hamil[4]+0.2651650429449552*jacobtot_inv[1]*cmag[3]*hamil[4]+0.1530931089239486*cmag[0]*jacobtot_inv[1]*hamil[4]+0.1530931089239486*jacobtot_inv[0]*cmag[1]*hamil[4])*rdvpar2)/m_; 
  alphaR[3] = ((1.026979795322186*cmag[5]*jacobtot_inv[5]*hamil[32]+0.592927061281571*cmag[1]*jacobtot_inv[5]*hamil[32]+0.592927061281571*jacobtot_inv[1]*cmag[5]*hamil[32]+1.026979795322186*cmag[3]*jacobtot_inv[3]*hamil[32]+0.592927061281571*cmag[0]*jacobtot_inv[3]*hamil[32]+0.592927061281571*jacobtot_inv[0]*cmag[3]*hamil[32]+0.3423265984407287*cmag[1]*jacobtot_inv[1]*hamil[32]+0.3423265984407287*cmag[0]*jacobtot_inv[0]*hamil[32])*rdvpar2)/m_; 
  alphaR[6] = ((1.026979795322186*cmag[3]*jacobtot_inv[5]*hamil[32]+0.592927061281571*cmag[0]*jacobtot_inv[5]*hamil[32]+1.026979795322186*jacobtot_inv[3]*cmag[5]*hamil[32]+0.592927061281571*jacobtot_inv[0]*cmag[5]*hamil[32]+0.592927061281571*cmag[1]*jacobtot_inv[3]*hamil[32]+0.592927061281571*jacobtot_inv[1]*cmag[3]*hamil[32]+0.3423265984407287*cmag[0]*jacobtot_inv[1]*hamil[32]+0.3423265984407287*jacobtot_inv[0]*cmag[1]*hamil[32])*rdvpar2)/m_; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.3354101966249678*alphaR[6]-0.3354101966249678*alphaR[3]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[0] = 1.0; 
  else  
    sgn_alpha_surfR[0] = -1.0; 
  
  if (0.25*alphaR[0]-0.25*alphaR[1] > 0.) 
    sgn_alpha_surfR[1] = 1.0; 
  else  
    sgn_alpha_surfR[1] = -1.0; 
  
  if (sgn_alpha_surfR[1] == sgn_alpha_surfR[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3354101966249678*alphaR[6])+0.3354101966249678*alphaR[3]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[2] = 1.0; 
  else  
    sgn_alpha_surfR[2] = -1.0; 
  
  if (sgn_alpha_surfR[2] == sgn_alpha_surfR[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3354101966249678*alphaR[6]-0.3354101966249678*alphaR[3]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[3] = 1.0; 
  else  
    sgn_alpha_surfR[3] = -1.0; 
  
  if (sgn_alpha_surfR[3] == sgn_alpha_surfR[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*alphaR[0]-0.25*alphaR[1] > 0.) 
    sgn_alpha_surfR[4] = 1.0; 
  else  
    sgn_alpha_surfR[4] = -1.0; 
  
  if (sgn_alpha_surfR[4] == sgn_alpha_surfR[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3354101966249678*alphaR[6])+0.3354101966249678*alphaR[3]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[5] = 1.0; 
  else  
    sgn_alpha_surfR[5] = -1.0; 
  
  if (sgn_alpha_surfR[5] == sgn_alpha_surfR[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3354101966249678*alphaR[6]-0.3354101966249678*alphaR[3]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[6] = 1.0; 
  else  
    sgn_alpha_surfR[6] = -1.0; 
  
  if (sgn_alpha_surfR[6] == sgn_alpha_surfR[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*alphaR[0]-0.25*alphaR[1] > 0.) 
    sgn_alpha_surfR[7] = 1.0; 
  else  
    sgn_alpha_surfR[7] = -1.0; 
  
  if (sgn_alpha_surfR[7] == sgn_alpha_surfR[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3354101966249678*alphaR[6])+0.3354101966249678*alphaR[3]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[8] = 1.0; 
  else  
    sgn_alpha_surfR[8] = -1.0; 
  
  if (sgn_alpha_surfR[8] == sgn_alpha_surfR[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3354101966249678*alphaR[6]-0.3354101966249678*alphaR[3]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[9] = 1.0; 
  else  
    sgn_alpha_surfR[9] = -1.0; 
  
  if (sgn_alpha_surfR[9] == sgn_alpha_surfR[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*alphaR[0]-0.25*alphaR[1] > 0.) 
    sgn_alpha_surfR[10] = 1.0; 
  else  
    sgn_alpha_surfR[10] = -1.0; 
  
  if (sgn_alpha_surfR[10] == sgn_alpha_surfR[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3354101966249678*alphaR[6])+0.3354101966249678*alphaR[3]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[11] = 1.0; 
  else  
    sgn_alpha_surfR[11] = -1.0; 
  
  if (sgn_alpha_surfR[11] == sgn_alpha_surfR[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaR[1]+alphaR[0])-0.3354101966249678*(alphaR[6]+alphaR[3]) > 0.) 
    sgn_alpha_surfR[12] = 1.0; 
  else  
    sgn_alpha_surfR[12] = -1.0; 
  
  if (sgn_alpha_surfR[12] == sgn_alpha_surfR[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[13] = 1.0; 
  else  
    sgn_alpha_surfR[13] = -1.0; 
  
  if (sgn_alpha_surfR[13] == sgn_alpha_surfR[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3354101966249678*(alphaR[6]+alphaR[3])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[14] = 1.0; 
  else  
    sgn_alpha_surfR[14] = -1.0; 
  
  if (sgn_alpha_surfR[14] == sgn_alpha_surfR[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaR[1]+alphaR[0])-0.3354101966249678*(alphaR[6]+alphaR[3]) > 0.) 
    sgn_alpha_surfR[15] = 1.0; 
  else  
    sgn_alpha_surfR[15] = -1.0; 
  
  if (sgn_alpha_surfR[15] == sgn_alpha_surfR[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[16] = 1.0; 
  else  
    sgn_alpha_surfR[16] = -1.0; 
  
  if (sgn_alpha_surfR[16] == sgn_alpha_surfR[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3354101966249678*(alphaR[6]+alphaR[3])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[17] = 1.0; 
  else  
    sgn_alpha_surfR[17] = -1.0; 
  
  if (sgn_alpha_surfR[17] == sgn_alpha_surfR[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaR[1]+alphaR[0])-0.3354101966249678*(alphaR[6]+alphaR[3]) > 0.) 
    sgn_alpha_surfR[18] = 1.0; 
  else  
    sgn_alpha_surfR[18] = -1.0; 
  
  if (sgn_alpha_surfR[18] == sgn_alpha_surfR[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[19] = 1.0; 
  else  
    sgn_alpha_surfR[19] = -1.0; 
  
  if (sgn_alpha_surfR[19] == sgn_alpha_surfR[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3354101966249678*(alphaR[6]+alphaR[3])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[20] = 1.0; 
  else  
    sgn_alpha_surfR[20] = -1.0; 
  
  if (sgn_alpha_surfR[20] == sgn_alpha_surfR[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaR[1]+alphaR[0])-0.3354101966249678*(alphaR[6]+alphaR[3]) > 0.) 
    sgn_alpha_surfR[21] = 1.0; 
  else  
    sgn_alpha_surfR[21] = -1.0; 
  
  if (sgn_alpha_surfR[21] == sgn_alpha_surfR[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[22] = 1.0; 
  else  
    sgn_alpha_surfR[22] = -1.0; 
  
  if (sgn_alpha_surfR[22] == sgn_alpha_surfR[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3354101966249678*(alphaR[6]+alphaR[3])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[23] = 1.0; 
  else  
    sgn_alpha_surfR[23] = -1.0; 
  
  if (sgn_alpha_surfR[23] == sgn_alpha_surfR[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
