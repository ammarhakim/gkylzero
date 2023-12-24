#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfz_3x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, 
  const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
  const double *phi, const double *Bstar_Bmag, 
  double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // jacobtot_inv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // cmag: coefficient multiplying parallel gradient.
  // b_i: covariant components of the field aligned unit vector.
  // phi: electrostatic potential.
  // Bstar_Bmag: Bstar/Bmag volume expansion, pre-computed time-independent part.
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

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wySq = w[1]*w[1];
  double rdy2Sq = rdy2*rdy2;
  double wzSq = w[2]*w[2];
  double rdz2Sq = rdz2*rdz2;
  double wvparSq = w[3]*w[3];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[4]*w[4];
  double rdmu2Sq = rdmu2*rdmu2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[8];
  const double *b_z = &b_i[16];

  const double *BstarZdBmag = &Bstar_Bmag[16]; 

  double hamil[48] = {0.}; 
  hamil[0] = 2.828427124746191*m_*wvparSq+2.0*bmag[0]*wmu+(0.9428090415820636*m_)/rdvpar2Sq+2.0*phi[0]*q_; 
  hamil[1] = 2.0*bmag[1]*wmu+2.0*phi[1]*q_; 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = 2.0*bmag[3]*wmu+2.0*phi[3]*q_; 
  hamil[4] = (3.265986323710906*m_*wvpar)/rdvpar2; 
  hamil[5] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[6] = 2.0*phi[4]*q_; 
  hamil[7] = 2.0*bmag[5]*wmu+2.0*phi[5]*q_; 
  hamil[8] = 2.0*phi[6]*q_; 
  hamil[12] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[14] = (1.154700538379252*bmag[3])/rdmu2; 
  hamil[16] = 2.0*phi[7]*q_; 
  hamil[21] = (1.154700538379252*bmag[5])/rdmu2; 
  hamil[32] = (0.8432740427115681*m_)/rdvpar2Sq; 

  double *alphaR2 = &alpha_surf[48];
  alphaR2[0] = (((0.4592793267718456*(b_x[0]*jacobtot_inv[5]+jacobtot_inv[0]*b_x[5]+b_x[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_x[3])+0.7954951288348656*(b_x[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_x[5])+0.2651650429449552*(b_x[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_x[1]))*hamil[16]+(0.4592793267718456*(b_x[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_x[5]+b_x[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_x[3])+0.7954951288348656*(b_x[5]*jacobtot_inv[5]+b_x[3]*jacobtot_inv[3])+0.2651650429449552*(b_x[1]*jacobtot_inv[1]+b_x[0]*jacobtot_inv[0]))*hamil[8]+(0.2651650429449552*(b_x[0]*jacobtot_inv[5]+jacobtot_inv[0]*b_x[5]+b_x[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_x[3])+0.4592793267718456*(b_x[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_x[5])+0.1530931089239486*(b_x[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_x[1]))*hamil[6]+hamil[2]*(0.2651650429449552*(b_x[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_x[5]+b_x[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_x[3])+0.4592793267718456*(b_x[5]*jacobtot_inv[5]+b_x[3]*jacobtot_inv[3])+0.1530931089239486*(b_x[1]*jacobtot_inv[1]+b_x[0]*jacobtot_inv[0])))*rdy2+(((-0.4592793267718456*(b_y[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_y[5]+b_y[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_y[3]))-0.7954951288348656*(b_y[5]*jacobtot_inv[5]+b_y[3]*jacobtot_inv[3])-0.2651650429449552*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[7]+hamil[1]*((-0.2651650429449552*(b_y[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_y[5]+b_y[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_y[3]))-0.4592793267718456*(b_y[5]*jacobtot_inv[5]+b_y[3]*jacobtot_inv[3])-0.1530931089239486*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0])))*rdx2)/q_+(((1.677050983124842*BstarZdBmag[6]+0.9682458365518543*BstarZdBmag[3])*hamil[32]+(0.75*BstarZdBmag[2]+0.4330127018922193*BstarZdBmag[0])*hamil[4])*rdvpar2)/m_; 
  alphaR2[1] = (((0.826702788189322*(b_x[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_x[5])+1.431891231902758*b_x[5]*jacobtot_inv[5]+0.4592793267718456*(b_x[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_x[3])+0.7954951288348656*b_x[3]*jacobtot_inv[3]+0.4772970773009194*b_x[1]*jacobtot_inv[1]+0.2651650429449552*b_x[0]*jacobtot_inv[0])*hamil[16]+(0.4592793267718456*(b_x[0]*jacobtot_inv[5]+jacobtot_inv[0]*b_x[5]+b_x[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_x[3])+0.7954951288348656*(b_x[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_x[5])+0.2651650429449552*(b_x[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_x[1]))*hamil[8]+(0.4772970773009194*(b_x[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_x[5])+0.826702788189322*b_x[5]*jacobtot_inv[5]+0.2651650429449552*(b_x[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_x[3])+0.4592793267718456*b_x[3]*jacobtot_inv[3]+0.2755675960631073*b_x[1]*jacobtot_inv[1]+0.1530931089239486*b_x[0]*jacobtot_inv[0])*hamil[6]+hamil[2]*(0.2651650429449552*(b_x[0]*jacobtot_inv[5]+jacobtot_inv[0]*b_x[5]+b_x[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_x[3])+0.4592793267718456*(b_x[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_x[5])+0.1530931089239486*(b_x[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_x[1])))*rdy2+(((-0.4592793267718456*(b_y[0]*jacobtot_inv[5]+jacobtot_inv[0]*b_y[5]+b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3]))-0.7954951288348656*(b_y[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_y[5])-0.2651650429449552*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[7]+hamil[1]*((-0.2651650429449552*(b_y[0]*jacobtot_inv[5]+jacobtot_inv[0]*b_y[5]+b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3]))-0.4592793267718456*(b_y[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_y[5])-0.1530931089239486*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1])))*rdx2)/q_+(((1.677050983124842*BstarZdBmag[7]+0.9682458365518543*BstarZdBmag[5])*hamil[32]+(0.75*BstarZdBmag[4]+0.4330127018922193*BstarZdBmag[1])*hamil[4])*rdvpar2)/m_; 
  alphaR2[2] = ((((-0.4592793267718456*(b_y[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_y[5]+b_y[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_y[3]))-0.7954951288348656*(b_y[5]*jacobtot_inv[5]+b_y[3]*jacobtot_inv[3])-0.2651650429449552*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[16]+((-0.2651650429449552*(b_y[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_y[5]+b_y[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_y[3]))-0.4592793267718456*(b_y[5]*jacobtot_inv[5]+b_y[3]*jacobtot_inv[3])-0.1530931089239486*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[6])*rdx2)/q_; 
  alphaR2[3] = (((1.677050983124842*BstarZdBmag[2]+0.9682458365518543*BstarZdBmag[0])*hamil[32]+hamil[4]*(0.75*BstarZdBmag[6]+0.4330127018922193*BstarZdBmag[3]))*rdvpar2)/m_; 
  alphaR2[4] = ((((-0.4592793267718456*(b_y[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_y[5]+b_y[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_y[3]))-0.7954951288348656*(b_y[5]*jacobtot_inv[5]+b_y[3]*jacobtot_inv[3])-0.2651650429449552*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[21]+((-0.2651650429449552*(b_y[1]*jacobtot_inv[5]+jacobtot_inv[1]*b_y[5]+b_y[0]*jacobtot_inv[3]+jacobtot_inv[0]*b_y[3]))-0.4592793267718456*(b_y[5]*jacobtot_inv[5]+b_y[3]*jacobtot_inv[3])-0.1530931089239486*(b_y[1]*jacobtot_inv[1]+b_y[0]*jacobtot_inv[0]))*hamil[12])*rdx2)/q_; 
  alphaR2[5] = ((((-0.4592793267718456*(b_y[0]*jacobtot_inv[5]+jacobtot_inv[0]*b_y[5]+b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3]))-0.7954951288348656*(b_y[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_y[5])-0.2651650429449552*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[16]+((-0.2651650429449552*(b_y[0]*jacobtot_inv[5]+jacobtot_inv[0]*b_y[5]+b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3]))-0.4592793267718456*(b_y[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_y[5])-0.1530931089239486*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[6])*rdx2)/q_; 
  alphaR2[6] = (((1.677050983124842*BstarZdBmag[4]+0.9682458365518543*BstarZdBmag[1])*hamil[32]+hamil[4]*(0.75*BstarZdBmag[7]+0.4330127018922193*BstarZdBmag[5]))*rdvpar2)/m_; 
  alphaR2[8] = ((((-0.4592793267718456*(b_y[0]*jacobtot_inv[5]+jacobtot_inv[0]*b_y[5]+b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3]))-0.7954951288348656*(b_y[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_y[5])-0.2651650429449552*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[21]+((-0.2651650429449552*(b_y[0]*jacobtot_inv[5]+jacobtot_inv[0]*b_y[5]+b_y[1]*jacobtot_inv[3]+jacobtot_inv[1]*b_y[3]))-0.4592793267718456*(b_y[3]*jacobtot_inv[5]+jacobtot_inv[3]*b_y[5])-0.1530931089239486*(b_y[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_y[1]))*hamil[12])*rdx2)/q_; 
  alphaR2[16] = ((1.5*BstarZdBmag[6]+0.8660254037844386*BstarZdBmag[3])*hamil[32]*rdvpar2)/m_; 
  alphaR2[17] = ((1.5*BstarZdBmag[7]+0.8660254037844387*BstarZdBmag[5])*hamil[32]*rdvpar2)/m_; 

  double *sgn_alpha_surf2 = &sgn_alpha_surf[48];
  int const_sgn_alpha_surf = 1; 
  
  if ((-0.2236067977499786*alphaR2[17])+0.2236067977499786*alphaR2[16]+0.25*alphaR2[8]+0.3354101966249678*alphaR2[6]+0.25*alphaR2[5]-0.25*alphaR2[4]-0.3354101966249678*alphaR2[3]-0.25*(alphaR2[2]+alphaR2[1])+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[0] = 1.0; 
  else  
    sgn_alpha_surf2[0] = -1.0; 
  
  if (0.2795084971874732*alphaR2[17]-0.2795084971874732*alphaR2[16]+0.25*alphaR2[8]+0.25*alphaR2[5]-0.25*alphaR2[4]-0.25*alphaR2[2]-0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[1] = 1.0; 
  else  
    sgn_alpha_surf2[1] = -1.0; 
  
  if (sgn_alpha_surf2[1] == sgn_alpha_surf2[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR2[17])+0.2236067977499786*alphaR2[16]+0.25*alphaR2[8]-0.3354101966249678*alphaR2[6]+0.25*alphaR2[5]-0.25*alphaR2[4]+0.3354101966249678*alphaR2[3]-0.25*alphaR2[2]-0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[2] = 1.0; 
  else  
    sgn_alpha_surf2[2] = -1.0; 
  
  if (sgn_alpha_surf2[2] == sgn_alpha_surf2[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR2[17])+0.2236067977499786*alphaR2[16]-0.25*alphaR2[8]+0.3354101966249678*alphaR2[6]+0.25*alphaR2[5]+0.25*alphaR2[4]-0.3354101966249678*alphaR2[3]-0.25*alphaR2[2]-0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[3] = 1.0; 
  else  
    sgn_alpha_surf2[3] = -1.0; 
  
  if (sgn_alpha_surf2[3] == sgn_alpha_surf2[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaR2[17]-0.2795084971874732*alphaR2[16]-0.25*alphaR2[8]+0.25*alphaR2[5]+0.25*alphaR2[4]-0.25*alphaR2[2]-0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[4] = 1.0; 
  else  
    sgn_alpha_surf2[4] = -1.0; 
  
  if (sgn_alpha_surf2[4] == sgn_alpha_surf2[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR2[17])+0.2236067977499786*alphaR2[16]-0.25*alphaR2[8]-0.3354101966249678*alphaR2[6]+0.25*alphaR2[5]+0.25*alphaR2[4]+0.3354101966249678*alphaR2[3]-0.25*alphaR2[2]-0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[5] = 1.0; 
  else  
    sgn_alpha_surf2[5] = -1.0; 
  
  if (sgn_alpha_surf2[5] == sgn_alpha_surf2[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR2[17])+0.2236067977499786*alphaR2[16]+0.25*alphaR2[8]+0.3354101966249678*alphaR2[6]-0.25*alphaR2[5]-0.25*alphaR2[4]-0.3354101966249678*alphaR2[3]+0.25*alphaR2[2]-0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[6] = 1.0; 
  else  
    sgn_alpha_surf2[6] = -1.0; 
  
  if (sgn_alpha_surf2[6] == sgn_alpha_surf2[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaR2[17]-0.2795084971874732*alphaR2[16]+0.25*alphaR2[8]-0.25*alphaR2[5]-0.25*alphaR2[4]+0.25*alphaR2[2]-0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[7] = 1.0; 
  else  
    sgn_alpha_surf2[7] = -1.0; 
  
  if (sgn_alpha_surf2[7] == sgn_alpha_surf2[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR2[17])+0.2236067977499786*alphaR2[16]+0.25*alphaR2[8]-0.3354101966249678*alphaR2[6]-0.25*alphaR2[5]-0.25*alphaR2[4]+0.3354101966249678*alphaR2[3]+0.25*alphaR2[2]-0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[8] = 1.0; 
  else  
    sgn_alpha_surf2[8] = -1.0; 
  
  if (sgn_alpha_surf2[8] == sgn_alpha_surf2[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR2[17])+0.2236067977499786*alphaR2[16]-0.25*alphaR2[8]+0.3354101966249678*alphaR2[6]-0.25*alphaR2[5]+0.25*alphaR2[4]-0.3354101966249678*alphaR2[3]+0.25*alphaR2[2]-0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[9] = 1.0; 
  else  
    sgn_alpha_surf2[9] = -1.0; 
  
  if (sgn_alpha_surf2[9] == sgn_alpha_surf2[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaR2[17]-0.2795084971874732*alphaR2[16]-0.25*alphaR2[8]-0.25*alphaR2[5]+0.25*alphaR2[4]+0.25*alphaR2[2]-0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[10] = 1.0; 
  else  
    sgn_alpha_surf2[10] = -1.0; 
  
  if (sgn_alpha_surf2[10] == sgn_alpha_surf2[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaR2[17])+0.2236067977499786*alphaR2[16]-0.25*alphaR2[8]-0.3354101966249678*alphaR2[6]-0.25*alphaR2[5]+0.25*alphaR2[4]+0.3354101966249678*alphaR2[3]+0.25*alphaR2[2]-0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[11] = 1.0; 
  else  
    sgn_alpha_surf2[11] = -1.0; 
  
  if (sgn_alpha_surf2[11] == sgn_alpha_surf2[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*alphaR2[17]+0.2236067977499786*alphaR2[16]-0.25*alphaR2[8]-0.3354101966249678*alphaR2[6]-0.25*alphaR2[5]-0.25*alphaR2[4]-0.3354101966249678*alphaR2[3]-0.25*alphaR2[2]+0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[12] = 1.0; 
  else  
    sgn_alpha_surf2[12] = -1.0; 
  
  if (sgn_alpha_surf2[12] == sgn_alpha_surf2[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*alphaR2[17])-0.2795084971874732*alphaR2[16]-0.25*alphaR2[8]-0.25*alphaR2[5]-0.25*alphaR2[4]-0.25*alphaR2[2]+0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[13] = 1.0; 
  else  
    sgn_alpha_surf2[13] = -1.0; 
  
  if (sgn_alpha_surf2[13] == sgn_alpha_surf2[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*alphaR2[17]+0.2236067977499786*alphaR2[16]-0.25*alphaR2[8]+0.3354101966249678*alphaR2[6]-0.25*alphaR2[5]-0.25*alphaR2[4]+0.3354101966249678*alphaR2[3]-0.25*alphaR2[2]+0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[14] = 1.0; 
  else  
    sgn_alpha_surf2[14] = -1.0; 
  
  if (sgn_alpha_surf2[14] == sgn_alpha_surf2[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*alphaR2[17]+0.2236067977499786*alphaR2[16]+0.25*alphaR2[8]-0.3354101966249678*alphaR2[6]-0.25*alphaR2[5]+0.25*alphaR2[4]-0.3354101966249678*alphaR2[3]-0.25*alphaR2[2]+0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[15] = 1.0; 
  else  
    sgn_alpha_surf2[15] = -1.0; 
  
  if (sgn_alpha_surf2[15] == sgn_alpha_surf2[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*alphaR2[17])-0.2795084971874732*alphaR2[16]+0.25*alphaR2[8]-0.25*alphaR2[5]+0.25*alphaR2[4]-0.25*alphaR2[2]+0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[16] = 1.0; 
  else  
    sgn_alpha_surf2[16] = -1.0; 
  
  if (sgn_alpha_surf2[16] == sgn_alpha_surf2[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*alphaR2[17]+0.2236067977499786*alphaR2[16]+0.25*alphaR2[8]+0.3354101966249678*alphaR2[6]-0.25*alphaR2[5]+0.25*alphaR2[4]+0.3354101966249678*alphaR2[3]-0.25*alphaR2[2]+0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[17] = 1.0; 
  else  
    sgn_alpha_surf2[17] = -1.0; 
  
  if (sgn_alpha_surf2[17] == sgn_alpha_surf2[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*alphaR2[17]+0.2236067977499786*alphaR2[16]-0.25*alphaR2[8]-0.3354101966249678*alphaR2[6]+0.25*alphaR2[5]-0.25*alphaR2[4]-0.3354101966249678*alphaR2[3]+0.25*alphaR2[2]+0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[18] = 1.0; 
  else  
    sgn_alpha_surf2[18] = -1.0; 
  
  if (sgn_alpha_surf2[18] == sgn_alpha_surf2[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*alphaR2[17])-0.2795084971874732*alphaR2[16]-0.25*alphaR2[8]+0.25*alphaR2[5]-0.25*alphaR2[4]+0.25*alphaR2[2]+0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[19] = 1.0; 
  else  
    sgn_alpha_surf2[19] = -1.0; 
  
  if (sgn_alpha_surf2[19] == sgn_alpha_surf2[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*alphaR2[17]+0.2236067977499786*alphaR2[16]-0.25*alphaR2[8]+0.3354101966249678*alphaR2[6]+0.25*alphaR2[5]-0.25*alphaR2[4]+0.3354101966249678*alphaR2[3]+0.25*alphaR2[2]+0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[20] = 1.0; 
  else  
    sgn_alpha_surf2[20] = -1.0; 
  
  if (sgn_alpha_surf2[20] == sgn_alpha_surf2[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*alphaR2[17]+0.2236067977499786*alphaR2[16]+0.25*alphaR2[8]-0.3354101966249678*alphaR2[6]+0.25*alphaR2[5]+0.25*alphaR2[4]-0.3354101966249678*alphaR2[3]+0.25*alphaR2[2]+0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[21] = 1.0; 
  else  
    sgn_alpha_surf2[21] = -1.0; 
  
  if (sgn_alpha_surf2[21] == sgn_alpha_surf2[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*alphaR2[17])-0.2795084971874732*alphaR2[16]+0.25*alphaR2[8]+0.25*alphaR2[5]+0.25*alphaR2[4]+0.25*alphaR2[2]+0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[22] = 1.0; 
  else  
    sgn_alpha_surf2[22] = -1.0; 
  
  if (sgn_alpha_surf2[22] == sgn_alpha_surf2[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*alphaR2[17]+0.2236067977499786*alphaR2[16]+0.25*alphaR2[8]+0.3354101966249678*alphaR2[6]+0.25*alphaR2[5]+0.25*alphaR2[4]+0.3354101966249678*alphaR2[3]+0.25*alphaR2[2]+0.25*alphaR2[1]+0.25*alphaR2[0] > 0.) 
    sgn_alpha_surf2[23] = 1.0; 
  else  
    sgn_alpha_surf2[23] = -1.0; 
  
  if (sgn_alpha_surf2[23] == sgn_alpha_surf2[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
