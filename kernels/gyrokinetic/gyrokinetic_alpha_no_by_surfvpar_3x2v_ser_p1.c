#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfvpar_3x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, 
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

  double *alphaL3 = &alpha_surf[72];
  alphaL3[0] = ((((-0.2651650429449552*b_x[3]*jacobtot_inv[5]*hamil[16])-0.2651650429449552*jacobtot_inv[3]*b_x[5]*hamil[16]-0.2651650429449552*b_x[5]*jacobtot_inv[5]*hamil[8]-0.2651650429449552*b_x[3]*jacobtot_inv[3]*hamil[8]-0.2651650429449552*jacobtot_inv[0]*b_x[5]*hamil[6]-0.2651650429449552*jacobtot_inv[1]*b_x[3]*hamil[6]-0.2651650429449552*jacobtot_inv[1]*hamil[2]*b_x[5]-0.2651650429449552*jacobtot_inv[0]*hamil[2]*b_x[3])*rdy2*rdz2+(0.2651650429449552*b_z[1]*jacobtot_inv[5]*hamil[16]+0.2651650429449552*jacobtot_inv[1]*b_z[5]*hamil[16]+0.2651650429449552*jacobtot_inv[0]*b_z[5]*hamil[8]+0.2651650429449552*b_z[1]*jacobtot_inv[3]*hamil[8]+0.2651650429449552*b_z[5]*jacobtot_inv[5]*hamil[6]+0.2651650429449552*b_z[1]*jacobtot_inv[1]*hamil[6]+0.2651650429449552*hamil[2]*jacobtot_inv[3]*b_z[5]+0.2651650429449552*jacobtot_inv[0]*b_z[1]*hamil[2])*rdx2*rdy2)*wvpar)/q_+(((0.2651650429449552*b_x[3]*jacobtot_inv[5]*hamil[16]+0.2651650429449552*jacobtot_inv[3]*b_x[5]*hamil[16]+0.2651650429449552*b_x[5]*jacobtot_inv[5]*hamil[8]+0.2651650429449552*b_x[3]*jacobtot_inv[3]*hamil[8]+0.2651650429449552*jacobtot_inv[0]*b_x[5]*hamil[6]+0.2651650429449552*jacobtot_inv[1]*b_x[3]*hamil[6]+0.2651650429449552*jacobtot_inv[1]*hamil[2]*b_x[5]+0.2651650429449552*jacobtot_inv[0]*hamil[2]*b_x[3])*rdy2*rdz2)/rdvpar2+(((-0.2651650429449552*b_z[1]*jacobtot_inv[5]*hamil[16])-0.2651650429449552*jacobtot_inv[1]*b_z[5]*hamil[16]-0.2651650429449552*jacobtot_inv[0]*b_z[5]*hamil[8]-0.2651650429449552*b_z[1]*jacobtot_inv[3]*hamil[8]-0.2651650429449552*b_z[5]*jacobtot_inv[5]*hamil[6]-0.2651650429449552*b_z[1]*jacobtot_inv[1]*hamil[6]-0.2651650429449552*hamil[2]*jacobtot_inv[3]*b_z[5]-0.2651650429449552*jacobtot_inv[0]*b_z[1]*hamil[2])*rdx2*rdy2)/rdvpar2)/q_+(((-0.1530931089239486*cmag[3]*jacobtot_inv[5]*hamil[7])-0.1530931089239486*jacobtot_inv[3]*cmag[5]*hamil[7]-0.1530931089239486*cmag[0]*jacobtot_inv[1]*hamil[7]-0.1530931089239486*jacobtot_inv[0]*cmag[1]*hamil[7]-0.1530931089239486*hamil[3]*cmag[5]*jacobtot_inv[5]-0.1530931089239486*cmag[3]*hamil[3]*jacobtot_inv[3]-0.1530931089239486*cmag[1]*jacobtot_inv[1]*hamil[3]-0.1530931089239486*cmag[0]*jacobtot_inv[0]*hamil[3])*rdz2)/m_; 
  alphaL3[1] = ((((-0.4772970773009194*b_x[5]*jacobtot_inv[5]*hamil[16])-0.2651650429449552*b_x[3]*jacobtot_inv[3]*hamil[16]-0.2651650429449552*b_x[3]*jacobtot_inv[5]*hamil[8]-0.2651650429449552*jacobtot_inv[3]*b_x[5]*hamil[8]-0.4772970773009194*jacobtot_inv[1]*b_x[5]*hamil[6]-0.2651650429449552*jacobtot_inv[0]*b_x[3]*hamil[6]-0.2651650429449552*jacobtot_inv[0]*hamil[2]*b_x[5]-0.2651650429449552*jacobtot_inv[1]*hamil[2]*b_x[3])*rdy2*rdz2+(0.2651650429449552*jacobtot_inv[0]*b_z[5]*hamil[16]+0.2651650429449552*b_z[1]*jacobtot_inv[3]*hamil[16]+0.2651650429449552*b_z[1]*jacobtot_inv[5]*hamil[8]+0.2651650429449552*jacobtot_inv[1]*b_z[5]*hamil[8]+0.2651650429449552*jacobtot_inv[3]*b_z[5]*hamil[6]+0.2651650429449552*jacobtot_inv[0]*b_z[1]*hamil[6]+0.2651650429449552*hamil[2]*b_z[5]*jacobtot_inv[5]+0.2651650429449552*b_z[1]*jacobtot_inv[1]*hamil[2])*rdx2*rdy2)*wvpar)/q_+(((0.4772970773009194*b_x[5]*jacobtot_inv[5]*hamil[16]+0.2651650429449552*b_x[3]*jacobtot_inv[3]*hamil[16]+0.2651650429449552*b_x[3]*jacobtot_inv[5]*hamil[8]+0.2651650429449552*jacobtot_inv[3]*b_x[5]*hamil[8]+0.4772970773009194*jacobtot_inv[1]*b_x[5]*hamil[6]+0.2651650429449552*jacobtot_inv[0]*b_x[3]*hamil[6]+0.2651650429449552*jacobtot_inv[0]*hamil[2]*b_x[5]+0.2651650429449552*jacobtot_inv[1]*hamil[2]*b_x[3])*rdy2*rdz2)/rdvpar2+(((-0.2651650429449552*jacobtot_inv[0]*b_z[5]*hamil[16])-0.2651650429449552*b_z[1]*jacobtot_inv[3]*hamil[16]-0.2651650429449552*b_z[1]*jacobtot_inv[5]*hamil[8]-0.2651650429449552*jacobtot_inv[1]*b_z[5]*hamil[8]-0.2651650429449552*jacobtot_inv[3]*b_z[5]*hamil[6]-0.2651650429449552*jacobtot_inv[0]*b_z[1]*hamil[6]-0.2651650429449552*hamil[2]*b_z[5]*jacobtot_inv[5]-0.2651650429449552*b_z[1]*jacobtot_inv[1]*hamil[2])*rdx2*rdy2)/rdvpar2)/q_+(((-0.2755675960631073*cmag[5]*jacobtot_inv[5]*hamil[7])-0.1530931089239486*cmag[3]*jacobtot_inv[3]*hamil[7]-0.2755675960631073*cmag[1]*jacobtot_inv[1]*hamil[7]-0.1530931089239486*cmag[0]*jacobtot_inv[0]*hamil[7]-0.1530931089239486*cmag[3]*hamil[3]*jacobtot_inv[5]-0.1530931089239486*hamil[3]*jacobtot_inv[3]*cmag[5]-0.1530931089239486*cmag[0]*jacobtot_inv[1]*hamil[3]-0.1530931089239486*jacobtot_inv[0]*cmag[1]*hamil[3])*rdz2)/m_; 
  alphaL3[2] = (((-0.1530931089239486*cmag[3]*jacobtot_inv[5]*hamil[16])-0.1530931089239486*jacobtot_inv[3]*cmag[5]*hamil[16]-0.1530931089239486*cmag[0]*jacobtot_inv[1]*hamil[16]-0.1530931089239486*jacobtot_inv[0]*cmag[1]*hamil[16]-0.1530931089239486*cmag[5]*jacobtot_inv[5]*hamil[8]-0.1530931089239486*cmag[3]*jacobtot_inv[3]*hamil[8]-0.1530931089239486*cmag[1]*jacobtot_inv[1]*hamil[8]-0.1530931089239486*cmag[0]*jacobtot_inv[0]*hamil[8])*rdz2)/m_; 
  alphaL3[3] = ((((-0.2651650429449552*jacobtot_inv[0]*b_x[5]*hamil[16])-0.2651650429449552*jacobtot_inv[1]*b_x[3]*hamil[16]-0.2651650429449552*jacobtot_inv[1]*b_x[5]*hamil[8]-0.2651650429449552*jacobtot_inv[0]*b_x[3]*hamil[8]-0.2651650429449552*b_x[3]*jacobtot_inv[5]*hamil[6]-0.2651650429449552*jacobtot_inv[3]*b_x[5]*hamil[6]-0.2651650429449552*hamil[2]*b_x[5]*jacobtot_inv[5]-0.2651650429449552*hamil[2]*b_x[3]*jacobtot_inv[3])*rdy2*rdz2+(0.4772970773009194*b_z[5]*jacobtot_inv[5]*hamil[16]+0.2651650429449552*b_z[1]*jacobtot_inv[1]*hamil[16]+0.4772970773009194*jacobtot_inv[3]*b_z[5]*hamil[8]+0.2651650429449552*jacobtot_inv[0]*b_z[1]*hamil[8]+0.2651650429449552*b_z[1]*jacobtot_inv[5]*hamil[6]+0.2651650429449552*jacobtot_inv[1]*b_z[5]*hamil[6]+0.2651650429449552*jacobtot_inv[0]*hamil[2]*b_z[5]+0.2651650429449552*b_z[1]*hamil[2]*jacobtot_inv[3])*rdx2*rdy2)*wvpar)/q_+(((0.2651650429449552*jacobtot_inv[0]*b_x[5]*hamil[16]+0.2651650429449552*jacobtot_inv[1]*b_x[3]*hamil[16]+0.2651650429449552*jacobtot_inv[1]*b_x[5]*hamil[8]+0.2651650429449552*jacobtot_inv[0]*b_x[3]*hamil[8]+0.2651650429449552*b_x[3]*jacobtot_inv[5]*hamil[6]+0.2651650429449552*jacobtot_inv[3]*b_x[5]*hamil[6]+0.2651650429449552*hamil[2]*b_x[5]*jacobtot_inv[5]+0.2651650429449552*hamil[2]*b_x[3]*jacobtot_inv[3])*rdy2*rdz2)/rdvpar2+(((-0.4772970773009194*b_z[5]*jacobtot_inv[5]*hamil[16])-0.2651650429449552*b_z[1]*jacobtot_inv[1]*hamil[16]-0.4772970773009194*jacobtot_inv[3]*b_z[5]*hamil[8]-0.2651650429449552*jacobtot_inv[0]*b_z[1]*hamil[8]-0.2651650429449552*b_z[1]*jacobtot_inv[5]*hamil[6]-0.2651650429449552*jacobtot_inv[1]*b_z[5]*hamil[6]-0.2651650429449552*jacobtot_inv[0]*hamil[2]*b_z[5]-0.2651650429449552*b_z[1]*hamil[2]*jacobtot_inv[3])*rdx2*rdy2)/rdvpar2)/q_+(((-0.1530931089239486*cmag[0]*jacobtot_inv[5]*hamil[7])-0.1530931089239486*jacobtot_inv[0]*cmag[5]*hamil[7]-0.1530931089239486*cmag[1]*jacobtot_inv[3]*hamil[7]-0.1530931089239486*jacobtot_inv[1]*cmag[3]*hamil[7]-0.1530931089239486*cmag[1]*hamil[3]*jacobtot_inv[5]-0.1530931089239486*jacobtot_inv[1]*hamil[3]*cmag[5]-0.1530931089239486*cmag[0]*hamil[3]*jacobtot_inv[3]-0.1530931089239486*jacobtot_inv[0]*cmag[3]*hamil[3])*rdz2)/m_; 
  alphaL3[4] = (((-0.1530931089239486*cmag[3]*jacobtot_inv[5]*hamil[21])-0.1530931089239486*jacobtot_inv[3]*cmag[5]*hamil[21]-0.1530931089239486*cmag[0]*jacobtot_inv[1]*hamil[21]-0.1530931089239486*jacobtot_inv[0]*cmag[1]*hamil[21]-0.1530931089239486*cmag[5]*jacobtot_inv[5]*hamil[14]-0.1530931089239486*cmag[3]*jacobtot_inv[3]*hamil[14]-0.1530931089239486*cmag[1]*jacobtot_inv[1]*hamil[14]-0.1530931089239486*cmag[0]*jacobtot_inv[0]*hamil[14])*rdz2)/m_; 
  alphaL3[5] = (((-0.2755675960631073*cmag[5]*jacobtot_inv[5]*hamil[16])-0.1530931089239486*cmag[3]*jacobtot_inv[3]*hamil[16]-0.2755675960631073*cmag[1]*jacobtot_inv[1]*hamil[16]-0.1530931089239486*cmag[0]*jacobtot_inv[0]*hamil[16]-0.1530931089239486*cmag[3]*jacobtot_inv[5]*hamil[8]-0.1530931089239486*jacobtot_inv[3]*cmag[5]*hamil[8]-0.1530931089239486*cmag[0]*jacobtot_inv[1]*hamil[8]-0.1530931089239486*jacobtot_inv[0]*cmag[1]*hamil[8])*rdz2)/m_; 
  alphaL3[6] = ((((-0.4772970773009194*jacobtot_inv[1]*b_x[5]*hamil[16])-0.2651650429449552*jacobtot_inv[0]*b_x[3]*hamil[16]-0.2651650429449552*jacobtot_inv[0]*b_x[5]*hamil[8]-0.2651650429449552*jacobtot_inv[1]*b_x[3]*hamil[8]-0.4772970773009194*b_x[5]*jacobtot_inv[5]*hamil[6]-0.2651650429449552*b_x[3]*jacobtot_inv[3]*hamil[6]-0.2651650429449552*hamil[2]*b_x[3]*jacobtot_inv[5]-0.2651650429449552*hamil[2]*jacobtot_inv[3]*b_x[5])*rdy2*rdz2+(0.4772970773009194*jacobtot_inv[3]*b_z[5]*hamil[16]+0.2651650429449552*jacobtot_inv[0]*b_z[1]*hamil[16]+0.4772970773009194*b_z[5]*jacobtot_inv[5]*hamil[8]+0.2651650429449552*b_z[1]*jacobtot_inv[1]*hamil[8]+0.2651650429449552*jacobtot_inv[0]*b_z[5]*hamil[6]+0.2651650429449552*b_z[1]*jacobtot_inv[3]*hamil[6]+0.2651650429449552*b_z[1]*hamil[2]*jacobtot_inv[5]+0.2651650429449552*jacobtot_inv[1]*hamil[2]*b_z[5])*rdx2*rdy2)*wvpar)/q_+(((0.4772970773009194*jacobtot_inv[1]*b_x[5]*hamil[16]+0.2651650429449552*jacobtot_inv[0]*b_x[3]*hamil[16]+0.2651650429449552*jacobtot_inv[0]*b_x[5]*hamil[8]+0.2651650429449552*jacobtot_inv[1]*b_x[3]*hamil[8]+0.4772970773009194*b_x[5]*jacobtot_inv[5]*hamil[6]+0.2651650429449552*b_x[3]*jacobtot_inv[3]*hamil[6]+0.2651650429449552*hamil[2]*b_x[3]*jacobtot_inv[5]+0.2651650429449552*hamil[2]*jacobtot_inv[3]*b_x[5])*rdy2*rdz2)/rdvpar2+(((-0.4772970773009194*jacobtot_inv[3]*b_z[5]*hamil[16])-0.2651650429449552*jacobtot_inv[0]*b_z[1]*hamil[16]-0.4772970773009194*b_z[5]*jacobtot_inv[5]*hamil[8]-0.2651650429449552*b_z[1]*jacobtot_inv[1]*hamil[8]-0.2651650429449552*jacobtot_inv[0]*b_z[5]*hamil[6]-0.2651650429449552*b_z[1]*jacobtot_inv[3]*hamil[6]-0.2651650429449552*b_z[1]*hamil[2]*jacobtot_inv[5]-0.2651650429449552*jacobtot_inv[1]*hamil[2]*b_z[5])*rdx2*rdy2)/rdvpar2)/q_+(((-0.2755675960631073*cmag[1]*jacobtot_inv[5]*hamil[7])-0.2755675960631073*jacobtot_inv[1]*cmag[5]*hamil[7]-0.1530931089239486*cmag[0]*jacobtot_inv[3]*hamil[7]-0.1530931089239486*jacobtot_inv[0]*cmag[3]*hamil[7]-0.1530931089239486*cmag[0]*hamil[3]*jacobtot_inv[5]-0.1530931089239486*jacobtot_inv[0]*hamil[3]*cmag[5]-0.1530931089239486*cmag[1]*hamil[3]*jacobtot_inv[3]-0.1530931089239486*jacobtot_inv[1]*cmag[3]*hamil[3])*rdz2)/m_; 
  alphaL3[7] = (((-0.1530931089239486*cmag[0]*jacobtot_inv[5]*hamil[16])-0.1530931089239486*jacobtot_inv[0]*cmag[5]*hamil[16]-0.1530931089239486*cmag[1]*jacobtot_inv[3]*hamil[16]-0.1530931089239486*jacobtot_inv[1]*cmag[3]*hamil[16]-0.1530931089239486*cmag[1]*jacobtot_inv[5]*hamil[8]-0.1530931089239486*jacobtot_inv[1]*cmag[5]*hamil[8]-0.1530931089239486*cmag[0]*jacobtot_inv[3]*hamil[8]-0.1530931089239486*jacobtot_inv[0]*cmag[3]*hamil[8])*rdz2)/m_; 
  alphaL3[8] = (((-0.2755675960631073*cmag[5]*jacobtot_inv[5]*hamil[21])-0.1530931089239486*cmag[3]*jacobtot_inv[3]*hamil[21]-0.2755675960631073*cmag[1]*jacobtot_inv[1]*hamil[21]-0.1530931089239486*cmag[0]*jacobtot_inv[0]*hamil[21]-0.1530931089239486*cmag[3]*jacobtot_inv[5]*hamil[14]-0.1530931089239486*jacobtot_inv[3]*cmag[5]*hamil[14]-0.1530931089239486*cmag[0]*jacobtot_inv[1]*hamil[14]-0.1530931089239486*jacobtot_inv[0]*cmag[1]*hamil[14])*rdz2)/m_; 
  alphaL3[10] = (((-0.1530931089239486*cmag[0]*jacobtot_inv[5]*hamil[21])-0.1530931089239486*jacobtot_inv[0]*cmag[5]*hamil[21]-0.1530931089239486*cmag[1]*jacobtot_inv[3]*hamil[21]-0.1530931089239486*jacobtot_inv[1]*cmag[3]*hamil[21]-0.1530931089239486*cmag[1]*jacobtot_inv[5]*hamil[14]-0.1530931089239486*jacobtot_inv[1]*cmag[5]*hamil[14]-0.1530931089239486*cmag[0]*jacobtot_inv[3]*hamil[14]-0.1530931089239486*jacobtot_inv[0]*cmag[3]*hamil[14])*rdz2)/m_; 
  alphaL3[11] = (((-0.2755675960631073*cmag[1]*jacobtot_inv[5]*hamil[16])-0.2755675960631073*jacobtot_inv[1]*cmag[5]*hamil[16]-0.1530931089239486*cmag[0]*jacobtot_inv[3]*hamil[16]-0.1530931089239486*jacobtot_inv[0]*cmag[3]*hamil[16]-0.1530931089239486*cmag[0]*jacobtot_inv[5]*hamil[8]-0.1530931089239486*jacobtot_inv[0]*cmag[5]*hamil[8]-0.1530931089239486*cmag[1]*jacobtot_inv[3]*hamil[8]-0.1530931089239486*jacobtot_inv[1]*cmag[3]*hamil[8])*rdz2)/m_; 
  alphaL3[13] = (((-0.2755675960631073*cmag[1]*jacobtot_inv[5]*hamil[21])-0.2755675960631073*jacobtot_inv[1]*cmag[5]*hamil[21]-0.1530931089239486*cmag[0]*jacobtot_inv[3]*hamil[21]-0.1530931089239486*jacobtot_inv[0]*cmag[3]*hamil[21]-0.1530931089239486*cmag[0]*jacobtot_inv[5]*hamil[14]-0.1530931089239486*jacobtot_inv[0]*cmag[5]*hamil[14]-0.1530931089239486*cmag[1]*jacobtot_inv[3]*hamil[14]-0.1530931089239486*jacobtot_inv[1]*cmag[3]*hamil[14])*rdz2)/m_; 

  double *sgn_alpha_surf3 = &sgn_alpha_surf[72];
  int const_sgn_alpha_surf = 1; 
  
  if ((-0.25*(alphaL3[13]+alphaL3[11]))+0.25*(alphaL3[10]+alphaL3[8]+alphaL3[7]+alphaL3[6]+alphaL3[5])-0.25*(alphaL3[4]+alphaL3[3]+alphaL3[2]+alphaL3[1])+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[0] = 1.0; 
  else  
    sgn_alpha_surf3[0] = -1.0; 
  
  if (0.25*alphaL3[13]-0.25*(alphaL3[11]+alphaL3[10]+alphaL3[8])+0.25*(alphaL3[7]+alphaL3[6]+alphaL3[5]+alphaL3[4])-0.25*(alphaL3[3]+alphaL3[2]+alphaL3[1])+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[1] = 1.0; 
  else  
    sgn_alpha_surf3[1] = -1.0; 
  
  if (sgn_alpha_surf3[1] == sgn_alpha_surf3[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL3[13]+alphaL3[11])-0.25*alphaL3[10]+0.25*alphaL3[8]-0.25*(alphaL3[7]+alphaL3[6])+0.25*alphaL3[5]-0.25*alphaL3[4]+0.25*alphaL3[3]-0.25*(alphaL3[2]+alphaL3[1])+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[2] = 1.0; 
  else  
    sgn_alpha_surf3[2] = -1.0; 
  
  if (sgn_alpha_surf3[2] == sgn_alpha_surf3[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.25*alphaL3[13])+0.25*(alphaL3[11]+alphaL3[10])-0.25*(alphaL3[8]+alphaL3[7]+alphaL3[6])+0.25*(alphaL3[5]+alphaL3[4]+alphaL3[3])-0.25*(alphaL3[2]+alphaL3[1])+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[3] = 1.0; 
  else  
    sgn_alpha_surf3[3] = -1.0; 
  
  if (sgn_alpha_surf3[3] == sgn_alpha_surf3[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.25*alphaL3[13])+0.25*(alphaL3[11]+alphaL3[10]+alphaL3[8])-0.25*alphaL3[7]+0.25*alphaL3[6]-0.25*(alphaL3[5]+alphaL3[4]+alphaL3[3])+0.25*alphaL3[2]-0.25*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[4] = 1.0; 
  else  
    sgn_alpha_surf3[4] = -1.0; 
  
  if (sgn_alpha_surf3[4] == sgn_alpha_surf3[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL3[13]+alphaL3[11])-0.25*(alphaL3[10]+alphaL3[8]+alphaL3[7])+0.25*alphaL3[6]-0.25*alphaL3[5]+0.25*alphaL3[4]-0.25*alphaL3[3]+0.25*alphaL3[2]-0.25*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[5] = 1.0; 
  else  
    sgn_alpha_surf3[5] = -1.0; 
  
  if (sgn_alpha_surf3[5] == sgn_alpha_surf3[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*alphaL3[13]-0.25*(alphaL3[11]+alphaL3[10])+0.25*(alphaL3[8]+alphaL3[7])-0.25*(alphaL3[6]+alphaL3[5]+alphaL3[4])+0.25*(alphaL3[3]+alphaL3[2])-0.25*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[6] = 1.0; 
  else  
    sgn_alpha_surf3[6] = -1.0; 
  
  if (sgn_alpha_surf3[6] == sgn_alpha_surf3[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.25*(alphaL3[13]+alphaL3[11]))+0.25*alphaL3[10]-0.25*alphaL3[8]+0.25*alphaL3[7]-0.25*(alphaL3[6]+alphaL3[5])+0.25*(alphaL3[4]+alphaL3[3]+alphaL3[2])-0.25*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[7] = 1.0; 
  else  
    sgn_alpha_surf3[7] = -1.0; 
  
  if (sgn_alpha_surf3[7] == sgn_alpha_surf3[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL3[13]+alphaL3[11]+alphaL3[10])-0.25*alphaL3[8]+0.25*alphaL3[7]-0.25*(alphaL3[6]+alphaL3[5]+alphaL3[4]+alphaL3[3]+alphaL3[2])+0.25*(alphaL3[1]+alphaL3[0]) > 0.) 
    sgn_alpha_surf3[8] = 1.0; 
  else  
    sgn_alpha_surf3[8] = -1.0; 
  
  if (sgn_alpha_surf3[8] == sgn_alpha_surf3[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.25*alphaL3[13])+0.25*alphaL3[11]-0.25*alphaL3[10]+0.25*(alphaL3[8]+alphaL3[7])-0.25*(alphaL3[6]+alphaL3[5])+0.25*alphaL3[4]-0.25*(alphaL3[3]+alphaL3[2])+0.25*(alphaL3[1]+alphaL3[0]) > 0.) 
    sgn_alpha_surf3[9] = 1.0; 
  else  
    sgn_alpha_surf3[9] = -1.0; 
  
  if (sgn_alpha_surf3[9] == sgn_alpha_surf3[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.25*(alphaL3[13]+alphaL3[11]+alphaL3[10]+alphaL3[8]+alphaL3[7]))+0.25*alphaL3[6]-0.25*(alphaL3[5]+alphaL3[4])+0.25*alphaL3[3]-0.25*alphaL3[2]+0.25*(alphaL3[1]+alphaL3[0]) > 0.) 
    sgn_alpha_surf3[10] = 1.0; 
  else  
    sgn_alpha_surf3[10] = -1.0; 
  
  if (sgn_alpha_surf3[10] == sgn_alpha_surf3[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*alphaL3[13]-0.25*alphaL3[11]+0.25*(alphaL3[10]+alphaL3[8])-0.25*alphaL3[7]+0.25*alphaL3[6]-0.25*alphaL3[5]+0.25*(alphaL3[4]+alphaL3[3])-0.25*alphaL3[2]+0.25*(alphaL3[1]+alphaL3[0]) > 0.) 
    sgn_alpha_surf3[11] = 1.0; 
  else  
    sgn_alpha_surf3[11] = -1.0; 
  
  if (sgn_alpha_surf3[11] == sgn_alpha_surf3[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*alphaL3[13]-0.25*alphaL3[11]+0.25*alphaL3[10]-0.25*(alphaL3[8]+alphaL3[7]+alphaL3[6])+0.25*alphaL3[5]-0.25*(alphaL3[4]+alphaL3[3])+0.25*(alphaL3[2]+alphaL3[1]+alphaL3[0]) > 0.) 
    sgn_alpha_surf3[12] = 1.0; 
  else  
    sgn_alpha_surf3[12] = -1.0; 
  
  if (sgn_alpha_surf3[12] == sgn_alpha_surf3[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.25*(alphaL3[13]+alphaL3[11]+alphaL3[10]))+0.25*alphaL3[8]-0.25*(alphaL3[7]+alphaL3[6])+0.25*(alphaL3[5]+alphaL3[4])-0.25*alphaL3[3]+0.25*(alphaL3[2]+alphaL3[1]+alphaL3[0]) > 0.) 
    sgn_alpha_surf3[13] = 1.0; 
  else  
    sgn_alpha_surf3[13] = -1.0; 
  
  if (sgn_alpha_surf3[13] == sgn_alpha_surf3[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.25*alphaL3[13])+0.25*alphaL3[11]-0.25*(alphaL3[10]+alphaL3[8])+0.25*(alphaL3[7]+alphaL3[6]+alphaL3[5])-0.25*alphaL3[4]+0.25*(alphaL3[3]+alphaL3[2]+alphaL3[1]+alphaL3[0]) > 0.) 
    sgn_alpha_surf3[14] = 1.0; 
  else  
    sgn_alpha_surf3[14] = -1.0; 
  
  if (sgn_alpha_surf3[14] == sgn_alpha_surf3[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL3[13]+alphaL3[11]+alphaL3[10]+alphaL3[8]+alphaL3[7]+alphaL3[6]+alphaL3[5]+alphaL3[4]+alphaL3[3]+alphaL3[2]+alphaL3[1]+alphaL3[0]) > 0.) 
    sgn_alpha_surf3[15] = 1.0; 
  else  
    sgn_alpha_surf3[15] = -1.0; 
  
  if (sgn_alpha_surf3[15] == sgn_alpha_surf3[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 