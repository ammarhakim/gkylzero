#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_surfy_3x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, 
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

  double *alphaL1 = &alpha_surf[24];
  alphaL1[0] = (((0.2651650429449552*jacobtot_inv[1]*hamil[4]*b_x[5]+0.2651650429449552*jacobtot_inv[0]*b_x[3]*hamil[4])*rdvpar2*rdz2+((-0.2651650429449552*jacobtot_inv[3]*hamil[4]*b_z[5])-0.2651650429449552*jacobtot_inv[0]*b_z[1]*hamil[4])*rdvpar2*rdx2)*wvpar)/q_+((0.3423265984407287*jacobtot_inv[1]*b_x[5]*hamil[32]+0.3423265984407287*jacobtot_inv[0]*b_x[3]*hamil[32]+0.2651650429449552*b_x[3]*jacobtot_inv[5]*hamil[16]+0.2651650429449552*jacobtot_inv[3]*b_x[5]*hamil[16]+0.2651650429449552*b_x[0]*jacobtot_inv[1]*hamil[16]+0.2651650429449552*jacobtot_inv[0]*b_x[1]*hamil[16]+0.2651650429449552*b_x[5]*jacobtot_inv[5]*hamil[8]+0.2651650429449552*b_x[3]*jacobtot_inv[3]*hamil[8]+0.2651650429449552*b_x[1]*jacobtot_inv[1]*hamil[8]+0.2651650429449552*b_x[0]*jacobtot_inv[0]*hamil[8]-0.1530931089239486*b_x[3]*jacobtot_inv[5]*hamil[7]-0.1530931089239486*jacobtot_inv[3]*b_x[5]*hamil[7]-0.1530931089239486*b_x[0]*jacobtot_inv[1]*hamil[7]-0.1530931089239486*jacobtot_inv[0]*b_x[1]*hamil[7]-0.1530931089239486*hamil[3]*b_x[5]*jacobtot_inv[5]-0.1530931089239486*b_x[3]*hamil[3]*jacobtot_inv[3]-0.1530931089239486*b_x[1]*jacobtot_inv[1]*hamil[3]-0.1530931089239486*b_x[0]*jacobtot_inv[0]*hamil[3])*rdz2+((-0.3423265984407287*jacobtot_inv[3]*b_z[5]*hamil[32])-0.3423265984407287*jacobtot_inv[0]*b_z[1]*hamil[32]-0.2651650429449552*b_z[1]*jacobtot_inv[5]*hamil[16]-0.2651650429449552*jacobtot_inv[1]*b_z[5]*hamil[16]-0.2651650429449552*b_z[0]*jacobtot_inv[3]*hamil[16]-0.2651650429449552*jacobtot_inv[0]*b_z[3]*hamil[16]+0.1530931089239486*b_z[1]*jacobtot_inv[5]*hamil[7]+0.1530931089239486*jacobtot_inv[1]*b_z[5]*hamil[7]+0.1530931089239486*b_z[0]*jacobtot_inv[3]*hamil[7]+0.1530931089239486*jacobtot_inv[0]*b_z[3]*hamil[7]-0.2651650429449552*b_z[5]*jacobtot_inv[5]*hamil[6]-0.2651650429449552*b_z[3]*jacobtot_inv[3]*hamil[6]-0.2651650429449552*b_z[1]*jacobtot_inv[1]*hamil[6]-0.2651650429449552*b_z[0]*jacobtot_inv[0]*hamil[6]+0.1530931089239486*hamil[1]*b_z[5]*jacobtot_inv[5]+0.1530931089239486*hamil[1]*b_z[3]*jacobtot_inv[3]+0.1530931089239486*b_z[1]*hamil[1]*jacobtot_inv[1]+0.1530931089239486*b_z[0]*jacobtot_inv[0]*hamil[1])*rdx2)/q_; 
  alphaL1[1] = (((0.2651650429449552*jacobtot_inv[0]*hamil[4]*b_x[5]+0.2651650429449552*jacobtot_inv[1]*b_x[3]*hamil[4])*rdvpar2*rdz2+((-0.2651650429449552*hamil[4]*b_z[5]*jacobtot_inv[5])-0.2651650429449552*b_z[1]*jacobtot_inv[1]*hamil[4])*rdvpar2*rdx2)*wvpar)/q_+((0.3423265984407287*jacobtot_inv[0]*b_x[5]*hamil[32]+0.3423265984407287*jacobtot_inv[1]*b_x[3]*hamil[32]+0.4772970773009194*b_x[5]*jacobtot_inv[5]*hamil[16]+0.2651650429449552*b_x[3]*jacobtot_inv[3]*hamil[16]+0.4772970773009194*b_x[1]*jacobtot_inv[1]*hamil[16]+0.2651650429449552*b_x[0]*jacobtot_inv[0]*hamil[16]+0.2651650429449552*b_x[3]*jacobtot_inv[5]*hamil[8]+0.2651650429449552*jacobtot_inv[3]*b_x[5]*hamil[8]+0.2651650429449552*b_x[0]*jacobtot_inv[1]*hamil[8]+0.2651650429449552*jacobtot_inv[0]*b_x[1]*hamil[8]-0.2755675960631073*b_x[5]*jacobtot_inv[5]*hamil[7]-0.1530931089239486*b_x[3]*jacobtot_inv[3]*hamil[7]-0.2755675960631073*b_x[1]*jacobtot_inv[1]*hamil[7]-0.1530931089239486*b_x[0]*jacobtot_inv[0]*hamil[7]-0.1530931089239486*b_x[3]*hamil[3]*jacobtot_inv[5]-0.1530931089239486*hamil[3]*jacobtot_inv[3]*b_x[5]-0.1530931089239486*b_x[0]*jacobtot_inv[1]*hamil[3]-0.1530931089239486*jacobtot_inv[0]*b_x[1]*hamil[3])*rdz2+((-0.3423265984407287*b_z[5]*jacobtot_inv[5]*hamil[32])-0.3423265984407287*b_z[1]*jacobtot_inv[1]*hamil[32]-0.2651650429449552*b_z[0]*jacobtot_inv[5]*hamil[16]-0.2651650429449552*jacobtot_inv[0]*b_z[5]*hamil[16]-0.2651650429449552*b_z[1]*jacobtot_inv[3]*hamil[16]-0.2651650429449552*jacobtot_inv[1]*b_z[3]*hamil[16]+0.1530931089239486*b_z[0]*jacobtot_inv[5]*hamil[7]+0.1530931089239486*jacobtot_inv[0]*b_z[5]*hamil[7]+0.1530931089239486*b_z[1]*jacobtot_inv[3]*hamil[7]+0.1530931089239486*jacobtot_inv[1]*b_z[3]*hamil[7]-0.2651650429449552*b_z[3]*jacobtot_inv[5]*hamil[6]-0.2651650429449552*jacobtot_inv[3]*b_z[5]*hamil[6]-0.2651650429449552*b_z[0]*jacobtot_inv[1]*hamil[6]-0.2651650429449552*jacobtot_inv[0]*b_z[1]*hamil[6]+0.1530931089239486*hamil[1]*b_z[3]*jacobtot_inv[5]+0.1530931089239486*hamil[1]*jacobtot_inv[3]*b_z[5]+0.1530931089239486*b_z[0]*hamil[1]*jacobtot_inv[1]+0.1530931089239486*jacobtot_inv[0]*b_z[1]*hamil[1])*rdx2)/q_; 
  alphaL1[2] = (((0.2651650429449552*hamil[4]*b_x[5]*jacobtot_inv[5]+0.2651650429449552*b_x[3]*jacobtot_inv[3]*hamil[4])*rdvpar2*rdz2+((-0.2651650429449552*jacobtot_inv[0]*hamil[4]*b_z[5])-0.2651650429449552*b_z[1]*jacobtot_inv[3]*hamil[4])*rdvpar2*rdx2)*wvpar)/q_+((0.3423265984407287*b_x[5]*jacobtot_inv[5]*hamil[32]+0.3423265984407287*b_x[3]*jacobtot_inv[3]*hamil[32]+0.2651650429449552*b_x[0]*jacobtot_inv[5]*hamil[16]+0.2651650429449552*jacobtot_inv[0]*b_x[5]*hamil[16]+0.2651650429449552*b_x[1]*jacobtot_inv[3]*hamil[16]+0.2651650429449552*jacobtot_inv[1]*b_x[3]*hamil[16]+0.2651650429449552*b_x[1]*jacobtot_inv[5]*hamil[8]+0.2651650429449552*jacobtot_inv[1]*b_x[5]*hamil[8]+0.2651650429449552*b_x[0]*jacobtot_inv[3]*hamil[8]+0.2651650429449552*jacobtot_inv[0]*b_x[3]*hamil[8]-0.1530931089239486*b_x[0]*jacobtot_inv[5]*hamil[7]-0.1530931089239486*jacobtot_inv[0]*b_x[5]*hamil[7]-0.1530931089239486*b_x[1]*jacobtot_inv[3]*hamil[7]-0.1530931089239486*jacobtot_inv[1]*b_x[3]*hamil[7]-0.1530931089239486*b_x[1]*hamil[3]*jacobtot_inv[5]-0.1530931089239486*jacobtot_inv[1]*hamil[3]*b_x[5]-0.1530931089239486*b_x[0]*hamil[3]*jacobtot_inv[3]-0.1530931089239486*jacobtot_inv[0]*b_x[3]*hamil[3])*rdz2+((-0.3423265984407287*jacobtot_inv[0]*b_z[5]*hamil[32])-0.3423265984407287*b_z[1]*jacobtot_inv[3]*hamil[32]-0.4772970773009194*b_z[5]*jacobtot_inv[5]*hamil[16]-0.4772970773009194*b_z[3]*jacobtot_inv[3]*hamil[16]-0.2651650429449552*b_z[1]*jacobtot_inv[1]*hamil[16]-0.2651650429449552*b_z[0]*jacobtot_inv[0]*hamil[16]+0.2755675960631073*b_z[5]*jacobtot_inv[5]*hamil[7]+0.2755675960631073*b_z[3]*jacobtot_inv[3]*hamil[7]+0.1530931089239486*b_z[1]*jacobtot_inv[1]*hamil[7]+0.1530931089239486*b_z[0]*jacobtot_inv[0]*hamil[7]-0.2651650429449552*b_z[1]*jacobtot_inv[5]*hamil[6]-0.2651650429449552*jacobtot_inv[1]*b_z[5]*hamil[6]-0.2651650429449552*b_z[0]*jacobtot_inv[3]*hamil[6]-0.2651650429449552*jacobtot_inv[0]*b_z[3]*hamil[6]+0.1530931089239486*b_z[1]*hamil[1]*jacobtot_inv[5]+0.1530931089239486*hamil[1]*jacobtot_inv[1]*b_z[5]+0.1530931089239486*b_z[0]*hamil[1]*jacobtot_inv[3]+0.1530931089239486*jacobtot_inv[0]*hamil[1]*b_z[3])*rdx2)/q_; 
  alphaL1[3] = (((0.592927061281571*jacobtot_inv[1]*b_x[5]*hamil[32]+0.592927061281571*jacobtot_inv[0]*b_x[3]*hamil[32])*rdvpar2*rdz2+((-0.592927061281571*jacobtot_inv[3]*b_z[5]*hamil[32])-0.592927061281571*jacobtot_inv[0]*b_z[1]*hamil[32])*rdvpar2*rdx2)*wvpar)/q_+((0.1530931089239486*jacobtot_inv[1]*hamil[4]*b_x[5]+0.1530931089239486*jacobtot_inv[0]*b_x[3]*hamil[4])*rdz2+((-0.1530931089239486*jacobtot_inv[3]*hamil[4]*b_z[5])-0.1530931089239486*jacobtot_inv[0]*b_z[1]*hamil[4])*rdx2)/q_; 
  alphaL1[4] = (((-0.1530931089239486*b_x[3]*jacobtot_inv[5]*hamil[21])-0.1530931089239486*jacobtot_inv[3]*b_x[5]*hamil[21]-0.1530931089239486*b_x[0]*jacobtot_inv[1]*hamil[21]-0.1530931089239486*jacobtot_inv[0]*b_x[1]*hamil[21]-0.1530931089239486*b_x[5]*jacobtot_inv[5]*hamil[14]-0.1530931089239486*b_x[3]*jacobtot_inv[3]*hamil[14]-0.1530931089239486*b_x[1]*jacobtot_inv[1]*hamil[14]-0.1530931089239486*b_x[0]*jacobtot_inv[0]*hamil[14])*rdz2+(0.1530931089239486*b_z[1]*jacobtot_inv[5]*hamil[21]+0.1530931089239486*jacobtot_inv[1]*b_z[5]*hamil[21]+0.1530931089239486*b_z[0]*jacobtot_inv[3]*hamil[21]+0.1530931089239486*jacobtot_inv[0]*b_z[3]*hamil[21]+0.1530931089239486*b_z[5]*jacobtot_inv[5]*hamil[12]+0.1530931089239486*b_z[3]*jacobtot_inv[3]*hamil[12]+0.1530931089239486*b_z[1]*jacobtot_inv[1]*hamil[12]+0.1530931089239486*b_z[0]*jacobtot_inv[0]*hamil[12])*rdx2)/q_; 
  alphaL1[5] = (((0.2651650429449552*b_x[3]*hamil[4]*jacobtot_inv[5]+0.2651650429449552*jacobtot_inv[3]*hamil[4]*b_x[5])*rdvpar2*rdz2+((-0.2651650429449552*b_z[1]*hamil[4]*jacobtot_inv[5])-0.2651650429449552*jacobtot_inv[1]*hamil[4]*b_z[5])*rdvpar2*rdx2)*wvpar)/q_+((0.3423265984407287*b_x[3]*jacobtot_inv[5]*hamil[32]+0.3423265984407287*jacobtot_inv[3]*b_x[5]*hamil[32]+0.4772970773009194*b_x[1]*jacobtot_inv[5]*hamil[16]+0.4772970773009194*jacobtot_inv[1]*b_x[5]*hamil[16]+0.2651650429449552*b_x[0]*jacobtot_inv[3]*hamil[16]+0.2651650429449552*jacobtot_inv[0]*b_x[3]*hamil[16]+0.2651650429449552*b_x[0]*jacobtot_inv[5]*hamil[8]+0.2651650429449552*jacobtot_inv[0]*b_x[5]*hamil[8]+0.2651650429449552*b_x[1]*jacobtot_inv[3]*hamil[8]+0.2651650429449552*jacobtot_inv[1]*b_x[3]*hamil[8]-0.2755675960631073*b_x[1]*jacobtot_inv[5]*hamil[7]-0.2755675960631073*jacobtot_inv[1]*b_x[5]*hamil[7]-0.1530931089239486*b_x[0]*jacobtot_inv[3]*hamil[7]-0.1530931089239486*jacobtot_inv[0]*b_x[3]*hamil[7]-0.1530931089239486*b_x[0]*hamil[3]*jacobtot_inv[5]-0.1530931089239486*jacobtot_inv[0]*hamil[3]*b_x[5]-0.1530931089239486*b_x[1]*hamil[3]*jacobtot_inv[3]-0.1530931089239486*jacobtot_inv[1]*b_x[3]*hamil[3])*rdz2+((-0.3423265984407287*b_z[1]*jacobtot_inv[5]*hamil[32])-0.3423265984407287*jacobtot_inv[1]*b_z[5]*hamil[32]-0.4772970773009194*b_z[3]*jacobtot_inv[5]*hamil[16]-0.4772970773009194*jacobtot_inv[3]*b_z[5]*hamil[16]-0.2651650429449552*b_z[0]*jacobtot_inv[1]*hamil[16]-0.2651650429449552*jacobtot_inv[0]*b_z[1]*hamil[16]+0.2755675960631073*b_z[3]*jacobtot_inv[5]*hamil[7]+0.2755675960631073*jacobtot_inv[3]*b_z[5]*hamil[7]+0.1530931089239486*b_z[0]*jacobtot_inv[1]*hamil[7]+0.1530931089239486*jacobtot_inv[0]*b_z[1]*hamil[7]-0.2651650429449552*b_z[0]*jacobtot_inv[5]*hamil[6]-0.2651650429449552*jacobtot_inv[0]*b_z[5]*hamil[6]-0.2651650429449552*b_z[1]*jacobtot_inv[3]*hamil[6]-0.2651650429449552*jacobtot_inv[1]*b_z[3]*hamil[6]+0.1530931089239486*b_z[0]*hamil[1]*jacobtot_inv[5]+0.1530931089239486*jacobtot_inv[0]*hamil[1]*b_z[5]+0.1530931089239486*b_z[1]*hamil[1]*jacobtot_inv[3]+0.1530931089239486*hamil[1]*jacobtot_inv[1]*b_z[3])*rdx2)/q_; 
  alphaL1[6] = (((0.592927061281571*jacobtot_inv[0]*b_x[5]*hamil[32]+0.592927061281571*jacobtot_inv[1]*b_x[3]*hamil[32])*rdvpar2*rdz2+((-0.592927061281571*b_z[5]*jacobtot_inv[5]*hamil[32])-0.592927061281571*b_z[1]*jacobtot_inv[1]*hamil[32])*rdvpar2*rdx2)*wvpar)/q_+((0.1530931089239486*jacobtot_inv[0]*hamil[4]*b_x[5]+0.1530931089239486*jacobtot_inv[1]*b_x[3]*hamil[4])*rdz2+((-0.1530931089239486*hamil[4]*b_z[5]*jacobtot_inv[5])-0.1530931089239486*b_z[1]*jacobtot_inv[1]*hamil[4])*rdx2)/q_; 
  alphaL1[7] = (((0.592927061281571*b_x[5]*jacobtot_inv[5]*hamil[32]+0.592927061281571*b_x[3]*jacobtot_inv[3]*hamil[32])*rdvpar2*rdz2+((-0.592927061281571*jacobtot_inv[0]*b_z[5]*hamil[32])-0.592927061281571*b_z[1]*jacobtot_inv[3]*hamil[32])*rdvpar2*rdx2)*wvpar)/q_+((0.1530931089239486*hamil[4]*b_x[5]*jacobtot_inv[5]+0.1530931089239486*b_x[3]*jacobtot_inv[3]*hamil[4])*rdz2+((-0.1530931089239486*jacobtot_inv[0]*hamil[4]*b_z[5])-0.1530931089239486*b_z[1]*jacobtot_inv[3]*hamil[4])*rdx2)/q_; 
  alphaL1[8] = (((-0.2755675960631073*b_x[5]*jacobtot_inv[5]*hamil[21])-0.1530931089239486*b_x[3]*jacobtot_inv[3]*hamil[21]-0.2755675960631073*b_x[1]*jacobtot_inv[1]*hamil[21]-0.1530931089239486*b_x[0]*jacobtot_inv[0]*hamil[21]-0.1530931089239486*b_x[3]*jacobtot_inv[5]*hamil[14]-0.1530931089239486*jacobtot_inv[3]*b_x[5]*hamil[14]-0.1530931089239486*b_x[0]*jacobtot_inv[1]*hamil[14]-0.1530931089239486*jacobtot_inv[0]*b_x[1]*hamil[14])*rdz2+(0.1530931089239486*b_z[0]*jacobtot_inv[5]*hamil[21]+0.1530931089239486*jacobtot_inv[0]*b_z[5]*hamil[21]+0.1530931089239486*b_z[1]*jacobtot_inv[3]*hamil[21]+0.1530931089239486*jacobtot_inv[1]*b_z[3]*hamil[21]+0.1530931089239486*b_z[3]*jacobtot_inv[5]*hamil[12]+0.1530931089239486*jacobtot_inv[3]*b_z[5]*hamil[12]+0.1530931089239486*b_z[0]*jacobtot_inv[1]*hamil[12]+0.1530931089239486*jacobtot_inv[0]*b_z[1]*hamil[12])*rdx2)/q_; 
  alphaL1[9] = (((-0.1530931089239486*b_x[0]*jacobtot_inv[5]*hamil[21])-0.1530931089239486*jacobtot_inv[0]*b_x[5]*hamil[21]-0.1530931089239486*b_x[1]*jacobtot_inv[3]*hamil[21]-0.1530931089239486*jacobtot_inv[1]*b_x[3]*hamil[21]-0.1530931089239486*b_x[1]*jacobtot_inv[5]*hamil[14]-0.1530931089239486*jacobtot_inv[1]*b_x[5]*hamil[14]-0.1530931089239486*b_x[0]*jacobtot_inv[3]*hamil[14]-0.1530931089239486*jacobtot_inv[0]*b_x[3]*hamil[14])*rdz2+(0.2755675960631073*b_z[5]*jacobtot_inv[5]*hamil[21]+0.2755675960631073*b_z[3]*jacobtot_inv[3]*hamil[21]+0.1530931089239486*b_z[1]*jacobtot_inv[1]*hamil[21]+0.1530931089239486*b_z[0]*jacobtot_inv[0]*hamil[21]+0.1530931089239486*b_z[1]*jacobtot_inv[5]*hamil[12]+0.1530931089239486*jacobtot_inv[1]*b_z[5]*hamil[12]+0.1530931089239486*b_z[0]*jacobtot_inv[3]*hamil[12]+0.1530931089239486*jacobtot_inv[0]*b_z[3]*hamil[12])*rdx2)/q_; 
  alphaL1[11] = (((0.592927061281571*b_x[3]*jacobtot_inv[5]*hamil[32]+0.592927061281571*jacobtot_inv[3]*b_x[5]*hamil[32])*rdvpar2*rdz2+((-0.592927061281571*b_z[1]*jacobtot_inv[5]*hamil[32])-0.592927061281571*jacobtot_inv[1]*b_z[5]*hamil[32])*rdvpar2*rdx2)*wvpar)/q_+((0.1530931089239486*b_x[3]*hamil[4]*jacobtot_inv[5]+0.1530931089239486*jacobtot_inv[3]*hamil[4]*b_x[5])*rdz2+((-0.1530931089239486*b_z[1]*hamil[4]*jacobtot_inv[5])-0.1530931089239486*jacobtot_inv[1]*hamil[4]*b_z[5])*rdx2)/q_; 
  alphaL1[12] = (((-0.2755675960631073*b_x[1]*jacobtot_inv[5]*hamil[21])-0.2755675960631073*jacobtot_inv[1]*b_x[5]*hamil[21]-0.1530931089239486*b_x[0]*jacobtot_inv[3]*hamil[21]-0.1530931089239486*jacobtot_inv[0]*b_x[3]*hamil[21]-0.1530931089239486*b_x[0]*jacobtot_inv[5]*hamil[14]-0.1530931089239486*jacobtot_inv[0]*b_x[5]*hamil[14]-0.1530931089239486*b_x[1]*jacobtot_inv[3]*hamil[14]-0.1530931089239486*jacobtot_inv[1]*b_x[3]*hamil[14])*rdz2+(0.2755675960631073*b_z[3]*jacobtot_inv[5]*hamil[21]+0.2755675960631073*jacobtot_inv[3]*b_z[5]*hamil[21]+0.1530931089239486*b_z[0]*jacobtot_inv[1]*hamil[21]+0.1530931089239486*jacobtot_inv[0]*b_z[1]*hamil[21]+0.1530931089239486*b_z[0]*jacobtot_inv[5]*hamil[12]+0.1530931089239486*jacobtot_inv[0]*b_z[5]*hamil[12]+0.1530931089239486*b_z[1]*jacobtot_inv[3]*hamil[12]+0.1530931089239486*jacobtot_inv[1]*b_z[3]*hamil[12])*rdx2)/q_; 
  alphaL1[16] = ((0.3061862178478971*jacobtot_inv[1]*b_x[5]*hamil[32]+0.3061862178478971*jacobtot_inv[0]*b_x[3]*hamil[32])*rdz2+((-0.3061862178478971*jacobtot_inv[3]*b_z[5]*hamil[32])-0.3061862178478971*jacobtot_inv[0]*b_z[1]*hamil[32])*rdx2)/q_; 
  alphaL1[17] = ((0.3061862178478971*jacobtot_inv[0]*b_x[5]*hamil[32]+0.3061862178478971*jacobtot_inv[1]*b_x[3]*hamil[32])*rdz2+((-0.3061862178478971*b_z[5]*jacobtot_inv[5]*hamil[32])-0.3061862178478971*b_z[1]*jacobtot_inv[1]*hamil[32])*rdx2)/q_; 
  alphaL1[18] = ((0.3061862178478971*b_x[5]*jacobtot_inv[5]*hamil[32]+0.3061862178478971*b_x[3]*jacobtot_inv[3]*hamil[32])*rdz2+((-0.3061862178478971*jacobtot_inv[0]*b_z[5]*hamil[32])-0.3061862178478971*b_z[1]*jacobtot_inv[3]*hamil[32])*rdx2)/q_; 
  alphaL1[20] = ((0.3061862178478971*b_x[3]*jacobtot_inv[5]*hamil[32]+0.3061862178478971*jacobtot_inv[3]*b_x[5]*hamil[32])*rdz2+((-0.3061862178478971*b_z[1]*jacobtot_inv[5]*hamil[32])-0.3061862178478971*jacobtot_inv[1]*b_z[5]*hamil[32])*rdx2)/q_; 

  double *sgn_alpha_surf1 = &sgn_alpha_surf[24];
  int const_sgn_alpha_surf = 1; 
  
  if (0.2236067977499786*alphaL1[20]-0.2236067977499786*(alphaL1[18]+alphaL1[17])+0.2236067977499786*alphaL1[16]-0.25*alphaL1[12]-0.3354101966249678*alphaL1[11]+0.25*(alphaL1[9]+alphaL1[8])+0.3354101966249678*(alphaL1[7]+alphaL1[6])+0.25*alphaL1[5]-0.25*alphaL1[4]-0.3354101966249678*alphaL1[3]-0.25*(alphaL1[2]+alphaL1[1])+0.25*alphaL1[0] > 0.) 
    sgn_alpha_surf1[0] = 1.0; 
  else  
    sgn_alpha_surf1[0] = -1.0; 
  
  if ((-0.2795084971874732*alphaL1[20])+0.2795084971874732*(alphaL1[18]+alphaL1[17])-0.2795084971874732*alphaL1[16]-0.25*alphaL1[12]+0.25*(alphaL1[9]+alphaL1[8]+alphaL1[5])-0.25*(alphaL1[4]+alphaL1[2]+alphaL1[1])+0.25*alphaL1[0] > 0.) 
    sgn_alpha_surf1[1] = 1.0; 
  else  
    sgn_alpha_surf1[1] = -1.0; 
  
  if (sgn_alpha_surf1[1] == sgn_alpha_surf1[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*alphaL1[20]-0.2236067977499786*(alphaL1[18]+alphaL1[17])+0.2236067977499786*alphaL1[16]-0.25*alphaL1[12]+0.3354101966249678*alphaL1[11]+0.25*(alphaL1[9]+alphaL1[8])-0.3354101966249678*(alphaL1[7]+alphaL1[6])+0.25*alphaL1[5]-0.25*alphaL1[4]+0.3354101966249678*alphaL1[3]-0.25*(alphaL1[2]+alphaL1[1])+0.25*alphaL1[0] > 0.) 
    sgn_alpha_surf1[2] = 1.0; 
  else  
    sgn_alpha_surf1[2] = -1.0; 
  
  if (sgn_alpha_surf1[2] == sgn_alpha_surf1[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*alphaL1[20]-0.2236067977499786*(alphaL1[18]+alphaL1[17])+0.2236067977499786*alphaL1[16]+0.25*alphaL1[12]-0.3354101966249678*alphaL1[11]-0.25*(alphaL1[9]+alphaL1[8])+0.3354101966249678*(alphaL1[7]+alphaL1[6])+0.25*(alphaL1[5]+alphaL1[4])-0.3354101966249678*alphaL1[3]-0.25*(alphaL1[2]+alphaL1[1])+0.25*alphaL1[0] > 0.) 
    sgn_alpha_surf1[3] = 1.0; 
  else  
    sgn_alpha_surf1[3] = -1.0; 
  
  if (sgn_alpha_surf1[3] == sgn_alpha_surf1[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*alphaL1[20])+0.2795084971874732*(alphaL1[18]+alphaL1[17])-0.2795084971874732*alphaL1[16]+0.25*alphaL1[12]-0.25*(alphaL1[9]+alphaL1[8])+0.25*(alphaL1[5]+alphaL1[4])-0.25*(alphaL1[2]+alphaL1[1])+0.25*alphaL1[0] > 0.) 
    sgn_alpha_surf1[4] = 1.0; 
  else  
    sgn_alpha_surf1[4] = -1.0; 
  
  if (sgn_alpha_surf1[4] == sgn_alpha_surf1[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*alphaL1[20]-0.2236067977499786*(alphaL1[18]+alphaL1[17])+0.2236067977499786*alphaL1[16]+0.25*alphaL1[12]+0.3354101966249678*alphaL1[11]-0.25*(alphaL1[9]+alphaL1[8])-0.3354101966249678*(alphaL1[7]+alphaL1[6])+0.25*(alphaL1[5]+alphaL1[4])+0.3354101966249678*alphaL1[3]-0.25*(alphaL1[2]+alphaL1[1])+0.25*alphaL1[0] > 0.) 
    sgn_alpha_surf1[5] = 1.0; 
  else  
    sgn_alpha_surf1[5] = -1.0; 
  
  if (sgn_alpha_surf1[5] == sgn_alpha_surf1[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL1[20])+0.2236067977499786*alphaL1[18]-0.2236067977499786*alphaL1[17]+0.2236067977499786*alphaL1[16]+0.25*alphaL1[12]+0.3354101966249678*alphaL1[11]-0.25*alphaL1[9]+0.25*alphaL1[8]-0.3354101966249678*alphaL1[7]+0.3354101966249678*alphaL1[6]-0.25*(alphaL1[5]+alphaL1[4])-0.3354101966249678*alphaL1[3]+0.25*alphaL1[2]-0.25*alphaL1[1]+0.25*alphaL1[0] > 0.) 
    sgn_alpha_surf1[6] = 1.0; 
  else  
    sgn_alpha_surf1[6] = -1.0; 
  
  if (sgn_alpha_surf1[6] == sgn_alpha_surf1[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL1[20]-0.2795084971874732*alphaL1[18]+0.2795084971874732*alphaL1[17]-0.2795084971874732*alphaL1[16]+0.25*alphaL1[12]-0.25*alphaL1[9]+0.25*alphaL1[8]-0.25*(alphaL1[5]+alphaL1[4])+0.25*alphaL1[2]-0.25*alphaL1[1]+0.25*alphaL1[0] > 0.) 
    sgn_alpha_surf1[7] = 1.0; 
  else  
    sgn_alpha_surf1[7] = -1.0; 
  
  if (sgn_alpha_surf1[7] == sgn_alpha_surf1[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL1[20])+0.2236067977499786*alphaL1[18]-0.2236067977499786*alphaL1[17]+0.2236067977499786*alphaL1[16]+0.25*alphaL1[12]-0.3354101966249678*alphaL1[11]-0.25*alphaL1[9]+0.25*alphaL1[8]+0.3354101966249678*alphaL1[7]-0.3354101966249678*alphaL1[6]-0.25*(alphaL1[5]+alphaL1[4])+0.3354101966249678*alphaL1[3]+0.25*alphaL1[2]-0.25*alphaL1[1]+0.25*alphaL1[0] > 0.) 
    sgn_alpha_surf1[8] = 1.0; 
  else  
    sgn_alpha_surf1[8] = -1.0; 
  
  if (sgn_alpha_surf1[8] == sgn_alpha_surf1[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL1[20])+0.2236067977499786*alphaL1[18]-0.2236067977499786*alphaL1[17]+0.2236067977499786*alphaL1[16]-0.25*alphaL1[12]+0.3354101966249678*alphaL1[11]+0.25*alphaL1[9]-0.25*alphaL1[8]-0.3354101966249678*alphaL1[7]+0.3354101966249678*alphaL1[6]-0.25*alphaL1[5]+0.25*alphaL1[4]-0.3354101966249678*alphaL1[3]+0.25*alphaL1[2]-0.25*alphaL1[1]+0.25*alphaL1[0] > 0.) 
    sgn_alpha_surf1[9] = 1.0; 
  else  
    sgn_alpha_surf1[9] = -1.0; 
  
  if (sgn_alpha_surf1[9] == sgn_alpha_surf1[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL1[20]-0.2795084971874732*alphaL1[18]+0.2795084971874732*alphaL1[17]-0.2795084971874732*alphaL1[16]-0.25*alphaL1[12]+0.25*alphaL1[9]-0.25*(alphaL1[8]+alphaL1[5])+0.25*(alphaL1[4]+alphaL1[2])-0.25*alphaL1[1]+0.25*alphaL1[0] > 0.) 
    sgn_alpha_surf1[10] = 1.0; 
  else  
    sgn_alpha_surf1[10] = -1.0; 
  
  if (sgn_alpha_surf1[10] == sgn_alpha_surf1[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*alphaL1[20])+0.2236067977499786*alphaL1[18]-0.2236067977499786*alphaL1[17]+0.2236067977499786*alphaL1[16]-0.25*alphaL1[12]-0.3354101966249678*alphaL1[11]+0.25*alphaL1[9]-0.25*alphaL1[8]+0.3354101966249678*alphaL1[7]-0.3354101966249678*alphaL1[6]-0.25*alphaL1[5]+0.25*alphaL1[4]+0.3354101966249678*alphaL1[3]+0.25*alphaL1[2]-0.25*alphaL1[1]+0.25*alphaL1[0] > 0.) 
    sgn_alpha_surf1[11] = 1.0; 
  else  
    sgn_alpha_surf1[11] = -1.0; 
  
  if (sgn_alpha_surf1[11] == sgn_alpha_surf1[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*(alphaL1[20]+alphaL1[18]))+0.2236067977499786*(alphaL1[17]+alphaL1[16])+0.25*alphaL1[12]+0.3354101966249678*alphaL1[11]+0.25*alphaL1[9]-0.25*alphaL1[8]+0.3354101966249678*alphaL1[7]-0.3354101966249678*alphaL1[6]-0.25*(alphaL1[5]+alphaL1[4])-0.3354101966249678*alphaL1[3]-0.25*alphaL1[2]+0.25*(alphaL1[1]+alphaL1[0]) > 0.) 
    sgn_alpha_surf1[12] = 1.0; 
  else  
    sgn_alpha_surf1[12] = -1.0; 
  
  if (sgn_alpha_surf1[12] == sgn_alpha_surf1[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*(alphaL1[20]+alphaL1[18])-0.2795084971874732*(alphaL1[17]+alphaL1[16])+0.25*(alphaL1[12]+alphaL1[9])-0.25*(alphaL1[8]+alphaL1[5]+alphaL1[4]+alphaL1[2])+0.25*(alphaL1[1]+alphaL1[0]) > 0.) 
    sgn_alpha_surf1[13] = 1.0; 
  else  
    sgn_alpha_surf1[13] = -1.0; 
  
  if (sgn_alpha_surf1[13] == sgn_alpha_surf1[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*(alphaL1[20]+alphaL1[18]))+0.2236067977499786*(alphaL1[17]+alphaL1[16])+0.25*alphaL1[12]-0.3354101966249678*alphaL1[11]+0.25*alphaL1[9]-0.25*alphaL1[8]-0.3354101966249678*alphaL1[7]+0.3354101966249678*alphaL1[6]-0.25*(alphaL1[5]+alphaL1[4])+0.3354101966249678*alphaL1[3]-0.25*alphaL1[2]+0.25*(alphaL1[1]+alphaL1[0]) > 0.) 
    sgn_alpha_surf1[14] = 1.0; 
  else  
    sgn_alpha_surf1[14] = -1.0; 
  
  if (sgn_alpha_surf1[14] == sgn_alpha_surf1[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*(alphaL1[20]+alphaL1[18]))+0.2236067977499786*(alphaL1[17]+alphaL1[16])-0.25*alphaL1[12]+0.3354101966249678*alphaL1[11]-0.25*alphaL1[9]+0.25*alphaL1[8]+0.3354101966249678*alphaL1[7]-0.3354101966249678*alphaL1[6]-0.25*alphaL1[5]+0.25*alphaL1[4]-0.3354101966249678*alphaL1[3]-0.25*alphaL1[2]+0.25*(alphaL1[1]+alphaL1[0]) > 0.) 
    sgn_alpha_surf1[15] = 1.0; 
  else  
    sgn_alpha_surf1[15] = -1.0; 
  
  if (sgn_alpha_surf1[15] == sgn_alpha_surf1[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*(alphaL1[20]+alphaL1[18])-0.2795084971874732*(alphaL1[17]+alphaL1[16])-0.25*(alphaL1[12]+alphaL1[9])+0.25*alphaL1[8]-0.25*alphaL1[5]+0.25*alphaL1[4]-0.25*alphaL1[2]+0.25*(alphaL1[1]+alphaL1[0]) > 0.) 
    sgn_alpha_surf1[16] = 1.0; 
  else  
    sgn_alpha_surf1[16] = -1.0; 
  
  if (sgn_alpha_surf1[16] == sgn_alpha_surf1[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2236067977499786*(alphaL1[20]+alphaL1[18]))+0.2236067977499786*(alphaL1[17]+alphaL1[16])-0.25*alphaL1[12]-0.3354101966249678*alphaL1[11]-0.25*alphaL1[9]+0.25*alphaL1[8]-0.3354101966249678*alphaL1[7]+0.3354101966249678*alphaL1[6]-0.25*alphaL1[5]+0.25*alphaL1[4]+0.3354101966249678*alphaL1[3]-0.25*alphaL1[2]+0.25*(alphaL1[1]+alphaL1[0]) > 0.) 
    sgn_alpha_surf1[17] = 1.0; 
  else  
    sgn_alpha_surf1[17] = -1.0; 
  
  if (sgn_alpha_surf1[17] == sgn_alpha_surf1[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL1[20]+alphaL1[18]+alphaL1[17]+alphaL1[16])-0.25*alphaL1[12]-0.3354101966249678*alphaL1[11]-0.25*(alphaL1[9]+alphaL1[8])-0.3354101966249678*(alphaL1[7]+alphaL1[6])+0.25*alphaL1[5]-0.25*alphaL1[4]-0.3354101966249678*alphaL1[3]+0.25*(alphaL1[2]+alphaL1[1]+alphaL1[0]) > 0.) 
    sgn_alpha_surf1[18] = 1.0; 
  else  
    sgn_alpha_surf1[18] = -1.0; 
  
  if (sgn_alpha_surf1[18] == sgn_alpha_surf1[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*(alphaL1[20]+alphaL1[18]+alphaL1[17]+alphaL1[16]))-0.25*(alphaL1[12]+alphaL1[9]+alphaL1[8])+0.25*alphaL1[5]-0.25*alphaL1[4]+0.25*(alphaL1[2]+alphaL1[1]+alphaL1[0]) > 0.) 
    sgn_alpha_surf1[19] = 1.0; 
  else  
    sgn_alpha_surf1[19] = -1.0; 
  
  if (sgn_alpha_surf1[19] == sgn_alpha_surf1[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL1[20]+alphaL1[18]+alphaL1[17]+alphaL1[16])-0.25*alphaL1[12]+0.3354101966249678*alphaL1[11]-0.25*(alphaL1[9]+alphaL1[8])+0.3354101966249678*(alphaL1[7]+alphaL1[6])+0.25*alphaL1[5]-0.25*alphaL1[4]+0.3354101966249678*alphaL1[3]+0.25*(alphaL1[2]+alphaL1[1]+alphaL1[0]) > 0.) 
    sgn_alpha_surf1[20] = 1.0; 
  else  
    sgn_alpha_surf1[20] = -1.0; 
  
  if (sgn_alpha_surf1[20] == sgn_alpha_surf1[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL1[20]+alphaL1[18]+alphaL1[17]+alphaL1[16])+0.25*alphaL1[12]-0.3354101966249678*alphaL1[11]+0.25*(alphaL1[9]+alphaL1[8])-0.3354101966249678*(alphaL1[7]+alphaL1[6])+0.25*(alphaL1[5]+alphaL1[4])-0.3354101966249678*alphaL1[3]+0.25*(alphaL1[2]+alphaL1[1]+alphaL1[0]) > 0.) 
    sgn_alpha_surf1[21] = 1.0; 
  else  
    sgn_alpha_surf1[21] = -1.0; 
  
  if (sgn_alpha_surf1[21] == sgn_alpha_surf1[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL1[12]+alphaL1[9]+alphaL1[8]+alphaL1[5]+alphaL1[4]+alphaL1[2]+alphaL1[1]+alphaL1[0])-0.2795084971874732*(alphaL1[20]+alphaL1[18]+alphaL1[17]+alphaL1[16]) > 0.) 
    sgn_alpha_surf1[22] = 1.0; 
  else  
    sgn_alpha_surf1[22] = -1.0; 
  
  if (sgn_alpha_surf1[22] == sgn_alpha_surf1[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2236067977499786*(alphaL1[20]+alphaL1[18]+alphaL1[17]+alphaL1[16])+0.25*alphaL1[12]+0.3354101966249678*alphaL1[11]+0.25*(alphaL1[9]+alphaL1[8])+0.3354101966249678*(alphaL1[7]+alphaL1[6])+0.25*(alphaL1[5]+alphaL1[4])+0.3354101966249678*alphaL1[3]+0.25*(alphaL1[2]+alphaL1[1]+alphaL1[0]) > 0.) 
    sgn_alpha_surf1[23] = 1.0; 
  else  
    sgn_alpha_surf1[23] = -1.0; 
  
  if (sgn_alpha_surf1[23] == sgn_alpha_surf1[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
