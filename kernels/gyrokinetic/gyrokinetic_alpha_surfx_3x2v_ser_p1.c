#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_alpha_surfx_3x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *Bstar_Bmag, double* GKYL_RESTRICT alpha_surf) 
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
  //             Note: Each cell owns their *lower* edge surface evaluation (evaluated at -1).

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

  const double *BstarXdBmag = &Bstar_Bmag[0]; 

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

  double *alphaL0 = &alpha_surf[0];
  alphaL0[0] = ((((0.4592793267718456*b_y[3]-0.7954951288348656*b_y[5])*jacobtot_inv[5]+jacobtot_inv[3]*(0.4592793267718456*b_y[5]-0.2651650429449552*b_y[3])+(0.4592793267718456*b_y[0]-0.7954951288348656*b_y[1])*jacobtot_inv[1]+jacobtot_inv[0]*(0.4592793267718456*b_y[1]-0.2651650429449552*b_y[0]))*hamil[7]+hamil[3]*((0.4592793267718456*b_y[5]-0.2651650429449552*b_y[3])*jacobtot_inv[5]+jacobtot_inv[3]*(0.1530931089239486*b_y[3]-0.2651650429449552*b_y[5])+(0.4592793267718456*b_y[1]-0.2651650429449552*b_y[0])*jacobtot_inv[1]+jacobtot_inv[0]*(0.1530931089239486*b_y[0]-0.2651650429449552*b_y[1])))*rdz2+(((0.7954951288348656*b_z[1]-0.4592793267718456*b_z[0])*jacobtot_inv[5]+(0.7954951288348656*jacobtot_inv[1]-0.4592793267718456*jacobtot_inv[0])*b_z[5]+(0.2651650429449552*b_z[0]-0.4592793267718456*b_z[1])*jacobtot_inv[3]+(0.2651650429449552*jacobtot_inv[0]-0.4592793267718456*jacobtot_inv[1])*b_z[3])*hamil[16]+((0.2651650429449552*b_z[0]-0.4592793267718456*b_z[1])*jacobtot_inv[5]+(0.2651650429449552*jacobtot_inv[0]-0.4592793267718456*jacobtot_inv[1])*b_z[5]+(0.2651650429449552*b_z[1]-0.1530931089239486*b_z[0])*jacobtot_inv[3]+(0.2651650429449552*jacobtot_inv[1]-0.1530931089239486*jacobtot_inv[0])*b_z[3])*hamil[8]+((0.7954951288348656*b_z[5]-0.4592793267718456*b_z[3])*jacobtot_inv[5]+jacobtot_inv[3]*(0.2651650429449552*b_z[3]-0.4592793267718456*b_z[5])+(0.7954951288348656*b_z[1]-0.4592793267718456*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(0.2651650429449552*b_z[0]-0.4592793267718456*b_z[1]))*hamil[6]+hamil[2]*((0.2651650429449552*b_z[3]-0.4592793267718456*b_z[5])*jacobtot_inv[5]+jacobtot_inv[3]*(0.2651650429449552*b_z[5]-0.1530931089239486*b_z[3])+(0.2651650429449552*b_z[0]-0.4592793267718456*b_z[1])*jacobtot_inv[1]+jacobtot_inv[0]*(0.2651650429449552*b_z[1]-0.1530931089239486*b_z[0])))*rdy2)/q_+(((0.9682458365518543*BstarXdBmag[3]-1.677050983124842*BstarXdBmag[5])*hamil[32]+(0.4330127018922193*BstarXdBmag[0]-0.75*BstarXdBmag[1])*hamil[4])*rdvpar2)/m_; 
  alphaL0[1] = ((((0.4592793267718456*b_y[3]-0.7954951288348656*b_y[5])*jacobtot_inv[5]+jacobtot_inv[3]*(0.4592793267718456*b_y[5]-0.2651650429449552*b_y[3])+(0.4592793267718456*b_y[0]-0.7954951288348656*b_y[1])*jacobtot_inv[1]+jacobtot_inv[0]*(0.4592793267718456*b_y[1]-0.2651650429449552*b_y[0]))*hamil[16]+((0.4592793267718456*b_y[5]-0.2651650429449552*b_y[3])*jacobtot_inv[5]+jacobtot_inv[3]*(0.1530931089239486*b_y[3]-0.2651650429449552*b_y[5])+(0.4592793267718456*b_y[1]-0.2651650429449552*b_y[0])*jacobtot_inv[1]+jacobtot_inv[0]*(0.1530931089239486*b_y[0]-0.2651650429449552*b_y[1]))*hamil[8])*rdz2)/q_; 
  alphaL0[2] = ((((0.4592793267718456*b_y[0]-0.7954951288348656*b_y[1])*jacobtot_inv[5]+(0.4592793267718456*jacobtot_inv[0]-0.7954951288348656*jacobtot_inv[1])*b_y[5]+(0.4592793267718456*b_y[1]-0.2651650429449552*b_y[0])*jacobtot_inv[3]+(0.4592793267718456*jacobtot_inv[1]-0.2651650429449552*jacobtot_inv[0])*b_y[3])*hamil[7]+hamil[3]*((0.4592793267718456*b_y[1]-0.2651650429449552*b_y[0])*jacobtot_inv[5]+(0.4592793267718456*jacobtot_inv[1]-0.2651650429449552*jacobtot_inv[0])*b_y[5]+(0.1530931089239486*b_y[0]-0.2651650429449552*b_y[1])*jacobtot_inv[3]+(0.1530931089239486*jacobtot_inv[0]-0.2651650429449552*jacobtot_inv[1])*b_y[3]))*rdz2+(((1.431891231902758*b_z[5]-0.826702788189322*b_z[3])*jacobtot_inv[5]+jacobtot_inv[3]*(0.4772970773009194*b_z[3]-0.826702788189322*b_z[5])+(0.7954951288348656*b_z[1]-0.4592793267718456*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*(0.2651650429449552*b_z[0]-0.4592793267718456*b_z[1]))*hamil[16]+((0.4772970773009194*b_z[3]-0.826702788189322*b_z[5])*jacobtot_inv[5]+jacobtot_inv[3]*(0.4772970773009194*b_z[5]-0.2755675960631073*b_z[3])+(0.2651650429449552*b_z[0]-0.4592793267718456*b_z[1])*jacobtot_inv[1]+jacobtot_inv[0]*(0.2651650429449552*b_z[1]-0.1530931089239486*b_z[0]))*hamil[8]+((0.7954951288348656*b_z[1]-0.4592793267718456*b_z[0])*jacobtot_inv[5]+(0.7954951288348656*jacobtot_inv[1]-0.4592793267718456*jacobtot_inv[0])*b_z[5]+(0.2651650429449552*b_z[0]-0.4592793267718456*b_z[1])*jacobtot_inv[3]+(0.2651650429449552*jacobtot_inv[0]-0.4592793267718456*jacobtot_inv[1])*b_z[3])*hamil[6]+hamil[2]*((0.2651650429449552*b_z[0]-0.4592793267718456*b_z[1])*jacobtot_inv[5]+(0.2651650429449552*jacobtot_inv[0]-0.4592793267718456*jacobtot_inv[1])*b_z[5]+(0.2651650429449552*b_z[1]-0.1530931089239486*b_z[0])*jacobtot_inv[3]+(0.2651650429449552*jacobtot_inv[1]-0.1530931089239486*jacobtot_inv[0])*b_z[3]))*rdy2)/q_+(((0.9682458365518543*BstarXdBmag[6]-1.677050983124842*BstarXdBmag[7])*hamil[32]+(0.4330127018922193*BstarXdBmag[2]-0.75*BstarXdBmag[4])*hamil[4])*rdvpar2)/m_; 
  alphaL0[3] = (((0.9682458365518543*BstarXdBmag[0]-1.677050983124842*BstarXdBmag[1])*hamil[32]+hamil[4]*(0.4330127018922193*BstarXdBmag[3]-0.75*BstarXdBmag[5]))*rdvpar2)/m_; 
  alphaL0[4] = ((((0.4592793267718456*b_y[3]-0.7954951288348656*b_y[5])*jacobtot_inv[5]+jacobtot_inv[3]*(0.4592793267718456*b_y[5]-0.2651650429449552*b_y[3])+(0.4592793267718456*b_y[0]-0.7954951288348656*b_y[1])*jacobtot_inv[1]+jacobtot_inv[0]*(0.4592793267718456*b_y[1]-0.2651650429449552*b_y[0]))*hamil[21]+((0.4592793267718456*b_y[5]-0.2651650429449552*b_y[3])*jacobtot_inv[5]+jacobtot_inv[3]*(0.1530931089239486*b_y[3]-0.2651650429449552*b_y[5])+(0.4592793267718456*b_y[1]-0.2651650429449552*b_y[0])*jacobtot_inv[1]+jacobtot_inv[0]*(0.1530931089239486*b_y[0]-0.2651650429449552*b_y[1]))*hamil[14])*rdz2)/q_; 
  alphaL0[5] = ((((0.4592793267718456*b_y[0]-0.7954951288348656*b_y[1])*jacobtot_inv[5]+(0.4592793267718456*jacobtot_inv[0]-0.7954951288348656*jacobtot_inv[1])*b_y[5]+(0.4592793267718456*b_y[1]-0.2651650429449552*b_y[0])*jacobtot_inv[3]+(0.4592793267718456*jacobtot_inv[1]-0.2651650429449552*jacobtot_inv[0])*b_y[3])*hamil[16]+((0.4592793267718456*b_y[1]-0.2651650429449552*b_y[0])*jacobtot_inv[5]+(0.4592793267718456*jacobtot_inv[1]-0.2651650429449552*jacobtot_inv[0])*b_y[5]+(0.1530931089239486*b_y[0]-0.2651650429449552*b_y[1])*jacobtot_inv[3]+(0.1530931089239486*jacobtot_inv[0]-0.2651650429449552*jacobtot_inv[1])*b_y[3])*hamil[8])*rdz2)/q_; 
  alphaL0[7] = (((0.9682458365518543*BstarXdBmag[2]-1.677050983124842*BstarXdBmag[4])*hamil[32]+hamil[4]*(0.4330127018922193*BstarXdBmag[6]-0.75*BstarXdBmag[7]))*rdvpar2)/m_; 
  alphaL0[9] = ((((0.4592793267718456*b_y[0]-0.7954951288348656*b_y[1])*jacobtot_inv[5]+(0.4592793267718456*jacobtot_inv[0]-0.7954951288348656*jacobtot_inv[1])*b_y[5]+(0.4592793267718456*b_y[1]-0.2651650429449552*b_y[0])*jacobtot_inv[3]+(0.4592793267718456*jacobtot_inv[1]-0.2651650429449552*jacobtot_inv[0])*b_y[3])*hamil[21]+((0.4592793267718456*b_y[1]-0.2651650429449552*b_y[0])*jacobtot_inv[5]+(0.4592793267718456*jacobtot_inv[1]-0.2651650429449552*jacobtot_inv[0])*b_y[5]+(0.1530931089239486*b_y[0]-0.2651650429449552*b_y[1])*jacobtot_inv[3]+(0.1530931089239486*jacobtot_inv[0]-0.2651650429449552*jacobtot_inv[1])*b_y[3])*hamil[14])*rdz2)/q_; 
  alphaL0[16] = ((0.8660254037844386*BstarXdBmag[3]-1.5*BstarXdBmag[5])*hamil[32]*rdvpar2)/m_; 
  alphaL0[18] = ((0.8660254037844387*BstarXdBmag[6]-1.5*BstarXdBmag[7])*hamil[32]*rdvpar2)/m_; 

} 