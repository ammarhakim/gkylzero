#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_surfvpar_3x2v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, 
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
  const double *b_y = &b_i[20];
  const double *b_z = &b_i[40];

  const double *BstarXdBmag = &Bstar_Bmag[0]; 
  const double *BstarYdBmag = &Bstar_Bmag[20]; 
  const double *BstarZdBmag = &Bstar_Bmag[40]; 

  double hamil[112] = {0.}; 
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
  hamil[16] = 2.0*bmag[7]*wmu+2.0*phi[7]*q_; 
  hamil[17] = 2.0*phi[8]*q_; 
  hamil[18] = 2.0*bmag[9]*wmu+2.0*phi[9]*q_; 
  hamil[19] = (0.8432740427115681*m_)/rdvpar2Sq; 
  hamil[21] = 2.0*phi[10]*q_; 
  hamil[26] = (1.154700538379252*bmag[5])/rdmu2; 
  hamil[31] = 2.0*phi[11]*q_; 
  hamil[32] = 2.0*phi[12]*q_; 
  hamil[33] = 2.0*bmag[13]*wmu+2.0*phi[13]*q_; 
  hamil[34] = 2.0*phi[14]*q_; 
  hamil[35] = 2.0*bmag[15]*wmu+2.0*phi[15]*q_; 
  hamil[36] = 2.0*phi[16]*q_; 
  hamil[43] = (1.154700538379251*bmag[7])/rdmu2; 
  hamil[45] = (1.154700538379251*bmag[9])/rdmu2; 
  hamil[56] = 2.0*phi[17]*q_; 
  hamil[57] = 2.0*phi[18]*q_; 
  hamil[58] = 2.0*phi[19]*q_; 
  hamil[70] = (1.154700538379251*bmag[13])/rdmu2; 
  hamil[72] = (1.154700538379251*bmag[15])/rdmu2; 

  double *alphaL3 = &alpha_surf[144];
  alphaL3[0] = (((1.677050983124842*BstarZdBmag[10]-0.9682458365518543*BstarZdBmag[4])*hamil[35]+(0.75*BstarZdBmag[13]-0.4330127018922194*BstarZdBmag[7])*hamil[33]+(1.677050983124842*BstarZdBmag[6]-0.9682458365518543*BstarZdBmag[2])*hamil[18]+(0.75*BstarZdBmag[5]-0.4330127018922193*BstarZdBmag[1])*hamil[7]+(0.75*BstarZdBmag[3]-0.4330127018922193*BstarZdBmag[0])*hamil[3])*rdz2+((0.75*BstarYdBmag[18]-0.4330127018922194*BstarYdBmag[12])*hamil[58]+(0.75*BstarYdBmag[17]-0.4330127018922194*BstarYdBmag[11])*hamil[56]+(0.75*BstarYdBmag[14]-0.4330127018922194*BstarYdBmag[8])*hamil[36]+(0.75*BstarYdBmag[13]-0.4330127018922194*BstarYdBmag[7])*hamil[31]+(0.75*BstarYdBmag[10]-0.4330127018922193*BstarYdBmag[4])*hamil[21]+(0.75*BstarYdBmag[6]-0.4330127018922193*BstarYdBmag[2])*hamil[8]+(0.75*BstarYdBmag[5]-0.4330127018922193*BstarYdBmag[1])*hamil[6]+hamil[2]*(0.75*BstarYdBmag[3]-0.4330127018922193*BstarYdBmag[0]))*rdy2+((0.75*BstarXdBmag[14]-0.4330127018922194*BstarXdBmag[8])*hamil[35]+(1.677050983124842*BstarXdBmag[10]-0.9682458365518543*BstarXdBmag[4])*hamil[33]+(1.677050983124842*BstarXdBmag[5]-0.9682458365518543*BstarXdBmag[1])*hamil[16]+(0.75*BstarXdBmag[6]-0.4330127018922193*BstarXdBmag[2])*hamil[7]+hamil[1]*(0.75*BstarXdBmag[3]-0.4330127018922193*BstarXdBmag[0]))*rdx2)/m_; 
  alphaL3[1] = (((1.5*BstarZdBmag[17]-0.8660254037844386*BstarZdBmag[11]+1.677050983124842*BstarZdBmag[6]-0.9682458365518543*BstarZdBmag[2])*hamil[35]+(0.6708203932499369*BstarZdBmag[5]-0.3872983346207417*BstarZdBmag[1])*hamil[33]+(1.677050983124842*BstarZdBmag[10]-0.9682458365518543*BstarZdBmag[4])*hamil[18]+hamil[7]*(0.6708203932499369*BstarZdBmag[13]-0.3872983346207416*BstarZdBmag[7]+0.75*BstarZdBmag[3]-0.4330127018922193*BstarZdBmag[0])+hamil[3]*(0.75*BstarZdBmag[5]-0.4330127018922193*BstarZdBmag[1]))*rdz2+((0.75*BstarYdBmag[14]-0.4330127018922193*BstarYdBmag[8])*hamil[58]+(0.6708203932499369*BstarYdBmag[10]-0.3872983346207416*BstarYdBmag[4])*hamil[56]+(0.75*BstarYdBmag[18]-0.4330127018922193*BstarYdBmag[12])*hamil[36]+(0.6708203932499369*BstarYdBmag[5]-0.3872983346207417*BstarYdBmag[1])*hamil[31]+(0.6708203932499369*BstarYdBmag[17]-0.3872983346207417*BstarYdBmag[11]+0.75*BstarYdBmag[6]-0.4330127018922193*BstarYdBmag[2])*hamil[21]+0.6708203932499369*hamil[6]*BstarYdBmag[13]+hamil[8]*(0.75*BstarYdBmag[10]-0.4330127018922193*BstarYdBmag[4])+hamil[6]*((-0.3872983346207416*BstarYdBmag[7])+0.75*BstarYdBmag[3]-0.4330127018922193*BstarYdBmag[0])+hamil[2]*(0.75*BstarYdBmag[5]-0.4330127018922193*BstarYdBmag[1]))*rdy2+((0.75*BstarXdBmag[18]-0.4330127018922193*BstarXdBmag[12])*hamil[35]+(1.5*BstarXdBmag[17]-0.8660254037844386*BstarXdBmag[11]+1.677050983124842*BstarXdBmag[6]-0.9682458365518543*BstarXdBmag[2])*hamil[33]+(1.5*BstarXdBmag[13]-0.8660254037844386*BstarXdBmag[7]+1.677050983124842*BstarXdBmag[3]-0.9682458365518543*BstarXdBmag[0])*hamil[16]+hamil[7]*(0.75*BstarXdBmag[10]-0.4330127018922193*BstarXdBmag[4])+hamil[1]*(0.75*BstarXdBmag[5]-0.4330127018922193*BstarXdBmag[1]))*rdx2)/m_; 
  alphaL3[2] = (((1.677050983124842*BstarZdBmag[10]-0.9682458365518543*BstarZdBmag[4])*hamil[58]+(0.75*BstarZdBmag[13]-0.4330127018922193*BstarZdBmag[7])*hamil[56]+(1.677050983124842*BstarZdBmag[6]-0.9682458365518543*BstarZdBmag[2])*hamil[36]+(0.75*BstarZdBmag[5]-0.4330127018922193*BstarZdBmag[1])*hamil[21]+(0.75*BstarZdBmag[3]-0.4330127018922193*BstarZdBmag[0])*hamil[8])*rdz2+((1.677050983124842*BstarYdBmag[10]-0.9682458365518543*BstarYdBmag[4])*hamil[57]+(1.677050983124842*BstarYdBmag[6]-0.9682458365518543*BstarYdBmag[2])*hamil[34]+(1.677050983124842*BstarYdBmag[5]-0.9682458365518543*BstarYdBmag[1])*hamil[32]+(1.677050983124842*BstarYdBmag[3]-0.9682458365518543*BstarYdBmag[0])*hamil[17])*rdy2+((0.75*BstarXdBmag[14]-0.4330127018922193*BstarXdBmag[8])*hamil[58]+(1.677050983124842*BstarXdBmag[10]-0.9682458365518543*BstarXdBmag[4])*hamil[56]+(1.677050983124842*BstarXdBmag[5]-0.9682458365518543*BstarXdBmag[1])*hamil[31]+(0.75*BstarXdBmag[6]-0.4330127018922193*BstarXdBmag[2])*hamil[21]+(0.75*BstarXdBmag[3]-0.4330127018922193*BstarXdBmag[0])*hamil[6])*rdx2)/m_; 
  alphaL3[3] = (((1.5*BstarZdBmag[18]-0.8660254037844386*BstarZdBmag[12]+1.677050983124842*BstarZdBmag[5]-0.9682458365518543*BstarZdBmag[1])*hamil[35]+(0.75*BstarZdBmag[17]-0.4330127018922193*BstarZdBmag[11])*hamil[33]+(1.5*BstarZdBmag[14]-0.8660254037844386*BstarZdBmag[8]+1.677050983124842*BstarZdBmag[3]-0.9682458365518543*BstarZdBmag[0])*hamil[18]+hamil[7]*(0.75*BstarZdBmag[10]-0.4330127018922193*BstarZdBmag[4])+hamil[3]*(0.75*BstarZdBmag[6]-0.4330127018922193*BstarZdBmag[2]))*rdz2+((0.6708203932499369*BstarYdBmag[10]-0.3872983346207416*BstarYdBmag[4])*hamil[58]+(0.75*BstarYdBmag[13]-0.4330127018922193*BstarYdBmag[7])*hamil[56]+(0.6708203932499369*BstarYdBmag[6]-0.3872983346207417*BstarYdBmag[2])*hamil[36]+(0.75*BstarYdBmag[17]-0.4330127018922193*BstarYdBmag[11])*hamil[31]+(0.6708203932499369*BstarYdBmag[18]-0.3872983346207417*BstarYdBmag[12]+0.75*BstarYdBmag[5]-0.4330127018922193*BstarYdBmag[1])*hamil[21]+0.6708203932499369*hamil[8]*BstarYdBmag[14]+0.75*hamil[6]*BstarYdBmag[10]-0.4330127018922193*(BstarYdBmag[0]*hamil[8]+BstarYdBmag[4]*hamil[6])+(0.75*BstarYdBmag[3]-0.3872983346207416*BstarYdBmag[8])*hamil[8]+hamil[2]*(0.75*BstarYdBmag[6]-0.4330127018922193*BstarYdBmag[2]))*rdy2+((0.6708203932499369*BstarXdBmag[6]-0.3872983346207417*BstarXdBmag[2])*hamil[35]+(1.5*BstarXdBmag[18]-0.8660254037844386*BstarXdBmag[12]+1.677050983124842*BstarXdBmag[5]-0.9682458365518543*BstarXdBmag[1])*hamil[33]+(1.677050983124842*BstarXdBmag[10]-0.9682458365518543*BstarXdBmag[4])*hamil[16]+hamil[7]*(0.6708203932499369*BstarXdBmag[14]-0.3872983346207416*BstarXdBmag[8]+0.75*BstarXdBmag[3]-0.4330127018922193*BstarXdBmag[0])+hamil[1]*(0.75*BstarXdBmag[6]-0.4330127018922193*BstarXdBmag[2]))*rdx2)/m_; 
  alphaL3[4] = (((1.677050983124842*BstarZdBmag[10]-0.9682458365518543*BstarZdBmag[4])*hamil[72]+(0.75*BstarZdBmag[13]-0.4330127018922193*BstarZdBmag[7])*hamil[70]+(1.677050983124842*BstarZdBmag[6]-0.9682458365518543*BstarZdBmag[2])*hamil[45]+(0.75*BstarZdBmag[5]-0.4330127018922193*BstarZdBmag[1])*hamil[26]+(0.75*BstarZdBmag[3]-0.4330127018922193*BstarZdBmag[0])*hamil[14])*rdz2+((0.75*BstarXdBmag[14]-0.4330127018922193*BstarXdBmag[8])*hamil[72]+(1.677050983124842*BstarXdBmag[10]-0.9682458365518543*BstarXdBmag[4])*hamil[70]+(1.677050983124842*BstarXdBmag[5]-0.9682458365518543*BstarXdBmag[1])*hamil[43]+(0.75*BstarXdBmag[6]-0.4330127018922193*BstarXdBmag[2])*hamil[26]+(0.75*BstarXdBmag[3]-0.4330127018922193*BstarXdBmag[0])*hamil[12])*rdx2)/m_; 
  alphaL3[5] = (((1.5*BstarZdBmag[17]-0.8660254037844387*BstarZdBmag[11]+1.677050983124842*BstarZdBmag[6]-0.9682458365518543*BstarZdBmag[2])*hamil[58]+(0.6708203932499369*BstarZdBmag[5]-0.3872983346207416*BstarZdBmag[1])*hamil[56]+(1.677050983124842*BstarZdBmag[10]-0.9682458365518543*BstarZdBmag[4])*hamil[36]+(0.6708203932499369*BstarZdBmag[13]-0.3872983346207416*BstarZdBmag[7]+0.75*BstarZdBmag[3]-0.4330127018922193*BstarZdBmag[0])*hamil[21]+(0.75*BstarZdBmag[5]-0.4330127018922193*BstarZdBmag[1])*hamil[8])*rdz2+((1.5*BstarYdBmag[17]-0.8660254037844387*BstarYdBmag[11]+1.677050983124842*BstarYdBmag[6]-0.9682458365518543*BstarYdBmag[2])*hamil[57]+(1.677050983124842*BstarYdBmag[10]-0.9682458365518543*BstarYdBmag[4])*hamil[34]+(1.5*BstarYdBmag[13]-0.8660254037844387*BstarYdBmag[7]+1.677050983124842*BstarYdBmag[3]-0.9682458365518543*BstarYdBmag[0])*hamil[32]+(1.677050983124842*BstarYdBmag[5]-0.9682458365518543*BstarYdBmag[1])*hamil[17])*rdy2+((0.75*BstarXdBmag[18]-0.4330127018922194*BstarXdBmag[12])*hamil[58]+(1.5*BstarXdBmag[17]-0.8660254037844387*BstarXdBmag[11]+1.677050983124842*BstarXdBmag[6]-0.9682458365518543*BstarXdBmag[2])*hamil[56]+(1.5*BstarXdBmag[13]-0.8660254037844387*BstarXdBmag[7]+1.677050983124842*BstarXdBmag[3]-0.9682458365518543*BstarXdBmag[0])*hamil[31]+(0.75*BstarXdBmag[10]-0.4330127018922193*BstarXdBmag[4])*hamil[21]+(0.75*BstarXdBmag[5]-0.4330127018922193*BstarXdBmag[1])*hamil[6])*rdx2)/m_; 
  alphaL3[6] = (((1.5*(BstarZdBmag[14]+BstarZdBmag[13])-0.8660254037844387*(BstarZdBmag[8]+BstarZdBmag[7])+1.677050983124842*BstarZdBmag[3]-0.9682458365518543*BstarZdBmag[0])*hamil[35]+(0.6708203932499369*BstarZdBmag[10]-0.3872983346207417*BstarZdBmag[4])*hamil[33]+(1.5*BstarZdBmag[18]-0.8660254037844387*BstarZdBmag[12]+1.677050983124842*BstarZdBmag[5]-0.9682458365518543*BstarZdBmag[1])*hamil[18]+hamil[7]*(0.6708203932499369*BstarZdBmag[17]-0.3872983346207417*BstarZdBmag[11])+0.75*hamil[3]*BstarZdBmag[10]-0.4330127018922193*(BstarZdBmag[2]*hamil[7]+hamil[3]*BstarZdBmag[4])+0.75*BstarZdBmag[6]*hamil[7])*rdz2+((0.6*BstarYdBmag[17]-0.3464101615137754*BstarYdBmag[11]+0.6708203932499369*BstarYdBmag[6]-0.3872983346207416*BstarYdBmag[2])*hamil[58]+(0.6*BstarYdBmag[18]-0.3464101615137754*BstarYdBmag[12]+0.6708203932499369*BstarYdBmag[5]-0.3872983346207416*BstarYdBmag[1])*hamil[56]+(0.6708203932499369*BstarYdBmag[10]-0.3872983346207417*BstarYdBmag[4])*(hamil[36]+hamil[31])+(0.6708203932499369*(BstarYdBmag[14]+BstarYdBmag[13])-0.3872983346207416*(BstarYdBmag[8]+BstarYdBmag[7])+0.75*BstarYdBmag[3]-0.4330127018922193*BstarYdBmag[0])*hamil[21]+0.6708203932499369*(hamil[8]*BstarYdBmag[18]+hamil[6]*BstarYdBmag[17])-0.3872983346207417*(hamil[8]*BstarYdBmag[12]+hamil[6]*BstarYdBmag[11])+0.75*(hamil[2]*BstarYdBmag[10]+BstarYdBmag[5]*hamil[8]+BstarYdBmag[6]*hamil[6])-0.4330127018922193*(BstarYdBmag[1]*hamil[8]+BstarYdBmag[2]*hamil[6]+hamil[2]*BstarYdBmag[4]))*rdy2+((0.6708203932499369*BstarXdBmag[10]-0.3872983346207417*BstarXdBmag[4])*hamil[35]+(1.5*(BstarXdBmag[14]+BstarXdBmag[13])-0.8660254037844387*(BstarXdBmag[8]+BstarXdBmag[7])+1.677050983124842*BstarXdBmag[3]-0.9682458365518543*BstarXdBmag[0])*hamil[33]+0.6708203932499369*hamil[7]*BstarXdBmag[18]+hamil[16]*(1.5*BstarXdBmag[17]-0.8660254037844387*BstarXdBmag[11]+1.677050983124842*BstarXdBmag[6]-0.9682458365518543*BstarXdBmag[2])-0.3872983346207417*hamil[7]*BstarXdBmag[12]+0.75*hamil[1]*BstarXdBmag[10]-0.4330127018922193*(BstarXdBmag[1]*hamil[7]+hamil[1]*BstarXdBmag[4])+0.75*BstarXdBmag[5]*hamil[7])*rdx2)/m_; 
  alphaL3[7] = (((1.5*BstarZdBmag[18]-0.8660254037844387*BstarZdBmag[12]+1.677050983124842*BstarZdBmag[5]-0.9682458365518543*BstarZdBmag[1])*hamil[58]+(0.75*BstarZdBmag[17]-0.4330127018922194*BstarZdBmag[11])*hamil[56]+(1.5*BstarZdBmag[14]-0.8660254037844387*BstarZdBmag[8]+1.677050983124842*BstarZdBmag[3]-0.9682458365518543*BstarZdBmag[0])*hamil[36]+(0.75*BstarZdBmag[10]-0.4330127018922193*BstarZdBmag[4])*hamil[21]+(0.75*BstarZdBmag[6]-0.4330127018922193*BstarZdBmag[2])*hamil[8])*rdz2+((1.5*BstarYdBmag[18]-0.8660254037844387*BstarYdBmag[12]+1.677050983124842*BstarYdBmag[5]-0.9682458365518543*BstarYdBmag[1])*hamil[57]+(1.5*BstarYdBmag[14]-0.8660254037844387*BstarYdBmag[8]+1.677050983124842*BstarYdBmag[3]-0.9682458365518543*BstarYdBmag[0])*hamil[34]+(1.677050983124842*BstarYdBmag[10]-0.9682458365518543*BstarYdBmag[4])*hamil[32]+(1.677050983124842*BstarYdBmag[6]-0.9682458365518543*BstarYdBmag[2])*hamil[17])*rdy2+((0.6708203932499369*BstarXdBmag[6]-0.3872983346207416*BstarXdBmag[2])*hamil[58]+(1.5*BstarXdBmag[18]-0.8660254037844387*BstarXdBmag[12]+1.677050983124842*BstarXdBmag[5]-0.9682458365518543*BstarXdBmag[1])*hamil[56]+(1.677050983124842*BstarXdBmag[10]-0.9682458365518543*BstarXdBmag[4])*hamil[31]+(0.6708203932499369*BstarXdBmag[14]-0.3872983346207416*BstarXdBmag[8]+0.75*BstarXdBmag[3]-0.4330127018922193*BstarXdBmag[0])*hamil[21]+(0.75*BstarXdBmag[6]-0.4330127018922193*BstarXdBmag[2])*hamil[6])*rdx2)/m_; 
  alphaL3[8] = (((1.5*BstarZdBmag[17]-0.8660254037844387*BstarZdBmag[11]+1.677050983124842*BstarZdBmag[6]-0.9682458365518543*BstarZdBmag[2])*hamil[72]+(0.6708203932499369*BstarZdBmag[5]-0.3872983346207416*BstarZdBmag[1])*hamil[70]+(1.677050983124842*BstarZdBmag[10]-0.9682458365518543*BstarZdBmag[4])*hamil[45]+(0.6708203932499369*BstarZdBmag[13]-0.3872983346207416*BstarZdBmag[7]+0.75*BstarZdBmag[3]-0.4330127018922193*BstarZdBmag[0])*hamil[26]+(0.75*BstarZdBmag[5]-0.4330127018922193*BstarZdBmag[1])*hamil[14])*rdz2+((0.75*BstarXdBmag[18]-0.4330127018922194*BstarXdBmag[12])*hamil[72]+(1.5*BstarXdBmag[17]-0.8660254037844387*BstarXdBmag[11]+1.677050983124842*BstarXdBmag[6]-0.9682458365518543*BstarXdBmag[2])*hamil[70]+(1.5*BstarXdBmag[13]-0.8660254037844387*BstarXdBmag[7]+1.677050983124842*BstarXdBmag[3]-0.9682458365518543*BstarXdBmag[0])*hamil[43]+(0.75*BstarXdBmag[10]-0.4330127018922193*BstarXdBmag[4])*hamil[26]+(0.75*BstarXdBmag[5]-0.4330127018922193*BstarXdBmag[1])*hamil[12])*rdx2)/m_; 
  alphaL3[10] = (((1.5*BstarZdBmag[18]-0.8660254037844387*BstarZdBmag[12]+1.677050983124842*BstarZdBmag[5]-0.9682458365518543*BstarZdBmag[1])*hamil[72]+(0.75*BstarZdBmag[17]-0.4330127018922194*BstarZdBmag[11])*hamil[70]+(1.5*BstarZdBmag[14]-0.8660254037844387*BstarZdBmag[8]+1.677050983124842*BstarZdBmag[3]-0.9682458365518543*BstarZdBmag[0])*hamil[45]+(0.75*BstarZdBmag[10]-0.4330127018922193*BstarZdBmag[4])*hamil[26]+(0.75*BstarZdBmag[6]-0.4330127018922193*BstarZdBmag[2])*hamil[14])*rdz2+((0.6708203932499369*BstarXdBmag[6]-0.3872983346207416*BstarXdBmag[2])*hamil[72]+(1.5*BstarXdBmag[18]-0.8660254037844387*BstarXdBmag[12]+1.677050983124842*BstarXdBmag[5]-0.9682458365518543*BstarXdBmag[1])*hamil[70]+(1.677050983124842*BstarXdBmag[10]-0.9682458365518543*BstarXdBmag[4])*hamil[43]+(0.6708203932499369*BstarXdBmag[14]-0.3872983346207416*BstarXdBmag[8]+0.75*BstarXdBmag[3]-0.4330127018922193*BstarXdBmag[0])*hamil[26]+(0.75*BstarXdBmag[6]-0.4330127018922193*BstarXdBmag[2])*hamil[12])*rdx2)/m_; 
  alphaL3[11] = (((1.5*BstarZdBmag[10]-0.8660254037844387*BstarZdBmag[4])*hamil[35]+(0.479157423749955*BstarZdBmag[13]-0.276641667586244*BstarZdBmag[7]+0.75*BstarZdBmag[3]-0.4330127018922194*BstarZdBmag[0])*hamil[33]+(1.677050983124842*BstarZdBmag[17]-0.9682458365518543*BstarZdBmag[11])*hamil[18]+0.75*hamil[3]*BstarZdBmag[13]+(0.6708203932499369*BstarZdBmag[5]-0.3872983346207416*BstarZdBmag[1])*hamil[7]-0.4330127018922193*hamil[3]*BstarZdBmag[7])*rdz2+((0.6708203932499369*BstarYdBmag[18]-0.3872983346207417*BstarYdBmag[12])*hamil[58]+(0.479157423749955*BstarYdBmag[17]-0.276641667586244*BstarYdBmag[11]+0.75*BstarYdBmag[6]-0.4330127018922193*BstarYdBmag[2])*hamil[56]+(0.479157423749955*BstarYdBmag[13]-0.276641667586244*BstarYdBmag[7]+0.75*BstarYdBmag[3]-0.4330127018922194*BstarYdBmag[0])*hamil[31]+(0.6708203932499369*BstarYdBmag[10]-0.3872983346207416*BstarYdBmag[4])*hamil[21]+0.75*(hamil[8]*BstarYdBmag[17]+hamil[2]*BstarYdBmag[13])-0.4330127018922194*hamil[8]*BstarYdBmag[11]-0.4330127018922193*hamil[2]*BstarYdBmag[7]+(0.6708203932499369*BstarYdBmag[5]-0.3872983346207416*BstarYdBmag[1])*hamil[6])*rdy2+((1.5*BstarXdBmag[10]-0.8660254037844387*BstarXdBmag[4])*hamil[33]+0.75*hamil[7]*BstarXdBmag[17]+(1.5*BstarXdBmag[5]-0.8660254037844386*BstarXdBmag[1])*hamil[16]+0.75*hamil[1]*BstarXdBmag[13]-0.4330127018922194*hamil[7]*BstarXdBmag[11]-0.4330127018922193*hamil[1]*BstarXdBmag[7])*rdx2)/m_; 
  alphaL3[12] = (((0.75*BstarZdBmag[5]-0.4330127018922193*BstarZdBmag[1])*hamil[57]+(0.75*BstarZdBmag[3]-0.4330127018922194*BstarZdBmag[0])*hamil[34])*rdz2+((0.75*BstarXdBmag[6]-0.4330127018922193*BstarXdBmag[2])*hamil[57]+(0.75*BstarXdBmag[3]-0.4330127018922194*BstarXdBmag[0])*hamil[32])*rdx2)/m_; 
  alphaL3[13] = (((1.5*BstarZdBmag[10]-0.8660254037844387*BstarZdBmag[4])*hamil[35]+(1.5*BstarZdBmag[6]-0.8660254037844386*BstarZdBmag[2])*hamil[18]+0.75*(hamil[7]*BstarZdBmag[18]+hamil[3]*BstarZdBmag[14])-0.4330127018922194*hamil[7]*BstarZdBmag[12]-0.4330127018922193*hamil[3]*BstarZdBmag[8])*rdz2+((0.479157423749955*BstarYdBmag[18]-0.276641667586244*BstarYdBmag[12]+0.75*BstarYdBmag[5]-0.4330127018922193*BstarYdBmag[1])*hamil[58]+(0.6708203932499369*BstarYdBmag[17]-0.3872983346207417*BstarYdBmag[11])*hamil[56]+(0.479157423749955*BstarYdBmag[14]-0.276641667586244*BstarYdBmag[8]+0.75*BstarYdBmag[3]-0.4330127018922194*BstarYdBmag[0])*hamil[36]+(0.6708203932499369*BstarYdBmag[10]-0.3872983346207416*BstarYdBmag[4])*hamil[21]+0.75*(hamil[6]*BstarYdBmag[18]+hamil[2]*BstarYdBmag[14])-0.4330127018922194*hamil[6]*BstarYdBmag[12]+(0.6708203932499369*BstarYdBmag[6]-0.3872983346207416*BstarYdBmag[2])*hamil[8]-0.4330127018922193*hamil[2]*BstarYdBmag[8])*rdy2+((0.479157423749955*BstarXdBmag[14]-0.276641667586244*BstarXdBmag[8]+0.75*BstarXdBmag[3]-0.4330127018922194*BstarXdBmag[0])*hamil[35]+(1.5*BstarXdBmag[10]-0.8660254037844387*BstarXdBmag[4])*hamil[33]+hamil[16]*(1.677050983124842*BstarXdBmag[18]-0.9682458365518543*BstarXdBmag[12])+hamil[1]*(0.75*BstarXdBmag[14]-0.4330127018922193*BstarXdBmag[8])+(0.6708203932499369*BstarXdBmag[6]-0.3872983346207416*BstarXdBmag[2])*hamil[7])*rdx2)/m_; 
  alphaL3[15] = (((1.5*(BstarZdBmag[14]+BstarZdBmag[13])-0.8660254037844386*(BstarZdBmag[8]+BstarZdBmag[7])+1.677050983124842*BstarZdBmag[3]-0.9682458365518543*BstarZdBmag[0])*hamil[58]+(0.6708203932499369*BstarZdBmag[10]-0.3872983346207416*BstarZdBmag[4])*hamil[56]+(1.5*BstarZdBmag[18]-0.8660254037844386*BstarZdBmag[12]+1.677050983124842*BstarZdBmag[5]-0.9682458365518543*BstarZdBmag[1])*hamil[36]+(0.6708203932499369*BstarZdBmag[17]-0.3872983346207417*BstarZdBmag[11]+0.75*BstarZdBmag[6]-0.4330127018922193*BstarZdBmag[2])*hamil[21]+hamil[8]*(0.75*BstarZdBmag[10]-0.4330127018922193*BstarZdBmag[4]))*rdz2+((1.5*(BstarYdBmag[14]+BstarYdBmag[13])-0.8660254037844386*(BstarYdBmag[8]+BstarYdBmag[7])+1.677050983124842*BstarYdBmag[3]-0.9682458365518543*BstarYdBmag[0])*hamil[57]+(1.5*BstarYdBmag[18]-0.8660254037844386*BstarYdBmag[12]+1.677050983124842*BstarYdBmag[5]-0.9682458365518543*BstarYdBmag[1])*hamil[34]+(1.5*BstarYdBmag[17]-0.8660254037844386*BstarYdBmag[11]+1.677050983124842*BstarYdBmag[6]-0.9682458365518543*BstarYdBmag[2])*hamil[32]+(1.677050983124842*BstarYdBmag[10]-0.9682458365518543*BstarYdBmag[4])*hamil[17])*rdy2+((0.6708203932499369*BstarXdBmag[10]-0.3872983346207416*BstarXdBmag[4])*hamil[58]+(1.5*(BstarXdBmag[14]+BstarXdBmag[13])-0.8660254037844386*(BstarXdBmag[8]+BstarXdBmag[7])+1.677050983124842*BstarXdBmag[3]-0.9682458365518543*BstarXdBmag[0])*hamil[56]+(1.5*BstarXdBmag[17]-0.8660254037844386*BstarXdBmag[11]+1.677050983124842*BstarXdBmag[6]-0.9682458365518543*BstarXdBmag[2])*hamil[31]+(0.6708203932499369*BstarXdBmag[18]-0.3872983346207417*BstarXdBmag[12]+0.75*BstarXdBmag[5]-0.4330127018922193*BstarXdBmag[1])*hamil[21]+hamil[6]*(0.75*BstarXdBmag[10]-0.4330127018922193*BstarXdBmag[4]))*rdx2)/m_; 
  alphaL3[17] = (((1.5*(BstarZdBmag[14]+BstarZdBmag[13])-0.8660254037844386*(BstarZdBmag[8]+BstarZdBmag[7])+1.677050983124842*BstarZdBmag[3]-0.9682458365518543*BstarZdBmag[0])*hamil[72]+(0.6708203932499369*BstarZdBmag[10]-0.3872983346207416*BstarZdBmag[4])*hamil[70]+(1.5*BstarZdBmag[18]-0.8660254037844386*BstarZdBmag[12]+1.677050983124842*BstarZdBmag[5]-0.9682458365518543*BstarZdBmag[1])*hamil[45]+(0.6708203932499369*BstarZdBmag[17]-0.3872983346207417*BstarZdBmag[11]+0.75*BstarZdBmag[6]-0.4330127018922193*BstarZdBmag[2])*hamil[26]+(0.75*BstarZdBmag[10]-0.4330127018922193*BstarZdBmag[4])*hamil[14])*rdz2+((0.6708203932499369*BstarXdBmag[10]-0.3872983346207416*BstarXdBmag[4])*hamil[72]+(1.5*(BstarXdBmag[14]+BstarXdBmag[13])-0.8660254037844386*(BstarXdBmag[8]+BstarXdBmag[7])+1.677050983124842*BstarXdBmag[3]-0.9682458365518543*BstarXdBmag[0])*hamil[70]+(1.5*BstarXdBmag[17]-0.8660254037844386*BstarXdBmag[11]+1.677050983124842*BstarXdBmag[6]-0.9682458365518543*BstarXdBmag[2])*hamil[43]+(0.6708203932499369*BstarXdBmag[18]-0.3872983346207417*BstarXdBmag[12]+0.75*BstarXdBmag[5]-0.4330127018922193*BstarXdBmag[1])*hamil[26]+(0.75*BstarXdBmag[10]-0.4330127018922193*BstarXdBmag[4])*hamil[12])*rdx2)/m_; 
  alphaL3[19] = (((1.5*BstarZdBmag[10]-0.8660254037844387*BstarZdBmag[4])*hamil[58]+(0.479157423749955*BstarZdBmag[13]-0.276641667586244*BstarZdBmag[7]+0.75*BstarZdBmag[3]-0.4330127018922194*BstarZdBmag[0])*hamil[56]+(1.677050983124842*BstarZdBmag[17]-0.9682458365518543*BstarZdBmag[11])*hamil[36]+(0.6708203932499369*BstarZdBmag[5]-0.3872983346207417*BstarZdBmag[1])*hamil[21]+hamil[8]*(0.75*BstarZdBmag[13]-0.4330127018922194*BstarZdBmag[7]))*rdz2+((1.5*BstarYdBmag[10]-0.8660254037844387*BstarYdBmag[4])*hamil[57]+(1.677050983124842*BstarYdBmag[17]-0.9682458365518543*BstarYdBmag[11])*hamil[34]+(1.5*BstarYdBmag[5]-0.8660254037844386*BstarYdBmag[1])*hamil[32]+(1.677050983124842*BstarYdBmag[13]-0.9682458365518543*BstarYdBmag[7])*hamil[17])*rdy2+((1.5*BstarXdBmag[10]-0.8660254037844387*BstarXdBmag[4])*hamil[56]+(1.5*BstarXdBmag[5]-0.8660254037844386*BstarXdBmag[1])*hamil[31]+(0.75*BstarXdBmag[17]-0.4330127018922193*BstarXdBmag[11])*hamil[21]+hamil[6]*(0.75*BstarXdBmag[13]-0.4330127018922194*BstarXdBmag[7]))*rdx2)/m_; 
  alphaL3[20] = (((0.6708203932499369*BstarZdBmag[13]-0.3872983346207417*BstarZdBmag[7]+0.75*BstarZdBmag[3]-0.4330127018922194*BstarZdBmag[0])*hamil[57]+(0.75*BstarZdBmag[5]-0.4330127018922193*BstarZdBmag[1])*hamil[34])*rdz2+((0.75*BstarXdBmag[10]-0.4330127018922194*BstarXdBmag[4])*hamil[57]+(0.75*BstarXdBmag[5]-0.4330127018922193*BstarXdBmag[1])*hamil[32])*rdx2)/m_; 
  alphaL3[21] = (((1.341640786499874*BstarZdBmag[18]-0.7745966692414834*BstarZdBmag[12]+1.5*BstarZdBmag[5]-0.8660254037844386*BstarZdBmag[1])*hamil[35]+(0.479157423749955*BstarZdBmag[17]-0.276641667586244*BstarZdBmag[11]+0.75*BstarZdBmag[6]-0.4330127018922193*BstarZdBmag[2])*hamil[33]+(1.677050983124842*BstarZdBmag[13]-0.9682458365518543*BstarZdBmag[7])*hamil[18]+hamil[3]*(0.75*BstarZdBmag[17]-0.4330127018922193*BstarZdBmag[11])+hamil[7]*(0.6708203932499369*BstarZdBmag[10]-0.3872983346207417*BstarZdBmag[4]))*rdz2+((0.5999999999999999*BstarYdBmag[10]-0.3464101615137754*BstarYdBmag[4])*hamil[58]+(0.6708203932499369*BstarYdBmag[14]+0.479157423749955*BstarYdBmag[13]-0.3872983346207417*BstarYdBmag[8]-0.276641667586244*BstarYdBmag[7]+0.75*BstarYdBmag[3]-0.4330127018922194*BstarYdBmag[0])*hamil[56]+(0.6708203932499369*BstarYdBmag[17]-0.3872983346207417*BstarYdBmag[11])*hamil[36]+(0.479157423749955*BstarYdBmag[17]-0.276641667586244*BstarYdBmag[11]+0.75*BstarYdBmag[6]-0.4330127018922193*BstarYdBmag[2])*hamil[31]+(0.5999999999999999*BstarYdBmag[18]-0.3464101615137755*BstarYdBmag[12]+0.6708203932499369*BstarYdBmag[5]-0.3872983346207417*BstarYdBmag[1])*hamil[21]+0.75*(hamil[2]*BstarYdBmag[17]+hamil[8]*BstarYdBmag[13])-0.4330127018922193*hamil[2]*BstarYdBmag[11]+0.6708203932499369*hamil[6]*BstarYdBmag[10]-0.4330127018922194*BstarYdBmag[7]*hamil[8]-0.3872983346207417*BstarYdBmag[4]*hamil[6])*rdy2+((0.6708203932499369*BstarXdBmag[17]-0.3872983346207417*BstarXdBmag[11])*hamil[35]+(1.341640786499874*BstarXdBmag[18]-0.7745966692414834*BstarXdBmag[12]+1.5*BstarXdBmag[5]-0.8660254037844386*BstarXdBmag[1])*hamil[33]+0.75*hamil[1]*BstarXdBmag[17]+(1.5*BstarXdBmag[10]-0.8660254037844387*BstarXdBmag[4])*hamil[16]+0.75*hamil[7]*BstarXdBmag[13]-0.4330127018922193*hamil[1]*BstarXdBmag[11]-0.4330127018922194*BstarXdBmag[7]*hamil[7])*rdx2)/m_; 
  alphaL3[22] = (((0.75*BstarZdBmag[10]-0.4330127018922194*BstarZdBmag[4])*hamil[57]+(0.75*BstarZdBmag[6]-0.4330127018922193*BstarZdBmag[2])*hamil[34])*rdz2+((0.6708203932499369*BstarXdBmag[14]-0.3872983346207417*BstarXdBmag[8]+0.75*BstarXdBmag[3]-0.4330127018922194*BstarXdBmag[0])*hamil[57]+(0.75*BstarXdBmag[6]-0.4330127018922193*BstarXdBmag[2])*hamil[32])*rdx2)/m_; 
  alphaL3[23] = (((1.341640786499874*BstarZdBmag[17]-0.7745966692414834*BstarZdBmag[11]+1.5*BstarZdBmag[6]-0.8660254037844386*BstarZdBmag[2])*hamil[35]+(0.6708203932499369*BstarZdBmag[18]-0.3872983346207417*BstarZdBmag[12])*hamil[33]+(1.5*BstarZdBmag[10]-0.8660254037844387*BstarZdBmag[4])*hamil[18]+0.75*(hamil[3]*BstarZdBmag[18]+hamil[7]*BstarZdBmag[14])-0.4330127018922193*hamil[3]*BstarZdBmag[12]-0.4330127018922194*hamil[7]*BstarZdBmag[8])*rdz2+((0.479157423749955*BstarYdBmag[14]+0.6708203932499369*BstarYdBmag[13]-0.276641667586244*BstarYdBmag[8]-0.3872983346207417*BstarYdBmag[7]+0.75*BstarYdBmag[3]-0.4330127018922194*BstarYdBmag[0])*hamil[58]+(0.5999999999999999*BstarYdBmag[10]-0.3464101615137754*BstarYdBmag[4])*hamil[56]+(0.479157423749955*BstarYdBmag[18]-0.276641667586244*BstarYdBmag[12]+0.75*BstarYdBmag[5]-0.4330127018922193*BstarYdBmag[1])*hamil[36]+(0.6708203932499369*BstarYdBmag[18]-0.3872983346207417*BstarYdBmag[12])*hamil[31]+(0.5999999999999999*BstarYdBmag[17]-0.3464101615137755*BstarYdBmag[11]+0.6708203932499369*BstarYdBmag[6]-0.3872983346207417*BstarYdBmag[2])*hamil[21]+0.75*(hamil[2]*BstarYdBmag[18]+hamil[6]*BstarYdBmag[14])-0.4330127018922193*hamil[2]*BstarYdBmag[12]+hamil[8]*(0.6708203932499369*BstarYdBmag[10]-0.3872983346207417*BstarYdBmag[4])-0.4330127018922194*hamil[6]*BstarYdBmag[8])*rdy2+((0.479157423749955*BstarXdBmag[18]-0.276641667586244*BstarXdBmag[12]+0.75*BstarXdBmag[5]-0.4330127018922193*BstarXdBmag[1])*hamil[35]+(1.341640786499874*BstarXdBmag[17]-0.7745966692414834*BstarXdBmag[11]+1.5*BstarXdBmag[6]-0.8660254037844386*BstarXdBmag[2])*hamil[33]+0.75*hamil[1]*BstarXdBmag[18]+(1.677050983124842*BstarXdBmag[14]-0.9682458365518543*BstarXdBmag[8])*hamil[16]-0.4330127018922193*hamil[1]*BstarXdBmag[12]+hamil[7]*(0.6708203932499369*BstarXdBmag[10]-0.3872983346207417*BstarXdBmag[4]))*rdx2)/m_; 
  alphaL3[24] = (((1.5*BstarZdBmag[10]-0.8660254037844387*BstarZdBmag[4])*hamil[58]+(1.5*BstarZdBmag[6]-0.8660254037844386*BstarZdBmag[2])*hamil[36]+(0.75*BstarZdBmag[18]-0.4330127018922193*BstarZdBmag[12])*hamil[21]+hamil[8]*(0.75*BstarZdBmag[14]-0.4330127018922194*BstarZdBmag[8]))*rdz2+((1.5*BstarYdBmag[10]-0.8660254037844387*BstarYdBmag[4])*hamil[57]+(1.5*BstarYdBmag[6]-0.8660254037844386*BstarYdBmag[2])*hamil[34]+(1.677050983124842*BstarYdBmag[18]-0.9682458365518543*BstarYdBmag[12])*hamil[32]+(1.677050983124842*BstarYdBmag[14]-0.9682458365518543*BstarYdBmag[8])*hamil[17])*rdy2+((0.479157423749955*BstarXdBmag[14]-0.276641667586244*BstarXdBmag[8]+0.75*BstarXdBmag[3]-0.4330127018922194*BstarXdBmag[0])*hamil[58]+(1.5*BstarXdBmag[10]-0.8660254037844387*BstarXdBmag[4])*hamil[56]+(1.677050983124842*BstarXdBmag[18]-0.9682458365518543*BstarXdBmag[12])*hamil[31]+(0.6708203932499369*BstarXdBmag[6]-0.3872983346207417*BstarXdBmag[2])*hamil[21]+hamil[6]*(0.75*BstarXdBmag[14]-0.4330127018922194*BstarXdBmag[8]))*rdx2)/m_; 
  alphaL3[25] = (((1.5*BstarZdBmag[10]-0.8660254037844387*BstarZdBmag[4])*hamil[72]+(0.479157423749955*BstarZdBmag[13]-0.276641667586244*BstarZdBmag[7]+0.75*BstarZdBmag[3]-0.4330127018922194*BstarZdBmag[0])*hamil[70]+(1.677050983124842*BstarZdBmag[17]-0.9682458365518543*BstarZdBmag[11])*hamil[45]+(0.6708203932499369*BstarZdBmag[5]-0.3872983346207417*BstarZdBmag[1])*hamil[26]+(0.75*BstarZdBmag[13]-0.4330127018922194*BstarZdBmag[7])*hamil[14])*rdz2+((1.5*BstarXdBmag[10]-0.8660254037844387*BstarXdBmag[4])*hamil[70]+(1.5*BstarXdBmag[5]-0.8660254037844386*BstarXdBmag[1])*hamil[43]+(0.75*BstarXdBmag[17]-0.4330127018922193*BstarXdBmag[11])*hamil[26]+hamil[12]*(0.75*BstarXdBmag[13]-0.4330127018922194*BstarXdBmag[7]))*rdx2)/m_; 
  alphaL3[27] = (((1.5*BstarZdBmag[10]-0.8660254037844387*BstarZdBmag[4])*hamil[72]+(1.5*BstarZdBmag[6]-0.8660254037844386*BstarZdBmag[2])*hamil[45]+(0.75*BstarZdBmag[18]-0.4330127018922193*BstarZdBmag[12])*hamil[26]+(0.75*BstarZdBmag[14]-0.4330127018922194*BstarZdBmag[8])*hamil[14])*rdz2+((0.479157423749955*BstarXdBmag[14]-0.276641667586244*BstarXdBmag[8]+0.75*BstarXdBmag[3]-0.4330127018922194*BstarXdBmag[0])*hamil[72]+(1.5*BstarXdBmag[10]-0.8660254037844387*BstarXdBmag[4])*hamil[70]+(1.677050983124842*BstarXdBmag[18]-0.9682458365518543*BstarXdBmag[12])*hamil[43]+(0.6708203932499369*BstarXdBmag[6]-0.3872983346207417*BstarXdBmag[2])*hamil[26]+hamil[12]*(0.75*BstarXdBmag[14]-0.4330127018922194*BstarXdBmag[8]))*rdx2)/m_; 
  alphaL3[32] = (((1.341640786499874*BstarZdBmag[18]-0.7745966692414834*BstarZdBmag[12]+1.5*BstarZdBmag[5]-0.8660254037844386*BstarZdBmag[1])*hamil[58]+(0.479157423749955*BstarZdBmag[17]-0.276641667586244*BstarZdBmag[11]+0.75*BstarZdBmag[6]-0.4330127018922193*BstarZdBmag[2])*hamil[56]+(1.677050983124842*BstarZdBmag[13]-0.9682458365518543*BstarZdBmag[7])*hamil[36]+(0.6708203932499369*BstarZdBmag[10]-0.3872983346207416*BstarZdBmag[4])*hamil[21]+hamil[8]*(0.75*BstarZdBmag[17]-0.4330127018922194*BstarZdBmag[11]))*rdz2+((1.341640786499874*BstarYdBmag[18]-0.7745966692414834*BstarYdBmag[12]+1.5*BstarYdBmag[5]-0.8660254037844386*BstarYdBmag[1])*hamil[57]+(1.677050983124842*BstarYdBmag[13]-0.9682458365518543*BstarYdBmag[7])*hamil[34]+(1.5*BstarYdBmag[10]-0.8660254037844387*BstarYdBmag[4])*hamil[32]+(1.677050983124842*BstarYdBmag[17]-0.9682458365518543*BstarYdBmag[11])*hamil[17])*rdy2+((0.6708203932499369*BstarXdBmag[17]-0.3872983346207417*BstarXdBmag[11])*hamil[58]+(1.341640786499874*BstarXdBmag[18]-0.7745966692414834*BstarXdBmag[12]+1.5*BstarXdBmag[5]-0.8660254037844386*BstarXdBmag[1])*hamil[56]+(1.5*BstarXdBmag[10]-0.8660254037844387*BstarXdBmag[4])*hamil[31]+(0.75*BstarXdBmag[13]-0.4330127018922193*BstarXdBmag[7])*hamil[21]+hamil[6]*(0.75*BstarXdBmag[17]-0.4330127018922194*BstarXdBmag[11]))*rdx2)/m_; 
  alphaL3[33] = (((0.6708203932499369*BstarZdBmag[17]-0.3872983346207417*BstarZdBmag[11]+0.75*BstarZdBmag[6]-0.4330127018922193*BstarZdBmag[2])*hamil[57]+(0.75*BstarZdBmag[10]-0.4330127018922194*BstarZdBmag[4])*hamil[34])*rdz2+((0.6708203932499369*BstarXdBmag[18]-0.3872983346207417*BstarXdBmag[12]+0.75*BstarXdBmag[5]-0.4330127018922193*BstarXdBmag[1])*hamil[57]+(0.75*BstarXdBmag[10]-0.4330127018922194*BstarXdBmag[4])*hamil[32])*rdx2)/m_; 
  alphaL3[34] = (((1.341640786499874*BstarZdBmag[17]-0.7745966692414834*BstarZdBmag[11]+1.5*BstarZdBmag[6]-0.8660254037844386*BstarZdBmag[2])*hamil[58]+(0.6708203932499369*BstarZdBmag[18]-0.3872983346207417*BstarZdBmag[12])*hamil[56]+(1.5*BstarZdBmag[10]-0.8660254037844387*BstarZdBmag[4])*hamil[36]+(0.75*BstarZdBmag[14]-0.4330127018922193*BstarZdBmag[8])*hamil[21]+hamil[8]*(0.75*BstarZdBmag[18]-0.4330127018922194*BstarZdBmag[12]))*rdz2+((1.341640786499874*BstarYdBmag[17]-0.7745966692414834*BstarYdBmag[11]+1.5*BstarYdBmag[6]-0.8660254037844386*BstarYdBmag[2])*hamil[57]+(1.5*BstarYdBmag[10]-0.8660254037844387*BstarYdBmag[4])*hamil[34]+(1.677050983124842*BstarYdBmag[14]-0.9682458365518543*BstarYdBmag[8])*hamil[32]+hamil[17]*(1.677050983124842*BstarYdBmag[18]-0.9682458365518543*BstarYdBmag[12]))*rdy2+((0.479157423749955*BstarXdBmag[18]-0.276641667586244*BstarXdBmag[12]+0.75*BstarXdBmag[5]-0.4330127018922193*BstarXdBmag[1])*hamil[58]+(1.341640786499874*BstarXdBmag[17]-0.7745966692414834*BstarXdBmag[11]+1.5*BstarXdBmag[6]-0.8660254037844386*BstarXdBmag[2])*hamil[56]+(1.677050983124842*BstarXdBmag[14]-0.9682458365518543*BstarXdBmag[8])*hamil[31]+(0.6708203932499369*BstarXdBmag[10]-0.3872983346207416*BstarXdBmag[4])*hamil[21]+hamil[6]*(0.75*BstarXdBmag[18]-0.4330127018922194*BstarXdBmag[12]))*rdx2)/m_; 
  alphaL3[37] = (((1.341640786499874*BstarZdBmag[18]-0.7745966692414834*BstarZdBmag[12]+1.5*BstarZdBmag[5]-0.8660254037844386*BstarZdBmag[1])*hamil[72]+(0.479157423749955*BstarZdBmag[17]-0.276641667586244*BstarZdBmag[11]+0.75*BstarZdBmag[6]-0.4330127018922193*BstarZdBmag[2])*hamil[70]+(1.677050983124842*BstarZdBmag[13]-0.9682458365518543*BstarZdBmag[7])*hamil[45]+(0.6708203932499369*BstarZdBmag[10]-0.3872983346207416*BstarZdBmag[4])*hamil[26]+hamil[14]*(0.75*BstarZdBmag[17]-0.4330127018922194*BstarZdBmag[11]))*rdz2+((0.6708203932499369*BstarXdBmag[17]-0.3872983346207417*BstarXdBmag[11])*hamil[72]+(1.341640786499874*BstarXdBmag[18]-0.7745966692414834*BstarXdBmag[12]+1.5*BstarXdBmag[5]-0.8660254037844386*BstarXdBmag[1])*hamil[70]+(1.5*BstarXdBmag[10]-0.8660254037844387*BstarXdBmag[4])*hamil[43]+(0.75*BstarXdBmag[13]-0.4330127018922193*BstarXdBmag[7])*hamil[26]+hamil[12]*(0.75*BstarXdBmag[17]-0.4330127018922194*BstarXdBmag[11]))*rdx2)/m_; 
  alphaL3[39] = (((1.341640786499874*BstarZdBmag[17]-0.7745966692414834*BstarZdBmag[11]+1.5*BstarZdBmag[6]-0.8660254037844386*BstarZdBmag[2])*hamil[72]+(0.6708203932499369*BstarZdBmag[18]-0.3872983346207417*BstarZdBmag[12])*hamil[70]+(1.5*BstarZdBmag[10]-0.8660254037844387*BstarZdBmag[4])*hamil[45]+(0.75*BstarZdBmag[14]-0.4330127018922193*BstarZdBmag[8])*hamil[26]+hamil[14]*(0.75*BstarZdBmag[18]-0.4330127018922194*BstarZdBmag[12]))*rdz2+((0.479157423749955*BstarXdBmag[18]-0.276641667586244*BstarXdBmag[12]+0.75*BstarXdBmag[5]-0.4330127018922193*BstarXdBmag[1])*hamil[72]+(1.341640786499874*BstarXdBmag[17]-0.7745966692414834*BstarXdBmag[11]+1.5*BstarXdBmag[6]-0.8660254037844386*BstarXdBmag[2])*hamil[70]+(1.677050983124842*BstarXdBmag[14]-0.9682458365518543*BstarXdBmag[8])*hamil[43]+(0.6708203932499369*BstarXdBmag[10]-0.3872983346207416*BstarXdBmag[4])*hamil[26]+hamil[12]*(0.75*BstarXdBmag[18]-0.4330127018922194*BstarXdBmag[12]))*rdx2)/m_; 

  double *sgn_alpha_surf3 = &sgn_alpha_surf[243];
  int const_sgn_alpha_surf = 1; 
  
  if (0.4024922359499623*(alphaL3[39]+alphaL3[37]+alphaL3[34]+alphaL3[33]+alphaL3[32])-0.3*(alphaL3[27]+alphaL3[25]+alphaL3[24]+alphaL3[23]+alphaL3[22]+alphaL3[21]+alphaL3[20]+alphaL3[19])-0.603738353924943*(alphaL3[17]+alphaL3[15])+0.2236067977499786*(alphaL3[13]+alphaL3[12]+alphaL3[11])+0.45*(alphaL3[10]+alphaL3[8]+alphaL3[7]+alphaL3[6]+alphaL3[5])-0.3354101966249678*(alphaL3[4]+alphaL3[3]+alphaL3[2]+alphaL3[1])+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[0] = 1.0; 
  else  
    sgn_alpha_surf3[0] = -1.0; 
  
  if (0.4024922359499623*alphaL3[34]+0.4024922359499623*alphaL3[33]+0.4024922359499623*alphaL3[32]-0.3*alphaL3[24]-0.3*alphaL3[23]-0.3*alphaL3[22]-0.3*alphaL3[21]-0.3*alphaL3[20]-0.3*alphaL3[19]-0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[7]+0.45*alphaL3[6]+0.45*alphaL3[5]-0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[1] = 1.0; 
  else  
    sgn_alpha_surf3[1] = -1.0; 
  
  if (sgn_alpha_surf3[1] == sgn_alpha_surf3[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[39])-0.4024922359499623*alphaL3[37]+0.4024922359499623*alphaL3[34]+0.4024922359499623*alphaL3[33]+0.4024922359499623*alphaL3[32]+0.3*alphaL3[27]+0.3*alphaL3[25]-0.3*alphaL3[24]-0.3*alphaL3[23]-0.3*alphaL3[22]-0.3*alphaL3[21]-0.3*alphaL3[20]-0.3*alphaL3[19]+0.603738353924943*alphaL3[17]-0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[10]-0.45*alphaL3[8]+0.45*alphaL3[7]+0.45*alphaL3[6]+0.45*alphaL3[5]+0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[2] = 1.0; 
  else  
    sgn_alpha_surf3[2] = -1.0; 
  
  if (sgn_alpha_surf3[2] == sgn_alpha_surf3[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[39])-0.5031152949374518*alphaL3[34]+0.375*alphaL3[27]-0.3*alphaL3[25]+0.375*alphaL3[24]+0.375*alphaL3[23]-0.3*alphaL3[20]-0.3*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[8]+0.45*alphaL3[5]-0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[3] = 1.0; 
  else  
    sgn_alpha_surf3[3] = -1.0; 
  
  if (sgn_alpha_surf3[3] == sgn_alpha_surf3[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[34])+0.375*alphaL3[24]+0.375*alphaL3[23]-0.3*alphaL3[20]-0.3*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[5]-0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[4] = 1.0; 
  else  
    sgn_alpha_surf3[4] = -1.0; 
  
  if (sgn_alpha_surf3[4] == sgn_alpha_surf3[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[39]-0.5031152949374518*alphaL3[34]-0.375*alphaL3[27]+0.3*alphaL3[25]+0.375*alphaL3[24]+0.375*alphaL3[23]-0.3*alphaL3[20]-0.3*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[8]+0.45*alphaL3[5]+0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[5] = 1.0; 
  else  
    sgn_alpha_surf3[5] = -1.0; 
  
  if (sgn_alpha_surf3[5] == sgn_alpha_surf3[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL3[39]-0.4024922359499623*alphaL3[37]+0.4024922359499623*alphaL3[34]-0.4024922359499623*alphaL3[33]-0.4024922359499623*alphaL3[32]-0.3*alphaL3[27]-0.3*alphaL3[25]-0.3*alphaL3[24]-0.3*alphaL3[23]+0.3*alphaL3[22]+0.3*alphaL3[21]-0.3*alphaL3[20]-0.3*alphaL3[19]+0.603738353924943*alphaL3[17]+0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[10]+0.45*alphaL3[8]-0.45*alphaL3[7]-0.45*alphaL3[6]+0.45*alphaL3[5]-0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[6] = 1.0; 
  else  
    sgn_alpha_surf3[6] = -1.0; 
  
  if (sgn_alpha_surf3[6] == sgn_alpha_surf3[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL3[34]-0.4024922359499623*alphaL3[33]-0.4024922359499623*alphaL3[32]-0.3*alphaL3[24]-0.3*alphaL3[23]+0.3*alphaL3[22]+0.3*alphaL3[21]-0.3*alphaL3[20]-0.3*alphaL3[19]+0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[7]-0.45*alphaL3[6]+0.45*alphaL3[5]+0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[7] = 1.0; 
  else  
    sgn_alpha_surf3[7] = -1.0; 
  
  if (sgn_alpha_surf3[7] == sgn_alpha_surf3[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[39])+0.4024922359499623*alphaL3[37]+0.4024922359499623*alphaL3[34]-0.4024922359499623*alphaL3[33]-0.4024922359499623*alphaL3[32]+0.3*alphaL3[27]+0.3*alphaL3[25]-0.3*alphaL3[24]-0.3*alphaL3[23]+0.3*alphaL3[22]+0.3*alphaL3[21]-0.3*alphaL3[20]-0.3*alphaL3[19]-0.603738353924943*alphaL3[17]+0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[10]-0.45*alphaL3[8]-0.45*alphaL3[7]-0.45*alphaL3[6]+0.45*alphaL3[5]+0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[8] = 1.0; 
  else  
    sgn_alpha_surf3[8] = -1.0; 
  
  if (sgn_alpha_surf3[8] == sgn_alpha_surf3[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL3[39]+0.4024922359499623*alphaL3[37]-0.5031152949374518*alphaL3[33]-0.3*alphaL3[27]-0.3*alphaL3[25]-0.3*alphaL3[23]+0.375*alphaL3[22]-0.3*alphaL3[21]+0.375*alphaL3[20]-0.603738353924943*alphaL3[17]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[10]+0.45*alphaL3[8]+0.45*alphaL3[6]-0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[9] = 1.0; 
  else  
    sgn_alpha_surf3[9] = -1.0; 
  
  if (sgn_alpha_surf3[9] == sgn_alpha_surf3[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[33])-0.3*alphaL3[23]+0.375*alphaL3[22]-0.3*alphaL3[21]+0.375*alphaL3[20]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[6]-0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[10] = 1.0; 
  else  
    sgn_alpha_surf3[10] = -1.0; 
  
  if (sgn_alpha_surf3[10] == sgn_alpha_surf3[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[39])-0.4024922359499623*alphaL3[37]-0.5031152949374518*alphaL3[33]+0.3*alphaL3[27]+0.3*alphaL3[25]-0.3*alphaL3[23]+0.375*alphaL3[22]-0.3*alphaL3[21]+0.375*alphaL3[20]+0.603738353924943*alphaL3[17]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[10]-0.45*alphaL3[8]+0.45*alphaL3[6]+0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[11] = 1.0; 
  else  
    sgn_alpha_surf3[11] = -1.0; 
  
  if (sgn_alpha_surf3[11] == sgn_alpha_surf3[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[39])+0.375*alphaL3[27]-0.3*alphaL3[25]+0.375*alphaL3[23]+0.375*alphaL3[20]-0.2795084971874732*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[8]-0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[12] = 1.0; 
  else  
    sgn_alpha_surf3[12] = -1.0; 
  
  if (sgn_alpha_surf3[12] == sgn_alpha_surf3[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaL3[23]+0.375*alphaL3[20]-0.2795084971874732*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[13] = 1.0; 
  else  
    sgn_alpha_surf3[13] = -1.0; 
  
  if (sgn_alpha_surf3[13] == sgn_alpha_surf3[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[39]-0.375*alphaL3[27]+0.3*alphaL3[25]+0.375*alphaL3[23]+0.375*alphaL3[20]-0.2795084971874732*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[8]+0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[14] = 1.0; 
  else  
    sgn_alpha_surf3[14] = -1.0; 
  
  if (sgn_alpha_surf3[14] == sgn_alpha_surf3[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL3[39]-0.4024922359499623*alphaL3[37]+0.5031152949374518*alphaL3[33]-0.3*alphaL3[27]-0.3*alphaL3[25]-0.3*alphaL3[23]-0.375*alphaL3[22]+0.3*alphaL3[21]+0.375*alphaL3[20]+0.603738353924943*alphaL3[17]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[10]+0.45*alphaL3[8]-0.45*alphaL3[6]-0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[15] = 1.0; 
  else  
    sgn_alpha_surf3[15] = -1.0; 
  
  if (sgn_alpha_surf3[15] == sgn_alpha_surf3[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[33]-0.3*alphaL3[23]-0.375*alphaL3[22]+0.3*alphaL3[21]+0.375*alphaL3[20]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[6]+0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[16] = 1.0; 
  else  
    sgn_alpha_surf3[16] = -1.0; 
  
  if (sgn_alpha_surf3[16] == sgn_alpha_surf3[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[39])+0.4024922359499623*alphaL3[37]+0.5031152949374518*alphaL3[33]+0.3*alphaL3[27]+0.3*alphaL3[25]-0.3*alphaL3[23]-0.375*alphaL3[22]+0.3*alphaL3[21]+0.375*alphaL3[20]-0.603738353924943*alphaL3[17]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[10]-0.45*alphaL3[8]-0.45*alphaL3[6]+0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[17] = 1.0; 
  else  
    sgn_alpha_surf3[17] = -1.0; 
  
  if (sgn_alpha_surf3[17] == sgn_alpha_surf3[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL3[39]+0.4024922359499623*alphaL3[37]-0.4024922359499623*alphaL3[34]+0.4024922359499623*alphaL3[33]-0.4024922359499623*alphaL3[32]-0.3*alphaL3[27]-0.3*alphaL3[25]+0.3*alphaL3[24]-0.3*alphaL3[23]-0.3*alphaL3[22]-0.3*alphaL3[21]-0.3*alphaL3[20]+0.3*alphaL3[19]-0.603738353924943*alphaL3[17]+0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[10]+0.45*alphaL3[8]-0.45*alphaL3[7]+0.45*alphaL3[6]-0.45*alphaL3[5]-0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[18] = 1.0; 
  else  
    sgn_alpha_surf3[18] = -1.0; 
  
  if (sgn_alpha_surf3[18] == sgn_alpha_surf3[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[34])+0.4024922359499623*alphaL3[33]-0.4024922359499623*alphaL3[32]+0.3*alphaL3[24]-0.3*alphaL3[23]-0.3*alphaL3[22]-0.3*alphaL3[21]-0.3*alphaL3[20]+0.3*alphaL3[19]+0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[7]+0.45*alphaL3[6]-0.45*alphaL3[5]-0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[19] = 1.0; 
  else  
    sgn_alpha_surf3[19] = -1.0; 
  
  if (sgn_alpha_surf3[19] == sgn_alpha_surf3[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[39])-0.4024922359499623*alphaL3[37]-0.4024922359499623*alphaL3[34]+0.4024922359499623*alphaL3[33]-0.4024922359499623*alphaL3[32]+0.3*alphaL3[27]+0.3*alphaL3[25]+0.3*alphaL3[24]-0.3*alphaL3[23]-0.3*alphaL3[22]-0.3*alphaL3[21]-0.3*alphaL3[20]+0.3*alphaL3[19]+0.603738353924943*alphaL3[17]+0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[10]-0.45*alphaL3[8]-0.45*alphaL3[7]+0.45*alphaL3[6]-0.45*alphaL3[5]+0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[20] = 1.0; 
  else  
    sgn_alpha_surf3[20] = -1.0; 
  
  if (sgn_alpha_surf3[20] == sgn_alpha_surf3[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[39])+0.5031152949374518*alphaL3[34]+0.375*alphaL3[27]-0.3*alphaL3[25]-0.375*alphaL3[24]+0.375*alphaL3[23]-0.3*alphaL3[20]+0.3*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[8]-0.45*alphaL3[5]-0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[21] = 1.0; 
  else  
    sgn_alpha_surf3[21] = -1.0; 
  
  if (sgn_alpha_surf3[21] == sgn_alpha_surf3[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[34]-0.375*alphaL3[24]+0.375*alphaL3[23]-0.3*alphaL3[20]+0.3*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[5]+0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[22] = 1.0; 
  else  
    sgn_alpha_surf3[22] = -1.0; 
  
  if (sgn_alpha_surf3[22] == sgn_alpha_surf3[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[39]+0.5031152949374518*alphaL3[34]-0.375*alphaL3[27]+0.3*alphaL3[25]-0.375*alphaL3[24]+0.375*alphaL3[23]-0.3*alphaL3[20]+0.3*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[8]-0.45*alphaL3[5]+0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[23] = 1.0; 
  else  
    sgn_alpha_surf3[23] = -1.0; 
  
  if (sgn_alpha_surf3[23] == sgn_alpha_surf3[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL3[39]-0.4024922359499623*alphaL3[37]-0.4024922359499623*alphaL3[34]-0.4024922359499623*alphaL3[33]+0.4024922359499623*alphaL3[32]-0.3*alphaL3[27]-0.3*alphaL3[25]+0.3*alphaL3[24]-0.3*alphaL3[23]+0.3*alphaL3[22]+0.3*alphaL3[21]-0.3*alphaL3[20]+0.3*alphaL3[19]+0.603738353924943*alphaL3[17]-0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[10]+0.45*alphaL3[8]+0.45*alphaL3[7]-0.45*alphaL3[6]-0.45*alphaL3[5]-0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[24] = 1.0; 
  else  
    sgn_alpha_surf3[24] = -1.0; 
  
  if (sgn_alpha_surf3[24] == sgn_alpha_surf3[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[34])-0.4024922359499623*alphaL3[33]+0.4024922359499623*alphaL3[32]+0.3*alphaL3[24]-0.3*alphaL3[23]+0.3*alphaL3[22]+0.3*alphaL3[21]-0.3*alphaL3[20]+0.3*alphaL3[19]-0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[7]-0.45*alphaL3[6]-0.45*alphaL3[5]+0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[25] = 1.0; 
  else  
    sgn_alpha_surf3[25] = -1.0; 
  
  if (sgn_alpha_surf3[25] == sgn_alpha_surf3[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[39])+0.4024922359499623*alphaL3[37]-0.4024922359499623*alphaL3[34]-0.4024922359499623*alphaL3[33]+0.4024922359499623*alphaL3[32]+0.3*alphaL3[27]+0.3*alphaL3[25]+0.3*alphaL3[24]-0.3*alphaL3[23]+0.3*alphaL3[22]+0.3*alphaL3[21]-0.3*alphaL3[20]+0.3*alphaL3[19]-0.603738353924943*alphaL3[17]-0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[10]-0.45*alphaL3[8]+0.45*alphaL3[7]-0.45*alphaL3[6]-0.45*alphaL3[5]+0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]-0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[26] = 1.0; 
  else  
    sgn_alpha_surf3[26] = -1.0; 
  
  if (sgn_alpha_surf3[26] == sgn_alpha_surf3[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[37])-0.5031152949374518*alphaL3[32]-0.3*alphaL3[27]+0.375*alphaL3[25]-0.3*alphaL3[24]-0.3*alphaL3[22]+0.375*alphaL3[21]+0.375*alphaL3[19]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]+0.45*alphaL3[10]+0.45*alphaL3[7]-0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[27] = 1.0; 
  else  
    sgn_alpha_surf3[27] = -1.0; 
  
  if (sgn_alpha_surf3[27] == sgn_alpha_surf3[26]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[32])-0.3*alphaL3[24]-0.3*alphaL3[22]+0.375*alphaL3[21]+0.375*alphaL3[19]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]+0.45*alphaL3[7]-0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[28] = 1.0; 
  else  
    sgn_alpha_surf3[28] = -1.0; 
  
  if (sgn_alpha_surf3[28] == sgn_alpha_surf3[27]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[37]-0.5031152949374518*alphaL3[32]+0.3*alphaL3[27]-0.375*alphaL3[25]-0.3*alphaL3[24]-0.3*alphaL3[22]+0.375*alphaL3[21]+0.375*alphaL3[19]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]-0.45*alphaL3[10]+0.45*alphaL3[7]+0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[29] = 1.0; 
  else  
    sgn_alpha_surf3[29] = -1.0; 
  
  if (sgn_alpha_surf3[29] == sgn_alpha_surf3[28]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaL3[27]+0.375*alphaL3[25]+0.375*alphaL3[24]+0.375*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]-0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[30] = 1.0; 
  else  
    sgn_alpha_surf3[30] = -1.0; 
  
  if (sgn_alpha_surf3[30] == sgn_alpha_surf3[29]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaL3[24]+0.375*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]-0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[31] = 1.0; 
  else  
    sgn_alpha_surf3[31] = -1.0; 
  
  if (sgn_alpha_surf3[31] == sgn_alpha_surf3[30]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*alphaL3[27])-0.375*alphaL3[25]+0.375*alphaL3[24]+0.375*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]+0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[32] = 1.0; 
  else  
    sgn_alpha_surf3[32] = -1.0; 
  
  if (sgn_alpha_surf3[32] == sgn_alpha_surf3[31]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[37]+0.5031152949374518*alphaL3[32]-0.3*alphaL3[27]+0.375*alphaL3[25]-0.3*alphaL3[24]+0.3*alphaL3[22]-0.375*alphaL3[21]+0.375*alphaL3[19]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]-0.45*alphaL3[10]-0.45*alphaL3[7]-0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[33] = 1.0; 
  else  
    sgn_alpha_surf3[33] = -1.0; 
  
  if (sgn_alpha_surf3[33] == sgn_alpha_surf3[32]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[32]-0.3*alphaL3[24]+0.3*alphaL3[22]-0.375*alphaL3[21]+0.375*alphaL3[19]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]-0.45*alphaL3[7]+0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[34] = 1.0; 
  else  
    sgn_alpha_surf3[34] = -1.0; 
  
  if (sgn_alpha_surf3[34] == sgn_alpha_surf3[33]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[37])+0.5031152949374518*alphaL3[32]+0.3*alphaL3[27]-0.375*alphaL3[25]-0.3*alphaL3[24]+0.3*alphaL3[22]-0.375*alphaL3[21]+0.375*alphaL3[19]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]+0.45*alphaL3[10]-0.45*alphaL3[7]+0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[35] = 1.0; 
  else  
    sgn_alpha_surf3[35] = -1.0; 
  
  if (sgn_alpha_surf3[35] == sgn_alpha_surf3[34]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[37])-0.3*alphaL3[27]+0.375*alphaL3[25]+0.375*alphaL3[22]+0.375*alphaL3[21]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]-0.2795084971874732*alphaL3[11]+0.45*alphaL3[10]-0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[36] = 1.0; 
  else  
    sgn_alpha_surf3[36] = -1.0; 
  
  if (sgn_alpha_surf3[36] == sgn_alpha_surf3[35]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaL3[22]+0.375*alphaL3[21]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]-0.2795084971874732*alphaL3[11]-0.3354101966249678*alphaL3[3]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[37] = 1.0; 
  else  
    sgn_alpha_surf3[37] = -1.0; 
  
  if (sgn_alpha_surf3[37] == sgn_alpha_surf3[36]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[37]+0.3*alphaL3[27]-0.375*alphaL3[25]+0.375*alphaL3[22]+0.375*alphaL3[21]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]-0.2795084971874732*alphaL3[11]-0.45*alphaL3[10]+0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[38] = 1.0; 
  else  
    sgn_alpha_surf3[38] = -1.0; 
  
  if (sgn_alpha_surf3[38] == sgn_alpha_surf3[37]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaL3[27]+0.375*alphaL3[25]-0.2795084971874732*alphaL3[13]-0.2795084971874732*alphaL3[12]-0.2795084971874732*alphaL3[11]-0.3354101966249678*alphaL3[4]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[39] = 1.0; 
  else  
    sgn_alpha_surf3[39] = -1.0; 
  
  if (sgn_alpha_surf3[39] == sgn_alpha_surf3[38]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.2795084971874732*alphaL3[13])-0.2795084971874732*alphaL3[12]-0.2795084971874732*alphaL3[11]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[40] = 1.0; 
  else  
    sgn_alpha_surf3[40] = -1.0; 
  
  if (sgn_alpha_surf3[40] == sgn_alpha_surf3[39]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*alphaL3[27])-0.375*alphaL3[25]-0.2795084971874732*alphaL3[13]-0.2795084971874732*alphaL3[12]-0.2795084971874732*alphaL3[11]+0.3354101966249678*alphaL3[4]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[41] = 1.0; 
  else  
    sgn_alpha_surf3[41] = -1.0; 
  
  if (sgn_alpha_surf3[41] == sgn_alpha_surf3[40]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[37]-0.3*alphaL3[27]+0.375*alphaL3[25]-0.375*alphaL3[22]-0.375*alphaL3[21]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]-0.2795084971874732*alphaL3[11]-0.45*alphaL3[10]-0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[42] = 1.0; 
  else  
    sgn_alpha_surf3[42] = -1.0; 
  
  if (sgn_alpha_surf3[42] == sgn_alpha_surf3[41]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*alphaL3[22])-0.375*alphaL3[21]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]-0.2795084971874732*alphaL3[11]+0.3354101966249678*alphaL3[3]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[43] = 1.0; 
  else  
    sgn_alpha_surf3[43] = -1.0; 
  
  if (sgn_alpha_surf3[43] == sgn_alpha_surf3[42]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[37])+0.3*alphaL3[27]-0.375*alphaL3[25]-0.375*alphaL3[22]-0.375*alphaL3[21]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]-0.2795084971874732*alphaL3[11]+0.45*alphaL3[10]+0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[44] = 1.0; 
  else  
    sgn_alpha_surf3[44] = -1.0; 
  
  if (sgn_alpha_surf3[44] == sgn_alpha_surf3[43]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[37])+0.5031152949374518*alphaL3[32]-0.3*alphaL3[27]+0.375*alphaL3[25]+0.3*alphaL3[24]-0.3*alphaL3[22]+0.375*alphaL3[21]-0.375*alphaL3[19]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]+0.45*alphaL3[10]-0.45*alphaL3[7]-0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[45] = 1.0; 
  else  
    sgn_alpha_surf3[45] = -1.0; 
  
  if (sgn_alpha_surf3[45] == sgn_alpha_surf3[44]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[32]+0.3*alphaL3[24]-0.3*alphaL3[22]+0.375*alphaL3[21]-0.375*alphaL3[19]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]-0.45*alphaL3[7]-0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[46] = 1.0; 
  else  
    sgn_alpha_surf3[46] = -1.0; 
  
  if (sgn_alpha_surf3[46] == sgn_alpha_surf3[45]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[37]+0.5031152949374518*alphaL3[32]+0.3*alphaL3[27]-0.375*alphaL3[25]+0.3*alphaL3[24]-0.3*alphaL3[22]+0.375*alphaL3[21]-0.375*alphaL3[19]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]-0.45*alphaL3[10]-0.45*alphaL3[7]+0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[47] = 1.0; 
  else  
    sgn_alpha_surf3[47] = -1.0; 
  
  if (sgn_alpha_surf3[47] == sgn_alpha_surf3[46]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.375*alphaL3[27]+0.375*alphaL3[25]-0.375*alphaL3[24]-0.375*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]-0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[48] = 1.0; 
  else  
    sgn_alpha_surf3[48] = -1.0; 
  
  if (sgn_alpha_surf3[48] == sgn_alpha_surf3[47]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*alphaL3[24])-0.375*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]+0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[49] = 1.0; 
  else  
    sgn_alpha_surf3[49] = -1.0; 
  
  if (sgn_alpha_surf3[49] == sgn_alpha_surf3[48]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*alphaL3[27])-0.375*alphaL3[25]-0.375*alphaL3[24]-0.375*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]+0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[50] = 1.0; 
  else  
    sgn_alpha_surf3[50] = -1.0; 
  
  if (sgn_alpha_surf3[50] == sgn_alpha_surf3[49]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[37]-0.5031152949374518*alphaL3[32]-0.3*alphaL3[27]+0.375*alphaL3[25]+0.3*alphaL3[24]+0.3*alphaL3[22]-0.375*alphaL3[21]-0.375*alphaL3[19]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]-0.45*alphaL3[10]+0.45*alphaL3[7]-0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[51] = 1.0; 
  else  
    sgn_alpha_surf3[51] = -1.0; 
  
  if (sgn_alpha_surf3[51] == sgn_alpha_surf3[50]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[32])+0.3*alphaL3[24]+0.3*alphaL3[22]-0.375*alphaL3[21]-0.375*alphaL3[19]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]+0.45*alphaL3[7]+0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[52] = 1.0; 
  else  
    sgn_alpha_surf3[52] = -1.0; 
  
  if (sgn_alpha_surf3[52] == sgn_alpha_surf3[51]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[37])-0.5031152949374518*alphaL3[32]+0.3*alphaL3[27]-0.375*alphaL3[25]+0.3*alphaL3[24]+0.3*alphaL3[22]-0.375*alphaL3[21]-0.375*alphaL3[19]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]-0.2795084971874732*alphaL3[11]+0.45*alphaL3[10]+0.45*alphaL3[7]+0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[53] = 1.0; 
  else  
    sgn_alpha_surf3[53] = -1.0; 
  
  if (sgn_alpha_surf3[53] == sgn_alpha_surf3[52]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[39])+0.4024922359499623*alphaL3[37]-0.4024922359499623*alphaL3[34]-0.4024922359499623*alphaL3[33]+0.4024922359499623*alphaL3[32]-0.3*alphaL3[27]-0.3*alphaL3[25]-0.3*alphaL3[24]+0.3*alphaL3[23]-0.3*alphaL3[22]-0.3*alphaL3[21]+0.3*alphaL3[20]-0.3*alphaL3[19]+0.603738353924943*alphaL3[17]+0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[10]-0.45*alphaL3[8]+0.45*alphaL3[7]-0.45*alphaL3[6]-0.45*alphaL3[5]-0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[54] = 1.0; 
  else  
    sgn_alpha_surf3[54] = -1.0; 
  
  if (sgn_alpha_surf3[54] == sgn_alpha_surf3[53]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[34])-0.4024922359499623*alphaL3[33]+0.4024922359499623*alphaL3[32]-0.3*alphaL3[24]+0.3*alphaL3[23]-0.3*alphaL3[22]-0.3*alphaL3[21]+0.3*alphaL3[20]-0.3*alphaL3[19]+0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[7]-0.45*alphaL3[6]-0.45*alphaL3[5]-0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[55] = 1.0; 
  else  
    sgn_alpha_surf3[55] = -1.0; 
  
  if (sgn_alpha_surf3[55] == sgn_alpha_surf3[54]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL3[39]-0.4024922359499623*alphaL3[37]-0.4024922359499623*alphaL3[34]-0.4024922359499623*alphaL3[33]+0.4024922359499623*alphaL3[32]+0.3*alphaL3[27]+0.3*alphaL3[25]-0.3*alphaL3[24]+0.3*alphaL3[23]-0.3*alphaL3[22]-0.3*alphaL3[21]+0.3*alphaL3[20]-0.3*alphaL3[19]-0.603738353924943*alphaL3[17]+0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[10]+0.45*alphaL3[8]+0.45*alphaL3[7]-0.45*alphaL3[6]-0.45*alphaL3[5]+0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[56] = 1.0; 
  else  
    sgn_alpha_surf3[56] = -1.0; 
  
  if (sgn_alpha_surf3[56] == sgn_alpha_surf3[55]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[39]+0.5031152949374518*alphaL3[34]+0.375*alphaL3[27]-0.3*alphaL3[25]+0.375*alphaL3[24]-0.375*alphaL3[23]+0.3*alphaL3[20]-0.3*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[8]-0.45*alphaL3[5]-0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[57] = 1.0; 
  else  
    sgn_alpha_surf3[57] = -1.0; 
  
  if (sgn_alpha_surf3[57] == sgn_alpha_surf3[56]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[34]+0.375*alphaL3[24]-0.375*alphaL3[23]+0.3*alphaL3[20]-0.3*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[5]-0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[58] = 1.0; 
  else  
    sgn_alpha_surf3[58] = -1.0; 
  
  if (sgn_alpha_surf3[58] == sgn_alpha_surf3[57]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[39])+0.5031152949374518*alphaL3[34]-0.375*alphaL3[27]+0.3*alphaL3[25]+0.375*alphaL3[24]-0.375*alphaL3[23]+0.3*alphaL3[20]-0.3*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[8]-0.45*alphaL3[5]+0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[59] = 1.0; 
  else  
    sgn_alpha_surf3[59] = -1.0; 
  
  if (sgn_alpha_surf3[59] == sgn_alpha_surf3[58]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[39])-0.4024922359499623*alphaL3[37]-0.4024922359499623*alphaL3[34]+0.4024922359499623*alphaL3[33]-0.4024922359499623*alphaL3[32]-0.3*alphaL3[27]-0.3*alphaL3[25]-0.3*alphaL3[24]+0.3*alphaL3[23]+0.3*alphaL3[22]+0.3*alphaL3[21]+0.3*alphaL3[20]-0.3*alphaL3[19]-0.603738353924943*alphaL3[17]-0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[10]-0.45*alphaL3[8]-0.45*alphaL3[7]+0.45*alphaL3[6]-0.45*alphaL3[5]-0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[60] = 1.0; 
  else  
    sgn_alpha_surf3[60] = -1.0; 
  
  if (sgn_alpha_surf3[60] == sgn_alpha_surf3[59]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[34])+0.4024922359499623*alphaL3[33]-0.4024922359499623*alphaL3[32]-0.3*alphaL3[24]+0.3*alphaL3[23]+0.3*alphaL3[22]+0.3*alphaL3[21]+0.3*alphaL3[20]-0.3*alphaL3[19]-0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[7]+0.45*alphaL3[6]-0.45*alphaL3[5]+0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[61] = 1.0; 
  else  
    sgn_alpha_surf3[61] = -1.0; 
  
  if (sgn_alpha_surf3[61] == sgn_alpha_surf3[60]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL3[39]+0.4024922359499623*alphaL3[37]-0.4024922359499623*alphaL3[34]+0.4024922359499623*alphaL3[33]-0.4024922359499623*alphaL3[32]+0.3*alphaL3[27]+0.3*alphaL3[25]-0.3*alphaL3[24]+0.3*alphaL3[23]+0.3*alphaL3[22]+0.3*alphaL3[21]+0.3*alphaL3[20]-0.3*alphaL3[19]+0.603738353924943*alphaL3[17]-0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[10]+0.45*alphaL3[8]-0.45*alphaL3[7]+0.45*alphaL3[6]-0.45*alphaL3[5]+0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]-0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[62] = 1.0; 
  else  
    sgn_alpha_surf3[62] = -1.0; 
  
  if (sgn_alpha_surf3[62] == sgn_alpha_surf3[61]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[39])+0.4024922359499623*alphaL3[37]+0.5031152949374518*alphaL3[33]-0.3*alphaL3[27]-0.3*alphaL3[25]+0.3*alphaL3[23]+0.375*alphaL3[22]-0.3*alphaL3[21]-0.375*alphaL3[20]+0.603738353924943*alphaL3[17]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[10]-0.45*alphaL3[8]-0.45*alphaL3[6]-0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[63] = 1.0; 
  else  
    sgn_alpha_surf3[63] = -1.0; 
  
  if (sgn_alpha_surf3[63] == sgn_alpha_surf3[62]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[33]+0.3*alphaL3[23]+0.375*alphaL3[22]-0.3*alphaL3[21]-0.375*alphaL3[20]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[6]-0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[64] = 1.0; 
  else  
    sgn_alpha_surf3[64] = -1.0; 
  
  if (sgn_alpha_surf3[64] == sgn_alpha_surf3[63]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL3[39]-0.4024922359499623*alphaL3[37]+0.5031152949374518*alphaL3[33]+0.3*alphaL3[27]+0.3*alphaL3[25]+0.3*alphaL3[23]+0.375*alphaL3[22]-0.3*alphaL3[21]-0.375*alphaL3[20]-0.603738353924943*alphaL3[17]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[10]+0.45*alphaL3[8]-0.45*alphaL3[6]+0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[65] = 1.0; 
  else  
    sgn_alpha_surf3[65] = -1.0; 
  
  if (sgn_alpha_surf3[65] == sgn_alpha_surf3[64]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[39]+0.375*alphaL3[27]-0.3*alphaL3[25]-0.375*alphaL3[23]-0.375*alphaL3[20]-0.2795084971874732*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[8]-0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[66] = 1.0; 
  else  
    sgn_alpha_surf3[66] = -1.0; 
  
  if (sgn_alpha_surf3[66] == sgn_alpha_surf3[65]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.375*alphaL3[23])-0.375*alphaL3[20]-0.2795084971874732*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[67] = 1.0; 
  else  
    sgn_alpha_surf3[67] = -1.0; 
  
  if (sgn_alpha_surf3[67] == sgn_alpha_surf3[66]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[39])-0.375*alphaL3[27]+0.3*alphaL3[25]-0.375*alphaL3[23]-0.375*alphaL3[20]-0.2795084971874732*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[8]+0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[68] = 1.0; 
  else  
    sgn_alpha_surf3[68] = -1.0; 
  
  if (sgn_alpha_surf3[68] == sgn_alpha_surf3[67]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[39])-0.4024922359499623*alphaL3[37]-0.5031152949374518*alphaL3[33]-0.3*alphaL3[27]-0.3*alphaL3[25]+0.3*alphaL3[23]-0.375*alphaL3[22]+0.3*alphaL3[21]-0.375*alphaL3[20]-0.603738353924943*alphaL3[17]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[10]-0.45*alphaL3[8]+0.45*alphaL3[6]-0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[69] = 1.0; 
  else  
    sgn_alpha_surf3[69] = -1.0; 
  
  if (sgn_alpha_surf3[69] == sgn_alpha_surf3[68]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[33])+0.3*alphaL3[23]-0.375*alphaL3[22]+0.3*alphaL3[21]-0.375*alphaL3[20]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[6]+0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[70] = 1.0; 
  else  
    sgn_alpha_surf3[70] = -1.0; 
  
  if (sgn_alpha_surf3[70] == sgn_alpha_surf3[69]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL3[39]+0.4024922359499623*alphaL3[37]-0.5031152949374518*alphaL3[33]+0.3*alphaL3[27]+0.3*alphaL3[25]+0.3*alphaL3[23]-0.375*alphaL3[22]+0.3*alphaL3[21]-0.375*alphaL3[20]+0.603738353924943*alphaL3[17]+0.2236067977499786*alphaL3[13]-0.2795084971874732*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[10]+0.45*alphaL3[8]+0.45*alphaL3[6]+0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[71] = 1.0; 
  else  
    sgn_alpha_surf3[71] = -1.0; 
  
  if (sgn_alpha_surf3[71] == sgn_alpha_surf3[70]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[39])+0.4024922359499623*alphaL3[37]+0.4024922359499623*alphaL3[34]-0.4024922359499623*alphaL3[33]-0.4024922359499623*alphaL3[32]-0.3*alphaL3[27]-0.3*alphaL3[25]+0.3*alphaL3[24]+0.3*alphaL3[23]-0.3*alphaL3[22]-0.3*alphaL3[21]+0.3*alphaL3[20]+0.3*alphaL3[19]+0.603738353924943*alphaL3[17]-0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[10]-0.45*alphaL3[8]-0.45*alphaL3[7]-0.45*alphaL3[6]+0.45*alphaL3[5]-0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[72] = 1.0; 
  else  
    sgn_alpha_surf3[72] = -1.0; 
  
  if (sgn_alpha_surf3[72] == sgn_alpha_surf3[71]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL3[34]-0.4024922359499623*alphaL3[33]-0.4024922359499623*alphaL3[32]+0.3*alphaL3[24]+0.3*alphaL3[23]-0.3*alphaL3[22]-0.3*alphaL3[21]+0.3*alphaL3[20]+0.3*alphaL3[19]-0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[7]-0.45*alphaL3[6]+0.45*alphaL3[5]-0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[73] = 1.0; 
  else  
    sgn_alpha_surf3[73] = -1.0; 
  
  if (sgn_alpha_surf3[73] == sgn_alpha_surf3[72]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL3[39]-0.4024922359499623*alphaL3[37]+0.4024922359499623*alphaL3[34]-0.4024922359499623*alphaL3[33]-0.4024922359499623*alphaL3[32]+0.3*alphaL3[27]+0.3*alphaL3[25]+0.3*alphaL3[24]+0.3*alphaL3[23]-0.3*alphaL3[22]-0.3*alphaL3[21]+0.3*alphaL3[20]+0.3*alphaL3[19]-0.603738353924943*alphaL3[17]-0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[10]+0.45*alphaL3[8]-0.45*alphaL3[7]-0.45*alphaL3[6]+0.45*alphaL3[5]+0.3354101966249678*alphaL3[4]-0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[74] = 1.0; 
  else  
    sgn_alpha_surf3[74] = -1.0; 
  
  if (sgn_alpha_surf3[74] == sgn_alpha_surf3[73]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5031152949374518*alphaL3[39]-0.5031152949374518*alphaL3[34]+0.375*alphaL3[27]-0.3*alphaL3[25]-0.375*alphaL3[24]-0.375*alphaL3[23]+0.3*alphaL3[20]+0.3*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[8]+0.45*alphaL3[5]-0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[75] = 1.0; 
  else  
    sgn_alpha_surf3[75] = -1.0; 
  
  if (sgn_alpha_surf3[75] == sgn_alpha_surf3[74]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[34])-0.375*alphaL3[24]-0.375*alphaL3[23]+0.3*alphaL3[20]+0.3*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[5]+0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[76] = 1.0; 
  else  
    sgn_alpha_surf3[76] = -1.0; 
  
  if (sgn_alpha_surf3[76] == sgn_alpha_surf3[75]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5031152949374518*alphaL3[39])-0.5031152949374518*alphaL3[34]-0.375*alphaL3[27]+0.3*alphaL3[25]-0.375*alphaL3[24]-0.375*alphaL3[23]+0.3*alphaL3[20]+0.3*alphaL3[19]-0.2795084971874732*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[8]+0.45*alphaL3[5]+0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[77] = 1.0; 
  else  
    sgn_alpha_surf3[77] = -1.0; 
  
  if (sgn_alpha_surf3[77] == sgn_alpha_surf3[76]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4024922359499623*alphaL3[39])-0.4024922359499623*alphaL3[37]+0.4024922359499623*alphaL3[34]+0.4024922359499623*alphaL3[33]+0.4024922359499623*alphaL3[32]-0.3*alphaL3[27]-0.3*alphaL3[25]+0.3*alphaL3[24]+0.3*alphaL3[23]+0.3*alphaL3[22]+0.3*alphaL3[21]+0.3*alphaL3[20]+0.3*alphaL3[19]-0.603738353924943*alphaL3[17]+0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]-0.45*alphaL3[10]-0.45*alphaL3[8]+0.45*alphaL3[7]+0.45*alphaL3[6]+0.45*alphaL3[5]-0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[78] = 1.0; 
  else  
    sgn_alpha_surf3[78] = -1.0; 
  
  if (sgn_alpha_surf3[78] == sgn_alpha_surf3[77]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL3[34]+0.4024922359499623*alphaL3[33]+0.4024922359499623*alphaL3[32]+0.3*alphaL3[24]+0.3*alphaL3[23]+0.3*alphaL3[22]+0.3*alphaL3[21]+0.3*alphaL3[20]+0.3*alphaL3[19]+0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[7]+0.45*alphaL3[6]+0.45*alphaL3[5]+0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[79] = 1.0; 
  else  
    sgn_alpha_surf3[79] = -1.0; 
  
  if (sgn_alpha_surf3[79] == sgn_alpha_surf3[78]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4024922359499623*alphaL3[39]+0.4024922359499623*alphaL3[37]+0.4024922359499623*alphaL3[34]+0.4024922359499623*alphaL3[33]+0.4024922359499623*alphaL3[32]+0.3*alphaL3[27]+0.3*alphaL3[25]+0.3*alphaL3[24]+0.3*alphaL3[23]+0.3*alphaL3[22]+0.3*alphaL3[21]+0.3*alphaL3[20]+0.3*alphaL3[19]+0.603738353924943*alphaL3[17]+0.603738353924943*alphaL3[15]+0.2236067977499786*alphaL3[13]+0.2236067977499786*alphaL3[12]+0.2236067977499786*alphaL3[11]+0.45*alphaL3[10]+0.45*alphaL3[8]+0.45*alphaL3[7]+0.45*alphaL3[6]+0.45*alphaL3[5]+0.3354101966249678*alphaL3[4]+0.3354101966249678*alphaL3[3]+0.3354101966249678*alphaL3[2]+0.3354101966249678*alphaL3[1]+0.25*alphaL3[0] > 0.) 
    sgn_alpha_surf3[80] = 1.0; 
  else  
    sgn_alpha_surf3[80] = -1.0; 
  
  if (sgn_alpha_surf3[80] == sgn_alpha_surf3[79]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
