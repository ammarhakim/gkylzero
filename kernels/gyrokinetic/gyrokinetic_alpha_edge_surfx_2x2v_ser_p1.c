#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_alpha_edge_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *Bstar_Bmag, double* GKYL_RESTRICT alpha_surf) 
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
  const double *b_y = &b_i[4];
  const double *b_z = &b_i[8];

  const double *BstarXdBmag = &Bstar_Bmag[0]; 

  double hamil[24] = {0.}; 
  hamil[0] = 2.0*m_*wvparSq+2.0*bmag[0]*wmu+(0.6666666666666666*m_)/rdvpar2Sq+2.0*phi[0]*q_; 
  hamil[1] = 2.0*bmag[1]*wmu+2.0*phi[1]*q_; 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[16] = (0.5962847939999438*m_)/rdvpar2Sq; 

  double *alphaR0 = &alpha_surf[0];
  alphaR0[0] = (((((-1.590990257669731*b_z[1])-0.9185586535436913*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*((-0.9185586535436913*b_z[1])-0.5303300858899105*b_z[0]))*hamil[5]+(((-0.9185586535436913*b_z[1])-0.5303300858899105*b_z[0])*jacobtot_inv[1]+jacobtot_inv[0]*((-0.5303300858899105*b_z[1])-0.3061862178478971*b_z[0]))*hamil[2])*rdy2)/q_+(((2.371708245126284*BstarXdBmag[3]+1.369306393762915*BstarXdBmag[2])*hamil[16]+(1.060660171779821*BstarXdBmag[1]+0.6123724356957944*BstarXdBmag[0])*hamil[3])*rdvpar2)/m_; 
  alphaR0[2] = (((2.371708245126284*BstarXdBmag[1]+1.369306393762915*BstarXdBmag[0])*hamil[16]+(1.060660171779821*BstarXdBmag[3]+0.6123724356957944*BstarXdBmag[2])*hamil[3])*rdvpar2)/m_; 
  alphaR0[8] = ((2.121320343559642*BstarXdBmag[3]+1.224744871391589*BstarXdBmag[2])*hamil[16]*rdvpar2)/m_; 

} 
