#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_alpha_edge_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *Bstar_Bmag, double* GKYL_RESTRICT alpha_surf) 
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
  double wvpar = w[1];
  double rdvpar2 = 2.0/dxv[1];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wvparSq = w[1]*w[1];
  double rdvpar2Sq = rdvpar2*rdvpar2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[2];
  const double *b_z = &b_i[4];

  const double *BstarZdBmag = &Bstar_Bmag[0]; 

  double hamil[6] = {0.}; 
  hamil[0] = m_*wvparSq+(0.3333333333333333*m_)/rdvpar2Sq+1.414213562373095*phi[0]*q_; 
  hamil[1] = 1.414213562373095*phi[1]*q_; 
  hamil[2] = (1.154700538379252*m_*wvpar)/rdvpar2; 
  hamil[4] = (0.2981423969999719*m_)/rdvpar2Sq; 

  double *alphaR0 = &alpha_surf[0];
  alphaR0[0] = ((1.5*BstarZdBmag[1]+0.8660254037844386*BstarZdBmag[0])*hamil[2]*rdvpar2)/m_; 
  alphaR0[1] = ((3.354101966249685*BstarZdBmag[1]+1.936491673103709*BstarZdBmag[0])*hamil[4]*rdvpar2)/m_; 

} 
