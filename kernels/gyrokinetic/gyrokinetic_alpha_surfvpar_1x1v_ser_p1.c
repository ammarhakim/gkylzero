#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, 
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
  double wvpar = w[1];
  double rdvpar2 = 2.0/dxv[1];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wvparSq = w[1]*w[1];
  double rdvpar2Sq = rdvpar2*rdvpar2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[2];
  const double *b_z = &b_i[4];

  double hamil[6] = {0.}; 
  hamil[0] = m_*wvparSq+(0.3333333333333333*m_)/rdvpar2Sq+1.414213562373095*phi[0]*q_; 
  hamil[1] = 1.414213562373095*phi[1]*q_; 
  hamil[2] = (1.154700538379252*m_*wvpar)/rdvpar2; 
  hamil[4] = (0.2981423969999719*m_)/rdvpar2Sq; 

  double *alphaL1 = &alpha_surf[3];
  alphaL1[0] = (((-0.6123724356957944*cmag[1]*hamil[1]*jacobtot_inv[1])-0.6123724356957944*cmag[0]*jacobtot_inv[0]*hamil[1])*rdx2)/m_; 
  alphaL1[1] = (((-0.6123724356957944*cmag[0]*hamil[1]*jacobtot_inv[1])-0.6123724356957944*jacobtot_inv[0]*cmag[1]*hamil[1])*rdx2)/m_; 

  double *sgn_alpha_surf1 = &sgn_alpha_surf[3];
  int const_sgn_alpha_surf = 1; 
  
  if (0.7071067811865468*alphaL1[0]-0.7071067811865468*alphaL1[1] > 0.) 
    sgn_alpha_surf1[0] = 1.0; 
  else  
    sgn_alpha_surf1[0] = -1.0; 
  
  if (0.7071067811865468*(alphaL1[1]+alphaL1[0]) > 0.) 
    sgn_alpha_surf1[1] = 1.0; 
  else  
    sgn_alpha_surf1[1] = -1.0; 
  
  if (sgn_alpha_surf1[1] == sgn_alpha_surf1[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
