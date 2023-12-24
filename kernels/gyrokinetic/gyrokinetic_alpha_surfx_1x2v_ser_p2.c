#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_surfx_1x2v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, 
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
  double wvpar = w[1];
  double rdvpar2 = 2.0/dxv[1];
  double wmu = w[2];
  double rdmu2 = 2.0/dxv[2];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wvparSq = w[1]*w[1];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[2]*w[2];
  double rdmu2Sq = rdmu2*rdmu2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[3];
  const double *b_z = &b_i[6];

  const double *BstarZdBmag = &Bstar_Bmag[0]; 

  double hamil[20] = {0.}; 
  hamil[0] = 1.414213562373095*m_*wvparSq+2.0*bmag[0]*wmu+(0.4714045207910317*m_)/rdvpar2Sq+2.0*phi[0]*q_; 
  hamil[1] = 2.0*bmag[1]*wmu+2.0*phi[1]*q_; 
  hamil[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[7] = 2.0*bmag[2]*wmu+2.0*phi[2]*q_; 
  hamil[8] = (0.421637021355784*m_)/rdvpar2Sq; 
  hamil[13] = (1.154700538379251*bmag[2])/rdmu2; 

  double *alphaL0 = &alpha_surf[0];
  alphaL0[0] = ((1.936491673103709*BstarZdBmag[2]-1.5*BstarZdBmag[1]+0.8660254037844386*BstarZdBmag[0])*hamil[2]*rdvpar2)/m_; 
  alphaL0[1] = ((4.330127018922193*BstarZdBmag[2]-3.354101966249685*BstarZdBmag[1]+1.936491673103709*BstarZdBmag[0])*hamil[8]*rdvpar2)/m_; 

  double *sgn_alpha_surf0 = &sgn_alpha_surf[0];
  int const_sgn_alpha_surf = 1; 
  
  if (0.5*alphaL0[0]-0.6708203932499357*alphaL0[1] > 0.) 
    sgn_alpha_surf0[0] = 1.0; 
  else  
    sgn_alpha_surf0[0] = -1.0; 
  
  if (0.5*alphaL0[0]-0.6708203932499357*alphaL0[1] > 0.) 
    sgn_alpha_surf0[1] = 1.0; 
  else  
    sgn_alpha_surf0[1] = -1.0; 
  
  if (sgn_alpha_surf0[1] == sgn_alpha_surf0[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5*alphaL0[0]-0.6708203932499357*alphaL0[1] > 0.) 
    sgn_alpha_surf0[2] = 1.0; 
  else  
    sgn_alpha_surf0[2] = -1.0; 
  
  if (sgn_alpha_surf0[2] == sgn_alpha_surf0[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5*alphaL0[0] > 0.) 
    sgn_alpha_surf0[3] = 1.0; 
  else  
    sgn_alpha_surf0[3] = -1.0; 
  
  if (sgn_alpha_surf0[3] == sgn_alpha_surf0[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5*alphaL0[0] > 0.) 
    sgn_alpha_surf0[4] = 1.0; 
  else  
    sgn_alpha_surf0[4] = -1.0; 
  
  if (sgn_alpha_surf0[4] == sgn_alpha_surf0[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5*alphaL0[0] > 0.) 
    sgn_alpha_surf0[5] = 1.0; 
  else  
    sgn_alpha_surf0[5] = -1.0; 
  
  if (sgn_alpha_surf0[5] == sgn_alpha_surf0[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.6708203932499357*alphaL0[1]+0.5*alphaL0[0] > 0.) 
    sgn_alpha_surf0[6] = 1.0; 
  else  
    sgn_alpha_surf0[6] = -1.0; 
  
  if (sgn_alpha_surf0[6] == sgn_alpha_surf0[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.6708203932499357*alphaL0[1]+0.5*alphaL0[0] > 0.) 
    sgn_alpha_surf0[7] = 1.0; 
  else  
    sgn_alpha_surf0[7] = -1.0; 
  
  if (sgn_alpha_surf0[7] == sgn_alpha_surf0[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.6708203932499357*alphaL0[1]+0.5*alphaL0[0] > 0.) 
    sgn_alpha_surf0[8] = 1.0; 
  else  
    sgn_alpha_surf0[8] = -1.0; 
  
  if (sgn_alpha_surf0[8] == sgn_alpha_surf0[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
