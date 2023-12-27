#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double q_, const double m_, 
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
  double wmu = w[2];
  double rdmu2 = 2.0/dxv[2];

  double wvparSq = wvpar*wvpar;
  double rdvpar2Sq = rdvpar2*rdvpar2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[2];
  const double *b_z = &b_i[4];

  double hamil[12] = {0.}; 
  hamil[0] = 1.414213562373095*m_*wvparSq+2.0*bmag[0]*wmu+(0.4714045207910317*m_)/rdvpar2Sq+2.0*phi[0]*q_; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[8] = (0.421637021355784*m_)/rdvpar2Sq; 

  double *alphaR0 = &alpha_surf[0];
  alphaR0[0] = ((1.837117307087383*cmag[1]*jacobtot_inv[1]*hamil[2]+1.060660171779821*cmag[0]*jacobtot_inv[1]*hamil[2]+1.060660171779821*jacobtot_inv[0]*cmag[1]*hamil[2]+0.6123724356957944*cmag[0]*jacobtot_inv[0]*hamil[2])*rdvpar2)/m_; 
  alphaR0[1] = ((4.107919181288745*cmag[1]*jacobtot_inv[1]*hamil[8]+2.371708245126284*cmag[0]*jacobtot_inv[1]*hamil[8]+2.371708245126284*jacobtot_inv[0]*cmag[1]*hamil[8]+1.369306393762915*cmag[0]*jacobtot_inv[0]*hamil[8])*rdvpar2)/m_; 

  double *sgn_alpha_surf0 = &sgn_alpha_surf[0];
  int const_sgn_alpha_surf = 1; 
  
  if (0.5*alphaR0[0]-0.6708203932499357*alphaR0[1] > 0.) 
    sgn_alpha_surf0[0] = 1.0; 
  else  
    sgn_alpha_surf0[0] = -1.0; 
  
  if (0.5*alphaR0[0] > 0.) 
    sgn_alpha_surf0[1] = 1.0; 
  else  
    sgn_alpha_surf0[1] = -1.0; 
  
  if (sgn_alpha_surf0[1] == sgn_alpha_surf0[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.6708203932499357*alphaR0[1]+0.5*alphaR0[0] > 0.) 
    sgn_alpha_surf0[2] = 1.0; 
  else  
    sgn_alpha_surf0[2] = -1.0; 
  
  if (sgn_alpha_surf0[2] == sgn_alpha_surf0[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5*alphaR0[0]-0.6708203932499357*alphaR0[1] > 0.) 
    sgn_alpha_surf0[3] = 1.0; 
  else  
    sgn_alpha_surf0[3] = -1.0; 
  
  if (sgn_alpha_surf0[3] == sgn_alpha_surf0[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5*alphaR0[0] > 0.) 
    sgn_alpha_surf0[4] = 1.0; 
  else  
    sgn_alpha_surf0[4] = -1.0; 
  
  if (sgn_alpha_surf0[4] == sgn_alpha_surf0[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.6708203932499357*alphaR0[1]+0.5*alphaR0[0] > 0.) 
    sgn_alpha_surf0[5] = 1.0; 
  else  
    sgn_alpha_surf0[5] = -1.0; 
  
  if (sgn_alpha_surf0[5] == sgn_alpha_surf0[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
