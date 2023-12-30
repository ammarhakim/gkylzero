#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_surfvpar_1x2v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, 
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
  const double *b_y = &b_i[3];
  const double *b_z = &b_i[6];

  double hamil[20] = {0.}; 
  hamil[0] = 1.414213562373095*m_*wvparSq+2.0*bmag[0]*wmu+(0.4714045207910317*m_)/rdvpar2Sq+2.0*phi[0]*q_; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[7] = 2.0*(bmag[2]*wmu+phi[2]*q_); 
  hamil[8] = (0.421637021355784*m_)/rdvpar2Sq; 
  hamil[13] = (1.154700538379251*bmag[2])/rdmu2; 

  double *alphaL = &alpha_surf[8];
  double *sgn_alpha_surfL = &sgn_alpha_surf[9];
  alphaL[0] = (((-1.224744871391589*cmag[1]*jacobtot_inv[2]*hamil[7])-1.224744871391589*jacobtot_inv[1]*cmag[2]*hamil[7]-1.369306393762915*cmag[0]*jacobtot_inv[1]*hamil[7]-1.369306393762915*jacobtot_inv[0]*cmag[1]*hamil[7]-0.6123724356957944*hamil[1]*cmag[2]*jacobtot_inv[2]-0.6123724356957944*cmag[1]*hamil[1]*jacobtot_inv[1]-0.6123724356957944*cmag[0]*jacobtot_inv[0]*hamil[1])*rdx2)/m_; 
  alphaL[1] = (((-2.151767190198866*cmag[2]*jacobtot_inv[2]*hamil[7])-1.224744871391589*cmag[0]*jacobtot_inv[2]*hamil[7]-1.224744871391589*jacobtot_inv[0]*cmag[2]*hamil[7]-2.464751508773246*cmag[1]*jacobtot_inv[1]*hamil[7]-1.369306393762915*cmag[0]*jacobtot_inv[0]*hamil[7]-0.5477225575051661*cmag[1]*hamil[1]*jacobtot_inv[2]-0.5477225575051661*hamil[1]*jacobtot_inv[1]*cmag[2]-0.6123724356957944*cmag[0]*hamil[1]*jacobtot_inv[1]-0.6123724356957944*jacobtot_inv[0]*cmag[1]*hamil[1])*rdx2)/m_; 
  alphaL[2] = (((-1.224744871391589*cmag[1]*jacobtot_inv[2]*hamil[13])-1.224744871391589*jacobtot_inv[1]*cmag[2]*hamil[13]-1.369306393762915*cmag[0]*jacobtot_inv[1]*hamil[13]-1.369306393762915*jacobtot_inv[0]*cmag[1]*hamil[13]-0.6123724356957944*cmag[2]*jacobtot_inv[2]*hamil[5]-0.6123724356957944*cmag[1]*jacobtot_inv[1]*hamil[5]-0.6123724356957944*cmag[0]*jacobtot_inv[0]*hamil[5])*rdx2)/m_; 
  alphaL[3] = (((-2.151767190198866*cmag[2]*jacobtot_inv[2]*hamil[13])-1.224744871391589*cmag[0]*jacobtot_inv[2]*hamil[13]-1.224744871391589*jacobtot_inv[0]*cmag[2]*hamil[13]-2.464751508773247*cmag[1]*jacobtot_inv[1]*hamil[13]-1.369306393762915*cmag[0]*jacobtot_inv[0]*hamil[13]-0.5477225575051661*cmag[1]*jacobtot_inv[2]*hamil[5]-0.5477225575051661*jacobtot_inv[1]*cmag[2]*hamil[5]-0.6123724356957944*cmag[0]*jacobtot_inv[1]*hamil[5]-0.6123724356957944*jacobtot_inv[0]*cmag[1]*hamil[5])*rdx2)/m_; 
  alphaL[4] = (((-2.151767190198866*cmag[1]*jacobtot_inv[2]*hamil[7])-2.151767190198866*jacobtot_inv[1]*cmag[2]*hamil[7]-1.224744871391589*cmag[0]*jacobtot_inv[1]*hamil[7]-1.224744871391589*jacobtot_inv[0]*cmag[1]*hamil[7]-0.3912303982179757*hamil[1]*cmag[2]*jacobtot_inv[2]-0.6123724356957944*cmag[0]*hamil[1]*jacobtot_inv[2]-0.6123724356957944*jacobtot_inv[0]*hamil[1]*cmag[2]-0.5477225575051661*cmag[1]*hamil[1]*jacobtot_inv[1])*rdx2)/m_; 
  alphaL[6] = (((-2.151767190198866*cmag[1]*jacobtot_inv[2]*hamil[13])-2.151767190198866*jacobtot_inv[1]*cmag[2]*hamil[13]-1.224744871391589*cmag[0]*jacobtot_inv[1]*hamil[13]-1.224744871391589*jacobtot_inv[0]*cmag[1]*hamil[13]-0.3912303982179757*cmag[2]*jacobtot_inv[2]*hamil[5]-0.6123724356957944*cmag[0]*jacobtot_inv[2]*hamil[5]-0.6123724356957944*jacobtot_inv[0]*cmag[2]*hamil[5]-0.5477225575051661*cmag[1]*jacobtot_inv[1]*hamil[5])*rdx2)/m_; 

  int const_sgn_alpha_surf = 1;  
  
  if ((-0.6*alphaL[6])+0.4472135954999572*alphaL[4]+0.9*alphaL[3]-0.6708203932499357*(alphaL[2]+alphaL[1])+0.5*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if (0.4472135954999572*alphaL[4]-0.6708203932499357*alphaL[1]+0.5*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.6*alphaL[6]+0.4472135954999572*alphaL[4]-0.9*alphaL[3]+0.6708203932499357*alphaL[2]-0.6708203932499357*alphaL[1]+0.5*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.75*alphaL[6]-0.5590169943749465*alphaL[4]-0.6708203932499357*alphaL[2]+0.5*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5*alphaL[0]-0.5590169943749465*alphaL[4] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.75*alphaL[6])-0.5590169943749465*alphaL[4]+0.6708203932499357*alphaL[2]+0.5*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.6*alphaL[6])+0.4472135954999572*alphaL[4]-0.9*alphaL[3]-0.6708203932499357*alphaL[2]+0.6708203932499357*alphaL[1]+0.5*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4472135954999572*alphaL[4]+0.6708203932499357*alphaL[1]+0.5*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.6*alphaL[6]+0.4472135954999572*alphaL[4]+0.9*alphaL[3]+0.6708203932499357*(alphaL[2]+alphaL[1])+0.5*alphaL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
