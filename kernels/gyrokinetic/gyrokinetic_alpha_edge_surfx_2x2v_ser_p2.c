#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfx_2x2v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, 
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
  double wvpar = w[2];
  double rdvpar2 = 2.0/dxv[2];
  double wmu = w[3];
  double rdmu2 = 2.0/dxv[3];

  double wvparSq = wvpar*wvpar;
  double rdvpar2Sq = rdvpar2*rdvpar2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[8];
  const double *b_z = &b_i[16];

  double hamil[48] = {0.}; 
  hamil[0] = 2.0*(m_*wvparSq+bmag[0]*wmu)+(0.6666666666666666*m_)/rdvpar2Sq+2.0*phi[0]*q_; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[11] = 2.0*(bmag[4]*wmu+phi[4]*q_); 
  hamil[12] = 2.0*phi[5]*q_; 
  hamil[13] = (0.5962847939999438*m_)/rdvpar2Sq; 
  hamil[19] = 2.0*phi[6]*q_; 
  hamil[20] = 2.0*phi[7]*q_; 
  hamil[25] = (1.154700538379251*bmag[4])/rdmu2; 

  double *alphaR0 = &alpha_surf[0];
  alphaR0[0] = (((-3.423265984407287*b_z[4]*jacobtot_inv[4]*hamil[19])-2.651650429449552*b_z[1]*jacobtot_inv[4]*hamil[19]-1.530931089239486*b_z[0]*jacobtot_inv[4]*hamil[19]-2.651650429449552*jacobtot_inv[1]*b_z[4]*hamil[19]-1.530931089239486*jacobtot_inv[0]*b_z[4]*hamil[19]-2.053959590644372*b_z[1]*jacobtot_inv[1]*hamil[19]-1.185854122563142*b_z[0]*jacobtot_inv[1]*hamil[19]-1.185854122563142*jacobtot_inv[0]*b_z[1]*hamil[19]-0.6846531968814574*b_z[0]*jacobtot_inv[0]*hamil[19]-2.651650429449552*b_z[4]*jacobtot_inv[4]*hamil[5]-2.053959590644372*b_z[1]*jacobtot_inv[4]*hamil[5]-1.185854122563142*b_z[0]*jacobtot_inv[4]*hamil[5]-2.053959590644372*jacobtot_inv[1]*b_z[4]*hamil[5]-1.185854122563142*jacobtot_inv[0]*b_z[4]*hamil[5]-1.590990257669731*b_z[1]*jacobtot_inv[1]*hamil[5]-0.9185586535436913*b_z[0]*jacobtot_inv[1]*hamil[5]-0.9185586535436913*jacobtot_inv[0]*b_z[1]*hamil[5]-0.5303300858899105*b_z[0]*jacobtot_inv[0]*hamil[5]-1.530931089239486*hamil[2]*b_z[4]*jacobtot_inv[4]-1.185854122563142*b_z[1]*hamil[2]*jacobtot_inv[4]-0.6846531968814573*b_z[0]*hamil[2]*jacobtot_inv[4]-1.185854122563142*jacobtot_inv[1]*hamil[2]*b_z[4]-0.6846531968814573*jacobtot_inv[0]*hamil[2]*b_z[4]-0.9185586535436913*b_z[1]*jacobtot_inv[1]*hamil[2]-0.5303300858899105*b_z[0]*jacobtot_inv[1]*hamil[2]-0.5303300858899105*jacobtot_inv[0]*b_z[1]*hamil[2]-0.3061862178478971*b_z[0]*jacobtot_inv[0]*hamil[2])*rdy2)/q_; 
  alphaR0[1] = (((-5.929270612815709*b_z[4]*jacobtot_inv[4]*hamil[20])-4.592793267718458*b_z[1]*jacobtot_inv[4]*hamil[20]-2.651650429449552*b_z[0]*jacobtot_inv[4]*hamil[20]-4.592793267718458*jacobtot_inv[1]*b_z[4]*hamil[20]-2.651650429449552*jacobtot_inv[0]*b_z[4]*hamil[20]-3.557562367689425*b_z[1]*jacobtot_inv[1]*hamil[20]-2.053959590644372*b_z[0]*jacobtot_inv[1]*hamil[20]-2.053959590644372*jacobtot_inv[0]*b_z[1]*hamil[20]-1.185854122563142*b_z[0]*jacobtot_inv[0]*hamil[20]-3.423265984407287*b_z[4]*jacobtot_inv[4]*hamil[12]-2.651650429449552*b_z[1]*jacobtot_inv[4]*hamil[12]-1.530931089239486*b_z[0]*jacobtot_inv[4]*hamil[12]-2.651650429449552*jacobtot_inv[1]*b_z[4]*hamil[12]-1.530931089239486*jacobtot_inv[0]*b_z[4]*hamil[12]-2.053959590644372*b_z[1]*jacobtot_inv[1]*hamil[12]-1.185854122563142*b_z[0]*jacobtot_inv[1]*hamil[12]-1.185854122563142*jacobtot_inv[0]*b_z[1]*hamil[12]-0.6846531968814573*b_z[0]*jacobtot_inv[0]*hamil[12])*rdy2)/q_; 

  double *sgn_alpha_surf0 = &sgn_alpha_surf[0];
  int const_sgn_alpha_surf = 1; 
  
  if (0.3535533905932734*alphaR0[0]-0.4743416490252568*alphaR0[1] > 0.) 
    sgn_alpha_surf0[0] = 1.0; 
  else  
    sgn_alpha_surf0[0] = -1.0; 
  
  if (0.3535533905932734*alphaR0[0]-0.4743416490252568*alphaR0[1] > 0.) 
    sgn_alpha_surf0[1] = 1.0; 
  else  
    sgn_alpha_surf0[1] = -1.0; 
  
  if (sgn_alpha_surf0[1] == sgn_alpha_surf0[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0]-0.4743416490252568*alphaR0[1] > 0.) 
    sgn_alpha_surf0[2] = 1.0; 
  else  
    sgn_alpha_surf0[2] = -1.0; 
  
  if (sgn_alpha_surf0[2] == sgn_alpha_surf0[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0]-0.4743416490252568*alphaR0[1] > 0.) 
    sgn_alpha_surf0[3] = 1.0; 
  else  
    sgn_alpha_surf0[3] = -1.0; 
  
  if (sgn_alpha_surf0[3] == sgn_alpha_surf0[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0]-0.4743416490252568*alphaR0[1] > 0.) 
    sgn_alpha_surf0[4] = 1.0; 
  else  
    sgn_alpha_surf0[4] = -1.0; 
  
  if (sgn_alpha_surf0[4] == sgn_alpha_surf0[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0]-0.4743416490252568*alphaR0[1] > 0.) 
    sgn_alpha_surf0[5] = 1.0; 
  else  
    sgn_alpha_surf0[5] = -1.0; 
  
  if (sgn_alpha_surf0[5] == sgn_alpha_surf0[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0]-0.4743416490252568*alphaR0[1] > 0.) 
    sgn_alpha_surf0[6] = 1.0; 
  else  
    sgn_alpha_surf0[6] = -1.0; 
  
  if (sgn_alpha_surf0[6] == sgn_alpha_surf0[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0]-0.4743416490252568*alphaR0[1] > 0.) 
    sgn_alpha_surf0[7] = 1.0; 
  else  
    sgn_alpha_surf0[7] = -1.0; 
  
  if (sgn_alpha_surf0[7] == sgn_alpha_surf0[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0]-0.4743416490252568*alphaR0[1] > 0.) 
    sgn_alpha_surf0[8] = 1.0; 
  else  
    sgn_alpha_surf0[8] = -1.0; 
  
  if (sgn_alpha_surf0[8] == sgn_alpha_surf0[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[9] = 1.0; 
  else  
    sgn_alpha_surf0[9] = -1.0; 
  
  if (sgn_alpha_surf0[9] == sgn_alpha_surf0[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[10] = 1.0; 
  else  
    sgn_alpha_surf0[10] = -1.0; 
  
  if (sgn_alpha_surf0[10] == sgn_alpha_surf0[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[11] = 1.0; 
  else  
    sgn_alpha_surf0[11] = -1.0; 
  
  if (sgn_alpha_surf0[11] == sgn_alpha_surf0[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[12] = 1.0; 
  else  
    sgn_alpha_surf0[12] = -1.0; 
  
  if (sgn_alpha_surf0[12] == sgn_alpha_surf0[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[13] = 1.0; 
  else  
    sgn_alpha_surf0[13] = -1.0; 
  
  if (sgn_alpha_surf0[13] == sgn_alpha_surf0[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[14] = 1.0; 
  else  
    sgn_alpha_surf0[14] = -1.0; 
  
  if (sgn_alpha_surf0[14] == sgn_alpha_surf0[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[15] = 1.0; 
  else  
    sgn_alpha_surf0[15] = -1.0; 
  
  if (sgn_alpha_surf0[15] == sgn_alpha_surf0[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[16] = 1.0; 
  else  
    sgn_alpha_surf0[16] = -1.0; 
  
  if (sgn_alpha_surf0[16] == sgn_alpha_surf0[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[17] = 1.0; 
  else  
    sgn_alpha_surf0[17] = -1.0; 
  
  if (sgn_alpha_surf0[17] == sgn_alpha_surf0[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*alphaR0[1]+0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[18] = 1.0; 
  else  
    sgn_alpha_surf0[18] = -1.0; 
  
  if (sgn_alpha_surf0[18] == sgn_alpha_surf0[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*alphaR0[1]+0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[19] = 1.0; 
  else  
    sgn_alpha_surf0[19] = -1.0; 
  
  if (sgn_alpha_surf0[19] == sgn_alpha_surf0[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*alphaR0[1]+0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[20] = 1.0; 
  else  
    sgn_alpha_surf0[20] = -1.0; 
  
  if (sgn_alpha_surf0[20] == sgn_alpha_surf0[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*alphaR0[1]+0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[21] = 1.0; 
  else  
    sgn_alpha_surf0[21] = -1.0; 
  
  if (sgn_alpha_surf0[21] == sgn_alpha_surf0[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*alphaR0[1]+0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[22] = 1.0; 
  else  
    sgn_alpha_surf0[22] = -1.0; 
  
  if (sgn_alpha_surf0[22] == sgn_alpha_surf0[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*alphaR0[1]+0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[23] = 1.0; 
  else  
    sgn_alpha_surf0[23] = -1.0; 
  
  if (sgn_alpha_surf0[23] == sgn_alpha_surf0[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*alphaR0[1]+0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[24] = 1.0; 
  else  
    sgn_alpha_surf0[24] = -1.0; 
  
  if (sgn_alpha_surf0[24] == sgn_alpha_surf0[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*alphaR0[1]+0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[25] = 1.0; 
  else  
    sgn_alpha_surf0[25] = -1.0; 
  
  if (sgn_alpha_surf0[25] == sgn_alpha_surf0[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*alphaR0[1]+0.3535533905932734*alphaR0[0] > 0.) 
    sgn_alpha_surf0[26] = 1.0; 
  else  
    sgn_alpha_surf0[26] = -1.0; 
  
  if (sgn_alpha_surf0[26] == sgn_alpha_surf0[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
