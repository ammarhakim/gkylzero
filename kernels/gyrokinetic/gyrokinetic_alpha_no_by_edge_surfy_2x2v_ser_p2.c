#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, 
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
  double wz = w[1];
  double rdz2 = 2.0/dxv[1];
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
  hamil[2] = 2.0*(bmag[2]*wmu+phi[2]*q_); 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*(bmag[3]*wmu+phi[3]*q_); 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[9] = (1.154700538379252*bmag[2])/rdmu2; 
  hamil[11] = 2.0*(bmag[4]*wmu+phi[4]*q_); 
  hamil[12] = 2.0*(bmag[5]*wmu+phi[5]*q_); 
  hamil[13] = (0.5962847939999438*m_)/rdvpar2Sq; 
  hamil[16] = (1.154700538379252*bmag[3])/rdmu2; 
  hamil[19] = 2.0*(bmag[6]*wmu+phi[6]*q_); 
  hamil[20] = 2.0*(bmag[7]*wmu+phi[7]*q_); 
  hamil[25] = (1.154700538379251*bmag[4])/rdmu2; 
  hamil[26] = (1.154700538379251*bmag[5])/rdmu2; 
  hamil[35] = (1.154700538379251*bmag[6])/rdmu2; 
  hamil[36] = (1.154700538379251*bmag[7])/rdmu2; 

  double *alphaR = &alpha_surf[20];
  double *sgn_alpha_surfR = &sgn_alpha_surf[27];
  alphaR[0] = ((1.530931089239486*hamil[3]*cmag[7]*jacobtot_inv[7]+1.185854122563142*cmag[3]*hamil[3]*jacobtot_inv[7]+0.6846531968814574*cmag[1]*hamil[3]*jacobtot_inv[7]+1.185854122563142*hamil[3]*jacobtot_inv[3]*cmag[7]+0.6846531968814574*jacobtot_inv[1]*hamil[3]*cmag[7]+0.9185586535436913*hamil[3]*cmag[6]*jacobtot_inv[6]+0.5303300858899104*hamil[3]*cmag[4]*jacobtot_inv[6]+0.5303300858899104*hamil[3]*jacobtot_inv[4]*cmag[6]+1.530931089239486*hamil[3]*cmag[5]*jacobtot_inv[5]+1.185854122563142*cmag[2]*hamil[3]*jacobtot_inv[5]+0.6846531968814573*cmag[0]*hamil[3]*jacobtot_inv[5]+1.185854122563142*jacobtot_inv[2]*hamil[3]*cmag[5]+0.6846531968814573*jacobtot_inv[0]*hamil[3]*cmag[5]+0.3061862178478971*hamil[3]*cmag[4]*jacobtot_inv[4]+0.9185586535436913*cmag[3]*hamil[3]*jacobtot_inv[3]+0.5303300858899105*cmag[1]*hamil[3]*jacobtot_inv[3]+0.5303300858899105*jacobtot_inv[1]*cmag[3]*hamil[3]+0.9185586535436913*cmag[2]*jacobtot_inv[2]*hamil[3]+0.5303300858899105*cmag[0]*jacobtot_inv[2]*hamil[3]+0.5303300858899105*jacobtot_inv[0]*cmag[2]*hamil[3]+0.3061862178478971*cmag[1]*jacobtot_inv[1]*hamil[3]+0.3061862178478971*cmag[0]*jacobtot_inv[0]*hamil[3])*rdvpar2)/m_; 
  alphaR[1] = ((1.060660171779821*hamil[3]*cmag[6]*jacobtot_inv[7]+1.530931089239486*hamil[3]*cmag[5]*jacobtot_inv[7]+0.6123724356957944*hamil[3]*cmag[4]*jacobtot_inv[7]+1.185854122563142*cmag[2]*hamil[3]*jacobtot_inv[7]+0.6846531968814574*cmag[0]*hamil[3]*jacobtot_inv[7]+1.060660171779821*hamil[3]*jacobtot_inv[6]*cmag[7]+1.530931089239486*hamil[3]*jacobtot_inv[5]*cmag[7]+0.6123724356957944*hamil[3]*jacobtot_inv[4]*cmag[7]+1.185854122563142*jacobtot_inv[2]*hamil[3]*cmag[7]+0.6846531968814574*jacobtot_inv[0]*hamil[3]*cmag[7]+0.821583836257749*cmag[3]*hamil[3]*jacobtot_inv[6]+0.4743416490252568*cmag[1]*hamil[3]*jacobtot_inv[6]+0.821583836257749*hamil[3]*jacobtot_inv[3]*cmag[6]+0.4743416490252568*jacobtot_inv[1]*hamil[3]*cmag[6]+1.185854122563142*cmag[3]*hamil[3]*jacobtot_inv[5]+0.6846531968814573*cmag[1]*hamil[3]*jacobtot_inv[5]+1.185854122563142*hamil[3]*jacobtot_inv[3]*cmag[5]+0.6846531968814573*jacobtot_inv[1]*hamil[3]*cmag[5]+0.4743416490252568*cmag[3]*hamil[3]*jacobtot_inv[4]+0.273861278752583*cmag[1]*hamil[3]*jacobtot_inv[4]+0.4743416490252568*hamil[3]*jacobtot_inv[3]*cmag[4]+0.273861278752583*jacobtot_inv[1]*hamil[3]*cmag[4]+0.9185586535436913*cmag[2]*hamil[3]*jacobtot_inv[3]+0.5303300858899105*cmag[0]*hamil[3]*jacobtot_inv[3]+0.9185586535436913*jacobtot_inv[2]*cmag[3]*hamil[3]+0.5303300858899105*jacobtot_inv[0]*cmag[3]*hamil[3]+0.5303300858899105*cmag[1]*jacobtot_inv[2]*hamil[3]+0.5303300858899105*jacobtot_inv[1]*cmag[2]*hamil[3]+0.3061862178478971*cmag[0]*jacobtot_inv[1]*hamil[3]+0.3061862178478971*jacobtot_inv[0]*cmag[1]*hamil[3])*rdvpar2)/m_; 
  alphaR[2] = ((3.423265984407287*cmag[7]*jacobtot_inv[7]*hamil[13]+2.651650429449552*cmag[3]*jacobtot_inv[7]*hamil[13]+1.530931089239486*cmag[1]*jacobtot_inv[7]*hamil[13]+2.651650429449552*jacobtot_inv[3]*cmag[7]*hamil[13]+1.530931089239486*jacobtot_inv[1]*cmag[7]*hamil[13]+2.053959590644372*cmag[6]*jacobtot_inv[6]*hamil[13]+1.185854122563142*cmag[4]*jacobtot_inv[6]*hamil[13]+1.185854122563142*jacobtot_inv[4]*cmag[6]*hamil[13]+3.423265984407287*cmag[5]*jacobtot_inv[5]*hamil[13]+2.651650429449552*cmag[2]*jacobtot_inv[5]*hamil[13]+1.530931089239486*cmag[0]*jacobtot_inv[5]*hamil[13]+2.651650429449552*jacobtot_inv[2]*cmag[5]*hamil[13]+1.530931089239486*jacobtot_inv[0]*cmag[5]*hamil[13]+0.6846531968814573*cmag[4]*jacobtot_inv[4]*hamil[13]+2.053959590644372*cmag[3]*jacobtot_inv[3]*hamil[13]+1.185854122563142*cmag[1]*jacobtot_inv[3]*hamil[13]+1.185854122563142*jacobtot_inv[1]*cmag[3]*hamil[13]+2.053959590644372*cmag[2]*jacobtot_inv[2]*hamil[13]+1.185854122563142*cmag[0]*jacobtot_inv[2]*hamil[13]+1.185854122563142*jacobtot_inv[0]*cmag[2]*hamil[13]+0.6846531968814573*cmag[1]*jacobtot_inv[1]*hamil[13]+0.6846531968814573*cmag[0]*jacobtot_inv[0]*hamil[13])*rdvpar2)/m_; 
  alphaR[4] = ((2.371708245126284*cmag[6]*jacobtot_inv[7]*hamil[13]+3.423265984407287*cmag[5]*jacobtot_inv[7]*hamil[13]+1.369306393762915*cmag[4]*jacobtot_inv[7]*hamil[13]+2.651650429449552*cmag[2]*jacobtot_inv[7]*hamil[13]+1.530931089239486*cmag[0]*jacobtot_inv[7]*hamil[13]+2.371708245126284*jacobtot_inv[6]*cmag[7]*hamil[13]+3.423265984407287*jacobtot_inv[5]*cmag[7]*hamil[13]+1.369306393762915*jacobtot_inv[4]*cmag[7]*hamil[13]+2.651650429449552*jacobtot_inv[2]*cmag[7]*hamil[13]+1.530931089239486*jacobtot_inv[0]*cmag[7]*hamil[13]+1.837117307087383*cmag[3]*jacobtot_inv[6]*hamil[13]+1.060660171779821*cmag[1]*jacobtot_inv[6]*hamil[13]+1.837117307087383*jacobtot_inv[3]*cmag[6]*hamil[13]+1.060660171779821*jacobtot_inv[1]*cmag[6]*hamil[13]+2.651650429449552*cmag[3]*jacobtot_inv[5]*hamil[13]+1.530931089239486*cmag[1]*jacobtot_inv[5]*hamil[13]+2.651650429449552*jacobtot_inv[3]*cmag[5]*hamil[13]+1.530931089239486*jacobtot_inv[1]*cmag[5]*hamil[13]+1.060660171779821*cmag[3]*jacobtot_inv[4]*hamil[13]+0.6123724356957944*cmag[1]*jacobtot_inv[4]*hamil[13]+1.060660171779821*jacobtot_inv[3]*cmag[4]*hamil[13]+0.6123724356957944*jacobtot_inv[1]*cmag[4]*hamil[13]+2.053959590644372*cmag[2]*jacobtot_inv[3]*hamil[13]+1.185854122563142*cmag[0]*jacobtot_inv[3]*hamil[13]+2.053959590644372*jacobtot_inv[2]*cmag[3]*hamil[13]+1.185854122563142*jacobtot_inv[0]*cmag[3]*hamil[13]+1.185854122563142*cmag[1]*jacobtot_inv[2]*hamil[13]+1.185854122563142*jacobtot_inv[1]*cmag[2]*hamil[13]+0.6846531968814573*cmag[0]*jacobtot_inv[1]*hamil[13]+0.6846531968814573*jacobtot_inv[0]*cmag[1]*hamil[13])*rdvpar2)/m_; 
  alphaR[7] = ((1.369306393762915*hamil[3]*cmag[7]*jacobtot_inv[7]+1.060660171779821*cmag[3]*hamil[3]*jacobtot_inv[7]+0.6123724356957944*cmag[1]*hamil[3]*jacobtot_inv[7]+1.060660171779821*hamil[3]*jacobtot_inv[3]*cmag[7]+0.6123724356957944*jacobtot_inv[1]*hamil[3]*cmag[7]+0.5868455973269634*hamil[3]*cmag[6]*jacobtot_inv[6]+1.185854122563142*hamil[3]*cmag[5]*jacobtot_inv[6]+0.3388154635894691*hamil[3]*cmag[4]*jacobtot_inv[6]+0.9185586535436916*cmag[2]*hamil[3]*jacobtot_inv[6]+0.5303300858899104*cmag[0]*hamil[3]*jacobtot_inv[6]+1.185854122563142*hamil[3]*jacobtot_inv[5]*cmag[6]+0.3388154635894691*hamil[3]*jacobtot_inv[4]*cmag[6]+0.9185586535436916*jacobtot_inv[2]*hamil[3]*cmag[6]+0.5303300858899104*jacobtot_inv[0]*hamil[3]*cmag[6]+0.6846531968814573*hamil[3]*cmag[4]*jacobtot_inv[5]+0.6846531968814573*hamil[3]*jacobtot_inv[4]*cmag[5]+0.1956151991089878*hamil[3]*cmag[4]*jacobtot_inv[4]+0.5303300858899105*cmag[2]*hamil[3]*jacobtot_inv[4]+0.3061862178478971*cmag[0]*hamil[3]*jacobtot_inv[4]+0.5303300858899105*jacobtot_inv[2]*hamil[3]*cmag[4]+0.3061862178478971*jacobtot_inv[0]*hamil[3]*cmag[4]+0.8215838362577489*cmag[3]*hamil[3]*jacobtot_inv[3]+0.4743416490252568*cmag[1]*hamil[3]*jacobtot_inv[3]+0.4743416490252568*jacobtot_inv[1]*cmag[3]*hamil[3]+0.273861278752583*cmag[1]*jacobtot_inv[1]*hamil[3])*rdvpar2)/m_; 
  alphaR[11] = ((3.061862178478972*cmag[7]*jacobtot_inv[7]*hamil[13]+2.371708245126284*cmag[3]*jacobtot_inv[7]*hamil[13]+1.369306393762915*cmag[1]*jacobtot_inv[7]*hamil[13]+2.371708245126284*jacobtot_inv[3]*cmag[7]*hamil[13]+1.369306393762915*jacobtot_inv[1]*cmag[7]*hamil[13]+1.312226647919559*cmag[6]*jacobtot_inv[6]*hamil[13]+2.651650429449552*cmag[5]*jacobtot_inv[6]*hamil[13]+0.7576144084141578*cmag[4]*jacobtot_inv[6]*hamil[13]+2.053959590644372*cmag[2]*jacobtot_inv[6]*hamil[13]+1.185854122563142*cmag[0]*jacobtot_inv[6]*hamil[13]+2.651650429449552*jacobtot_inv[5]*cmag[6]*hamil[13]+0.7576144084141578*jacobtot_inv[4]*cmag[6]*hamil[13]+2.053959590644372*jacobtot_inv[2]*cmag[6]*hamil[13]+1.185854122563142*jacobtot_inv[0]*cmag[6]*hamil[13]+1.530931089239486*cmag[4]*jacobtot_inv[5]*hamil[13]+1.530931089239486*jacobtot_inv[4]*cmag[5]*hamil[13]+0.4374088826398531*cmag[4]*jacobtot_inv[4]*hamil[13]+1.185854122563142*cmag[2]*jacobtot_inv[4]*hamil[13]+0.6846531968814574*cmag[0]*jacobtot_inv[4]*hamil[13]+1.185854122563142*jacobtot_inv[2]*cmag[4]*hamil[13]+0.6846531968814574*jacobtot_inv[0]*cmag[4]*hamil[13]+1.837117307087383*cmag[3]*jacobtot_inv[3]*hamil[13]+1.060660171779821*cmag[1]*jacobtot_inv[3]*hamil[13]+1.060660171779821*jacobtot_inv[1]*cmag[3]*hamil[13]+0.6123724356957944*cmag[1]*jacobtot_inv[1]*hamil[13])*rdvpar2)/m_; 

  int const_sgn_alpha_surf = 1;  
  
  if ((-0.4242640687119286*alphaR[11])+0.3162277660168378*alphaR[7]+0.6363961030678927*alphaR[4]-0.4743416490252568*(alphaR[2]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[0] = 1.0; 
  else  
    sgn_alpha_surfR[0] = -1.0; 
  
  if ((-0.4242640687119286*alphaR[11])+0.3162277660168378*alphaR[7]+0.6363961030678927*alphaR[4]-0.4743416490252568*(alphaR[2]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[1] = 1.0; 
  else  
    sgn_alpha_surfR[1] = -1.0; 
  
  if (sgn_alpha_surfR[1] == sgn_alpha_surfR[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119286*alphaR[11])+0.3162277660168378*alphaR[7]+0.6363961030678927*alphaR[4]-0.4743416490252568*(alphaR[2]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[2] = 1.0; 
  else  
    sgn_alpha_surfR[2] = -1.0; 
  
  if (sgn_alpha_surfR[2] == sgn_alpha_surfR[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*alphaR[7]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[3] = 1.0; 
  else  
    sgn_alpha_surfR[3] = -1.0; 
  
  if (sgn_alpha_surfR[3] == sgn_alpha_surfR[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*alphaR[7]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[4] = 1.0; 
  else  
    sgn_alpha_surfR[4] = -1.0; 
  
  if (sgn_alpha_surfR[4] == sgn_alpha_surfR[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*alphaR[7]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[5] = 1.0; 
  else  
    sgn_alpha_surfR[5] = -1.0; 
  
  if (sgn_alpha_surfR[5] == sgn_alpha_surfR[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[7]-0.6363961030678927*alphaR[4]+0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[6] = 1.0; 
  else  
    sgn_alpha_surfR[6] = -1.0; 
  
  if (sgn_alpha_surfR[6] == sgn_alpha_surfR[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[7]-0.6363961030678927*alphaR[4]+0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[7] = 1.0; 
  else  
    sgn_alpha_surfR[7] = -1.0; 
  
  if (sgn_alpha_surfR[7] == sgn_alpha_surfR[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[7]-0.6363961030678927*alphaR[4]+0.4743416490252568*alphaR[2]-0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[8] = 1.0; 
  else  
    sgn_alpha_surfR[8] = -1.0; 
  
  if (sgn_alpha_surfR[8] == sgn_alpha_surfR[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaR[11]-0.3952847075210471*alphaR[7]-0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[9] = 1.0; 
  else  
    sgn_alpha_surfR[9] = -1.0; 
  
  if (sgn_alpha_surfR[9] == sgn_alpha_surfR[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaR[11]-0.3952847075210471*alphaR[7]-0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[10] = 1.0; 
  else  
    sgn_alpha_surfR[10] = -1.0; 
  
  if (sgn_alpha_surfR[10] == sgn_alpha_surfR[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.5303300858899102*alphaR[11]-0.3952847075210471*alphaR[7]-0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[11] = 1.0; 
  else  
    sgn_alpha_surfR[11] = -1.0; 
  
  if (sgn_alpha_surfR[11] == sgn_alpha_surfR[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR[0]-0.3952847075210471*alphaR[7] > 0.) 
    sgn_alpha_surfR[12] = 1.0; 
  else  
    sgn_alpha_surfR[12] = -1.0; 
  
  if (sgn_alpha_surfR[12] == sgn_alpha_surfR[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR[0]-0.3952847075210471*alphaR[7] > 0.) 
    sgn_alpha_surfR[13] = 1.0; 
  else  
    sgn_alpha_surfR[13] = -1.0; 
  
  if (sgn_alpha_surfR[13] == sgn_alpha_surfR[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*alphaR[0]-0.3952847075210471*alphaR[7] > 0.) 
    sgn_alpha_surfR[14] = 1.0; 
  else  
    sgn_alpha_surfR[14] = -1.0; 
  
  if (sgn_alpha_surfR[14] == sgn_alpha_surfR[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaR[11])-0.3952847075210471*alphaR[7]+0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[15] = 1.0; 
  else  
    sgn_alpha_surfR[15] = -1.0; 
  
  if (sgn_alpha_surfR[15] == sgn_alpha_surfR[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaR[11])-0.3952847075210471*alphaR[7]+0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[16] = 1.0; 
  else  
    sgn_alpha_surfR[16] = -1.0; 
  
  if (sgn_alpha_surfR[16] == sgn_alpha_surfR[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.5303300858899102*alphaR[11])-0.3952847075210471*alphaR[7]+0.4743416490252568*alphaR[2]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[17] = 1.0; 
  else  
    sgn_alpha_surfR[17] = -1.0; 
  
  if (sgn_alpha_surfR[17] == sgn_alpha_surfR[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119286*alphaR[11])+0.3162277660168378*alphaR[7]-0.6363961030678927*alphaR[4]-0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[18] = 1.0; 
  else  
    sgn_alpha_surfR[18] = -1.0; 
  
  if (sgn_alpha_surfR[18] == sgn_alpha_surfR[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119286*alphaR[11])+0.3162277660168378*alphaR[7]-0.6363961030678927*alphaR[4]-0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[19] = 1.0; 
  else  
    sgn_alpha_surfR[19] = -1.0; 
  
  if (sgn_alpha_surfR[19] == sgn_alpha_surfR[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.4242640687119286*alphaR[11])+0.3162277660168378*alphaR[7]-0.6363961030678927*alphaR[4]-0.4743416490252568*alphaR[2]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[20] = 1.0; 
  else  
    sgn_alpha_surfR[20] = -1.0; 
  
  if (sgn_alpha_surfR[20] == sgn_alpha_surfR[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*alphaR[7]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[21] = 1.0; 
  else  
    sgn_alpha_surfR[21] = -1.0; 
  
  if (sgn_alpha_surfR[21] == sgn_alpha_surfR[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*alphaR[7]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[22] = 1.0; 
  else  
    sgn_alpha_surfR[22] = -1.0; 
  
  if (sgn_alpha_surfR[22] == sgn_alpha_surfR[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*alphaR[7]+0.4743416490252568*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[23] = 1.0; 
  else  
    sgn_alpha_surfR[23] = -1.0; 
  
  if (sgn_alpha_surfR[23] == sgn_alpha_surfR[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[7]+0.6363961030678927*alphaR[4]+0.4743416490252568*(alphaR[2]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[24] = 1.0; 
  else  
    sgn_alpha_surfR[24] = -1.0; 
  
  if (sgn_alpha_surfR[24] == sgn_alpha_surfR[23]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[7]+0.6363961030678927*alphaR[4]+0.4743416490252568*(alphaR[2]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[25] = 1.0; 
  else  
    sgn_alpha_surfR[25] = -1.0; 
  
  if (sgn_alpha_surfR[25] == sgn_alpha_surfR[24]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4242640687119286*alphaR[11]+0.3162277660168378*alphaR[7]+0.6363961030678927*alphaR[4]+0.4743416490252568*(alphaR[2]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[26] = 1.0; 
  else  
    sgn_alpha_surfR[26] = -1.0; 
  
  if (sgn_alpha_surfR[26] == sgn_alpha_surfR[25]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
