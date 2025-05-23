#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_, const double *bmag, const double *jacobtot_inv,
    const double *cmag, const double *b_i, const double *phi, double* GKYL_RESTRICT alpha_surf,
    double* GKYL_RESTRICT sgn_alpha_surf) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // vmap: velocity space mapping.
  // vmapSq: velocity space mapping squared.
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

  double rdx2 = 2.0/dxv[0];
  double rdz2 = 2.0/dxv[1];
  double rdvpar2 = 2.0/dxv[2];

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[4];
  const double *b_z = &b_i[8];

  double hamil[24] = {0.}; 
  hamil[0] = 2.0*phi[0]*q_+1.414213562373095*(vmapSq[0]*m_+bmag[0]*vmap[2]); 
  hamil[1] = 2.0*phi[1]*q_+1.414213562373095*bmag[1]*vmap[2]; 
  hamil[2] = 2.0*phi[2]*q_+1.414213562373095*bmag[2]*vmap[2]; 
  hamil[3] = 1.414213562373095*vmapSq[1]*m_; 
  hamil[4] = 1.414213562373095*bmag[0]*vmap[3]; 
  hamil[5] = 2.0*phi[3]*q_+1.414213562373095*vmap[2]*bmag[3]; 
  hamil[8] = 1.414213562373095*bmag[1]*vmap[3]; 
  hamil[9] = 1.414213562373095*bmag[2]*vmap[3]; 
  hamil[12] = 1.414213562373095*bmag[3]*vmap[3]; 
  hamil[16] = 1.414213562373095*vmapSq[2]*m_; 

  double *alphaR = &alpha_surf[0];
  double *sgn_alpha_surfR = &sgn_alpha_surf[0];
  alphaR[0] = (((-2.053959590644372*jacobtot_inv[1]*b_y[3]*hamil[16])-1.185854122563142*jacobtot_inv[0]*b_y[3]*hamil[16]-1.185854122563142*jacobtot_inv[1]*b_y[2]*hamil[16]-0.6846531968814573*jacobtot_inv[0]*b_y[2]*hamil[16]+1.590990257669731*b_y[3]*jacobtot_inv[3]*hamil[5]+0.9185586535436913*b_y[2]*jacobtot_inv[3]*hamil[5]+0.9185586535436913*jacobtot_inv[2]*b_y[3]*hamil[5]+0.5303300858899105*b_y[2]*jacobtot_inv[2]*hamil[5]+1.590990257669731*b_y[1]*jacobtot_inv[1]*hamil[5]+0.9185586535436913*b_y[0]*jacobtot_inv[1]*hamil[5]+0.9185586535436913*jacobtot_inv[0]*b_y[1]*hamil[5]+0.5303300858899105*b_y[0]*jacobtot_inv[0]*hamil[5]+0.9185586535436913*hamil[2]*b_y[3]*jacobtot_inv[3]+0.5303300858899105*b_y[2]*hamil[2]*jacobtot_inv[3]+0.5303300858899105*hamil[2]*jacobtot_inv[2]*b_y[3]+0.3061862178478971*b_y[2]*hamil[2]*jacobtot_inv[2]+0.9185586535436913*b_y[1]*jacobtot_inv[1]*hamil[2]+0.5303300858899105*b_y[0]*jacobtot_inv[1]*hamil[2]+0.5303300858899105*jacobtot_inv[0]*b_y[1]*hamil[2]+0.3061862178478971*b_y[0]*jacobtot_inv[0]*hamil[2])*rdz2+(vmap[0]*((-0.9185586535436913*jacobtot_inv[1]*b_y[3]*hamil[3])-0.5303300858899105*jacobtot_inv[0]*b_y[3]*hamil[3]-0.5303300858899105*jacobtot_inv[1]*b_y[2]*hamil[3]-0.3061862178478971*jacobtot_inv[0]*b_y[2]*hamil[3])*rdz2)/vmap[1])/q_; 
  alphaR[1] = (((-2.053959590644372*b_y[3]*jacobtot_inv[3]*hamil[16])-1.185854122563142*b_y[2]*jacobtot_inv[3]*hamil[16]-1.185854122563142*jacobtot_inv[2]*b_y[3]*hamil[16]-0.6846531968814573*b_y[2]*jacobtot_inv[2]*hamil[16]+1.590990257669731*b_y[1]*jacobtot_inv[3]*hamil[5]+0.9185586535436913*b_y[0]*jacobtot_inv[3]*hamil[5]+1.590990257669731*jacobtot_inv[1]*b_y[3]*hamil[5]+0.9185586535436913*jacobtot_inv[0]*b_y[3]*hamil[5]+0.9185586535436913*b_y[1]*jacobtot_inv[2]*hamil[5]+0.5303300858899105*b_y[0]*jacobtot_inv[2]*hamil[5]+0.9185586535436913*jacobtot_inv[1]*b_y[2]*hamil[5]+0.5303300858899105*jacobtot_inv[0]*b_y[2]*hamil[5]+0.9185586535436913*b_y[1]*hamil[2]*jacobtot_inv[3]+0.5303300858899105*b_y[0]*hamil[2]*jacobtot_inv[3]+0.9185586535436913*jacobtot_inv[1]*hamil[2]*b_y[3]+0.5303300858899105*jacobtot_inv[0]*hamil[2]*b_y[3]+0.5303300858899105*b_y[1]*hamil[2]*jacobtot_inv[2]+0.3061862178478971*b_y[0]*hamil[2]*jacobtot_inv[2]+0.5303300858899105*jacobtot_inv[1]*b_y[2]*hamil[2]+0.3061862178478971*jacobtot_inv[0]*b_y[2]*hamil[2])*rdz2+(vmap[0]*((-0.9185586535436913*b_y[3]*hamil[3]*jacobtot_inv[3])-0.5303300858899105*b_y[2]*hamil[3]*jacobtot_inv[3]-0.5303300858899105*jacobtot_inv[2]*b_y[3]*hamil[3]-0.3061862178478971*b_y[2]*jacobtot_inv[2]*hamil[3])*rdz2)/vmap[1])/q_; 
  alphaR[2] = ((vmap[0]*((-2.053959590644372*jacobtot_inv[1]*b_y[3]*hamil[16])-1.185854122563142*jacobtot_inv[0]*b_y[3]*hamil[16]-1.185854122563142*jacobtot_inv[1]*b_y[2]*hamil[16]-0.6846531968814573*jacobtot_inv[0]*b_y[2]*hamil[16])*rdz2)/vmap[1]+((-0.9185586535436913*jacobtot_inv[1]*b_y[3]*hamil[3])-0.5303300858899105*jacobtot_inv[0]*b_y[3]*hamil[3]-0.5303300858899105*jacobtot_inv[1]*b_y[2]*hamil[3]-0.3061862178478971*jacobtot_inv[0]*b_y[2]*hamil[3])*rdz2)/q_; 
  alphaR[3] = ((1.590990257669731*b_y[3]*jacobtot_inv[3]*hamil[12]+0.9185586535436913*b_y[2]*jacobtot_inv[3]*hamil[12]+0.9185586535436913*jacobtot_inv[2]*b_y[3]*hamil[12]+0.5303300858899105*b_y[2]*jacobtot_inv[2]*hamil[12]+1.590990257669731*b_y[1]*jacobtot_inv[1]*hamil[12]+0.9185586535436913*b_y[0]*jacobtot_inv[1]*hamil[12]+0.9185586535436913*jacobtot_inv[0]*b_y[1]*hamil[12]+0.5303300858899105*b_y[0]*jacobtot_inv[0]*hamil[12]+0.9185586535436913*b_y[3]*jacobtot_inv[3]*hamil[9]+0.5303300858899105*b_y[2]*jacobtot_inv[3]*hamil[9]+0.5303300858899105*jacobtot_inv[2]*b_y[3]*hamil[9]+0.3061862178478971*b_y[2]*jacobtot_inv[2]*hamil[9]+0.9185586535436913*b_y[1]*jacobtot_inv[1]*hamil[9]+0.5303300858899105*b_y[0]*jacobtot_inv[1]*hamil[9]+0.5303300858899105*jacobtot_inv[0]*b_y[1]*hamil[9]+0.3061862178478971*b_y[0]*jacobtot_inv[0]*hamil[9])*rdz2)/q_; 
  alphaR[4] = ((vmap[0]*((-2.053959590644372*b_y[3]*jacobtot_inv[3]*hamil[16])-1.185854122563142*b_y[2]*jacobtot_inv[3]*hamil[16]-1.185854122563142*jacobtot_inv[2]*b_y[3]*hamil[16]-0.6846531968814573*b_y[2]*jacobtot_inv[2]*hamil[16])*rdz2)/vmap[1]+((-0.9185586535436913*b_y[3]*hamil[3]*jacobtot_inv[3])-0.5303300858899105*b_y[2]*hamil[3]*jacobtot_inv[3]-0.5303300858899105*jacobtot_inv[2]*b_y[3]*hamil[3]-0.3061862178478971*b_y[2]*jacobtot_inv[2]*hamil[3])*rdz2)/q_; 
  alphaR[5] = ((1.590990257669731*b_y[1]*jacobtot_inv[3]*hamil[12]+0.9185586535436913*b_y[0]*jacobtot_inv[3]*hamil[12]+1.590990257669731*jacobtot_inv[1]*b_y[3]*hamil[12]+0.9185586535436913*jacobtot_inv[0]*b_y[3]*hamil[12]+0.9185586535436913*b_y[1]*jacobtot_inv[2]*hamil[12]+0.5303300858899105*b_y[0]*jacobtot_inv[2]*hamil[12]+0.9185586535436913*jacobtot_inv[1]*b_y[2]*hamil[12]+0.5303300858899105*jacobtot_inv[0]*b_y[2]*hamil[12]+0.9185586535436913*b_y[1]*jacobtot_inv[3]*hamil[9]+0.5303300858899105*b_y[0]*jacobtot_inv[3]*hamil[9]+0.9185586535436913*jacobtot_inv[1]*b_y[3]*hamil[9]+0.5303300858899105*jacobtot_inv[0]*b_y[3]*hamil[9]+0.5303300858899105*b_y[1]*jacobtot_inv[2]*hamil[9]+0.3061862178478971*b_y[0]*jacobtot_inv[2]*hamil[9]+0.5303300858899105*jacobtot_inv[1]*b_y[2]*hamil[9]+0.3061862178478971*jacobtot_inv[0]*b_y[2]*hamil[9])*rdz2)/q_; 
  alphaR[8] = (((-1.837117307087383*jacobtot_inv[1]*b_y[3]*hamil[16])-1.060660171779821*jacobtot_inv[0]*b_y[3]*hamil[16]-1.060660171779821*jacobtot_inv[1]*b_y[2]*hamil[16]-0.6123724356957944*jacobtot_inv[0]*b_y[2]*hamil[16])*rdz2)/q_; 
  alphaR[9] = (((-1.837117307087383*b_y[3]*jacobtot_inv[3]*hamil[16])-1.060660171779821*b_y[2]*jacobtot_inv[3]*hamil[16]-1.060660171779821*jacobtot_inv[2]*b_y[3]*hamil[16]-0.6123724356957944*b_y[2]*jacobtot_inv[2]*hamil[16])*rdz2)/q_; 

  int const_sgn_alpha_surf = 1;  
  
  if ((-0.3162277660168378*alphaR[9])+0.3162277660168378*alphaR[8]+0.3535533905932734*alphaR[5]+0.4743416490252568*alphaR[4]-0.3535533905932734*alphaR[3]-0.4743416490252568*alphaR[2]-0.3535533905932734*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[0] = 1.0; 
  else  
    sgn_alpha_surfR[0] = -1.0; 
  
  if (0.3952847075210473*alphaR[9]-0.3952847075210471*alphaR[8]+0.3535533905932734*alphaR[5]-0.3535533905932734*(alphaR[3]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[1] = 1.0; 
  else  
    sgn_alpha_surfR[1] = -1.0; 
  
  if (sgn_alpha_surfR[1] == sgn_alpha_surfR[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3162277660168378*alphaR[9])+0.3162277660168378*alphaR[8]+0.3535533905932734*alphaR[5]-0.4743416490252568*alphaR[4]-0.3535533905932734*alphaR[3]+0.4743416490252568*alphaR[2]-0.3535533905932734*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[2] = 1.0; 
  else  
    sgn_alpha_surfR[2] = -1.0; 
  
  if (sgn_alpha_surfR[2] == sgn_alpha_surfR[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3162277660168378*alphaR[9])+0.3162277660168378*alphaR[8]-0.3535533905932734*alphaR[5]+0.4743416490252568*alphaR[4]+0.3535533905932734*alphaR[3]-0.4743416490252568*alphaR[2]-0.3535533905932734*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[3] = 1.0; 
  else  
    sgn_alpha_surfR[3] = -1.0; 
  
  if (sgn_alpha_surfR[3] == sgn_alpha_surfR[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3952847075210473*alphaR[9]-0.3952847075210471*alphaR[8]-0.3535533905932734*alphaR[5]+0.3535533905932734*alphaR[3]-0.3535533905932734*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[4] = 1.0; 
  else  
    sgn_alpha_surfR[4] = -1.0; 
  
  if (sgn_alpha_surfR[4] == sgn_alpha_surfR[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3162277660168378*alphaR[9])+0.3162277660168378*alphaR[8]-0.3535533905932734*alphaR[5]-0.4743416490252568*alphaR[4]+0.3535533905932734*alphaR[3]+0.4743416490252568*alphaR[2]-0.3535533905932734*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[5] = 1.0; 
  else  
    sgn_alpha_surfR[5] = -1.0; 
  
  if (sgn_alpha_surfR[5] == sgn_alpha_surfR[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*(alphaR[9]+alphaR[8])-0.3535533905932734*alphaR[5]-0.4743416490252568*alphaR[4]-0.3535533905932734*alphaR[3]-0.4743416490252568*alphaR[2]+0.3535533905932734*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[6] = 1.0; 
  else  
    sgn_alpha_surfR[6] = -1.0; 
  
  if (sgn_alpha_surfR[6] == sgn_alpha_surfR[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3952847075210473*alphaR[9])-0.3952847075210471*alphaR[8]-0.3535533905932734*(alphaR[5]+alphaR[3])+0.3535533905932734*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[7] = 1.0; 
  else  
    sgn_alpha_surfR[7] = -1.0; 
  
  if (sgn_alpha_surfR[7] == sgn_alpha_surfR[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*(alphaR[9]+alphaR[8])-0.3535533905932734*alphaR[5]+0.4743416490252568*alphaR[4]-0.3535533905932734*alphaR[3]+0.4743416490252568*alphaR[2]+0.3535533905932734*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[8] = 1.0; 
  else  
    sgn_alpha_surfR[8] = -1.0; 
  
  if (sgn_alpha_surfR[8] == sgn_alpha_surfR[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*(alphaR[9]+alphaR[8])+0.3535533905932734*alphaR[5]-0.4743416490252568*alphaR[4]+0.3535533905932734*alphaR[3]-0.4743416490252568*alphaR[2]+0.3535533905932734*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[9] = 1.0; 
  else  
    sgn_alpha_surfR[9] = -1.0; 
  
  if (sgn_alpha_surfR[9] == sgn_alpha_surfR[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-0.3952847075210473*alphaR[9])-0.3952847075210471*alphaR[8]+0.3535533905932734*(alphaR[5]+alphaR[3]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[10] = 1.0; 
  else  
    sgn_alpha_surfR[10] = -1.0; 
  
  if (sgn_alpha_surfR[10] == sgn_alpha_surfR[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3162277660168378*(alphaR[9]+alphaR[8])+0.3535533905932734*alphaR[5]+0.4743416490252568*alphaR[4]+0.3535533905932734*alphaR[3]+0.4743416490252568*alphaR[2]+0.3535533905932734*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[11] = 1.0; 
  else  
    sgn_alpha_surfR[11] = -1.0; 
  
  if (sgn_alpha_surfR[11] == sgn_alpha_surfR[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
