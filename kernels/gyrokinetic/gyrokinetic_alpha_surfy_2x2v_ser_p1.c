#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
    const double q_, const double m_,
    const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, 
    const double *bmag_surf, const double *jacobtot_inv_surf, const double *cmag_surf, const double *b_i_surf, 
    const double *phi, double* GKYL_RESTRICT alpha_surf, double* GKYL_RESTRICT sgn_alpha_surf) 
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
  // bmag_surf: bmag represented on the surface.
  // jacobtot_inv_surf: jacobtot_inv represented on the surface.
  // cmag_surf: cmag represented on the surface.
  // b_i_surf: b_i represented on the surface.
  // alpha_surf: output surface phase space flux in each direction (cdim + 1 components).
  //             Note: Each cell owns their *lower* edge surface evaluation.
  // sgn_alpha_surf: output sign(alpha_surf) in each direction at quadrature points (cdim + 1 components).
  //                 Note: Each cell owns their *lower* edge sign(alpha_surf).
  // returns int const_sgn_alpha (true if sign(alpha_surf) is only one sign, either +1 or -1).

  double rdx2 = 2.0/dxv[0];
  double rdz2 = 2.0/dxv[1];
  double rdvpar2 = 2.0/dxv[2];

  const double *b_x = &b_i[0];
  const double *b_x_surf = &b_i_surf[0];
  const double *b_y = &b_i[4];
  const double *b_y_surf = &b_i_surf[2];
  const double *b_z = &b_i[8];
  const double *b_z_surf = &b_i_surf[4];

  double hamil[12] = {0.}; 
  hamil[0] = (1.4142135623730951*phi[0]-2.4494897427831783*phi[2])*q_+vmapSq[0]*m_+1.4142135623730951*bmag_surf[0]*vmap[2]; 
  hamil[1] = 1.4142135623730951*(phi[1]*q_+bmag_surf[1]*vmap[2])-2.4494897427831783*phi[3]*q_; 
  hamil[2] = vmapSq[1]*m_; 
  hamil[3] = 1.4142135623730951*bmag_surf[0]*vmap[3]; 
  hamil[5] = 1.4142135623730951*bmag_surf[1]*vmap[3]; 
  hamil[8] = vmapSq[2]*m_; 

  double *alphaL = &alpha_surf[12];
  double *sgn_alpha_surfL = &sgn_alpha_surf[12];
  alphaL[0] = ((2.9047375096555625*jacobtot_inv[2]*b_y[3]*hamil[8]-1.6770509831248424*jacobtot_inv[0]*b_y[3]*hamil[8]-1.6770509831248424*b_y[1]*jacobtot_inv[2]*hamil[8]+0.9682458365518543*jacobtot_inv[0]*b_y[1]*hamil[8]-0.8660254037844386*b_y_surf[1]*hamil[1]*jacobtot_inv_surf[1]-0.8660254037844386*b_y_surf[0]*jacobtot_inv_surf[0]*hamil[1])*rdx2+(vmap[0]*(1.2990381056766578*hamil[2]*jacobtot_inv[2]*b_y[3]-0.75*jacobtot_inv[0]*hamil[2]*b_y[3]-0.75*b_y[1]*hamil[2]*jacobtot_inv[2]+0.4330127018922193*jacobtot_inv[0]*b_y[1]*hamil[2])*rdx2)/vmap[1])/q_+(1.060660171779821*hamil[2]*cmag[3]*jacobtot_inv[3]-0.6123724356957944*cmag[1]*hamil[2]*jacobtot_inv[3]-0.6123724356957944*jacobtot_inv[1]*hamil[2]*cmag[3]+1.060660171779821*cmag[2]*hamil[2]*jacobtot_inv[2]-0.6123724356957944*cmag[0]*hamil[2]*jacobtot_inv[2]-0.6123724356957944*jacobtot_inv[0]*cmag[2]*hamil[2]+0.3535533905932737*cmag[1]*jacobtot_inv[1]*hamil[2]+0.3535533905932737*cmag[0]*jacobtot_inv[0]*hamil[2])/(vmap[1]*m_); 
  alphaL[1] = ((2.9047375096555625*b_y[3]*jacobtot_inv[3]*hamil[8]-1.6770509831248424*b_y[1]*jacobtot_inv[3]*hamil[8]-1.6770509831248424*jacobtot_inv[1]*b_y[3]*hamil[8]+0.9682458365518543*b_y[1]*jacobtot_inv[1]*hamil[8]-0.8660254037844386*b_y_surf[0]*hamil[1]*jacobtot_inv_surf[1]-0.8660254037844386*jacobtot_inv_surf[0]*b_y_surf[1]*hamil[1])*rdx2+(vmap[0]*(1.2990381056766578*hamil[2]*b_y[3]*jacobtot_inv[3]-0.75*b_y[1]*hamil[2]*jacobtot_inv[3]-0.75*jacobtot_inv[1]*hamil[2]*b_y[3]+0.4330127018922193*b_y[1]*jacobtot_inv[1]*hamil[2])*rdx2)/vmap[1])/q_+(1.060660171779821*cmag[2]*hamil[2]*jacobtot_inv[3]-0.6123724356957944*cmag[0]*hamil[2]*jacobtot_inv[3]+1.060660171779821*hamil[2]*jacobtot_inv[2]*cmag[3]-0.6123724356957944*jacobtot_inv[0]*hamil[2]*cmag[3]-0.6123724356957944*cmag[1]*hamil[2]*jacobtot_inv[2]-0.6123724356957944*jacobtot_inv[1]*cmag[2]*hamil[2]+0.3535533905932737*cmag[0]*jacobtot_inv[1]*hamil[2]+0.3535533905932737*jacobtot_inv[0]*cmag[1]*hamil[2])/(vmap[1]*m_); 
  alphaL[2] = ((vmap[0]*(2.9047375096555625*jacobtot_inv[2]*b_y[3]*hamil[8]-1.6770509831248424*jacobtot_inv[0]*b_y[3]*hamil[8]-1.6770509831248424*b_y[1]*jacobtot_inv[2]*hamil[8]+0.9682458365518543*jacobtot_inv[0]*b_y[1]*hamil[8])*rdx2)/vmap[1]+(1.2990381056766578*hamil[2]*jacobtot_inv[2]*b_y[3]-0.75*jacobtot_inv[0]*hamil[2]*b_y[3]-0.75*b_y[1]*hamil[2]*jacobtot_inv[2]+0.4330127018922193*jacobtot_inv[0]*b_y[1]*hamil[2])*rdx2)/q_+(2.371708245126284*cmag[3]*jacobtot_inv[3]*hamil[8]-1.369306393762915*cmag[1]*jacobtot_inv[3]*hamil[8]-1.369306393762915*jacobtot_inv[1]*cmag[3]*hamil[8]+2.371708245126284*cmag[2]*jacobtot_inv[2]*hamil[8]-1.369306393762915*cmag[0]*jacobtot_inv[2]*hamil[8]-1.369306393762915*jacobtot_inv[0]*cmag[2]*hamil[8]+0.7905694150420947*cmag[1]*jacobtot_inv[1]*hamil[8]+0.7905694150420947*cmag[0]*jacobtot_inv[0]*hamil[8])/(vmap[1]*m_); 
  alphaL[3] = ((-(0.8660254037844386*b_y_surf[1]*jacobtot_inv_surf[1]*hamil[5])-0.8660254037844386*b_y_surf[0]*jacobtot_inv_surf[0]*hamil[5])*rdx2)/q_; 
  alphaL[4] = ((vmap[0]*(2.9047375096555625*b_y[3]*jacobtot_inv[3]*hamil[8]-1.6770509831248424*b_y[1]*jacobtot_inv[3]*hamil[8]-1.6770509831248424*jacobtot_inv[1]*b_y[3]*hamil[8]+0.9682458365518543*b_y[1]*jacobtot_inv[1]*hamil[8])*rdx2)/vmap[1]+(1.2990381056766578*hamil[2]*b_y[3]*jacobtot_inv[3]-0.75*b_y[1]*hamil[2]*jacobtot_inv[3]-0.75*jacobtot_inv[1]*hamil[2]*b_y[3]+0.4330127018922193*b_y[1]*jacobtot_inv[1]*hamil[2])*rdx2)/q_+(2.371708245126284*cmag[2]*jacobtot_inv[3]*hamil[8]-1.369306393762915*cmag[0]*jacobtot_inv[3]*hamil[8]+2.371708245126284*jacobtot_inv[2]*cmag[3]*hamil[8]-1.369306393762915*jacobtot_inv[0]*cmag[3]*hamil[8]-1.369306393762915*cmag[1]*jacobtot_inv[2]*hamil[8]-1.369306393762915*jacobtot_inv[1]*cmag[2]*hamil[8]+0.7905694150420947*cmag[0]*jacobtot_inv[1]*hamil[8]+0.7905694150420947*jacobtot_inv[0]*cmag[1]*hamil[8])/(vmap[1]*m_); 
  alphaL[5] = ((-(0.8660254037844386*b_y_surf[0]*jacobtot_inv_surf[1]*hamil[5])-0.8660254037844386*jacobtot_inv_surf[0]*b_y_surf[1]*hamil[5])*rdx2)/q_; 
  alphaL[8] = ((2.5980762113533156*jacobtot_inv[2]*b_y[3]*hamil[8]-1.5*jacobtot_inv[0]*b_y[3]*hamil[8]-1.5*b_y[1]*jacobtot_inv[2]*hamil[8]+0.8660254037844386*jacobtot_inv[0]*b_y[1]*hamil[8])*rdx2)/q_; 
  alphaL[9] = ((2.598076211353316*b_y[3]*jacobtot_inv[3]*hamil[8]-1.5*b_y[1]*jacobtot_inv[3]*hamil[8]-1.5*jacobtot_inv[1]*b_y[3]*hamil[8]+0.8660254037844387*b_y[1]*jacobtot_inv[1]*hamil[8])*rdx2)/q_; 

  int const_sgn_alpha_surf = 1;  
  
  if (-(0.31622776601683783*alphaL[9])+0.31622776601683783*alphaL[8]+0.3535533905932734*alphaL[5]+0.4743416490252568*alphaL[4]-0.3535533905932734*alphaL[3]-0.4743416490252568*alphaL[2]-0.3535533905932734*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if (0.3952847075210473*alphaL[9]-0.3952847075210471*alphaL[8]+0.3535533905932734*alphaL[5]-0.3535533905932734*(alphaL[3]+alphaL[1])+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.31622776601683783*alphaL[9])+0.31622776601683783*alphaL[8]+0.3535533905932734*alphaL[5]-0.4743416490252568*alphaL[4]-0.3535533905932734*alphaL[3]+0.4743416490252568*alphaL[2]-0.3535533905932734*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.31622776601683783*alphaL[9])+0.31622776601683783*alphaL[8]-0.3535533905932734*alphaL[5]+0.4743416490252568*alphaL[4]+0.3535533905932734*alphaL[3]-0.4743416490252568*alphaL[2]-0.3535533905932734*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3952847075210473*alphaL[9]-0.3952847075210471*alphaL[8]-0.3535533905932734*alphaL[5]+0.3535533905932734*alphaL[3]-0.3535533905932734*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.31622776601683783*alphaL[9])+0.31622776601683783*alphaL[8]-0.3535533905932734*alphaL[5]-0.4743416490252568*alphaL[4]+0.3535533905932734*alphaL[3]+0.4743416490252568*alphaL[2]-0.3535533905932734*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.31622776601683783*(alphaL[9]+alphaL[8])-0.3535533905932734*alphaL[5]-0.4743416490252568*alphaL[4]-0.3535533905932734*alphaL[3]-0.4743416490252568*alphaL[2]+0.3535533905932734*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.3952847075210473*alphaL[9])-0.3952847075210471*alphaL[8]-0.3535533905932734*(alphaL[5]+alphaL[3])+0.3535533905932734*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.31622776601683783*(alphaL[9]+alphaL[8])-0.3535533905932734*alphaL[5]+0.4743416490252568*alphaL[4]-0.3535533905932734*alphaL[3]+0.4743416490252568*alphaL[2]+0.3535533905932734*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.31622776601683783*(alphaL[9]+alphaL[8])+0.3535533905932734*alphaL[5]-0.4743416490252568*alphaL[4]+0.3535533905932734*alphaL[3]-0.4743416490252568*alphaL[2]+0.3535533905932734*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.3952847075210473*alphaL[9])-0.3952847075210471*alphaL[8]+0.3535533905932734*(alphaL[5]+alphaL[3]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.31622776601683783*(alphaL[9]+alphaL[8])+0.3535533905932734*alphaL[5]+0.4743416490252568*alphaL[4]+0.3535533905932734*alphaL[3]+0.4743416490252568*alphaL[2]+0.3535533905932734*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
