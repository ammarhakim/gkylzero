#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfy_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
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
  alphaL[0] = (0.7071067811865475*cmag_surf[1]*jacobtot_inv_surf[1]*hamil[2]+0.7071067811865475*cmag_surf[0]*jacobtot_inv_surf[0]*hamil[2])/(vmap[1]*m_); 
  alphaL[1] = (0.7071067811865475*cmag_surf[0]*jacobtot_inv_surf[1]*hamil[2]+0.7071067811865475*jacobtot_inv_surf[0]*cmag_surf[1]*hamil[2])/(vmap[1]*m_); 
  alphaL[2] = (1.5811388300841895*cmag_surf[1]*jacobtot_inv_surf[1]*hamil[8]+1.5811388300841895*cmag_surf[0]*jacobtot_inv_surf[0]*hamil[8])/(vmap[1]*m_); 
  alphaL[4] = (1.5811388300841895*cmag_surf[0]*jacobtot_inv_surf[1]*hamil[8]+1.5811388300841895*jacobtot_inv_surf[0]*cmag_surf[1]*hamil[8])/(vmap[1]*m_); 

  int const_sgn_alpha_surf = 1;  
  
  if (0.4743416490252568*alphaL[4]-0.4743416490252568*alphaL[2]-0.3535533905932734*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if (0.3535533905932734*alphaL[0]-0.3535533905932734*alphaL[1] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.4743416490252568*alphaL[4])+0.4743416490252568*alphaL[2]-0.3535533905932734*alphaL[1]+0.3535533905932734*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*(alphaL[1]+alphaL[0])-0.4743416490252568*(alphaL[4]+alphaL[2]) > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.4743416490252568*(alphaL[4]+alphaL[2])+0.3535533905932734*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
