#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
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
  const double *b_y_surf = &b_i_surf[4];
  const double *b_z = &b_i[8];
  const double *b_z_surf = &b_i_surf[8];

  double hamil[12] = {0.}; 
  hamil[0] = (2.4494897427831783*phi[1]+1.4142135623730951*phi[0])*q_+vmapSq[0]*m_+0.7071067811865475*(bmag_surf[3]*vmap[3]+bmag_surf[0]*vmap[2]); 
  hamil[1] = (2.4494897427831783*phi[3]+1.4142135623730951*phi[2])*q_+0.7071067811865475*(vmap[3]*bmag_surf[5]+bmag_surf[1]*vmap[2]); 
  hamil[2] = vmapSq[1]*m_+0.7071067811865475*(vmap[3]*bmag_surf[6]+bmag_surf[2]*vmap[2]); 
  hamil[3] = 0.7071067811865475*(bmag_surf[0]*vmap[3]+vmap[2]*bmag_surf[3]); 
  hamil[4] = 0.7071067811865475*(vmap[3]*bmag_surf[7]+vmap[2]*bmag_surf[4]); 
  hamil[5] = 0.7071067811865475*(vmap[2]*bmag_surf[5]+bmag_surf[1]*vmap[3]); 
  hamil[6] = 0.7071067811865475*(vmap[2]*bmag_surf[6]+bmag_surf[2]*vmap[3]); 
  hamil[7] = 0.7071067811865475*(vmap[2]*bmag_surf[7]+vmap[3]*bmag_surf[4]); 
  hamil[8] = vmapSq[2]*m_+0.7071067811865475*(vmap[3]*bmag_surf[10]+vmap[2]*bmag_surf[8]); 
  hamil[9] = 0.7071067811865475*(vmap[3]*bmag_surf[11]+vmap[2]*bmag_surf[9]); 
  hamil[10] = 0.7071067811865475*(vmap[2]*bmag_surf[10]+vmap[3]*bmag_surf[8]); 
  hamil[11] = 0.7071067811865475*(vmap[2]*bmag_surf[11]+vmap[3]*bmag_surf[9]); 

  double *alphaR = &alpha_surf[0];
  double *sgn_alpha_surfR = &sgn_alpha_surf[0];

  int const_sgn_alpha_surf = 1;  
  
  if (0.0 > 0.) 
    sgn_alpha_surfR[0] = 1.0; 
  else  
    sgn_alpha_surfR[0] = -1.0; 
  
  if (0.0 > 0.) 
    sgn_alpha_surfR[1] = 1.0; 
  else  
    sgn_alpha_surfR[1] = -1.0; 
  
  if (sgn_alpha_surfR[1] == sgn_alpha_surfR[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
