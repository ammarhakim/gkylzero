#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_no_by_edge_surfz_3x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
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
  double rdy2 = 2.0/dxv[1];
  double rdz2 = 2.0/dxv[2];
  double rdvpar2 = 2.0/dxv[3];

  const double *b_x = &b_i[0];
  const double *b_x_surf = &b_i_surf[0];
  const double *b_y = &b_i[8];
  const double *b_y_surf = &b_i_surf[8];
  const double *b_z = &b_i[16];
  const double *b_z_surf = &b_i_surf[16];

  double hamil[24] = {0.}; 
  hamil[0] = 2.4494897427831783*phi[3]*q_+1.4142135623730951*(phi[0]*q_+vmapSq[0]*m_)+0.7071067811865475*(vmap[3]*bmag_surf[4]+bmag_surf[0]*vmap[2]); 
  hamil[1] = (2.4494897427831783*phi[5]+1.4142135623730951*phi[1])*q_+0.7071067811865475*(vmap[3]*bmag_surf[8]+bmag_surf[1]*vmap[2]); 
  hamil[2] = (2.4494897427831783*phi[6]+1.4142135623730951*phi[2])*q_; 
  hamil[3] = 1.4142135623730951*vmapSq[1]*m_+0.7071067811865475*(vmap[3]*bmag_surf[10]+vmap[2]*bmag_surf[3]); 
  hamil[4] = 0.7071067811865475*(vmap[2]*bmag_surf[4]+bmag_surf[0]*vmap[3]); 
  hamil[5] = (2.4494897427831783*phi[7]+1.4142135623730951*phi[4])*q_; 
  hamil[6] = 0.7071067811865475*(vmap[3]*bmag_surf[13]+vmap[2]*bmag_surf[6]); 
  hamil[8] = 0.7071067811865475*(vmap[2]*bmag_surf[8]+bmag_surf[1]*vmap[3]); 
  hamil[10] = 0.7071067811865475*(vmap[2]*bmag_surf[10]+bmag_surf[3]*vmap[3]); 
  hamil[13] = 0.7071067811865475*(vmap[2]*bmag_surf[13]+vmap[3]*bmag_surf[6]); 
  hamil[16] = 1.4142135623730951*vmapSq[2]*m_+0.7071067811865475*(vmap[3]*bmag_surf[19]+vmap[2]*bmag_surf[16]); 
  hamil[17] = 0.7071067811865475*(vmap[3]*bmag_surf[21]+vmap[2]*bmag_surf[17]); 
  hamil[19] = 0.7071067811865475*(vmap[2]*bmag_surf[19]+vmap[3]*bmag_surf[16]); 
  hamil[21] = 0.7071067811865475*(vmap[2]*bmag_surf[21]+vmap[3]*bmag_surf[17]); 

  double *alphaR = &alpha_surf[48];
  double *sgn_alpha_surfR = &sgn_alpha_surf[48];
  alphaR[0] = (1.7320508075688772*hamil[2]*inFlds_e[7]*inFlds_e[10]*rdy2)/q_+(((0.6495190528383289*cmag[3]*jacobtot_inv[5]*hamil[6])/inFlds_e[13][1]+(0.375*cmag[0]*jacobtot_inv[5]*hamil[6])/inFlds_e[13][1]+(0.6495190528383289*jacobtot_inv[3]*cmag[5]*hamil[6])/inFlds_e[13][1]+(0.375*jacobtot_inv[0]*cmag[5]*hamil[6])/inFlds_e[13][1]+(0.375*cmag[1]*jacobtot_inv[3]*hamil[6])/inFlds_e[13][1]+(0.375*jacobtot_inv[1]*cmag[3]*hamil[6])/inFlds_e[13][1]+(0.21650635094610965*cmag[0]*jacobtot_inv[1]*hamil[6])/inFlds_e[13][1]+(0.21650635094610965*jacobtot_inv[0]*cmag[1]*hamil[6])/inFlds_e[13][1]+(0.6495190528383289*hamil[3]*cmag[5]*jacobtot_inv[5])/inFlds_e[13][1]+(0.375*cmag[1]*hamil[3]*jacobtot_inv[5])/inFlds_e[13][1]+(0.375*jacobtot_inv[1]*hamil[3]*cmag[5])/inFlds_e[13][1]+(0.6495190528383289*cmag[3]*hamil[3]*jacobtot_inv[3])/inFlds_e[13][1]+(0.375*cmag[0]*hamil[3]*jacobtot_inv[3])/inFlds_e[13][1]+(0.375*jacobtot_inv[0]*cmag[3]*hamil[3])/inFlds_e[13][1]+(0.21650635094610965*cmag[1]*jacobtot_inv[1]*hamil[3])/inFlds_e[13][1]+(0.21650635094610965*cmag[0]*jacobtot_inv[0]*hamil[3])/inFlds_e[13][1])*rdvpar2)/m_; 
  alphaR[1] = (1.7320508075688772*hamil[5]*inFlds_e[7]*inFlds_e[10]*rdy2)/q_+(((1.1691342951089918*cmag[5]*jacobtot_inv[5]*hamil[6])/inFlds_e[13][1]+(0.675*cmag[1]*jacobtot_inv[5]*hamil[6])/inFlds_e[13][1]+(0.675*jacobtot_inv[1]*cmag[5]*hamil[6])/inFlds_e[13][1]+(0.6495190528383289*cmag[3]*jacobtot_inv[3]*hamil[6])/inFlds_e[13][1]+(0.375*cmag[0]*jacobtot_inv[3]*hamil[6])/inFlds_e[13][1]+(0.375*jacobtot_inv[0]*cmag[3]*hamil[6])/inFlds_e[13][1]+(0.38971143170299727*cmag[1]*jacobtot_inv[1]*hamil[6])/inFlds_e[13][1]+(0.21650635094610965*cmag[0]*jacobtot_inv[0]*hamil[6])/inFlds_e[13][1]+(0.6495190528383289*cmag[3]*hamil[3]*jacobtot_inv[5])/inFlds_e[13][1]+(0.375*cmag[0]*hamil[3]*jacobtot_inv[5])/inFlds_e[13][1]+(0.6495190528383289*hamil[3]*jacobtot_inv[3]*cmag[5])/inFlds_e[13][1]+(0.375*jacobtot_inv[0]*hamil[3]*cmag[5])/inFlds_e[13][1]+(0.375*cmag[1]*hamil[3]*jacobtot_inv[3])/inFlds_e[13][1]+(0.375*jacobtot_inv[1]*cmag[3]*hamil[3])/inFlds_e[13][1]+(0.21650635094610965*cmag[0]*jacobtot_inv[1]*hamil[3])/inFlds_e[13][1]+(0.21650635094610965*jacobtot_inv[0]*cmag[1]*hamil[3])/inFlds_e[13][1])*rdvpar2)/m_; 
  alphaR[3] = (((1.4523687548277815*cmag[3]*jacobtot_inv[5]*hamil[17])/inFlds_e[13][1]+(0.8385254915624211*cmag[0]*jacobtot_inv[5]*hamil[17])/inFlds_e[13][1]+(1.4523687548277815*jacobtot_inv[3]*cmag[5]*hamil[17])/inFlds_e[13][1]+(0.8385254915624211*jacobtot_inv[0]*cmag[5]*hamil[17])/inFlds_e[13][1]+(0.8385254915624211*cmag[1]*jacobtot_inv[3]*hamil[17])/inFlds_e[13][1]+(0.8385254915624211*jacobtot_inv[1]*cmag[3]*hamil[17])/inFlds_e[13][1]+(0.4841229182759271*cmag[0]*jacobtot_inv[1]*hamil[17])/inFlds_e[13][1]+(0.4841229182759271*jacobtot_inv[0]*cmag[1]*hamil[17])/inFlds_e[13][1]+(1.4523687548277813*cmag[5]*jacobtot_inv[5]*hamil[16])/inFlds_e[13][1]+(0.8385254915624212*cmag[1]*jacobtot_inv[5]*hamil[16])/inFlds_e[13][1]+(0.8385254915624212*jacobtot_inv[1]*cmag[5]*hamil[16])/inFlds_e[13][1]+(1.4523687548277813*cmag[3]*jacobtot_inv[3]*hamil[16])/inFlds_e[13][1]+(0.8385254915624212*cmag[0]*jacobtot_inv[3]*hamil[16])/inFlds_e[13][1]+(0.8385254915624212*jacobtot_inv[0]*cmag[3]*hamil[16])/inFlds_e[13][1]+(0.4841229182759271*cmag[1]*jacobtot_inv[1]*hamil[16])/inFlds_e[13][1]+(0.4841229182759271*cmag[0]*jacobtot_inv[0]*hamil[16])/inFlds_e[13][1])*rdvpar2)/m_; 
  alphaR[4] = (((0.6495190528383289*cmag[3]*jacobtot_inv[5]*hamil[13])/inFlds_e[13][1]+(0.375*cmag[0]*jacobtot_inv[5]*hamil[13])/inFlds_e[13][1]+(0.6495190528383289*jacobtot_inv[3]*cmag[5]*hamil[13])/inFlds_e[13][1]+(0.375*jacobtot_inv[0]*cmag[5]*hamil[13])/inFlds_e[13][1]+(0.375*cmag[1]*jacobtot_inv[3]*hamil[13])/inFlds_e[13][1]+(0.375*jacobtot_inv[1]*cmag[3]*hamil[13])/inFlds_e[13][1]+(0.21650635094610965*cmag[0]*jacobtot_inv[1]*hamil[13])/inFlds_e[13][1]+(0.21650635094610965*jacobtot_inv[0]*cmag[1]*hamil[13])/inFlds_e[13][1]+(0.6495190528383289*cmag[5]*jacobtot_inv[5]*hamil[10])/inFlds_e[13][1]+(0.375*cmag[1]*jacobtot_inv[5]*hamil[10])/inFlds_e[13][1]+(0.375*jacobtot_inv[1]*cmag[5]*hamil[10])/inFlds_e[13][1]+(0.6495190528383289*cmag[3]*jacobtot_inv[3]*hamil[10])/inFlds_e[13][1]+(0.375*cmag[0]*jacobtot_inv[3]*hamil[10])/inFlds_e[13][1]+(0.375*jacobtot_inv[0]*cmag[3]*hamil[10])/inFlds_e[13][1]+(0.21650635094610965*cmag[1]*jacobtot_inv[1]*hamil[10])/inFlds_e[13][1]+(0.21650635094610965*cmag[0]*jacobtot_inv[0]*hamil[10])/inFlds_e[13][1])*rdvpar2)/m_; 
  alphaR[6] = (((2.6142637586900066*cmag[5]*jacobtot_inv[5]*hamil[17])/inFlds_e[13][1]+(1.5093458848123575*cmag[1]*jacobtot_inv[5]*hamil[17])/inFlds_e[13][1]+(1.5093458848123575*jacobtot_inv[1]*cmag[5]*hamil[17])/inFlds_e[13][1]+(1.4523687548277815*cmag[3]*jacobtot_inv[3]*hamil[17])/inFlds_e[13][1]+(0.8385254915624211*cmag[0]*jacobtot_inv[3]*hamil[17])/inFlds_e[13][1]+(0.8385254915624211*jacobtot_inv[0]*cmag[3]*hamil[17])/inFlds_e[13][1]+(0.8714212528966688*cmag[1]*jacobtot_inv[1]*hamil[17])/inFlds_e[13][1]+(0.4841229182759271*cmag[0]*jacobtot_inv[0]*hamil[17])/inFlds_e[13][1]+(1.4523687548277813*cmag[3]*jacobtot_inv[5]*hamil[16])/inFlds_e[13][1]+(0.8385254915624212*cmag[0]*jacobtot_inv[5]*hamil[16])/inFlds_e[13][1]+(1.4523687548277813*jacobtot_inv[3]*cmag[5]*hamil[16])/inFlds_e[13][1]+(0.8385254915624212*jacobtot_inv[0]*cmag[5]*hamil[16])/inFlds_e[13][1]+(0.8385254915624212*cmag[1]*jacobtot_inv[3]*hamil[16])/inFlds_e[13][1]+(0.8385254915624212*jacobtot_inv[1]*cmag[3]*hamil[16])/inFlds_e[13][1]+(0.4841229182759271*cmag[0]*jacobtot_inv[1]*hamil[16])/inFlds_e[13][1]+(0.4841229182759271*jacobtot_inv[0]*cmag[1]*hamil[16])/inFlds_e[13][1])*rdvpar2)/m_; 
  alphaR[8] = (((1.1691342951089918*cmag[5]*jacobtot_inv[5]*hamil[13])/inFlds_e[13][1]+(0.675*cmag[1]*jacobtot_inv[5]*hamil[13])/inFlds_e[13][1]+(0.675*jacobtot_inv[1]*cmag[5]*hamil[13])/inFlds_e[13][1]+(0.6495190528383289*cmag[3]*jacobtot_inv[3]*hamil[13])/inFlds_e[13][1]+(0.375*cmag[0]*jacobtot_inv[3]*hamil[13])/inFlds_e[13][1]+(0.375*jacobtot_inv[0]*cmag[3]*hamil[13])/inFlds_e[13][1]+(0.38971143170299727*cmag[1]*jacobtot_inv[1]*hamil[13])/inFlds_e[13][1]+(0.21650635094610965*cmag[0]*jacobtot_inv[0]*hamil[13])/inFlds_e[13][1]+(0.6495190528383289*cmag[3]*jacobtot_inv[5]*hamil[10])/inFlds_e[13][1]+(0.375*cmag[0]*jacobtot_inv[5]*hamil[10])/inFlds_e[13][1]+(0.6495190528383289*jacobtot_inv[3]*cmag[5]*hamil[10])/inFlds_e[13][1]+(0.375*jacobtot_inv[0]*cmag[5]*hamil[10])/inFlds_e[13][1]+(0.375*cmag[1]*jacobtot_inv[3]*hamil[10])/inFlds_e[13][1]+(0.375*jacobtot_inv[1]*cmag[3]*hamil[10])/inFlds_e[13][1]+(0.21650635094610965*cmag[0]*jacobtot_inv[1]*hamil[10])/inFlds_e[13][1]+(0.21650635094610965*jacobtot_inv[0]*cmag[1]*hamil[10])/inFlds_e[13][1])*rdvpar2)/m_; 
  alphaR[10] = (((1.4523687548277813*cmag[3]*jacobtot_inv[5]*hamil[21])/inFlds_e[13][1]+(0.8385254915624212*cmag[0]*jacobtot_inv[5]*hamil[21])/inFlds_e[13][1]+(1.4523687548277813*jacobtot_inv[3]*cmag[5]*hamil[21])/inFlds_e[13][1]+(0.8385254915624212*jacobtot_inv[0]*cmag[5]*hamil[21])/inFlds_e[13][1]+(0.8385254915624212*cmag[1]*jacobtot_inv[3]*hamil[21])/inFlds_e[13][1]+(0.8385254915624212*jacobtot_inv[1]*cmag[3]*hamil[21])/inFlds_e[13][1]+(0.4841229182759271*cmag[0]*jacobtot_inv[1]*hamil[21])/inFlds_e[13][1]+(0.4841229182759271*jacobtot_inv[0]*cmag[1]*hamil[21])/inFlds_e[13][1]+(1.4523687548277815*cmag[5]*jacobtot_inv[5]*hamil[19])/inFlds_e[13][1]+(0.8385254915624211*cmag[1]*jacobtot_inv[5]*hamil[19])/inFlds_e[13][1]+(0.8385254915624211*jacobtot_inv[1]*cmag[5]*hamil[19])/inFlds_e[13][1]+(1.4523687548277815*cmag[3]*jacobtot_inv[3]*hamil[19])/inFlds_e[13][1]+(0.8385254915624211*cmag[0]*jacobtot_inv[3]*hamil[19])/inFlds_e[13][1]+(0.8385254915624211*jacobtot_inv[0]*cmag[3]*hamil[19])/inFlds_e[13][1]+(0.4841229182759271*cmag[1]*jacobtot_inv[1]*hamil[19])/inFlds_e[13][1]+(0.4841229182759271*cmag[0]*jacobtot_inv[0]*hamil[19])/inFlds_e[13][1])*rdvpar2)/m_; 
  alphaR[13] = (((2.6142637586900053*cmag[5]*jacobtot_inv[5]*hamil[21])/inFlds_e[13][1]+(1.509345884812358*cmag[1]*jacobtot_inv[5]*hamil[21])/inFlds_e[13][1]+(1.509345884812358*jacobtot_inv[1]*cmag[5]*hamil[21])/inFlds_e[13][1]+(1.4523687548277813*cmag[3]*jacobtot_inv[3]*hamil[21])/inFlds_e[13][1]+(0.8385254915624212*cmag[0]*jacobtot_inv[3]*hamil[21])/inFlds_e[13][1]+(0.8385254915624212*jacobtot_inv[0]*cmag[3]*hamil[21])/inFlds_e[13][1]+(0.8714212528966685*cmag[1]*jacobtot_inv[1]*hamil[21])/inFlds_e[13][1]+(0.4841229182759271*cmag[0]*jacobtot_inv[0]*hamil[21])/inFlds_e[13][1]+(1.4523687548277815*cmag[3]*jacobtot_inv[5]*hamil[19])/inFlds_e[13][1]+(0.8385254915624211*cmag[0]*jacobtot_inv[5]*hamil[19])/inFlds_e[13][1]+(1.4523687548277815*jacobtot_inv[3]*cmag[5]*hamil[19])/inFlds_e[13][1]+(0.8385254915624211*jacobtot_inv[0]*cmag[5]*hamil[19])/inFlds_e[13][1]+(0.8385254915624211*cmag[1]*jacobtot_inv[3]*hamil[19])/inFlds_e[13][1]+(0.8385254915624211*jacobtot_inv[1]*cmag[3]*hamil[19])/inFlds_e[13][1]+(0.4841229182759271*cmag[0]*jacobtot_inv[1]*hamil[19])/inFlds_e[13][1]+(0.4841229182759271*jacobtot_inv[0]*cmag[1]*hamil[19])/inFlds_e[13][1])*rdvpar2)/m_; 

  int const_sgn_alpha_surf = 1;  
  
  if ((0.5809475019311127*alphaR[13]-0.5809475019311127*alphaR[10]-0.43301270189221924*alphaR[8]+0.43301270189221924*alphaR[4])*mu+0.33541019662496785*alphaR[6]-0.33541019662496785*alphaR[3]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[0] = 1.0; 
  else  
    sgn_alpha_surfR[0] = -1.0; 
  
  if ((0.43301270189221924*alphaR[4]-0.43301270189221924*alphaR[8])*mu-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[1] = 1.0; 
  else  
    sgn_alpha_surfR[1] = -1.0; 
  
  if (sgn_alpha_surfR[1] == sgn_alpha_surfR[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((-(0.5809475019311127*alphaR[13])+0.5809475019311127*alphaR[10]-0.43301270189221924*alphaR[8]+0.43301270189221924*alphaR[4])*mu-0.33541019662496785*alphaR[6]+0.33541019662496785*alphaR[3]-0.25*alphaR[1]+0.25*alphaR[0] > 0.) 
    sgn_alpha_surfR[2] = 1.0; 
  else  
    sgn_alpha_surfR[2] = -1.0; 
  
  if (sgn_alpha_surfR[2] == sgn_alpha_surfR[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((0.43301270189221924*(alphaR[8]+alphaR[4])-0.5809475019311127*(alphaR[13]+alphaR[10]))*mu-0.33541019662496785*(alphaR[6]+alphaR[3])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[3] = 1.0; 
  else  
    sgn_alpha_surfR[3] = -1.0; 
  
  if (sgn_alpha_surfR[3] == sgn_alpha_surfR[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.43301270189221924*(alphaR[8]+alphaR[4])*mu+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[4] = 1.0; 
  else  
    sgn_alpha_surfR[4] = -1.0; 
  
  if (sgn_alpha_surfR[4] == sgn_alpha_surfR[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if ((0.5809475019311127*(alphaR[13]+alphaR[10])+0.43301270189221924*(alphaR[8]+alphaR[4]))*mu+0.33541019662496785*(alphaR[6]+alphaR[3])+0.25*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[5] = 1.0; 
  else  
    sgn_alpha_surfR[5] = -1.0; 
  
  if (sgn_alpha_surfR[5] == sgn_alpha_surfR[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
