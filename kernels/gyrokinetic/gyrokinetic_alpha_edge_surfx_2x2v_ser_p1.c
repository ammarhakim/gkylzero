#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_edge_surfx_2x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
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
  alphaR[0] = (vmap[1]*(-((3.557562367689425*b_y[3]*jacobtot_inv[3]*hamil[9])/inFlds_e[13][1])-(2.0539595906443724*b_y[2]*jacobtot_inv[3]*hamil[9])/inFlds_e[13][1]-(2.0539595906443724*jacobtot_inv[2]*b_y[3]*hamil[9])/inFlds_e[13][1]-(1.1858541225631418*b_y[2]*jacobtot_inv[2]*hamil[9])/inFlds_e[13][1]-(3.5575623676894255*jacobtot_inv[1]*b_y[3]*hamil[8])/inFlds_e[13][1]-(2.053959590644372*jacobtot_inv[0]*b_y[3]*hamil[8])/inFlds_e[13][1]-(2.053959590644372*jacobtot_inv[1]*b_y[2]*hamil[8])/inFlds_e[13][1]-(1.185854122563142*jacobtot_inv[0]*b_y[2]*hamil[8])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[0]*(-((1.5909902576697312*b_y[3]*jacobtot_inv[3]*hamil[4])/inFlds_e[13][1])-(0.9185586535436913*b_y[2]*jacobtot_inv[3]*hamil[4])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[2]*b_y[3]*hamil[4])/inFlds_e[13][1]-(0.5303300858899105*b_y[2]*jacobtot_inv[2]*hamil[4])/inFlds_e[13][1]-(1.5909902576697312*jacobtot_inv[1]*hamil[2]*b_y[3])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[0]*hamil[2]*b_y[3])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[1]*b_y[2]*hamil[2])/inFlds_e[13][1]-(0.5303300858899105*jacobtot_inv[0]*b_y[2]*hamil[2])/inFlds_e[13][1])*rdvpar2*rdz2+1.7320508075688772*hamil[1]*inFlds_e[8]*inFlds_e[10]*rdz2)/q_; 
  alphaR[1] = (vmap[1]*(-((3.557562367689425*jacobtot_inv[1]*b_y[3]*hamil[9])/inFlds_e[13][1])-(2.0539595906443724*jacobtot_inv[0]*b_y[3]*hamil[9])/inFlds_e[13][1]-(2.0539595906443724*jacobtot_inv[1]*b_y[2]*hamil[9])/inFlds_e[13][1]-(1.1858541225631418*jacobtot_inv[0]*b_y[2]*hamil[9])/inFlds_e[13][1]-(3.5575623676894255*b_y[3]*jacobtot_inv[3]*hamil[8])/inFlds_e[13][1]-(2.053959590644372*b_y[2]*jacobtot_inv[3]*hamil[8])/inFlds_e[13][1]-(2.053959590644372*jacobtot_inv[2]*b_y[3]*hamil[8])/inFlds_e[13][1]-(1.185854122563142*b_y[2]*jacobtot_inv[2]*hamil[8])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[0]*(-((1.5909902576697312*jacobtot_inv[1]*b_y[3]*hamil[4])/inFlds_e[13][1])-(0.9185586535436913*jacobtot_inv[0]*b_y[3]*hamil[4])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[1]*b_y[2]*hamil[4])/inFlds_e[13][1]-(0.5303300858899105*jacobtot_inv[0]*b_y[2]*hamil[4])/inFlds_e[13][1]-(1.5909902576697312*hamil[2]*b_y[3]*jacobtot_inv[3])/inFlds_e[13][1]-(0.9185586535436913*b_y[2]*hamil[2]*jacobtot_inv[3])/inFlds_e[13][1]-(0.9185586535436913*hamil[2]*jacobtot_inv[2]*b_y[3])/inFlds_e[13][1]-(0.5303300858899105*b_y[2]*hamil[2]*jacobtot_inv[2])/inFlds_e[13][1])*rdvpar2*rdz2)/q_; 
  alphaR[2] = (vmap[0]*(-((3.557562367689425*b_y[3]*jacobtot_inv[3]*hamil[9])/inFlds_e[13][1])-(2.0539595906443724*b_y[2]*jacobtot_inv[3]*hamil[9])/inFlds_e[13][1]-(2.0539595906443724*jacobtot_inv[2]*b_y[3]*hamil[9])/inFlds_e[13][1]-(1.1858541225631418*b_y[2]*jacobtot_inv[2]*hamil[9])/inFlds_e[13][1]-(3.5575623676894255*jacobtot_inv[1]*b_y[3]*hamil[8])/inFlds_e[13][1]-(2.053959590644372*jacobtot_inv[0]*b_y[3]*hamil[8])/inFlds_e[13][1]-(2.053959590644372*jacobtot_inv[1]*b_y[2]*hamil[8])/inFlds_e[13][1]-(1.185854122563142*jacobtot_inv[0]*b_y[2]*hamil[8])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[1]*(-((1.5909902576697312*b_y[3]*jacobtot_inv[3]*hamil[4])/inFlds_e[13][1])-(0.9185586535436913*b_y[2]*jacobtot_inv[3]*hamil[4])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[2]*b_y[3]*hamil[4])/inFlds_e[13][1]-(0.5303300858899105*b_y[2]*jacobtot_inv[2]*hamil[4])/inFlds_e[13][1]-(1.5909902576697312*jacobtot_inv[1]*hamil[2]*b_y[3])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[0]*hamil[2]*b_y[3])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[1]*b_y[2]*hamil[2])/inFlds_e[13][1]-(0.5303300858899105*jacobtot_inv[0]*b_y[2]*hamil[2])/inFlds_e[13][1])*rdvpar2*rdz2+1.7320508075688772*hamil[4]*inFlds_e[8]*inFlds_e[10]*rdz2)/q_; 
  alphaR[3] = (vmap[1]*(-((3.5575623676894255*b_y[3]*jacobtot_inv[3]*hamil[11])/inFlds_e[13][1])-(2.053959590644372*b_y[2]*jacobtot_inv[3]*hamil[11])/inFlds_e[13][1]-(2.053959590644372*jacobtot_inv[2]*b_y[3]*hamil[11])/inFlds_e[13][1]-(1.185854122563142*b_y[2]*jacobtot_inv[2]*hamil[11])/inFlds_e[13][1]-(3.557562367689425*jacobtot_inv[1]*b_y[3]*hamil[10])/inFlds_e[13][1]-(2.0539595906443724*jacobtot_inv[0]*b_y[3]*hamil[10])/inFlds_e[13][1]-(2.0539595906443724*jacobtot_inv[1]*b_y[2]*hamil[10])/inFlds_e[13][1]-(1.1858541225631418*jacobtot_inv[0]*b_y[2]*hamil[10])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[0]*(-((1.5909902576697312*b_y[3]*jacobtot_inv[3]*hamil[7])/inFlds_e[13][1])-(0.9185586535436913*b_y[2]*jacobtot_inv[3]*hamil[7])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[2]*b_y[3]*hamil[7])/inFlds_e[13][1]-(0.5303300858899105*b_y[2]*jacobtot_inv[2]*hamil[7])/inFlds_e[13][1]-(1.5909902576697312*jacobtot_inv[1]*b_y[3]*hamil[6])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[0]*b_y[3]*hamil[6])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[1]*b_y[2]*hamil[6])/inFlds_e[13][1]-(0.5303300858899105*jacobtot_inv[0]*b_y[2]*hamil[6])/inFlds_e[13][1])*rdvpar2*rdz2+1.7320508075688772*hamil[5]*inFlds_e[8]*inFlds_e[10]*rdz2)/q_; 
  alphaR[4] = (vmap[0]*(-((3.557562367689425*jacobtot_inv[1]*b_y[3]*hamil[9])/inFlds_e[13][1])-(2.0539595906443724*jacobtot_inv[0]*b_y[3]*hamil[9])/inFlds_e[13][1]-(2.0539595906443724*jacobtot_inv[1]*b_y[2]*hamil[9])/inFlds_e[13][1]-(1.1858541225631418*jacobtot_inv[0]*b_y[2]*hamil[9])/inFlds_e[13][1]-(3.5575623676894255*b_y[3]*jacobtot_inv[3]*hamil[8])/inFlds_e[13][1]-(2.053959590644372*b_y[2]*jacobtot_inv[3]*hamil[8])/inFlds_e[13][1]-(2.053959590644372*jacobtot_inv[2]*b_y[3]*hamil[8])/inFlds_e[13][1]-(1.185854122563142*b_y[2]*jacobtot_inv[2]*hamil[8])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[1]*(-((1.5909902576697312*jacobtot_inv[1]*b_y[3]*hamil[4])/inFlds_e[13][1])-(0.9185586535436913*jacobtot_inv[0]*b_y[3]*hamil[4])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[1]*b_y[2]*hamil[4])/inFlds_e[13][1]-(0.5303300858899105*jacobtot_inv[0]*b_y[2]*hamil[4])/inFlds_e[13][1]-(1.5909902576697312*hamil[2]*b_y[3]*jacobtot_inv[3])/inFlds_e[13][1]-(0.9185586535436913*b_y[2]*hamil[2]*jacobtot_inv[3])/inFlds_e[13][1]-(0.9185586535436913*hamil[2]*jacobtot_inv[2]*b_y[3])/inFlds_e[13][1]-(0.5303300858899105*b_y[2]*hamil[2]*jacobtot_inv[2])/inFlds_e[13][1])*rdvpar2*rdz2)/q_; 
  alphaR[5] = (vmap[1]*(-((3.5575623676894255*jacobtot_inv[1]*b_y[3]*hamil[11])/inFlds_e[13][1])-(2.053959590644372*jacobtot_inv[0]*b_y[3]*hamil[11])/inFlds_e[13][1]-(2.053959590644372*jacobtot_inv[1]*b_y[2]*hamil[11])/inFlds_e[13][1]-(1.185854122563142*jacobtot_inv[0]*b_y[2]*hamil[11])/inFlds_e[13][1]-(3.557562367689425*b_y[3]*jacobtot_inv[3]*hamil[10])/inFlds_e[13][1]-(2.0539595906443724*b_y[2]*jacobtot_inv[3]*hamil[10])/inFlds_e[13][1]-(2.0539595906443724*jacobtot_inv[2]*b_y[3]*hamil[10])/inFlds_e[13][1]-(1.1858541225631418*b_y[2]*jacobtot_inv[2]*hamil[10])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[0]*(-((1.5909902576697312*jacobtot_inv[1]*b_y[3]*hamil[7])/inFlds_e[13][1])-(0.9185586535436913*jacobtot_inv[0]*b_y[3]*hamil[7])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[1]*b_y[2]*hamil[7])/inFlds_e[13][1]-(0.5303300858899105*jacobtot_inv[0]*b_y[2]*hamil[7])/inFlds_e[13][1]-(1.5909902576697312*b_y[3]*jacobtot_inv[3]*hamil[6])/inFlds_e[13][1]-(0.9185586535436913*b_y[2]*jacobtot_inv[3]*hamil[6])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[2]*b_y[3]*hamil[6])/inFlds_e[13][1]-(0.5303300858899105*b_y[2]*jacobtot_inv[2]*hamil[6])/inFlds_e[13][1])*rdvpar2*rdz2)/q_; 
  alphaR[6] = (vmap[0]*(-((3.5575623676894255*b_y[3]*jacobtot_inv[3]*hamil[11])/inFlds_e[13][1])-(2.053959590644372*b_y[2]*jacobtot_inv[3]*hamil[11])/inFlds_e[13][1]-(2.053959590644372*jacobtot_inv[2]*b_y[3]*hamil[11])/inFlds_e[13][1]-(1.185854122563142*b_y[2]*jacobtot_inv[2]*hamil[11])/inFlds_e[13][1]-(3.557562367689425*jacobtot_inv[1]*b_y[3]*hamil[10])/inFlds_e[13][1]-(2.0539595906443724*jacobtot_inv[0]*b_y[3]*hamil[10])/inFlds_e[13][1]-(2.0539595906443724*jacobtot_inv[1]*b_y[2]*hamil[10])/inFlds_e[13][1]-(1.1858541225631418*jacobtot_inv[0]*b_y[2]*hamil[10])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[1]*(-((1.5909902576697312*b_y[3]*jacobtot_inv[3]*hamil[7])/inFlds_e[13][1])-(0.9185586535436913*b_y[2]*jacobtot_inv[3]*hamil[7])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[2]*b_y[3]*hamil[7])/inFlds_e[13][1]-(0.5303300858899105*b_y[2]*jacobtot_inv[2]*hamil[7])/inFlds_e[13][1]-(1.5909902576697312*jacobtot_inv[1]*b_y[3]*hamil[6])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[0]*b_y[3]*hamil[6])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[1]*b_y[2]*hamil[6])/inFlds_e[13][1]-(0.5303300858899105*jacobtot_inv[0]*b_y[2]*hamil[6])/inFlds_e[13][1])*rdvpar2*rdz2+1.7320508075688772*hamil[7]*inFlds_e[8]*inFlds_e[10]*rdz2)/q_; 
  alphaR[7] = (vmap[0]*(-((3.5575623676894255*jacobtot_inv[1]*b_y[3]*hamil[11])/inFlds_e[13][1])-(2.053959590644372*jacobtot_inv[0]*b_y[3]*hamil[11])/inFlds_e[13][1]-(2.053959590644372*jacobtot_inv[1]*b_y[2]*hamil[11])/inFlds_e[13][1]-(1.185854122563142*jacobtot_inv[0]*b_y[2]*hamil[11])/inFlds_e[13][1]-(3.557562367689425*b_y[3]*jacobtot_inv[3]*hamil[10])/inFlds_e[13][1]-(2.0539595906443724*b_y[2]*jacobtot_inv[3]*hamil[10])/inFlds_e[13][1]-(2.0539595906443724*jacobtot_inv[2]*b_y[3]*hamil[10])/inFlds_e[13][1]-(1.1858541225631418*b_y[2]*jacobtot_inv[2]*hamil[10])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[1]*(-((1.5909902576697312*jacobtot_inv[1]*b_y[3]*hamil[7])/inFlds_e[13][1])-(0.9185586535436913*jacobtot_inv[0]*b_y[3]*hamil[7])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[1]*b_y[2]*hamil[7])/inFlds_e[13][1]-(0.5303300858899105*jacobtot_inv[0]*b_y[2]*hamil[7])/inFlds_e[13][1]-(1.5909902576697312*b_y[3]*jacobtot_inv[3]*hamil[6])/inFlds_e[13][1]-(0.9185586535436913*b_y[2]*jacobtot_inv[3]*hamil[6])/inFlds_e[13][1]-(0.9185586535436913*jacobtot_inv[2]*b_y[3]*hamil[6])/inFlds_e[13][1]-(0.5303300858899105*b_y[2]*jacobtot_inv[2]*hamil[6])/inFlds_e[13][1])*rdvpar2*rdz2)/q_; 
  alphaR[8] = (vmap[1]*(-((3.181980515339463*b_y[3]*jacobtot_inv[3]*hamil[9])/inFlds_e[13][1])-(1.8371173070873832*b_y[2]*jacobtot_inv[3]*hamil[9])/inFlds_e[13][1]-(1.8371173070873832*jacobtot_inv[2]*b_y[3]*hamil[9])/inFlds_e[13][1]-(1.060660171779821*b_y[2]*jacobtot_inv[2]*hamil[9])/inFlds_e[13][1]-(3.181980515339463*jacobtot_inv[1]*b_y[3]*hamil[8])/inFlds_e[13][1]-(1.837117307087383*jacobtot_inv[0]*b_y[3]*hamil[8])/inFlds_e[13][1]-(1.837117307087383*jacobtot_inv[1]*b_y[2]*hamil[8])/inFlds_e[13][1]-(1.060660171779821*jacobtot_inv[0]*b_y[2]*hamil[8])/inFlds_e[13][1])*rdvpar2*rdz2+1.7320508075688774*inFlds_e[8]*hamil[9]*inFlds_e[10]*rdz2)/q_; 
  alphaR[9] = (vmap[1]*(-((3.181980515339463*jacobtot_inv[1]*b_y[3]*hamil[9])/inFlds_e[13][1])-(1.837117307087383*jacobtot_inv[0]*b_y[3]*hamil[9])/inFlds_e[13][1]-(1.837117307087383*jacobtot_inv[1]*b_y[2]*hamil[9])/inFlds_e[13][1]-(1.060660171779821*jacobtot_inv[0]*b_y[2]*hamil[9])/inFlds_e[13][1]-(3.181980515339463*b_y[3]*jacobtot_inv[3]*hamil[8])/inFlds_e[13][1]-(1.8371173070873832*b_y[2]*jacobtot_inv[3]*hamil[8])/inFlds_e[13][1]-(1.8371173070873832*jacobtot_inv[2]*b_y[3]*hamil[8])/inFlds_e[13][1]-(1.060660171779821*b_y[2]*jacobtot_inv[2]*hamil[8])/inFlds_e[13][1])*rdvpar2*rdz2)/q_; 
  alphaR[10] = (vmap[1]*(-((3.181980515339463*b_y[3]*jacobtot_inv[3]*hamil[11])/inFlds_e[13][1])-(1.8371173070873832*b_y[2]*jacobtot_inv[3]*hamil[11])/inFlds_e[13][1]-(1.8371173070873832*jacobtot_inv[2]*b_y[3]*hamil[11])/inFlds_e[13][1]-(1.060660171779821*b_y[2]*jacobtot_inv[2]*hamil[11])/inFlds_e[13][1]-(3.181980515339463*jacobtot_inv[1]*b_y[3]*hamil[10])/inFlds_e[13][1]-(1.837117307087383*jacobtot_inv[0]*b_y[3]*hamil[10])/inFlds_e[13][1]-(1.837117307087383*jacobtot_inv[1]*b_y[2]*hamil[10])/inFlds_e[13][1]-(1.060660171779821*jacobtot_inv[0]*b_y[2]*hamil[10])/inFlds_e[13][1])*rdvpar2*rdz2+1.7320508075688774*inFlds_e[8]*inFlds_e[10]*hamil[11]*rdz2)/q_; 
  alphaR[11] = (vmap[1]*(-((3.181980515339463*jacobtot_inv[1]*b_y[3]*hamil[11])/inFlds_e[13][1])-(1.837117307087383*jacobtot_inv[0]*b_y[3]*hamil[11])/inFlds_e[13][1]-(1.837117307087383*jacobtot_inv[1]*b_y[2]*hamil[11])/inFlds_e[13][1]-(1.060660171779821*jacobtot_inv[0]*b_y[2]*hamil[11])/inFlds_e[13][1]-(3.181980515339463*b_y[3]*jacobtot_inv[3]*hamil[10])/inFlds_e[13][1]-(1.8371173070873832*b_y[2]*jacobtot_inv[3]*hamil[10])/inFlds_e[13][1]-(1.8371173070873832*jacobtot_inv[2]*b_y[3]*hamil[10])/inFlds_e[13][1]-(1.060660171779821*b_y[2]*jacobtot_inv[2]*hamil[10])/inFlds_e[13][1])*rdvpar2*rdz2)/q_; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.31622776601683783*alphaR[11]-0.31622776601683783*(alphaR[10]+alphaR[9])+0.31622776601683783*alphaR[8]-0.4743416490252568*alphaR[7]+0.4743416490252568*alphaR[6]+0.3535533905932734*alphaR[5]+0.4743416490252568*alphaR[4]-0.3535533905932734*alphaR[3]-0.4743416490252568*alphaR[2]-0.3535533905932734*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[0] = 1.0; 
  else  
    sgn_alpha_surfR[0] = -1.0; 
  
  if (-(0.3952847075210471*alphaR[11])+0.3952847075210473*(alphaR[10]+alphaR[9])-0.3952847075210471*alphaR[8]+0.3535533905932734*alphaR[5]-0.3535533905932734*(alphaR[3]+alphaR[1])+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[1] = 1.0; 
  else  
    sgn_alpha_surfR[1] = -1.0; 
  
  if (sgn_alpha_surfR[1] == sgn_alpha_surfR[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.31622776601683783*alphaR[11]-0.31622776601683783*(alphaR[10]+alphaR[9])+0.31622776601683783*alphaR[8]+0.4743416490252568*alphaR[7]-0.4743416490252568*alphaR[6]+0.3535533905932734*alphaR[5]-0.4743416490252568*alphaR[4]-0.3535533905932734*alphaR[3]+0.4743416490252568*alphaR[2]-0.3535533905932734*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[2] = 1.0; 
  else  
    sgn_alpha_surfR[2] = -1.0; 
  
  if (sgn_alpha_surfR[2] == sgn_alpha_surfR[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.31622776601683783*alphaR[11])+0.31622776601683783*alphaR[10]-0.31622776601683783*alphaR[9]+0.31622776601683783*alphaR[8]+0.4743416490252568*alphaR[7]-0.4743416490252568*alphaR[6]-0.3535533905932734*alphaR[5]+0.4743416490252568*alphaR[4]+0.3535533905932734*alphaR[3]-0.4743416490252568*alphaR[2]-0.3535533905932734*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[3] = 1.0; 
  else  
    sgn_alpha_surfR[3] = -1.0; 
  
  if (sgn_alpha_surfR[3] == sgn_alpha_surfR[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3952847075210471*alphaR[11]-0.3952847075210473*alphaR[10]+0.3952847075210473*alphaR[9]-0.3952847075210471*alphaR[8]-0.3535533905932734*alphaR[5]+0.3535533905932734*alphaR[3]-0.3535533905932734*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[4] = 1.0; 
  else  
    sgn_alpha_surfR[4] = -1.0; 
  
  if (sgn_alpha_surfR[4] == sgn_alpha_surfR[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.31622776601683783*alphaR[11])+0.31622776601683783*alphaR[10]-0.31622776601683783*alphaR[9]+0.31622776601683783*alphaR[8]-0.4743416490252568*alphaR[7]+0.4743416490252568*alphaR[6]-0.3535533905932734*alphaR[5]-0.4743416490252568*alphaR[4]+0.3535533905932734*alphaR[3]+0.4743416490252568*alphaR[2]-0.3535533905932734*alphaR[1]+0.3535533905932734*alphaR[0] > 0.) 
    sgn_alpha_surfR[5] = 1.0; 
  else  
    sgn_alpha_surfR[5] = -1.0; 
  
  if (sgn_alpha_surfR[5] == sgn_alpha_surfR[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.31622776601683783*(alphaR[11]+alphaR[10]))+0.31622776601683783*(alphaR[9]+alphaR[8])+0.4743416490252568*(alphaR[7]+alphaR[6])-0.3535533905932734*alphaR[5]-0.4743416490252568*alphaR[4]-0.3535533905932734*alphaR[3]-0.4743416490252568*alphaR[2]+0.3535533905932734*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[6] = 1.0; 
  else  
    sgn_alpha_surfR[6] = -1.0; 
  
  if (sgn_alpha_surfR[6] == sgn_alpha_surfR[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3952847075210471*alphaR[11]+0.3952847075210473*alphaR[10]-0.3952847075210473*alphaR[9]-0.3952847075210471*alphaR[8]-0.3535533905932734*(alphaR[5]+alphaR[3])+0.3535533905932734*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[7] = 1.0; 
  else  
    sgn_alpha_surfR[7] = -1.0; 
  
  if (sgn_alpha_surfR[7] == sgn_alpha_surfR[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.31622776601683783*(alphaR[11]+alphaR[10]))+0.31622776601683783*(alphaR[9]+alphaR[8])-0.4743416490252568*(alphaR[7]+alphaR[6])-0.3535533905932734*alphaR[5]+0.4743416490252568*alphaR[4]-0.3535533905932734*alphaR[3]+0.4743416490252568*alphaR[2]+0.3535533905932734*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[8] = 1.0; 
  else  
    sgn_alpha_surfR[8] = -1.0; 
  
  if (sgn_alpha_surfR[8] == sgn_alpha_surfR[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.31622776601683783*(alphaR[11]+alphaR[10]+alphaR[9]+alphaR[8])-0.4743416490252568*(alphaR[7]+alphaR[6])+0.3535533905932734*alphaR[5]-0.4743416490252568*alphaR[4]+0.3535533905932734*alphaR[3]-0.4743416490252568*alphaR[2]+0.3535533905932734*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[9] = 1.0; 
  else  
    sgn_alpha_surfR[9] = -1.0; 
  
  if (sgn_alpha_surfR[9] == sgn_alpha_surfR[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.3952847075210471*alphaR[11])-0.3952847075210473*(alphaR[10]+alphaR[9])-0.3952847075210471*alphaR[8]+0.3535533905932734*(alphaR[5]+alphaR[3]+alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[10] = 1.0; 
  else  
    sgn_alpha_surfR[10] = -1.0; 
  
  if (sgn_alpha_surfR[10] == sgn_alpha_surfR[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.31622776601683783*(alphaR[11]+alphaR[10]+alphaR[9]+alphaR[8])+0.4743416490252568*(alphaR[7]+alphaR[6])+0.3535533905932734*alphaR[5]+0.4743416490252568*alphaR[4]+0.3535533905932734*alphaR[3]+0.4743416490252568*alphaR[2]+0.3535533905932734*(alphaR[1]+alphaR[0]) > 0.) 
    sgn_alpha_surfR[11] = 1.0; 
  else  
    sgn_alpha_surfR[11] = -1.0; 
  
  if (sgn_alpha_surfR[11] == sgn_alpha_surfR[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
