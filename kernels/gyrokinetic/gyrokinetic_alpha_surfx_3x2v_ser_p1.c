#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_surfx_3x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
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
  hamil[0] = -(2.4494897427831783*phi[1]*q_)+1.4142135623730951*(phi[0]*q_+vmapSq[0]*m_)+0.7071067811865475*(vmap[3]*bmag_surf[4]+bmag_surf[0]*vmap[2]); 
  hamil[1] = (1.4142135623730951*phi[2]-2.4494897427831783*phi[4])*q_; 
  hamil[2] = (1.4142135623730951*phi[3]-2.4494897427831783*phi[5])*q_+0.7071067811865475*(vmap[3]*bmag_surf[9]+bmag_surf[2]*vmap[2]); 
  hamil[3] = 1.4142135623730951*vmapSq[1]*m_+0.7071067811865475*(vmap[3]*bmag_surf[10]+vmap[2]*bmag_surf[3]); 
  hamil[4] = 0.7071067811865475*(vmap[2]*bmag_surf[4]+bmag_surf[0]*vmap[3]); 
  hamil[5] = (1.4142135623730951*phi[6]-2.4494897427831783*phi[7])*q_; 
  hamil[7] = 0.7071067811865475*(vmap[3]*bmag_surf[14]+vmap[2]*bmag_surf[7]); 
  hamil[9] = 0.7071067811865475*(vmap[2]*bmag_surf[9]+bmag_surf[2]*vmap[3]); 
  hamil[10] = 0.7071067811865475*(vmap[2]*bmag_surf[10]+bmag_surf[3]*vmap[3]); 
  hamil[14] = 0.7071067811865475*(vmap[2]*bmag_surf[14]+vmap[3]*bmag_surf[7]); 
  hamil[16] = 1.4142135623730951*vmapSq[2]*m_+0.7071067811865475*(vmap[3]*bmag_surf[19]+vmap[2]*bmag_surf[16]); 
  hamil[18] = 0.7071067811865475*(vmap[3]*bmag_surf[22]+vmap[2]*bmag_surf[18]); 
  hamil[19] = 0.7071067811865475*(vmap[2]*bmag_surf[19]+vmap[3]*bmag_surf[16]); 
  hamil[22] = 0.7071067811865475*(vmap[2]*bmag_surf[22]+vmap[3]*bmag_surf[18]); 

  double *alphaL = &alpha_surf[0];
  double *sgn_alpha_surfL = &sgn_alpha_surf[0];
  alphaL[0] = (vmap[1]*(-((1.7787811838447125*b_y[5]*jacobtot_inv[5]*hamil[18])/inFlds_e[13][1])+(1.0269797953221862*b_y[3]*jacobtot_inv[5]*hamil[18])/inFlds_e[13][1]+(1.0269797953221862*jacobtot_inv[3]*b_y[5]*hamil[18])/inFlds_e[13][1]-(0.5929270612815709*b_y[3]*jacobtot_inv[3]*hamil[18])/inFlds_e[13][1]-(1.7787811838447127*jacobtot_inv[1]*b_y[5]*hamil[16])/inFlds_e[13][1]+(1.026979795322186*jacobtot_inv[0]*b_y[5]*hamil[16])/inFlds_e[13][1]+(1.026979795322186*jacobtot_inv[1]*b_y[3]*hamil[16])/inFlds_e[13][1]-(0.592927061281571*jacobtot_inv[0]*b_y[3]*hamil[16])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[0]*(-((0.7954951288348656*b_y[5]*jacobtot_inv[5]*hamil[7])/inFlds_e[13][1])+(0.45927932677184563*b_y[3]*jacobtot_inv[5]*hamil[7])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[3]*b_y[5]*hamil[7])/inFlds_e[13][1]-(0.26516504294495524*b_y[3]*jacobtot_inv[3]*hamil[7])/inFlds_e[13][1]-(0.7954951288348656*jacobtot_inv[1]*hamil[3]*b_y[5])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[0]*hamil[3]*b_y[5])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[1]*b_y[3]*hamil[3])/inFlds_e[13][1]-(0.26516504294495524*jacobtot_inv[0]*b_y[3]*hamil[3])/inFlds_e[13][1])*rdvpar2*rdz2+1.7320508075688772*hamil[2]*inFlds_e[8]*inFlds_e[10]*rdz2-1.7320508075688772*hamil[1]*inFlds_e[9]*inFlds_e[10]*rdy2)/q_; 
  alphaL[1] = (1.7320508075688772*hamil[5]*inFlds_e[8]*inFlds_e[10]*rdz2)/q_; 
  alphaL[2] = (vmap[1]*(-((1.7787811838447125*jacobtot_inv[1]*b_y[5]*hamil[18])/inFlds_e[13][1])+(1.0269797953221862*jacobtot_inv[0]*b_y[5]*hamil[18])/inFlds_e[13][1]+(1.0269797953221862*jacobtot_inv[1]*b_y[3]*hamil[18])/inFlds_e[13][1]-(0.5929270612815709*jacobtot_inv[0]*b_y[3]*hamil[18])/inFlds_e[13][1]-(1.7787811838447127*b_y[5]*jacobtot_inv[5]*hamil[16])/inFlds_e[13][1]+(1.026979795322186*b_y[3]*jacobtot_inv[5]*hamil[16])/inFlds_e[13][1]+(1.026979795322186*jacobtot_inv[3]*b_y[5]*hamil[16])/inFlds_e[13][1]-(0.592927061281571*b_y[3]*jacobtot_inv[3]*hamil[16])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[0]*(-((0.7954951288348656*jacobtot_inv[1]*b_y[5]*hamil[7])/inFlds_e[13][1])+(0.45927932677184563*jacobtot_inv[0]*b_y[5]*hamil[7])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[1]*b_y[3]*hamil[7])/inFlds_e[13][1]-(0.26516504294495524*jacobtot_inv[0]*b_y[3]*hamil[7])/inFlds_e[13][1]-(0.7954951288348656*hamil[3]*b_y[5]*jacobtot_inv[5])/inFlds_e[13][1]+(0.45927932677184563*b_y[3]*hamil[3]*jacobtot_inv[5])/inFlds_e[13][1]+(0.45927932677184563*hamil[3]*jacobtot_inv[3]*b_y[5])/inFlds_e[13][1]-(0.26516504294495524*b_y[3]*hamil[3]*jacobtot_inv[3])/inFlds_e[13][1])*rdvpar2*rdz2-1.7320508075688772*hamil[5]*inFlds_e[9]*inFlds_e[10]*rdy2)/q_; 
  alphaL[3] = (vmap[0]*(-((1.7787811838447125*b_y[5]*jacobtot_inv[5]*hamil[18])/inFlds_e[13][1])+(1.0269797953221862*b_y[3]*jacobtot_inv[5]*hamil[18])/inFlds_e[13][1]+(1.0269797953221862*jacobtot_inv[3]*b_y[5]*hamil[18])/inFlds_e[13][1]-(0.5929270612815709*b_y[3]*jacobtot_inv[3]*hamil[18])/inFlds_e[13][1]-(1.7787811838447127*jacobtot_inv[1]*b_y[5]*hamil[16])/inFlds_e[13][1]+(1.026979795322186*jacobtot_inv[0]*b_y[5]*hamil[16])/inFlds_e[13][1]+(1.026979795322186*jacobtot_inv[1]*b_y[3]*hamil[16])/inFlds_e[13][1]-(0.592927061281571*jacobtot_inv[0]*b_y[3]*hamil[16])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[1]*(-((0.7954951288348656*b_y[5]*jacobtot_inv[5]*hamil[7])/inFlds_e[13][1])+(0.45927932677184563*b_y[3]*jacobtot_inv[5]*hamil[7])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[3]*b_y[5]*hamil[7])/inFlds_e[13][1]-(0.26516504294495524*b_y[3]*jacobtot_inv[3]*hamil[7])/inFlds_e[13][1]-(0.7954951288348656*jacobtot_inv[1]*hamil[3]*b_y[5])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[0]*hamil[3]*b_y[5])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[1]*b_y[3]*hamil[3])/inFlds_e[13][1]-(0.26516504294495524*jacobtot_inv[0]*b_y[3]*hamil[3])/inFlds_e[13][1])*rdvpar2*rdz2+1.7320508075688772*hamil[7]*inFlds_e[8]*inFlds_e[10]*rdz2)/q_; 
  alphaL[4] = (vmap[1]*(-((1.7787811838447127*b_y[5]*jacobtot_inv[5]*hamil[22])/inFlds_e[13][1])+(1.026979795322186*b_y[3]*jacobtot_inv[5]*hamil[22])/inFlds_e[13][1]+(1.026979795322186*jacobtot_inv[3]*b_y[5]*hamil[22])/inFlds_e[13][1]-(0.592927061281571*b_y[3]*jacobtot_inv[3]*hamil[22])/inFlds_e[13][1]-(1.7787811838447125*jacobtot_inv[1]*b_y[5]*hamil[19])/inFlds_e[13][1]+(1.0269797953221862*jacobtot_inv[0]*b_y[5]*hamil[19])/inFlds_e[13][1]+(1.0269797953221862*jacobtot_inv[1]*b_y[3]*hamil[19])/inFlds_e[13][1]-(0.5929270612815709*jacobtot_inv[0]*b_y[3]*hamil[19])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[0]*(-((0.7954951288348656*b_y[5]*jacobtot_inv[5]*hamil[14])/inFlds_e[13][1])+(0.45927932677184563*b_y[3]*jacobtot_inv[5]*hamil[14])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[3]*b_y[5]*hamil[14])/inFlds_e[13][1]-(0.26516504294495524*b_y[3]*jacobtot_inv[3]*hamil[14])/inFlds_e[13][1]-(0.7954951288348656*jacobtot_inv[1]*b_y[5]*hamil[10])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[0]*b_y[5]*hamil[10])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[1]*b_y[3]*hamil[10])/inFlds_e[13][1]-(0.26516504294495524*jacobtot_inv[0]*b_y[3]*hamil[10])/inFlds_e[13][1])*rdvpar2*rdz2+1.7320508075688772*inFlds_e[8]*hamil[9]*inFlds_e[10]*rdz2)/q_; 
  alphaL[7] = (vmap[0]*(-((1.7787811838447125*jacobtot_inv[1]*b_y[5]*hamil[18])/inFlds_e[13][1])+(1.0269797953221862*jacobtot_inv[0]*b_y[5]*hamil[18])/inFlds_e[13][1]+(1.0269797953221862*jacobtot_inv[1]*b_y[3]*hamil[18])/inFlds_e[13][1]-(0.5929270612815709*jacobtot_inv[0]*b_y[3]*hamil[18])/inFlds_e[13][1]-(1.7787811838447127*b_y[5]*jacobtot_inv[5]*hamil[16])/inFlds_e[13][1]+(1.026979795322186*b_y[3]*jacobtot_inv[5]*hamil[16])/inFlds_e[13][1]+(1.026979795322186*jacobtot_inv[3]*b_y[5]*hamil[16])/inFlds_e[13][1]-(0.592927061281571*b_y[3]*jacobtot_inv[3]*hamil[16])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[1]*(-((0.7954951288348656*jacobtot_inv[1]*b_y[5]*hamil[7])/inFlds_e[13][1])+(0.45927932677184563*jacobtot_inv[0]*b_y[5]*hamil[7])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[1]*b_y[3]*hamil[7])/inFlds_e[13][1]-(0.26516504294495524*jacobtot_inv[0]*b_y[3]*hamil[7])/inFlds_e[13][1]-(0.7954951288348656*hamil[3]*b_y[5]*jacobtot_inv[5])/inFlds_e[13][1]+(0.45927932677184563*b_y[3]*hamil[3]*jacobtot_inv[5])/inFlds_e[13][1]+(0.45927932677184563*hamil[3]*jacobtot_inv[3]*b_y[5])/inFlds_e[13][1]-(0.26516504294495524*b_y[3]*hamil[3]*jacobtot_inv[3])/inFlds_e[13][1])*rdvpar2*rdz2)/q_; 
  alphaL[9] = (vmap[1]*(-((1.7787811838447127*jacobtot_inv[1]*b_y[5]*hamil[22])/inFlds_e[13][1])+(1.026979795322186*jacobtot_inv[0]*b_y[5]*hamil[22])/inFlds_e[13][1]+(1.026979795322186*jacobtot_inv[1]*b_y[3]*hamil[22])/inFlds_e[13][1]-(0.592927061281571*jacobtot_inv[0]*b_y[3]*hamil[22])/inFlds_e[13][1]-(1.7787811838447125*b_y[5]*jacobtot_inv[5]*hamil[19])/inFlds_e[13][1]+(1.0269797953221862*b_y[3]*jacobtot_inv[5]*hamil[19])/inFlds_e[13][1]+(1.0269797953221862*jacobtot_inv[3]*b_y[5]*hamil[19])/inFlds_e[13][1]-(0.5929270612815709*b_y[3]*jacobtot_inv[3]*hamil[19])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[0]*(-((0.7954951288348656*jacobtot_inv[1]*b_y[5]*hamil[14])/inFlds_e[13][1])+(0.45927932677184563*jacobtot_inv[0]*b_y[5]*hamil[14])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[1]*b_y[3]*hamil[14])/inFlds_e[13][1]-(0.26516504294495524*jacobtot_inv[0]*b_y[3]*hamil[14])/inFlds_e[13][1]-(0.7954951288348656*b_y[5]*jacobtot_inv[5]*hamil[10])/inFlds_e[13][1]+(0.45927932677184563*b_y[3]*jacobtot_inv[5]*hamil[10])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[3]*b_y[5]*hamil[10])/inFlds_e[13][1]-(0.26516504294495524*b_y[3]*jacobtot_inv[3]*hamil[10])/inFlds_e[13][1])*rdvpar2*rdz2)/q_; 
  alphaL[10] = (vmap[0]*(-((1.7787811838447127*b_y[5]*jacobtot_inv[5]*hamil[22])/inFlds_e[13][1])+(1.026979795322186*b_y[3]*jacobtot_inv[5]*hamil[22])/inFlds_e[13][1]+(1.026979795322186*jacobtot_inv[3]*b_y[5]*hamil[22])/inFlds_e[13][1]-(0.592927061281571*b_y[3]*jacobtot_inv[3]*hamil[22])/inFlds_e[13][1]-(1.7787811838447125*jacobtot_inv[1]*b_y[5]*hamil[19])/inFlds_e[13][1]+(1.0269797953221862*jacobtot_inv[0]*b_y[5]*hamil[19])/inFlds_e[13][1]+(1.0269797953221862*jacobtot_inv[1]*b_y[3]*hamil[19])/inFlds_e[13][1]-(0.5929270612815709*jacobtot_inv[0]*b_y[3]*hamil[19])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[1]*(-((0.7954951288348656*b_y[5]*jacobtot_inv[5]*hamil[14])/inFlds_e[13][1])+(0.45927932677184563*b_y[3]*jacobtot_inv[5]*hamil[14])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[3]*b_y[5]*hamil[14])/inFlds_e[13][1]-(0.26516504294495524*b_y[3]*jacobtot_inv[3]*hamil[14])/inFlds_e[13][1]-(0.7954951288348656*jacobtot_inv[1]*b_y[5]*hamil[10])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[0]*b_y[5]*hamil[10])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[1]*b_y[3]*hamil[10])/inFlds_e[13][1]-(0.26516504294495524*jacobtot_inv[0]*b_y[3]*hamil[10])/inFlds_e[13][1])*rdvpar2*rdz2+1.7320508075688772*inFlds_e[8]*inFlds_e[10]*hamil[14]*rdz2)/q_; 
  alphaL[14] = (vmap[0]*(-((1.7787811838447127*jacobtot_inv[1]*b_y[5]*hamil[22])/inFlds_e[13][1])+(1.026979795322186*jacobtot_inv[0]*b_y[5]*hamil[22])/inFlds_e[13][1]+(1.026979795322186*jacobtot_inv[1]*b_y[3]*hamil[22])/inFlds_e[13][1]-(0.592927061281571*jacobtot_inv[0]*b_y[3]*hamil[22])/inFlds_e[13][1]-(1.7787811838447125*b_y[5]*jacobtot_inv[5]*hamil[19])/inFlds_e[13][1]+(1.0269797953221862*b_y[3]*jacobtot_inv[5]*hamil[19])/inFlds_e[13][1]+(1.0269797953221862*jacobtot_inv[3]*b_y[5]*hamil[19])/inFlds_e[13][1]-(0.5929270612815709*b_y[3]*jacobtot_inv[3]*hamil[19])/inFlds_e[13][1])*rdvpar2*rdz2+vmap[1]*(-((0.7954951288348656*jacobtot_inv[1]*b_y[5]*hamil[14])/inFlds_e[13][1])+(0.45927932677184563*jacobtot_inv[0]*b_y[5]*hamil[14])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[1]*b_y[3]*hamil[14])/inFlds_e[13][1]-(0.26516504294495524*jacobtot_inv[0]*b_y[3]*hamil[14])/inFlds_e[13][1]-(0.7954951288348656*b_y[5]*jacobtot_inv[5]*hamil[10])/inFlds_e[13][1]+(0.45927932677184563*b_y[3]*jacobtot_inv[5]*hamil[10])/inFlds_e[13][1]+(0.45927932677184563*jacobtot_inv[3]*b_y[5]*hamil[10])/inFlds_e[13][1]-(0.26516504294495524*b_y[3]*jacobtot_inv[3]*hamil[10])/inFlds_e[13][1])*rdvpar2*rdz2)/q_; 
  alphaL[16] = (vmap[1]*(-((1.590990257669731*b_y[5]*jacobtot_inv[5]*hamil[18])/inFlds_e[13][1])+(0.9185586535436916*b_y[3]*jacobtot_inv[5]*hamil[18])/inFlds_e[13][1]+(0.9185586535436916*jacobtot_inv[3]*b_y[5]*hamil[18])/inFlds_e[13][1]-(0.5303300858899104*b_y[3]*jacobtot_inv[3]*hamil[18])/inFlds_e[13][1]-(1.5909902576697312*jacobtot_inv[1]*b_y[5]*hamil[16])/inFlds_e[13][1]+(0.9185586535436913*jacobtot_inv[0]*b_y[5]*hamil[16])/inFlds_e[13][1]+(0.9185586535436913*jacobtot_inv[1]*b_y[3]*hamil[16])/inFlds_e[13][1]-(0.5303300858899105*jacobtot_inv[0]*b_y[3]*hamil[16])/inFlds_e[13][1])*rdvpar2*rdz2+1.7320508075688774*inFlds_e[8]*inFlds_e[10]*hamil[18]*rdz2)/q_; 
  alphaL[18] = (vmap[1]*(-((1.5909902576697312*jacobtot_inv[1]*b_y[5]*hamil[18])/inFlds_e[13][1])+(0.9185586535436913*jacobtot_inv[0]*b_y[5]*hamil[18])/inFlds_e[13][1]+(0.9185586535436913*jacobtot_inv[1]*b_y[3]*hamil[18])/inFlds_e[13][1]-(0.5303300858899105*jacobtot_inv[0]*b_y[3]*hamil[18])/inFlds_e[13][1]-(1.590990257669731*b_y[5]*jacobtot_inv[5]*hamil[16])/inFlds_e[13][1]+(0.9185586535436916*b_y[3]*jacobtot_inv[5]*hamil[16])/inFlds_e[13][1]+(0.9185586535436916*jacobtot_inv[3]*b_y[5]*hamil[16])/inFlds_e[13][1]-(0.5303300858899104*b_y[3]*jacobtot_inv[3]*hamil[16])/inFlds_e[13][1])*rdvpar2*rdz2)/q_; 
  alphaL[19] = (vmap[1]*(-((1.590990257669731*b_y[5]*jacobtot_inv[5]*hamil[22])/inFlds_e[13][1])+(0.9185586535436916*b_y[3]*jacobtot_inv[5]*hamil[22])/inFlds_e[13][1]+(0.9185586535436916*jacobtot_inv[3]*b_y[5]*hamil[22])/inFlds_e[13][1]-(0.5303300858899104*b_y[3]*jacobtot_inv[3]*hamil[22])/inFlds_e[13][1]-(1.5909902576697312*jacobtot_inv[1]*b_y[5]*hamil[19])/inFlds_e[13][1]+(0.9185586535436913*jacobtot_inv[0]*b_y[5]*hamil[19])/inFlds_e[13][1]+(0.9185586535436913*jacobtot_inv[1]*b_y[3]*hamil[19])/inFlds_e[13][1]-(0.5303300858899105*jacobtot_inv[0]*b_y[3]*hamil[19])/inFlds_e[13][1])*rdvpar2*rdz2+1.7320508075688774*inFlds_e[8]*inFlds_e[10]*hamil[22]*rdz2)/q_; 
  alphaL[22] = (vmap[1]*(-((1.5909902576697312*jacobtot_inv[1]*b_y[5]*hamil[22])/inFlds_e[13][1])+(0.9185586535436913*jacobtot_inv[0]*b_y[5]*hamil[22])/inFlds_e[13][1]+(0.9185586535436913*jacobtot_inv[1]*b_y[3]*hamil[22])/inFlds_e[13][1]-(0.5303300858899105*jacobtot_inv[0]*b_y[3]*hamil[22])/inFlds_e[13][1]-(1.590990257669731*b_y[5]*jacobtot_inv[5]*hamil[19])/inFlds_e[13][1]+(0.9185586535436916*b_y[3]*jacobtot_inv[5]*hamil[19])/inFlds_e[13][1]+(0.9185586535436916*jacobtot_inv[3]*b_y[5]*hamil[19])/inFlds_e[13][1]-(0.5303300858899104*b_y[3]*jacobtot_inv[3]*hamil[19])/inFlds_e[13][1])*rdvpar2*rdz2)/q_; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.22360679774997858*alphaL[22]-0.22360679774997858*(alphaL[19]+alphaL[18])+0.22360679774997858*alphaL[16]-0.33541019662496785*alphaL[14]+0.33541019662496785*alphaL[10]+0.25*alphaL[9]+0.33541019662496785*alphaL[7]-0.25*alphaL[4]-0.33541019662496785*alphaL[3]-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if (-(0.2795084971874732*alphaL[22])+0.2795084971874732*(alphaL[19]+alphaL[18])-0.2795084971874732*alphaL[16]+0.25*alphaL[9]-0.25*(alphaL[4]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.22360679774997858*alphaL[22]-0.22360679774997858*(alphaL[19]+alphaL[18])+0.22360679774997858*alphaL[16]+0.33541019662496785*alphaL[14]-0.33541019662496785*alphaL[10]+0.25*alphaL[9]-0.33541019662496785*alphaL[7]-0.25*alphaL[4]+0.33541019662496785*alphaL[3]-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*alphaL[22])+0.22360679774997858*alphaL[19]-0.22360679774997858*alphaL[18]+0.22360679774997858*alphaL[16]+0.33541019662496785*alphaL[14]-0.33541019662496785*alphaL[10]-0.25*alphaL[9]+0.33541019662496785*alphaL[7]+0.25*alphaL[4]-0.33541019662496785*alphaL[3]-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL[22]-0.2795084971874732*alphaL[19]+0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]-0.25*alphaL[9]+0.25*alphaL[4]-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*alphaL[22])+0.22360679774997858*alphaL[19]-0.22360679774997858*alphaL[18]+0.22360679774997858*alphaL[16]-0.33541019662496785*alphaL[14]+0.33541019662496785*alphaL[10]-0.25*alphaL[9]-0.33541019662496785*alphaL[7]+0.25*alphaL[4]+0.33541019662496785*alphaL[3]-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*(alphaL[22]+alphaL[19]))+0.22360679774997858*(alphaL[18]+alphaL[16])+0.33541019662496785*(alphaL[14]+alphaL[10])-0.25*alphaL[9]-0.33541019662496785*alphaL[7]-0.25*alphaL[4]-0.33541019662496785*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*(alphaL[22]+alphaL[19])-0.2795084971874732*(alphaL[18]+alphaL[16])-0.25*(alphaL[9]+alphaL[4])+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*(alphaL[22]+alphaL[19]))+0.22360679774997858*(alphaL[18]+alphaL[16])-0.33541019662496785*(alphaL[14]+alphaL[10])-0.25*alphaL[9]+0.33541019662496785*alphaL[7]-0.25*alphaL[4]+0.33541019662496785*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.22360679774997858*(alphaL[22]+alphaL[19]+alphaL[18]+alphaL[16])-0.33541019662496785*(alphaL[14]+alphaL[10])+0.25*alphaL[9]-0.33541019662496785*alphaL[7]+0.25*alphaL[4]-0.33541019662496785*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.2795084971874732*(alphaL[22]+alphaL[19]+alphaL[18]+alphaL[16]))+0.25*(alphaL[9]+alphaL[4]+alphaL[2])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.22360679774997858*(alphaL[22]+alphaL[19]+alphaL[18]+alphaL[16])+0.33541019662496785*(alphaL[14]+alphaL[10])+0.25*alphaL[9]+0.33541019662496785*alphaL[7]+0.25*alphaL[4]+0.33541019662496785*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.22360679774997858*alphaL[22]-0.22360679774997858*(alphaL[19]+alphaL[18])+0.22360679774997858*alphaL[16]-0.33541019662496785*alphaL[14]+0.33541019662496785*alphaL[10]+0.25*alphaL[9]+0.33541019662496785*alphaL[7]-0.25*alphaL[4]-0.33541019662496785*alphaL[3]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.2795084971874732*alphaL[22])+0.2795084971874732*(alphaL[19]+alphaL[18])-0.2795084971874732*alphaL[16]+0.25*alphaL[9]-0.25*(alphaL[4]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.22360679774997858*alphaL[22]-0.22360679774997858*(alphaL[19]+alphaL[18])+0.22360679774997858*alphaL[16]+0.33541019662496785*alphaL[14]-0.33541019662496785*alphaL[10]+0.25*alphaL[9]-0.33541019662496785*alphaL[7]-0.25*alphaL[4]+0.33541019662496785*alphaL[3]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*alphaL[22])+0.22360679774997858*alphaL[19]-0.22360679774997858*alphaL[18]+0.22360679774997858*alphaL[16]+0.33541019662496785*alphaL[14]-0.33541019662496785*alphaL[10]-0.25*alphaL[9]+0.33541019662496785*alphaL[7]+0.25*alphaL[4]-0.33541019662496785*alphaL[3]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL[22]-0.2795084971874732*alphaL[19]+0.2795084971874732*alphaL[18]-0.2795084971874732*alphaL[16]-0.25*alphaL[9]+0.25*alphaL[4]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[16] = 1.0; 
  else  
    sgn_alpha_surfL[16] = -1.0; 
  
  if (sgn_alpha_surfL[16] == sgn_alpha_surfL[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*alphaL[22])+0.22360679774997858*alphaL[19]-0.22360679774997858*alphaL[18]+0.22360679774997858*alphaL[16]-0.33541019662496785*alphaL[14]+0.33541019662496785*alphaL[10]-0.25*alphaL[9]-0.33541019662496785*alphaL[7]+0.25*alphaL[4]+0.33541019662496785*alphaL[3]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[17] = 1.0; 
  else  
    sgn_alpha_surfL[17] = -1.0; 
  
  if (sgn_alpha_surfL[17] == sgn_alpha_surfL[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*(alphaL[22]+alphaL[19]))+0.22360679774997858*(alphaL[18]+alphaL[16])+0.33541019662496785*(alphaL[14]+alphaL[10])-0.25*alphaL[9]-0.33541019662496785*alphaL[7]-0.25*alphaL[4]-0.33541019662496785*alphaL[3]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[18] = 1.0; 
  else  
    sgn_alpha_surfL[18] = -1.0; 
  
  if (sgn_alpha_surfL[18] == sgn_alpha_surfL[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*(alphaL[22]+alphaL[19])-0.2795084971874732*(alphaL[18]+alphaL[16])-0.25*(alphaL[9]+alphaL[4])+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[19] = 1.0; 
  else  
    sgn_alpha_surfL[19] = -1.0; 
  
  if (sgn_alpha_surfL[19] == sgn_alpha_surfL[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*(alphaL[22]+alphaL[19]))+0.22360679774997858*(alphaL[18]+alphaL[16])-0.33541019662496785*(alphaL[14]+alphaL[10])-0.25*alphaL[9]+0.33541019662496785*alphaL[7]-0.25*alphaL[4]+0.33541019662496785*alphaL[3]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[20] = 1.0; 
  else  
    sgn_alpha_surfL[20] = -1.0; 
  
  if (sgn_alpha_surfL[20] == sgn_alpha_surfL[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.22360679774997858*(alphaL[22]+alphaL[19]+alphaL[18]+alphaL[16])-0.33541019662496785*(alphaL[14]+alphaL[10])+0.25*alphaL[9]-0.33541019662496785*alphaL[7]+0.25*alphaL[4]-0.33541019662496785*alphaL[3]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[21] = 1.0; 
  else  
    sgn_alpha_surfL[21] = -1.0; 
  
  if (sgn_alpha_surfL[21] == sgn_alpha_surfL[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL[9]+alphaL[4]+alphaL[2]+alphaL[1]+alphaL[0])-0.2795084971874732*(alphaL[22]+alphaL[19]+alphaL[18]+alphaL[16]) > 0.) 
    sgn_alpha_surfL[22] = 1.0; 
  else  
    sgn_alpha_surfL[22] = -1.0; 
  
  if (sgn_alpha_surfL[22] == sgn_alpha_surfL[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.22360679774997858*(alphaL[22]+alphaL[19]+alphaL[18]+alphaL[16])+0.33541019662496785*(alphaL[14]+alphaL[10])+0.25*alphaL[9]+0.33541019662496785*alphaL[7]+0.25*alphaL[4]+0.33541019662496785*alphaL[3]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[23] = 1.0; 
  else  
    sgn_alpha_surfL[23] = -1.0; 
  
  if (sgn_alpha_surfL[23] == sgn_alpha_surfL[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
