#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH int gyrokinetic_alpha_no_by_surfy_3x2v_ser_p1(const double *w, const double *dxv, const double *vmap, const double *vmapSq,
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
  const double *b_y_surf = &b_i_surf[4];
  const double *b_z = &b_i[16];
  const double *b_z_surf = &b_i_surf[8];

  double hamil[24] = {0.}; 
  hamil[0] = 1.4142135623730951*(phi[0]*q_+vmapSq[0]*m_+bmag_surf[0]*vmap[2])-2.4494897427831783*phi[2]*q_; 
  hamil[1] = 1.4142135623730951*(phi[1]*q_+bmag_surf[1]*vmap[2])-2.4494897427831783*phi[4]*q_; 
  hamil[2] = 1.4142135623730951*(phi[3]*q_+bmag_surf[2]*vmap[2])-2.4494897427831783*phi[6]*q_; 
  hamil[3] = 1.4142135623730951*vmapSq[1]*m_; 
  hamil[4] = 1.4142135623730951*bmag_surf[0]*vmap[3]; 
  hamil[5] = 1.4142135623730951*(phi[5]*q_+vmap[2]*bmag_surf[3])-2.4494897427831783*phi[7]*q_; 
  hamil[8] = 1.4142135623730951*bmag_surf[1]*vmap[3]; 
  hamil[9] = 1.4142135623730951*bmag_surf[2]*vmap[3]; 
  hamil[12] = 1.4142135623730951*bmag_surf[3]*vmap[3]; 
  hamil[16] = 1.4142135623730951*vmapSq[2]*m_; 

  double *alphaL = &alpha_surf[24];
  double *sgn_alpha_surfL = &sgn_alpha_surf[24];
  alphaL[0] = ((0.4841229182759271*jacobtot_inv[1]*b_x[5]*hamil[16]+0.4841229182759271*jacobtot_inv[0]*b_x[3]*hamil[16]-0.4330127018922193*b_x_surf[2]*jacobtot_inv_surf[3]*hamil[5]-0.4330127018922193*jacobtot_inv_surf[2]*b_x_surf[3]*hamil[5]-0.4330127018922193*b_x_surf[0]*jacobtot_inv_surf[1]*hamil[5]-0.4330127018922193*jacobtot_inv_surf[0]*b_x_surf[1]*hamil[5]-0.4330127018922193*hamil[2]*b_x_surf[3]*jacobtot_inv_surf[3]-0.4330127018922193*b_x_surf[2]*hamil[2]*jacobtot_inv_surf[2]-0.4330127018922193*b_x_surf[1]*jacobtot_inv_surf[1]*hamil[2]-0.4330127018922193*b_x_surf[0]*jacobtot_inv_surf[0]*hamil[2])*rdz2+(vmap[0]*(0.21650635094610965*jacobtot_inv[1]*hamil[3]*b_x[5]+0.21650635094610965*jacobtot_inv[0]*b_x[3]*hamil[3])*rdz2)/vmap[1]+(-(0.4841229182759271*jacobtot_inv[3]*b_z[5]*hamil[16])-0.4841229182759271*jacobtot_inv[0]*b_z[1]*hamil[16]+0.4330127018922193*b_z_surf[1]*jacobtot_inv_surf[3]*hamil[5]+0.4330127018922193*jacobtot_inv_surf[1]*b_z_surf[3]*hamil[5]+0.4330127018922193*b_z_surf[0]*jacobtot_inv_surf[2]*hamil[5]+0.4330127018922193*jacobtot_inv_surf[0]*b_z_surf[2]*hamil[5]+0.4330127018922193*hamil[1]*b_z_surf[3]*jacobtot_inv_surf[3]+0.4330127018922193*hamil[1]*b_z_surf[2]*jacobtot_inv_surf[2]+0.4330127018922193*b_z_surf[1]*hamil[1]*jacobtot_inv_surf[1]+0.4330127018922193*b_z_surf[0]*jacobtot_inv_surf[0]*hamil[1])*rdx2+(vmap[0]*(-(0.21650635094610965*hamil[3]*jacobtot_inv[3]*b_z[5])-0.21650635094610965*jacobtot_inv[0]*b_z[1]*hamil[3])*rdx2)/vmap[1])/q_; 
  alphaL[1] = ((0.4841229182759271*jacobtot_inv[0]*b_x[5]*hamil[16]+0.4841229182759271*jacobtot_inv[1]*b_x[3]*hamil[16]-0.7794228634059945*b_x_surf[3]*jacobtot_inv_surf[3]*hamil[5]-0.4330127018922193*b_x_surf[2]*jacobtot_inv_surf[2]*hamil[5]-0.7794228634059945*b_x_surf[1]*jacobtot_inv_surf[1]*hamil[5]-0.4330127018922193*b_x_surf[0]*jacobtot_inv_surf[0]*hamil[5]-0.4330127018922193*b_x_surf[2]*hamil[2]*jacobtot_inv_surf[3]-0.4330127018922193*hamil[2]*jacobtot_inv_surf[2]*b_x_surf[3]-0.4330127018922193*b_x_surf[0]*jacobtot_inv_surf[1]*hamil[2]-0.4330127018922193*jacobtot_inv_surf[0]*b_x_surf[1]*hamil[2])*rdz2+(vmap[0]*(0.21650635094610965*jacobtot_inv[0]*hamil[3]*b_x[5]+0.21650635094610965*jacobtot_inv[1]*b_x[3]*hamil[3])*rdz2)/vmap[1]+(-(0.4841229182759271*b_z[5]*jacobtot_inv[5]*hamil[16])-0.4841229182759271*b_z[1]*jacobtot_inv[1]*hamil[16]+0.4330127018922193*b_z_surf[0]*jacobtot_inv_surf[3]*hamil[5]+0.4330127018922193*jacobtot_inv_surf[0]*b_z_surf[3]*hamil[5]+0.4330127018922193*b_z_surf[1]*jacobtot_inv_surf[2]*hamil[5]+0.4330127018922193*jacobtot_inv_surf[1]*b_z_surf[2]*hamil[5]+0.4330127018922193*hamil[1]*b_z_surf[2]*jacobtot_inv_surf[3]+0.4330127018922193*hamil[1]*jacobtot_inv_surf[2]*b_z_surf[3]+0.4330127018922193*b_z_surf[0]*hamil[1]*jacobtot_inv_surf[1]+0.4330127018922193*jacobtot_inv_surf[0]*b_z_surf[1]*hamil[1])*rdx2+(vmap[0]*(-(0.21650635094610965*hamil[3]*b_z[5]*jacobtot_inv[5])-0.21650635094610965*b_z[1]*jacobtot_inv[1]*hamil[3])*rdx2)/vmap[1])/q_; 
  alphaL[2] = ((0.4841229182759271*b_x[5]*jacobtot_inv[5]*hamil[16]+0.4841229182759271*b_x[3]*jacobtot_inv[3]*hamil[16]-0.4330127018922193*b_x_surf[0]*jacobtot_inv_surf[3]*hamil[5]-0.4330127018922193*jacobtot_inv_surf[0]*b_x_surf[3]*hamil[5]-0.4330127018922193*b_x_surf[1]*jacobtot_inv_surf[2]*hamil[5]-0.4330127018922193*jacobtot_inv_surf[1]*b_x_surf[2]*hamil[5]-0.4330127018922193*b_x_surf[1]*hamil[2]*jacobtot_inv_surf[3]-0.4330127018922193*jacobtot_inv_surf[1]*hamil[2]*b_x_surf[3]-0.4330127018922193*b_x_surf[0]*hamil[2]*jacobtot_inv_surf[2]-0.4330127018922193*jacobtot_inv_surf[0]*b_x_surf[2]*hamil[2])*rdz2+(vmap[0]*(0.21650635094610965*hamil[3]*b_x[5]*jacobtot_inv[5]+0.21650635094610965*b_x[3]*hamil[3]*jacobtot_inv[3])*rdz2)/vmap[1]+(-(0.4841229182759271*jacobtot_inv[0]*b_z[5]*hamil[16])-0.4841229182759271*b_z[1]*jacobtot_inv[3]*hamil[16]+0.7794228634059945*b_z_surf[3]*jacobtot_inv_surf[3]*hamil[5]+0.7794228634059945*b_z_surf[2]*jacobtot_inv_surf[2]*hamil[5]+0.4330127018922193*b_z_surf[1]*jacobtot_inv_surf[1]*hamil[5]+0.4330127018922193*b_z_surf[0]*jacobtot_inv_surf[0]*hamil[5]+0.4330127018922193*b_z_surf[1]*hamil[1]*jacobtot_inv_surf[3]+0.4330127018922193*hamil[1]*jacobtot_inv_surf[1]*b_z_surf[3]+0.4330127018922193*b_z_surf[0]*hamil[1]*jacobtot_inv_surf[2]+0.4330127018922193*jacobtot_inv_surf[0]*hamil[1]*b_z_surf[2])*rdx2+(vmap[0]*(-(0.21650635094610965*jacobtot_inv[0]*hamil[3]*b_z[5])-0.21650635094610965*b_z[1]*hamil[3]*jacobtot_inv[3])*rdx2)/vmap[1])/q_; 
  alphaL[3] = ((vmap[0]*(0.4841229182759271*jacobtot_inv[1]*b_x[5]*hamil[16]+0.4841229182759271*jacobtot_inv[0]*b_x[3]*hamil[16])*rdz2)/vmap[1]+(0.21650635094610965*jacobtot_inv[1]*hamil[3]*b_x[5]+0.21650635094610965*jacobtot_inv[0]*b_x[3]*hamil[3])*rdz2+(vmap[0]*(-(0.4841229182759271*jacobtot_inv[3]*b_z[5]*hamil[16])-0.4841229182759271*jacobtot_inv[0]*b_z[1]*hamil[16])*rdx2)/vmap[1]+(-(0.21650635094610965*hamil[3]*jacobtot_inv[3]*b_z[5])-0.21650635094610965*jacobtot_inv[0]*b_z[1]*hamil[3])*rdx2)/q_; 
  alphaL[4] = ((-(0.4330127018922193*b_x_surf[2]*jacobtot_inv_surf[3]*hamil[12])-0.4330127018922193*jacobtot_inv_surf[2]*b_x_surf[3]*hamil[12]-0.4330127018922193*b_x_surf[0]*jacobtot_inv_surf[1]*hamil[12]-0.4330127018922193*jacobtot_inv_surf[0]*b_x_surf[1]*hamil[12]-0.4330127018922193*b_x_surf[3]*jacobtot_inv_surf[3]*hamil[9]-0.4330127018922193*b_x_surf[2]*jacobtot_inv_surf[2]*hamil[9]-0.4330127018922193*b_x_surf[1]*jacobtot_inv_surf[1]*hamil[9]-0.4330127018922193*b_x_surf[0]*jacobtot_inv_surf[0]*hamil[9])*rdz2+(0.4330127018922193*b_z_surf[1]*jacobtot_inv_surf[3]*hamil[12]+0.4330127018922193*jacobtot_inv_surf[1]*b_z_surf[3]*hamil[12]+0.4330127018922193*b_z_surf[0]*jacobtot_inv_surf[2]*hamil[12]+0.4330127018922193*jacobtot_inv_surf[0]*b_z_surf[2]*hamil[12]+0.4330127018922193*b_z_surf[3]*jacobtot_inv_surf[3]*hamil[8]+0.4330127018922193*b_z_surf[2]*jacobtot_inv_surf[2]*hamil[8]+0.4330127018922193*b_z_surf[1]*jacobtot_inv_surf[1]*hamil[8]+0.4330127018922193*b_z_surf[0]*jacobtot_inv_surf[0]*hamil[8])*rdx2)/q_; 
  alphaL[5] = ((0.4841229182759271*b_x[3]*jacobtot_inv[5]*hamil[16]+0.4841229182759271*jacobtot_inv[3]*b_x[5]*hamil[16]-0.7794228634059945*b_x_surf[1]*jacobtot_inv_surf[3]*hamil[5]-0.7794228634059945*jacobtot_inv_surf[1]*b_x_surf[3]*hamil[5]-0.4330127018922193*b_x_surf[0]*jacobtot_inv_surf[2]*hamil[5]-0.4330127018922193*jacobtot_inv_surf[0]*b_x_surf[2]*hamil[5]-0.4330127018922193*b_x_surf[0]*hamil[2]*jacobtot_inv_surf[3]-0.4330127018922193*jacobtot_inv_surf[0]*hamil[2]*b_x_surf[3]-0.4330127018922193*b_x_surf[1]*hamil[2]*jacobtot_inv_surf[2]-0.4330127018922193*jacobtot_inv_surf[1]*b_x_surf[2]*hamil[2])*rdz2+(vmap[0]*(0.21650635094610965*b_x[3]*hamil[3]*jacobtot_inv[5]+0.21650635094610965*hamil[3]*jacobtot_inv[3]*b_x[5])*rdz2)/vmap[1]+(-(0.4841229182759271*b_z[1]*jacobtot_inv[5]*hamil[16])-0.4841229182759271*jacobtot_inv[1]*b_z[5]*hamil[16]+0.7794228634059945*b_z_surf[2]*jacobtot_inv_surf[3]*hamil[5]+0.7794228634059945*jacobtot_inv_surf[2]*b_z_surf[3]*hamil[5]+0.4330127018922193*b_z_surf[0]*jacobtot_inv_surf[1]*hamil[5]+0.4330127018922193*jacobtot_inv_surf[0]*b_z_surf[1]*hamil[5]+0.4330127018922193*b_z_surf[0]*hamil[1]*jacobtot_inv_surf[3]+0.4330127018922193*jacobtot_inv_surf[0]*hamil[1]*b_z_surf[3]+0.4330127018922193*b_z_surf[1]*hamil[1]*jacobtot_inv_surf[2]+0.4330127018922193*hamil[1]*jacobtot_inv_surf[1]*b_z_surf[2])*rdx2+(vmap[0]*(-(0.21650635094610965*b_z[1]*hamil[3]*jacobtot_inv[5])-0.21650635094610965*jacobtot_inv[1]*hamil[3]*b_z[5])*rdx2)/vmap[1])/q_; 
  alphaL[6] = ((vmap[0]*(0.4841229182759271*jacobtot_inv[0]*b_x[5]*hamil[16]+0.4841229182759271*jacobtot_inv[1]*b_x[3]*hamil[16])*rdz2)/vmap[1]+(0.21650635094610965*jacobtot_inv[0]*hamil[3]*b_x[5]+0.21650635094610965*jacobtot_inv[1]*b_x[3]*hamil[3])*rdz2+(vmap[0]*(-(0.4841229182759271*b_z[5]*jacobtot_inv[5]*hamil[16])-0.4841229182759271*b_z[1]*jacobtot_inv[1]*hamil[16])*rdx2)/vmap[1]+(-(0.21650635094610965*hamil[3]*b_z[5]*jacobtot_inv[5])-0.21650635094610965*b_z[1]*jacobtot_inv[1]*hamil[3])*rdx2)/q_; 
  alphaL[7] = ((vmap[0]*(0.4841229182759271*b_x[5]*jacobtot_inv[5]*hamil[16]+0.4841229182759271*b_x[3]*jacobtot_inv[3]*hamil[16])*rdz2)/vmap[1]+(0.21650635094610965*hamil[3]*b_x[5]*jacobtot_inv[5]+0.21650635094610965*b_x[3]*hamil[3]*jacobtot_inv[3])*rdz2+(vmap[0]*(-(0.4841229182759271*jacobtot_inv[0]*b_z[5]*hamil[16])-0.4841229182759271*b_z[1]*jacobtot_inv[3]*hamil[16])*rdx2)/vmap[1]+(-(0.21650635094610965*jacobtot_inv[0]*hamil[3]*b_z[5])-0.21650635094610965*b_z[1]*hamil[3]*jacobtot_inv[3])*rdx2)/q_; 
  alphaL[8] = ((-(0.7794228634059945*b_x_surf[3]*jacobtot_inv_surf[3]*hamil[12])-0.4330127018922193*b_x_surf[2]*jacobtot_inv_surf[2]*hamil[12]-0.7794228634059945*b_x_surf[1]*jacobtot_inv_surf[1]*hamil[12]-0.4330127018922193*b_x_surf[0]*jacobtot_inv_surf[0]*hamil[12]-0.4330127018922193*b_x_surf[2]*jacobtot_inv_surf[3]*hamil[9]-0.4330127018922193*jacobtot_inv_surf[2]*b_x_surf[3]*hamil[9]-0.4330127018922193*b_x_surf[0]*jacobtot_inv_surf[1]*hamil[9]-0.4330127018922193*jacobtot_inv_surf[0]*b_x_surf[1]*hamil[9])*rdz2+(0.4330127018922193*b_z_surf[0]*jacobtot_inv_surf[3]*hamil[12]+0.4330127018922193*jacobtot_inv_surf[0]*b_z_surf[3]*hamil[12]+0.4330127018922193*b_z_surf[1]*jacobtot_inv_surf[2]*hamil[12]+0.4330127018922193*jacobtot_inv_surf[1]*b_z_surf[2]*hamil[12]+0.4330127018922193*b_z_surf[2]*jacobtot_inv_surf[3]*hamil[8]+0.4330127018922193*jacobtot_inv_surf[2]*b_z_surf[3]*hamil[8]+0.4330127018922193*b_z_surf[0]*jacobtot_inv_surf[1]*hamil[8]+0.4330127018922193*jacobtot_inv_surf[0]*b_z_surf[1]*hamil[8])*rdx2)/q_; 
  alphaL[9] = ((-(0.4330127018922193*b_x_surf[0]*jacobtot_inv_surf[3]*hamil[12])-0.4330127018922193*jacobtot_inv_surf[0]*b_x_surf[3]*hamil[12]-0.4330127018922193*b_x_surf[1]*jacobtot_inv_surf[2]*hamil[12]-0.4330127018922193*jacobtot_inv_surf[1]*b_x_surf[2]*hamil[12]-0.4330127018922193*b_x_surf[1]*jacobtot_inv_surf[3]*hamil[9]-0.4330127018922193*jacobtot_inv_surf[1]*b_x_surf[3]*hamil[9]-0.4330127018922193*b_x_surf[0]*jacobtot_inv_surf[2]*hamil[9]-0.4330127018922193*jacobtot_inv_surf[0]*b_x_surf[2]*hamil[9])*rdz2+(0.7794228634059945*b_z_surf[3]*jacobtot_inv_surf[3]*hamil[12]+0.7794228634059945*b_z_surf[2]*jacobtot_inv_surf[2]*hamil[12]+0.4330127018922193*b_z_surf[1]*jacobtot_inv_surf[1]*hamil[12]+0.4330127018922193*b_z_surf[0]*jacobtot_inv_surf[0]*hamil[12]+0.4330127018922193*b_z_surf[1]*jacobtot_inv_surf[3]*hamil[8]+0.4330127018922193*jacobtot_inv_surf[1]*b_z_surf[3]*hamil[8]+0.4330127018922193*b_z_surf[0]*jacobtot_inv_surf[2]*hamil[8]+0.4330127018922193*jacobtot_inv_surf[0]*b_z_surf[2]*hamil[8])*rdx2)/q_; 
  alphaL[11] = ((vmap[0]*(0.4841229182759271*b_x[3]*jacobtot_inv[5]*hamil[16]+0.4841229182759271*jacobtot_inv[3]*b_x[5]*hamil[16])*rdz2)/vmap[1]+(0.21650635094610965*b_x[3]*hamil[3]*jacobtot_inv[5]+0.21650635094610965*hamil[3]*jacobtot_inv[3]*b_x[5])*rdz2+(vmap[0]*(-(0.4841229182759271*b_z[1]*jacobtot_inv[5]*hamil[16])-0.4841229182759271*jacobtot_inv[1]*b_z[5]*hamil[16])*rdx2)/vmap[1]+(-(0.21650635094610965*b_z[1]*hamil[3]*jacobtot_inv[5])-0.21650635094610965*jacobtot_inv[1]*hamil[3]*b_z[5])*rdx2)/q_; 
  alphaL[12] = ((-(0.7794228634059945*b_x_surf[1]*jacobtot_inv_surf[3]*hamil[12])-0.7794228634059945*jacobtot_inv_surf[1]*b_x_surf[3]*hamil[12]-0.4330127018922193*b_x_surf[0]*jacobtot_inv_surf[2]*hamil[12]-0.4330127018922193*jacobtot_inv_surf[0]*b_x_surf[2]*hamil[12]-0.4330127018922193*b_x_surf[0]*jacobtot_inv_surf[3]*hamil[9]-0.4330127018922193*jacobtot_inv_surf[0]*b_x_surf[3]*hamil[9]-0.4330127018922193*b_x_surf[1]*jacobtot_inv_surf[2]*hamil[9]-0.4330127018922193*jacobtot_inv_surf[1]*b_x_surf[2]*hamil[9])*rdz2+(0.7794228634059945*b_z_surf[2]*jacobtot_inv_surf[3]*hamil[12]+0.7794228634059945*jacobtot_inv_surf[2]*b_z_surf[3]*hamil[12]+0.4330127018922193*b_z_surf[0]*jacobtot_inv_surf[1]*hamil[12]+0.4330127018922193*jacobtot_inv_surf[0]*b_z_surf[1]*hamil[12]+0.4330127018922193*b_z_surf[0]*jacobtot_inv_surf[3]*hamil[8]+0.4330127018922193*jacobtot_inv_surf[0]*b_z_surf[3]*hamil[8]+0.4330127018922193*b_z_surf[1]*jacobtot_inv_surf[2]*hamil[8]+0.4330127018922193*jacobtot_inv_surf[1]*b_z_surf[2]*hamil[8])*rdx2)/q_; 
  alphaL[16] = ((0.4330127018922193*jacobtot_inv[1]*b_x[5]*hamil[16]+0.4330127018922193*jacobtot_inv[0]*b_x[3]*hamil[16])*rdz2+(-(0.4330127018922193*jacobtot_inv[3]*b_z[5]*hamil[16])-0.4330127018922193*jacobtot_inv[0]*b_z[1]*hamil[16])*rdx2)/q_; 
  alphaL[17] = ((0.43301270189221935*jacobtot_inv[0]*b_x[5]*hamil[16]+0.43301270189221935*jacobtot_inv[1]*b_x[3]*hamil[16])*rdz2+(-(0.43301270189221935*b_z[5]*jacobtot_inv[5]*hamil[16])-0.43301270189221935*b_z[1]*jacobtot_inv[1]*hamil[16])*rdx2)/q_; 
  alphaL[18] = ((0.43301270189221935*b_x[5]*jacobtot_inv[5]*hamil[16]+0.43301270189221935*b_x[3]*jacobtot_inv[3]*hamil[16])*rdz2+(-(0.43301270189221935*jacobtot_inv[0]*b_z[5]*hamil[16])-0.43301270189221935*b_z[1]*jacobtot_inv[3]*hamil[16])*rdx2)/q_; 
  alphaL[20] = ((0.4330127018922193*b_x[3]*jacobtot_inv[5]*hamil[16]+0.4330127018922193*jacobtot_inv[3]*b_x[5]*hamil[16])*rdz2+(-(0.4330127018922193*b_z[1]*jacobtot_inv[5]*hamil[16])-0.4330127018922193*jacobtot_inv[1]*b_z[5]*hamil[16])*rdx2)/q_; 

  int const_sgn_alpha_surf = 1;  
  
  if (0.22360679774997858*alphaL[20]-0.22360679774997858*(alphaL[18]+alphaL[17])+0.22360679774997858*alphaL[16]-0.25*alphaL[12]-0.33541019662496785*alphaL[11]+0.25*(alphaL[9]+alphaL[8])+0.33541019662496785*(alphaL[7]+alphaL[6])+0.25*alphaL[5]-0.25*alphaL[4]-0.33541019662496785*alphaL[3]-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[0] = 1.0; 
  else  
    sgn_alpha_surfL[0] = -1.0; 
  
  if (-(0.2795084971874732*alphaL[20])+0.2795084971874732*(alphaL[18]+alphaL[17])-0.2795084971874732*alphaL[16]-0.25*alphaL[12]+0.25*(alphaL[9]+alphaL[8]+alphaL[5])-0.25*(alphaL[4]+alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[1] = 1.0; 
  else  
    sgn_alpha_surfL[1] = -1.0; 
  
  if (sgn_alpha_surfL[1] == sgn_alpha_surfL[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.22360679774997858*alphaL[20]-0.22360679774997858*(alphaL[18]+alphaL[17])+0.22360679774997858*alphaL[16]-0.25*alphaL[12]+0.33541019662496785*alphaL[11]+0.25*(alphaL[9]+alphaL[8])-0.33541019662496785*(alphaL[7]+alphaL[6])+0.25*alphaL[5]-0.25*alphaL[4]+0.33541019662496785*alphaL[3]-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[2] = 1.0; 
  else  
    sgn_alpha_surfL[2] = -1.0; 
  
  if (sgn_alpha_surfL[2] == sgn_alpha_surfL[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.22360679774997858*alphaL[20]-0.22360679774997858*(alphaL[18]+alphaL[17])+0.22360679774997858*alphaL[16]+0.25*alphaL[12]-0.33541019662496785*alphaL[11]-0.25*(alphaL[9]+alphaL[8])+0.33541019662496785*(alphaL[7]+alphaL[6])+0.25*(alphaL[5]+alphaL[4])-0.33541019662496785*alphaL[3]-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[3] = 1.0; 
  else  
    sgn_alpha_surfL[3] = -1.0; 
  
  if (sgn_alpha_surfL[3] == sgn_alpha_surfL[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.2795084971874732*alphaL[20])+0.2795084971874732*(alphaL[18]+alphaL[17])-0.2795084971874732*alphaL[16]+0.25*alphaL[12]-0.25*(alphaL[9]+alphaL[8])+0.25*(alphaL[5]+alphaL[4])-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[4] = 1.0; 
  else  
    sgn_alpha_surfL[4] = -1.0; 
  
  if (sgn_alpha_surfL[4] == sgn_alpha_surfL[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.22360679774997858*alphaL[20]-0.22360679774997858*(alphaL[18]+alphaL[17])+0.22360679774997858*alphaL[16]+0.25*alphaL[12]+0.33541019662496785*alphaL[11]-0.25*(alphaL[9]+alphaL[8])-0.33541019662496785*(alphaL[7]+alphaL[6])+0.25*(alphaL[5]+alphaL[4])+0.33541019662496785*alphaL[3]-0.25*(alphaL[2]+alphaL[1])+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[5] = 1.0; 
  else  
    sgn_alpha_surfL[5] = -1.0; 
  
  if (sgn_alpha_surfL[5] == sgn_alpha_surfL[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*alphaL[20])+0.22360679774997858*alphaL[18]-0.22360679774997858*alphaL[17]+0.22360679774997858*alphaL[16]+0.25*alphaL[12]+0.33541019662496785*alphaL[11]-0.25*alphaL[9]+0.25*alphaL[8]-0.33541019662496785*alphaL[7]+0.33541019662496785*alphaL[6]-0.25*(alphaL[5]+alphaL[4])-0.33541019662496785*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[6] = 1.0; 
  else  
    sgn_alpha_surfL[6] = -1.0; 
  
  if (sgn_alpha_surfL[6] == sgn_alpha_surfL[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL[20]-0.2795084971874732*alphaL[18]+0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]+0.25*alphaL[12]-0.25*alphaL[9]+0.25*alphaL[8]-0.25*(alphaL[5]+alphaL[4])+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[7] = 1.0; 
  else  
    sgn_alpha_surfL[7] = -1.0; 
  
  if (sgn_alpha_surfL[7] == sgn_alpha_surfL[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*alphaL[20])+0.22360679774997858*alphaL[18]-0.22360679774997858*alphaL[17]+0.22360679774997858*alphaL[16]+0.25*alphaL[12]-0.33541019662496785*alphaL[11]-0.25*alphaL[9]+0.25*alphaL[8]+0.33541019662496785*alphaL[7]-0.33541019662496785*alphaL[6]-0.25*(alphaL[5]+alphaL[4])+0.33541019662496785*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[8] = 1.0; 
  else  
    sgn_alpha_surfL[8] = -1.0; 
  
  if (sgn_alpha_surfL[8] == sgn_alpha_surfL[7]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*alphaL[20])+0.22360679774997858*alphaL[18]-0.22360679774997858*alphaL[17]+0.22360679774997858*alphaL[16]-0.25*alphaL[12]+0.33541019662496785*alphaL[11]+0.25*alphaL[9]-0.25*alphaL[8]-0.33541019662496785*alphaL[7]+0.33541019662496785*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]-0.33541019662496785*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[9] = 1.0; 
  else  
    sgn_alpha_surfL[9] = -1.0; 
  
  if (sgn_alpha_surfL[9] == sgn_alpha_surfL[8]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*alphaL[20]-0.2795084971874732*alphaL[18]+0.2795084971874732*alphaL[17]-0.2795084971874732*alphaL[16]-0.25*alphaL[12]+0.25*alphaL[9]-0.25*(alphaL[8]+alphaL[5])+0.25*(alphaL[4]+alphaL[2])-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[10] = 1.0; 
  else  
    sgn_alpha_surfL[10] = -1.0; 
  
  if (sgn_alpha_surfL[10] == sgn_alpha_surfL[9]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*alphaL[20])+0.22360679774997858*alphaL[18]-0.22360679774997858*alphaL[17]+0.22360679774997858*alphaL[16]-0.25*alphaL[12]-0.33541019662496785*alphaL[11]+0.25*alphaL[9]-0.25*alphaL[8]+0.33541019662496785*alphaL[7]-0.33541019662496785*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]+0.33541019662496785*alphaL[3]+0.25*alphaL[2]-0.25*alphaL[1]+0.25*alphaL[0] > 0.) 
    sgn_alpha_surfL[11] = 1.0; 
  else  
    sgn_alpha_surfL[11] = -1.0; 
  
  if (sgn_alpha_surfL[11] == sgn_alpha_surfL[10]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*(alphaL[20]+alphaL[18]))+0.22360679774997858*(alphaL[17]+alphaL[16])+0.25*alphaL[12]+0.33541019662496785*alphaL[11]+0.25*alphaL[9]-0.25*alphaL[8]+0.33541019662496785*alphaL[7]-0.33541019662496785*alphaL[6]-0.25*(alphaL[5]+alphaL[4])-0.33541019662496785*alphaL[3]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[12] = 1.0; 
  else  
    sgn_alpha_surfL[12] = -1.0; 
  
  if (sgn_alpha_surfL[12] == sgn_alpha_surfL[11]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*(alphaL[20]+alphaL[18])-0.2795084971874732*(alphaL[17]+alphaL[16])+0.25*(alphaL[12]+alphaL[9])-0.25*(alphaL[8]+alphaL[5]+alphaL[4]+alphaL[2])+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[13] = 1.0; 
  else  
    sgn_alpha_surfL[13] = -1.0; 
  
  if (sgn_alpha_surfL[13] == sgn_alpha_surfL[12]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*(alphaL[20]+alphaL[18]))+0.22360679774997858*(alphaL[17]+alphaL[16])+0.25*alphaL[12]-0.33541019662496785*alphaL[11]+0.25*alphaL[9]-0.25*alphaL[8]-0.33541019662496785*alphaL[7]+0.33541019662496785*alphaL[6]-0.25*(alphaL[5]+alphaL[4])+0.33541019662496785*alphaL[3]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[14] = 1.0; 
  else  
    sgn_alpha_surfL[14] = -1.0; 
  
  if (sgn_alpha_surfL[14] == sgn_alpha_surfL[13]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*(alphaL[20]+alphaL[18]))+0.22360679774997858*(alphaL[17]+alphaL[16])-0.25*alphaL[12]+0.33541019662496785*alphaL[11]-0.25*alphaL[9]+0.25*alphaL[8]+0.33541019662496785*alphaL[7]-0.33541019662496785*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]-0.33541019662496785*alphaL[3]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[15] = 1.0; 
  else  
    sgn_alpha_surfL[15] = -1.0; 
  
  if (sgn_alpha_surfL[15] == sgn_alpha_surfL[14]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.2795084971874732*(alphaL[20]+alphaL[18])-0.2795084971874732*(alphaL[17]+alphaL[16])-0.25*(alphaL[12]+alphaL[9])+0.25*alphaL[8]-0.25*alphaL[5]+0.25*alphaL[4]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[16] = 1.0; 
  else  
    sgn_alpha_surfL[16] = -1.0; 
  
  if (sgn_alpha_surfL[16] == sgn_alpha_surfL[15]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.22360679774997858*(alphaL[20]+alphaL[18]))+0.22360679774997858*(alphaL[17]+alphaL[16])-0.25*alphaL[12]-0.33541019662496785*alphaL[11]-0.25*alphaL[9]+0.25*alphaL[8]-0.33541019662496785*alphaL[7]+0.33541019662496785*alphaL[6]-0.25*alphaL[5]+0.25*alphaL[4]+0.33541019662496785*alphaL[3]-0.25*alphaL[2]+0.25*(alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[17] = 1.0; 
  else  
    sgn_alpha_surfL[17] = -1.0; 
  
  if (sgn_alpha_surfL[17] == sgn_alpha_surfL[16]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.22360679774997858*(alphaL[20]+alphaL[18]+alphaL[17]+alphaL[16])-0.25*alphaL[12]-0.33541019662496785*alphaL[11]-0.25*(alphaL[9]+alphaL[8])-0.33541019662496785*(alphaL[7]+alphaL[6])+0.25*alphaL[5]-0.25*alphaL[4]-0.33541019662496785*alphaL[3]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[18] = 1.0; 
  else  
    sgn_alpha_surfL[18] = -1.0; 
  
  if (sgn_alpha_surfL[18] == sgn_alpha_surfL[17]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.2795084971874732*(alphaL[20]+alphaL[18]+alphaL[17]+alphaL[16]))-0.25*(alphaL[12]+alphaL[9]+alphaL[8])+0.25*alphaL[5]-0.25*alphaL[4]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[19] = 1.0; 
  else  
    sgn_alpha_surfL[19] = -1.0; 
  
  if (sgn_alpha_surfL[19] == sgn_alpha_surfL[18]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.22360679774997858*(alphaL[20]+alphaL[18]+alphaL[17]+alphaL[16])-0.25*alphaL[12]+0.33541019662496785*alphaL[11]-0.25*(alphaL[9]+alphaL[8])+0.33541019662496785*(alphaL[7]+alphaL[6])+0.25*alphaL[5]-0.25*alphaL[4]+0.33541019662496785*alphaL[3]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[20] = 1.0; 
  else  
    sgn_alpha_surfL[20] = -1.0; 
  
  if (sgn_alpha_surfL[20] == sgn_alpha_surfL[19]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.22360679774997858*(alphaL[20]+alphaL[18]+alphaL[17]+alphaL[16])+0.25*alphaL[12]-0.33541019662496785*alphaL[11]+0.25*(alphaL[9]+alphaL[8])-0.33541019662496785*(alphaL[7]+alphaL[6])+0.25*(alphaL[5]+alphaL[4])-0.33541019662496785*alphaL[3]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[21] = 1.0; 
  else  
    sgn_alpha_surfL[21] = -1.0; 
  
  if (sgn_alpha_surfL[21] == sgn_alpha_surfL[20]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.25*(alphaL[12]+alphaL[9]+alphaL[8]+alphaL[5]+alphaL[4]+alphaL[2]+alphaL[1]+alphaL[0])-0.2795084971874732*(alphaL[20]+alphaL[18]+alphaL[17]+alphaL[16]) > 0.) 
    sgn_alpha_surfL[22] = 1.0; 
  else  
    sgn_alpha_surfL[22] = -1.0; 
  
  if (sgn_alpha_surfL[22] == sgn_alpha_surfL[21]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.22360679774997858*(alphaL[20]+alphaL[18]+alphaL[17]+alphaL[16])+0.25*alphaL[12]+0.33541019662496785*alphaL[11]+0.25*(alphaL[9]+alphaL[8])+0.33541019662496785*(alphaL[7]+alphaL[6])+0.25*(alphaL[5]+alphaL[4])+0.33541019662496785*alphaL[3]+0.25*(alphaL[2]+alphaL[1]+alphaL[0]) > 0.) 
    sgn_alpha_surfL[23] = 1.0; 
  else  
    sgn_alpha_surfL[23] = -1.0; 
  
  if (sgn_alpha_surfL[23] == sgn_alpha_surfL[22]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 

} 
