#include <gkyl_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_alpha_surfy_2x2v_ser_p2(const double *w, const double *dxv, const double q_, const double m_, const double *bmag, const double *jacobtot_inv, const double *cmag, const double *b_i, const double *phi, const double *Bstar_Bmag, double* GKYL_RESTRICT alpha_surf) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // q_,m_: species charge and mass.
  // bmag: magnetic field amplitude.
  // jacobtot_inv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // cmag: coefficient multiplying parallel gradient.
  // b_i: covariant components of the field aligned unit vector.
  // phi: electrostatic potential.
  // Bstar_Bmag: Bstar/Bmag volume expansion, pre-computed time-independent part.
  // alpha_surf: output surface phase space flux in each direction (cdim + 1 components).
  //             Note: Each cell owns their *lower* edge surface evaluation (evaluated at -1).

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wy = w[1];
  double rdy2 = 2.0/dxv[1];
  double wvpar = w[2];
  double rdvpar2 = 2.0/dxv[2];
  double wmu = w[3];
  double rdmu2 = 2.0/dxv[3];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wySq = w[1]*w[1];
  double rdy2Sq = rdy2*rdy2;
  double wvparSq = w[2]*w[2];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[3]*w[3];
  double rdmu2Sq = rdmu2*rdmu2;

  const double *b_x = &b_i[0];
  const double *b_y = &b_i[8];
  const double *b_z = &b_i[16];

  const double *BstarYdBmag = &Bstar_Bmag[8]; 

  double hamil[48] = {0.}; 
  hamil[0] = 2.0*m_*wvparSq+2.0*bmag[0]*wmu+(0.6666666666666666*m_)/rdvpar2Sq+2.0*phi[0]*q_; 
  hamil[1] = 2.0*bmag[1]*wmu+2.0*phi[1]*q_; 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[11] = 2.0*bmag[4]*wmu+2.0*phi[4]*q_; 
  hamil[12] = 2.0*phi[5]*q_; 
  hamil[13] = (0.5962847939999438*m_)/rdvpar2Sq; 
  hamil[19] = 2.0*phi[6]*q_; 
  hamil[20] = 2.0*phi[7]*q_; 
  hamil[25] = (1.154700538379251*bmag[4])/rdmu2; 

  double *alphaL1 = &alpha_surf[20];
  alphaL1[0] = ((0.6846531968814574*(b_z[4]*jacobtot_inv[4]+b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*hamil[20]+((-1.060660171779821*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4]))-1.185854122563142*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*hamil[19]+(0.6123724356957944*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+0.6846531968814573*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*hamil[11]+(b_z[4]*jacobtot_inv[4]+b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*(0.3061862178478971*hamil[1]-0.5303300858899105*hamil[5]))*rdx2)/q_+((1.369306393762915*BstarYdBmag[2]*hamil[13]+0.6123724356957944*BstarYdBmag[0]*hamil[3])*rdvpar2)/m_; 
  alphaL1[1] = (((0.6123724356957944*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+0.6846531968814574*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*hamil[20]+((-1.060660171779821*(b_z[0]*jacobtot_inv[4]+jacobtot_inv[0]*b_z[4]))-1.86348504974208*b_z[4]*jacobtot_inv[4]-2.134537420613654*b_z[1]*jacobtot_inv[1]-1.185854122563142*b_z[0]*jacobtot_inv[0])*hamil[19]+(0.6123724356957944*(b_z[0]*jacobtot_inv[4]+jacobtot_inv[0]*b_z[4])+1.075883595099433*b_z[4]*jacobtot_inv[4]+1.232375754386623*b_z[1]*jacobtot_inv[1]+0.6846531968814573*b_z[0]*jacobtot_inv[0])*hamil[11]+((-0.4743416490252568*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4]))-0.5303300858899105*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*hamil[5]+hamil[1]*(0.273861278752583*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+0.3061862178478971*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1])))*rdx2)/q_+((1.369306393762915*BstarYdBmag[3]*hamil[13]+0.6123724356957944*BstarYdBmag[1]*hamil[3])*rdvpar2)/m_; 
  alphaL1[2] = ((1.369306393762915*BstarYdBmag[0]*hamil[13]+0.6123724356957944*BstarYdBmag[2]*hamil[3])*rdvpar2)/m_; 
  alphaL1[3] = (((0.6123724356957944*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+0.6846531968814574*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*hamil[25]+0.3061862178478971*(b_z[4]*jacobtot_inv[4]+b_z[1]*jacobtot_inv[1]+b_z[0]*jacobtot_inv[0])*hamil[8])*rdx2)/q_; 
  alphaL1[4] = ((1.369306393762915*BstarYdBmag[1]*hamil[13]+0.6123724356957944*BstarYdBmag[3]*hamil[3])*rdvpar2)/m_; 
  alphaL1[5] = (((0.6123724356957944*(b_z[0]*jacobtot_inv[4]+jacobtot_inv[0]*b_z[4])+1.075883595099433*b_z[4]*jacobtot_inv[4]+1.232375754386623*b_z[1]*jacobtot_inv[1]+0.6846531968814574*b_z[0]*jacobtot_inv[0])*hamil[25]+(0.273861278752583*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+0.3061862178478971*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*hamil[8])*rdx2)/q_; 
  alphaL1[7] = (((0.6846531968814574*(b_z[0]*jacobtot_inv[4]+jacobtot_inv[0]*b_z[4])+0.4374088826398531*b_z[4]*jacobtot_inv[4]+0.6123724356957944*b_z[1]*jacobtot_inv[1])*hamil[20]+((-1.86348504974208*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4]))-1.060660171779821*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*hamil[19]+(1.075883595099433*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+0.6123724356957944*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*hamil[11]+((-0.5303300858899105*(b_z[0]*jacobtot_inv[4]+jacobtot_inv[0]*b_z[4]))-0.3388154635894691*b_z[4]*jacobtot_inv[4]-0.4743416490252568*b_z[1]*jacobtot_inv[1])*hamil[5]+hamil[1]*(0.3061862178478971*(b_z[0]*jacobtot_inv[4]+jacobtot_inv[0]*b_z[4])+0.1956151991089878*b_z[4]*jacobtot_inv[4]+0.273861278752583*b_z[1]*jacobtot_inv[1]))*rdx2)/q_+((1.369306393762915*BstarYdBmag[6]*hamil[13]+0.6123724356957944*hamil[3]*BstarYdBmag[4])*rdvpar2)/m_; 
  alphaL1[8] = (1.224744871391589*BstarYdBmag[2]*hamil[13]*rdvpar2)/m_; 
  alphaL1[11] = ((1.369306393762915*BstarYdBmag[4]*hamil[13]+0.6123724356957944*hamil[3]*BstarYdBmag[6])*rdvpar2)/m_; 
  alphaL1[12] = (1.224744871391589*BstarYdBmag[3]*hamil[13]*rdvpar2)/m_; 
  alphaL1[13] = (((1.075883595099433*(b_z[1]*jacobtot_inv[4]+jacobtot_inv[1]*b_z[4])+0.6123724356957944*(b_z[0]*jacobtot_inv[1]+jacobtot_inv[0]*b_z[1]))*hamil[25]+(0.3061862178478971*(b_z[0]*jacobtot_inv[4]+jacobtot_inv[0]*b_z[4])+0.1956151991089878*b_z[4]*jacobtot_inv[4]+0.273861278752583*b_z[1]*jacobtot_inv[1])*hamil[8])*rdx2)/q_; 

} 